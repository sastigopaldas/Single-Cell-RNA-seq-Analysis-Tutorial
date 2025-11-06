# Single-Cell-RNA-seq-Analysis-Tutorial
This repository demonstrates a complete single-cell RNA-seq analysis workflow, from environment setup and reference genome preparation to quantification and R-based downstream analysis with Seurat

## 1 Environment Setup
```
# Create and activate Conda environment
conda create -n Single_Cell_RNA_seq_simpleaf python=3.10 -y
conda activate Single_Cell_RNA_seq_simpleaf

# Install necessary tools
conda install -c bioconda -c conda-forge simpleaf alevin-fry salmon gffread -y

# Check installation
simpleaf --version || true
alevin-fry --version || true
salmon --version || true
gffread --version || true
```
## Reference Preparation

#### Download reference genome and annotation (GENCODE v49)
```
wget -O gencode.v49.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz

wget -O GRCh38.primary_assembly.genome.fa.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz

gunzip -f gencode.v49.annotation.gtf.gz
gunzip -f GRCh38.primary_assembly.genome.fa.gz

# Generate transcriptome FASTA
gffread -w gencode.v49.transcriptome.fa \
  -g GRCh38.primary_assembly.genome.fa gencode.v49.annotation.gtf

# Inspect reference files
head gencode.v49.transcriptome.fa
head GRCh38.primary_assembly.genome.fa
head gencode.v49.annotation.gtf
```
#### Create transcript-to-gene mapping file
```
awk '{
  if ($3 == "transcript") {
    match($0, /gene_id "([^"]+)"/, g);
    match($0, /transcript_id "([^"]+)"/, t);
    if (t[1] && g[1]) print t[1] "\t" g[1];
  }
}' gencode.v49.annotation.gtf > t2g.tsv

head -5 t2g.tsv
```
## Index Building
```
# Build Salmon index for transcriptome
salmon index -t ref/gencode.v49.transcriptome.fa -i grch38_idx -p 4
```

## Quantification with Alevin-fry
```
# Run Salmon Alevin quantification (10x Chromium V3)

salmon alevin -i grch38_idx -p 4 -l IU --chromiumV3 --sketch \
  -1 fastq/sampel_1.fastq.gz \
  -2 fastq/sampel_2.fastq.gz \
  -o sampel_name_map \
  --tgMap ref/t2g.tsv

# Generate barcode permit list
alevin-fry generate-permit-list -d fw -k -i sampel_name_map -o sampel_name_quant

# Collate mapping output
alevin-fry collate -t 4 -i sampel_name_quant -r sampel_name_map

# quantification
alevin-fry quant -t 4 -i sampel_name_quant -o sampel_name_quant \
--tg-map ref/t2g.tsv --resolution cr-like --use-mtx

# Generate barcode permit list
alevin-fry generate-permit-list -d fw -k -i sampel_name_map -o sampel_name_quant

# Collate mapping output
alevin-fry collate -t 4 -i sampel_name_quant -r sampel_name_map

# Quantify gene/cell counts
alevin-fry quant -t 4 -i sampel_name_quant -o sampel_name_quant \
  --tg-map ref/t2g.tsv --resolution cr-like --use-mtx
```


## Output Inspection
```
# Quick checks for matrix and metadata files
head -10 sampel_name_quant/alevin/quants_mat.mtx
head -10 sampel_name_quant/alevin/quants_mat_rows.txt
head -10 sampel_name_quant/alevin/quants_mat_cols.txt

wc -l sampel_name_quant/alevin/quants_mat_rows.txt
wc -l sampel_name_quant/alevin/quants_mat_cols.txt

conda deactivate
```



## R Analysis (Seurat)
### Load required R packages:
```
library(Matrix)
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

# (Optional: install packages if missing)
required_packages <- c("Matrix", "Seurat", "dplyr", "stringr", "ggplot2")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}
```


### Read matrix data and create Seurat object:
```
counts <- ReadMtx(
  mtx = "sampel_name_quant/alevin/quants_mat.mtx",
  features = "sampel_name_quant/alevin/quants_mat_rows.txt",
  cells = "sampel_name_quant/alevin/quants_mat_cols.txt",
  feature.column = 1
)
counts <- t(counts)
seurat_obj <- CreateSeuratObject(counts = counts, project = "pbmc10k")
cat("Initial dimensions - Genes:", nrow(seurat_obj), "Cells:", ncol(seurat_obj), "\n")

```


## Gene Name Conversion & QC
```
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat_obj <- subset(seurat_obj, subset = 
                     nFeature_RNA > 200 & nFeature_RNA < 5000 & 
                     nCount_RNA > 300 & percent.mt < 15)

```


## Normalization to Clustering
```
# Normalize, feature selection & scaling
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(seurat_obj)
top10 <- head(VariableFeatures(seurat_obj), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# Dimensionality reduction
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
ElbowPlot(seurat_obj)

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
```

## Marker Identification
```
# Find marker genes per cluster
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "wilcox")
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

FeaturePlot(seurat_obj, features = top_markers$gene[1:9])
```

## Final Output
```
saveRDS(seurat_obj, file = "pbmc10k_processed_seurat.rds")
write.csv(cluster_markers, "pbmc10k_cluster_markers.csv", row.names = FALSE)
```
## References
```


    GENCODE Human Reference

    Salmon & Alevin-fry Documentation

    Seurat R package

```
