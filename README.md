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
[![Tutorial PDF Preview](images/pdf_preview.png)](ElbowPlot_600dpi.pdf)
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
seurat_obj <- CreateSeuratObject(counts = counts, project = "Tutorial")
cat("Initial dimensions - Genes:", nrow(seurat_obj), "Cells:", ncol(seurat_obj), "\n")

```


## Gene Name Conversion and Duplicate Handling

```
### Reads GTF annotation file
gtf <- read.table("ref/gencode.v49.annotation.gtf", 
                  sep = "\t", comment.char = "#", quote = "", stringsAsFactors = FALSE)

### Extracts clean gene information from GTF
gene_info <- gtf %>%
  filter(V3 == "gene") %>%
  mutate(
    gene_id = str_extract(V9, 'gene_id "[^"]+"') %>% str_remove_all('gene_id |"'),
    gene_name = str_extract(V9, 'gene_name "[^"]+"') %>% str_remove_all('gene_name |"')
  ) %>%
  select(gene_id, gene_name) %>%
  distinct()

### Removes Ensembl gene version numbers (ENSG00000000003.14 → ENSG00000000003)
rownames(seurat_obj) <- gsub("\\..*$", "", rownames(seurat_obj))
gene_info$gene_id <- gsub("\\..*$", "", gene_info$gene_id)

### Maps gene symbols to Seurat object
gene_symbols <- gene_info$gene_name[match(rownames(seurat_obj), gene_info$gene_id)]
valid <- !is.na(gene_symbols)
rownames(seurat_obj)[valid] <- gene_symbols[valid]

### Handles duplicate gene names
counts_matrix <- GetAssayData(seurat_obj, slot = "counts")
unique_genes <- unique(rownames(seurat_obj))
keep_indices_list <- list()
counter <- 1

for(gene in unique_genes) {
  indices <- which(rownames(seurat_obj) == gene)
  if(length(indices) == 1) {
    keep_indices_list[[counter]] <- indices
  } else {
    expr_sums <- Matrix::rowSums(counts_matrix[indices, ])
    keep_indices_list[[counter]] <- indices[which.max(expr_sums)]
  }
  counter <- counter + 1
}

keep_indices <- unlist(keep_indices_list)
seurat_obj <- seurat_obj[sort(keep_indices), ]
```


## Quality Control
```
### Calculates mitochondrial and ribosomal percentages
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

### Visualizes QC metric distributions
p <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Save as high-resolution PDF
pdf("VlnPlot_qc_600dpi.pdf", width = 10, height = 8, paper = "special")
print(p)
dev.off()

### Filters low-quality cells
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 300 & nFeature_RNA < 2500 &
           nCount_RNA > 500 & nCount_RNA < 10000 &
           percent.mt < 15
)

mtx_q75 <- quantile(seurat_obj$percent.mt, 0.75)
mtx_q25 <- quantile(seurat_obj$percent.mt, 0.25)
mtx_threshold <- mtx_q75 + 1.5 * (mtx_q75 - mtx_q25)
cat("Data-driven mt% threshold:", mtx_threshold, "\n")

```

## Duplicate Gene Removal (Post-QC)
```
gene_counts <- table(rownames(seurat_obj))
duplicated_genes <- names(gene_counts[gene_counts > 1])
seurat_obj <- seurat_obj[!rownames(seurat_obj) %in% duplicated_genes, ]
```

## Normalization and Feature Selection
```
### Normalizes counts using log transformation
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

### Identifies 2,000 most variable genes
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

### Visualizes variable features with top 10 labeled

# Generate the variable feature plot
plot1 <- VariableFeaturePlot(seurat_obj)

# Extract top 10 variable features
top10 <- head(VariableFeatures(seurat_obj), 10)

# Label top 10 variable genes
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

pdf("VariableFeaturePlot_Top10.pdf", width = 6, height = 5, paper = "special")
print(plot2)
dev.off()


### Scales data (zero-mean, unit-variance)

seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
```

## Dimensionality Reduction
```
### Runs PCA on variable features only
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

### Prints top genes for first 5 PCs
print(seurat_obj[["pca"]], dims = 1:10, nfeatures = 10)

# Generate the Elbow Plot
p <- ElbowPlot(seurat_obj)
pdf("ElbowPlot_600dpi.pdf", width = 6, height = 5, paper = "special")
print(p)
dev.off()
```
## Clustering and UMAP
```
### Clusters cells using KNN graph and Louvain algorithm
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

### Generates 2D UMAP visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Generate the UMAP plot with cluster labels
p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
pdf("UMAP_DimPlot_600dpi.pdf", width = 6, height = 5, paper = "special")
print(p)
dev.off()
```
## Marker Gene Identification
```
### Finds cluster-specific markers
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.5,  # More stringent
                                   test.use = "wilcox")  # Explicitly use Wilcoxon

### Extracts top 10 genes per cluster
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

### Generate FeaturePlot for top 9 marker genes
p <- FeaturePlot(seurat_obj, features = top_markers$gene[1:9])
pdf("FeaturePlot_Top9_600dpi.pdf", width = 9, height = 9, paper = "special")
print(p)
dev.off()
```

## Save Results
```
saveRDS(seurat_obj, file = "Tutorial_processed_seurat.rds")
write.csv(cluster_markers, "Tutorial_cluster_markers.csv", row.names = FALSE)
```
## References
```


    GENCODE Human Reference

    Salmon & Alevin-fry Documentation

    Seurat R package

```
##### sastigopaldas05@gmail.com
