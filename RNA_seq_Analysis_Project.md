
# Systems biological assessment of immunity COVID-19 infections (GSE152418)

✔ Organism

Homo sapiens

✔ Tissue / Cell Type

Peripheral Blood Mononuclear Cells (PBMCs)

10 COVID-19 samples (case)
``` text
SRR12007823  
SRR12007824  
SRR12007825  
SRR12007826  
SRR12007827  
SRR12007828  
SRR12007829  
SRR12007830  
SRR12007831  
SRR12007832  
```

10 healthy samples (control)
``` text
SRR12007855  
SRR12007856  
SRR12007857  
SRR12007858  
SRR12007859  
SRR12007860  
SRR12007861  
SRR12007862  
SRR12007863  
SRR12007864  
```


✔ Total Used in Analysis

N = 20 samples


📌 What You Achieved So Far

You have completed:
``` text
✔ FASTQ download
✔ QC
✔ Trimming
✔ HISAT2 alignment
✔ Sorting / MarkDuplicates
✔ rRNA removal
✔ MT removal
✔ Final cleaned BAM
✔ featureCounts
✔ DESeq2 differential expression
✔ Volcano plots
✔ Heatmaps
✔ Exported gene lists for pathway analysis
```
This is a full RNA-seq pipeline.

# COVID-19 Bulk RNA-Seq Pipeline 

This repository contains a full, reproducible RNA-seq pipeline used to analyze
COVID-19 PBMC RNA-seq data from GEO/SRA using:

- HISAT2 (alignment)
- SAMtools + Picard (sorting, deduplication, QC)
- rRNA + MT removal
- featureCounts (gene quantification)
- DESeq2 (differential expression)
- Volcano plots & heatmaps

---

## 📁 Repository Structure

``` text
covid_rnaseq_pipeline/
│
├── README.md
│
├── data/
│   ├── raw_fastq/            # SRR FASTQ files
│   ├── trim_fastp/           # trimmed FASTQs
│   ├── hisat2_bam/           # alignments
│   ├── sorted_bam/
│   ├── markdup_bam/
│   ├── clean_noMT/
│   ├── clean_no_rRNA/
│   ├── final_bam_clean/
│   ├── counts/
│   └── annotation/           # GTF, rRNA.bed, genome index paths
│
├── scripts/
│   ├── 01_download_fastq.sh
│   ├── 02_fastp_trim.sh
│   ├── 03_hisat2_align.sh
│   ├── 04_sort_markdup.sh
│   ├── 05_remove_MT_rRNA.sh
│   ├── 06_featureCounts.sh
│   ├── 07_deseq2_analysis.R
│   ├── 08_volcano_heatmap_plots.R
│   └── utils/                # small helper functions
│
├── results/
│   ├── alignment_summary.tsv
│   ├── gene_counts/
│   ├── deseq2/
│   ├── plots/
│ 
│
└── environment/
    ├── conda_env.yaml
    ├── sessionInfo.txt
    └── software_versions.txt

```



---

## 🔬 Processing Workflow

1. **Download FASTQ files from SRA**
2. **Trim adapters using fastp**
3. **Align reads using HISAT2 (GRCh38)**
4. **Sort + mark duplicates**
5. **Remove mitochondrial and rRNA reads**
6. **Index final BAMs**
7. **Generate gene-level counts using featureCounts**
8. **Run DESeq2**
9. **Visualize results (volcano plots, heatmaps)**
10. **Export DEGs for pathway analysis**

---

## 🧬 Requirements

- HISAT2
- SAMtools
- Picard Tools
- fastp
- subread (featureCounts)
- R + DESeq2 + tidyverse + pheatmap + ggrepel + rtracklayer

Install via:

```bash
conda env create -f environment/conda_env.yaml

📊 Final Outputs

Differential expression tables (Up / Down / All)

Volcano plots

Heatmaps

Gene lists for pathway analysis (GSEA, Enrichr, etc.)
```

### **01_download_fastq.sh**
```bash
#!/bin/bash
while read SRR; do
    fasterq-dump --split-files --threads 8 --outdir data/raw_fastq $SRR
done < SRR_Acc_List.txt
```

02_fastp_trim.sh
```bash
#!/bin/bash
mkdir -p data/trim_fastp data/qc_fastp
for f in data/raw_fastq/*.fastq; do
    base=$(basename "$f" .fastq)
    fastp -i $f -o data/trim_fastp/${base}.trimmed.fastq \
          -h data/qc_fastp/${base}.html \
          -j data/qc_fastp/${base}.json \
          -w 8
done
```

03_hisat2_align.sh
```
#!/bin/bash
mkdir -p data/hisat2_bam data/hisat2_logs
for f in data/trim_fastp/*.trimmed.fastq; do
    r=$(basename $f .trimmed.fastq)
    hisat2 -p 16 --dta -x data/annotation/grch38/genome \
           -U $f \
           -S data/hisat2_bam/${r}.sam \
           --summary-file data/hisat2_logs/${r}.summary.txt
done
```

04_sort_markdup.sh
```
#!/bin/bash
mkdir -p data/sorted_bam data/markdup_bam
for f in data/hisat2_bam/*.sam; do
    base=$(basename $f .sam)
    samtools sort -@ 8 -o data/sorted_bam/${base}.sorted.bam $f
    picard MarkDuplicates I=data/sorted_bam/${base}.sorted.bam \
                          O=data/markdup_bam/${base}.markdup.bam \
                          M=data/markdup_bam/${base}.metrics.txt \
                          REMOVE_DUPLICATES=true CREATE_INDEX=true
done
```

05_remove_MT_rRNA.sh

```
#!/bin/bash
mkdir -p data/clean_noMT data/clean_no_rRNA

# Remove MT reads
for f in data/markdup_bam/*.markdup.bam; do
    b=$(basename $f .markdup.bam)
    samtools view -h $f | grep -v -w "MT" | \
        samtools view -b -o data/clean_noMT/${b}.noMT.bam
    samtools index data/clean_noMT/${b}.noMT.bam
done

# Remove rRNA reads
for f in data/clean_noMT/*.noMT.bam; do
    b=$(basename $f .noMT.bam)
    samtools view -h -L data/annotation/grch38_rRNA.bed \
                  -U data/clean_no_rRNA/${b}.noMT.noRrna.bam $f > /dev/null
    samtools index data/clean_no_rRNA/${b}.noMT.noRrna.bam
done

```


06_featureCounts.sh

```
#!/bin/bash
featureCounts -T 12 -s 0 -t gene -g gene_id \
    -a data/annotation/Homo_sapiens.GRCh38.110.gtf \
    -o data/counts/caRNA_gene_counts.txt \
    data/final_bam_clean/*.clean.bam
```

07_deseq2_analysis.R — clean, simplified

```bash
############################################################
#E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/

# Load libraries
############################################################
############################################################
## 0. Libraries
############################################################
library(DESeq2)
library(tidyverse)
library(rtracklayer)

############################################################
## 1. Load featureCounts matrix
############################################################
# Path to your counts matrix (change if needed)
counts_file <- "E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/caRNA_gene_counts.matrix.csv"

raw <- read.csv(counts_file, header = TRUE, check.names = FALSE)

# Separate annotation columns from counts
gene_annot_counts <- raw[, c("Geneid","Chr","Start","End","Strand","Length")]

count_data <- raw[, !(names(raw) %in% c("Chr","Start","End","Strand","Length"))]

# Set rownames as Geneid
rownames(count_data) <- count_data$Geneid
count_data$Geneid <- NULL

############################################################
## 2. Clean sample names (strip paths & .clean.bam)
############################################################
clean_names <- basename(colnames(count_data))             # remove folder path
clean_names <- sub("\\.clean\\.bam$", "", clean_names)    # remove .clean.bam
colnames(count_data) <- clean_names

print(colnames(count_data))

############################################################
## 3. Define case/control samples (your exact lists)
############################################################
case_samples <- c(
  "SRR12007823","SRR12007824","SRR12007825","SRR12007826","SRR12007827",
  "SRR12007828","SRR12007829","SRR12007830","SRR12007831","SRR12007832"
)

control_samples <- c(
  "SRR12007855","SRR12007856","SRR12007857","SRR12007858","SRR12007859",
  "SRR12007860","SRR12007861","SRR12007862","SRR12007863","SRR12007864"
)

samples <- colnames(count_data)

condition <- ifelse(samples %in% case_samples, "COVID",
                    ifelse(samples %in% control_samples, "Control", NA))

# Quick sanity check
print(condition)

coldata <- data.frame(
  row.names = samples,
  condition = factor(condition, levels = c("Control", "COVID"))
)

############################################################
## 4. Build DESeq2 object
############################################################
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_data)),
  colData   = coldata,
  design    = ~ condition
)

# Filter out very low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

############################################################
## 5. Run DESeq2
############################################################
dds <- DESeq(dds)

############################################################
## 6. Get results (COVID vs Control) + shrinkage (no apeglm)
############################################################
# Basic results
res <- results(dds, contrast = c("condition", "COVID", "Control"))

# Optional shrinkage using DESeq2's internal method (no apeglm)
res_shrunk <- lfcShrink(dds,
                        coef = "condition_COVID_vs_Control",
                        type = "normal")

# We'll use the shrunken results for export/plots
res_use <- res_shrunk

# Order by adjusted p-value
res_ordered <- res_use[order(res_use$padj), ]

# Convert to data.frame and add Geneid column
res_df <- as.data.frame(res_ordered)
res_df$Geneid <- rownames(res_df)

############################################################
## 7. Add gene annotation from counts + gene names from GTF
############################################################
# 7a. Merge with featureCounts annotation (Chr, Start, End, Strand, Length)
res_merged <- gene_annot_counts %>%
  inner_join(res_df, by = "Geneid")

# 7b. Load GTF and extract gene_id / gene_name
gtf_path <- "E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/Homo_sapiens.GRCh38.110.gtf"  # adjust path if needed
gtf <- import(gtf_path)
gtf_df <- as.data.frame(gtf)

gene_annot_gtf <- gtf_df %>%
  dplyr::filter(type == "gene") %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct()

# 7c. Add gene_name
res_with_names <- res_merged %>%
  dplyr::rename(gene_id = Geneid) %>%
  left_join(gene_annot_gtf, by = "gene_id") %>%
  relocate(gene_id, gene_name, .before = Chr)

# Check
head(res_with_names)

############################################################
## 8. Flag DEGs (padj < 0.05 & |log2FC| > 1)
############################################################
res_with_names <- res_with_names %>%
  mutate(
    is_sig = !is.na(padj) & padj < 0.05,
    DEG_direction = case_when(
      is_sig & log2FoldChange >  1  ~ "Up",
      is_sig & log2FoldChange < -1  ~ "Down",
      TRUE                          ~ "NS"
    )
  )

table(res_with_names$DEG_direction)

# All DEGs
deg_all <- res_with_names %>%
  filter(DEG_direction %in% c("Up", "Down")) %>%
  arrange(padj)

deg_up <- deg_all %>% filter(DEG_direction == "Up")
deg_down <- deg_all %>% filter(DEG_direction == "Down")

############################################################
## 9. Save tables
############################################################
write.csv(res_with_names,
          file = "E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/DESeq2_COVID_vs_Control_with_gene_names_all.csv",
          row.names = FALSE)

write.csv(deg_all,
          file = "E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/DESeq2_COVID_vs_Control_DEGs_all.csv",
          row.names = FALSE)

write.csv(deg_up,
          file = "E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/DESeq2_COVID_vs_Control_DEGs_up.csv",
          row.names = FALSE)

write.csv(deg_down,
          file = "E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/DESeq2_COVID_vs_Control_DEGs_down.csv",
          row.names = FALSE)

############################################################
## 10. Plots (MA & PCA)
############################################################
# MA plot from shrunken results
pdf("E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/MA_plot_COVID_vs_Control.pdf")
plotMA(res_use, ylim = c(-5, 5))
dev.off()

# PCA plot
vsd <- vst(dds, blind = FALSE)
pdf("E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/PCA_plot_COVID_vs_Control.pdf")
plotPCA(vsd, intgroup = "condition")
dev.off()

############################################################
## Done
############################################################

library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)


volc_df <- res_with_names %>%
  mutate(
    negLog10Padj = -log10(padj),
    status = case_when(
      DEG_direction == "Up"   ~ "Up",
      DEG_direction == "Down" ~ "Down",
      TRUE                    ~ "NS"
    )
  )

ggplot(volc_df, aes(x = log2FoldChange, y = negLog10Padj, color = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) +
  labs(x = "log2 Fold Change (COVID vs Control)",
       y = "-log10 adj p-value",
       color = "DEG") +
  theme_minimal()










#Label only the top 20 most significant DEGs (to keep it readable).

# same data frame as above
volc_df <- volc_df %>%
  mutate(label = ifelse(DEG_direction %in% c("Up","Down"), gene_name, NA))

# pick top 20 by padj among significant
top_labels <- volc_df %>%
  filter(DEG_direction %in% c("Up","Down")) %>%
  arrange(padj) %>%
  slice(1:20)

ggplot(volc_df, aes(x = log2FoldChange, y = negLog10Padj, color = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = top_labels,
    aes(label = gene_name),
    size = 3,
    max.overlaps = 100
  ) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) +
  labs(x = "log2 Fold Change (COVID vs Control)",
       y = "-log10 adj p-value",
       color = "DEG") +
  theme_minimal()












library(pheatmap)
library(SummarizedExperiment)

# pick DEGs, ordered by padj
deg_all <- res_with_names %>%
  filter(DEG_direction %in% c("Up","Down")) %>%
  arrange(padj)

# choose top 50 genes
top_deg <- deg_all %>% slice(1:50)

# extract matrix from VST-normalized counts
mat <- assay(vsd)[top_deg$gene_id, ]

# replace row names with gene names
rownames(mat) <- top_deg$gene_name

# sample annotation
ann_col <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])

# heatmap
pheatmap(
  mat,
  annotation_col = ann_col,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  fontsize_row = 6
)


















#🔧 Fix the error and make the heatmap
library(pheatmap)
library(SummarizedExperiment)

# 1) Make sure we are working with a plain data.frame
res_df <- as.data.frame(res_with_names)

# 2) Pick DEGs (Up/Down) ordered by padj
deg_all <- res_df[res_df$DEG_direction %in% c("Up","Down"), ]
deg_all <- deg_all[order(deg_all$padj), ]

# 3) Choose top 50 genes (or fewer if you have less)
n_top <- min(50, nrow(deg_all))
top_deg <- deg_all[1:n_top, ]

# 4) Extract VST-normalized matrix for these genes
mat <- assay(vsd)[top_deg$gene_id, ]

# 5) Use gene names as rownames
rownames(mat) <- top_deg$gene_name

# 6) Sample annotation (COVID vs Control)
ann_col <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])

# 7) Heatmap
pheatmap(
  mat,
  annotation_col = ann_col,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  fontsize_row = 6
)















##############1. Convert results to a plain data.frame
#R code to export ONLY the gene symbols (gene names) for Up and Down DEGs.

res_df <- as.data.frame(res_with_names)


#2. Extract Up / Down DEGs (padj < 0.05 and |log2FC| > 1)
deg_up <- res_df %>%
  filter(DEG_direction == "Up" & !is.na(gene_name)) %>%
  pull(gene_name)

deg_down <- res_df %>%
  filter(DEG_direction == "Down" & !is.na(gene_name)) %>%
  pull(gene_name)



#Save gene symbol lists (for GSEA/Enrichr)
write.table(deg_up,
            file = "E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/DEG_Up_gene_symbols.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(deg_down,
            file = "E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/DEG_Down_gene_symbols.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)





#Optional: export all significant DEGs (mixed up/down)
deg_all <- res_df %>%
  filter(DEG_direction %in% c("Up","Down") & !is.na(gene_name)) %>%
  pull(gene_name)

write.table(deg_all,
            file = "E:/Purdue Fall 2025/Course/BIOL - 59500Q Fall 2025 Introduction to Bioinformatics/Project/RNAseq Project/DEG_All_gene_symbols.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)




```

08_volcano_heatmap_plots.R

```
Volcano plot (Up/Down/NS)

Volcano plot with gene labels (top 20 DEGs)

Heatmap of top 50 DEGs (rows = gene names, columns = samples)
```

volcano_heatmap_plots.R
```
############################################################
## 08_volcano_heatmap_plots.R
## Volcano plots + heatmap for COVID vs Control DEGs
############################################################

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(SummarizedExperiment)

############################################################
## 1. Load DESeq2 results with gene names
############################################################

res_file <- "results/deseq2/DESeq2_results_with_gene_names.csv"
res_with_names <- read.csv(res_file, header = TRUE, check.names = FALSE)

# Make sure we have a plain data.frame
res_df <- as.data.frame(res_with_names)

############################################################
## 2. Rebuild dds + vsd from counts (for heatmap)
############################################################

counts_file <- "data/counts/caRNA_gene_counts.matrix.csv"
raw_counts <- read.csv(counts_file, header = TRUE, check.names = FALSE)

gene_annot_counts <- raw_counts[, c("Geneid","Chr","Start","End","Strand","Length")]
count_data <- raw_counts[, !(names(raw_counts) %in% c("Chr","Start","End","Strand","Length"))]

rownames(count_data) <- count_data$Geneid
count_data$Geneid <- NULL

# Clean sample names
clean_names <- basename(colnames(count_data))
clean_names <- sub("\\.clean\\.bam$", "", clean_names)
colnames(count_data) <- clean_names

# Define case/control (same as Script 07)
case_samples <- c(
  "SRR12007823","SRR12007824","SRR12007825","SRR12007826","SRR12007827",
  "SRR12007828","SRR12007829","SRR12007830","SRR12007831","SRR12007832"
)

control_samples <- c(
  "SRR12007855","SRR12007856","SRR12007857","SRR12007858","SRR12007859",
  "SRR12007860","SRR12007861","SRR12007862","SRR12007863","SRR12007864"
)

samples <- colnames(count_data)

condition <- ifelse(samples %in% case_samples, "COVID",
             ifelse(samples %in% control_samples, "Control", NA))

coldata <- data.frame(
  row.names = samples,
  condition = factor(condition, levels = c("Control", "COVID"))
)

# DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_data)),
  colData   = coldata,
  design    = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

# VST transform for heatmap
vsd <- vst(dds, blind = FALSE)

############################################################
## 3. Volcano Plot – basic (Up/Down/NS)
############################################################

dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

volc_df <- res_df %>%
  mutate(
    negLog10Padj = -log10(padj),
    status = case_when(
      DEG_direction == "Up"   ~ "Up",
      DEG_direction == "Down" ~ "Down",
      TRUE                    ~ "NS"
    )
  )

p_volcano_basic <- ggplot(volc_df, aes(x = log2FoldChange, y = negLog10Padj, color = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) +
  labs(x = "log2 Fold Change (COVID vs Control)",
       y = "-log10 adjusted p-value",
       color = "DEG") +
  theme_minimal()

ggsave("results/plots/volcano_basic.png", p_volcano_basic, width = 6, height = 5, dpi = 300)

############################################################
## 4. Volcano Plot – with gene labels (top 20 DEGs)
############################################################

# Choose top 20 significant DEGs (Up or Down) by padj
top_labels <- volc_df %>%
  filter(DEG_direction %in% c("Up","Down")) %>%
  arrange(padj) %>%
  head(20)

p_volcano_labeled <- ggplot(volc_df, aes(x = log2FoldChange, y = negLog10Padj, color = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = top_labels,
    aes(label = gene_name),
    size = 3,
    max.overlaps = 100
  ) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) +
  labs(x = "log2 Fold Change (COVID vs Control)",
       y = "-log10 adjusted p-value",
       color = "DEG") +
  theme_minimal()

ggsave("results/plots/volcano_labeled_top20.png", p_volcano_labeled, width = 7, height = 6, dpi = 300)

############################################################
## 5. Heatmap – top 50 DEGs (by padj)
############################################################

# Work with plain data.frame
deg_all <- res_df[res_df$DEG_direction %in% c("Up","Down"), ]
deg_all <- deg_all[order(deg_all$padj), ]

n_top <- min(50, nrow(deg_all))
top_deg <- deg_all[1:n_top, ]

# Extract VST-normalized expression
mat <- assay(vsd)[top_deg$gene_id, ]

# Use gene names as rownames
rownames(mat) <- top_deg$gene_name

# Sample annotation
ann_col <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])

pheatmap(
  mat,
  annotation_col = ann_col,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  fontsize_row = 6,
  filename = "results/plots/heatmap_top50_DEGs.png",
  width = 7,
  height = 9
)

############################################################
## Done
############################################################

message("Plots saved to results/plots/:",
        "\n - volcano_basic.png",
        "\n - volcano_labeled_top20.png",
        "\n - heatmap_top50_DEGs.png")

```























<img width="2100" height="1500" alt="volcano_up_down_significant" src="https://github.com/user-attachments/assets/7a0e95a2-cc96-410b-9f2e-728c6e7b1ad7" />















