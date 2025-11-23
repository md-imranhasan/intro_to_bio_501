
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



