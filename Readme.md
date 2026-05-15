# NRF1 Target Gene GO Enrichment Analysis

A computational pipeline to identify **NRF1 (Nuclear Respiratory Factor 1) target genes** from human promoter sequences and perform **Gene Ontology (GO) enrichment analysis** on the identified targets.

---

## Overview

NRF1 binds a palindromic GC-rich motif (`GCGCNNNGCGC`) in gene promoters and is a master regulator of mitochondrial biogenesis and cell cycle progression. This pipeline:

1. Extracts TSS (Transcription Start Site) coordinates from a human gene annotation
2. Scans upstream promoter sequences for the NRF1 consensus motif
3. Runs GO Biological Process enrichment on the hit gene list
4. Produces publication-quality visualizations

---

## Repository Structure

```
.
├── make_bed.py          # Step 1 — Build TSS BED file from gene annotation
├── find_nrf1.py         # Step 2 — Scan promoter FASTA for NRF1 motif hits
├── go_enrichment.R      # Step 3 — GO enrichment + plots (clusterProfiler)
├── go_results.csv       # Output — Full enrichment results table
├── go_barplot.png       # Output — Bar plot (top 20 GO terms)
├── go_dotplot.png       # Output — Dot plot (top 20 GO terms)
└── go_cnetplot.png      # Output — Gene-concept network (top 10 GO terms)
```

---

## Requirements

### Python (≥ 3.8)
Standard library only — no extra packages needed.

### R (≥ 4.0)
```r
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
install.packages("ggplot2")
```

### External Tools
- [`bedtools getfasta`](https://bedtools.readthedocs.io/) — to extract promoter sequences from the genome

---

## Usage

### Step 1 — Build TSS BED file

Provide a gzipped TSV from [Ensembl BioMart](https://www.ensembl.org/biomart) with columns including chromosome, strand, TSS, and gene name.

```bash
python make_bed.py
# Input:  human_gene_annotation.tsv.gz
# Output: tss.bed
```

Expected TSV columns (0-indexed): `4=chrom`, `5=strand`, `6=gene_name`, `7=TSS`

### Step 2 — Extract Promoter Sequences

Using `bedtools`, extend each TSS upstream by 500 bp (adjust as needed):

```bash
bedtools slop -i tss.bed -g hg38.chrom.sizes -l 500 -r 0 -s > promoters.bed
bedtools getfasta -fi hg38.fa -bed promoters.bed -s -name > promoter_sequences.fa
```

### Step 3 — Find NRF1 Motif Hits

```bash
python find_nrf1.py
# Input:  promoter_sequences.fa
# Output: nrf1_genes.txt
```

The script scans for the regex `GCGC[ACGT]{2}GCGC` (NRF1 core consensus) and writes one gene symbol per line.

### Step 4 — GO Enrichment Analysis

```bash
Rscript go_enrichment.R
# Input:  nrf1_genes.txt
# Output: go_results.csv, go_dotplot.png, go_barplot.png, go_cnetplot.png
```

---

## Results

The top enriched GO Biological Process terms for NRF1 target genes include:

| Category | Representative Terms |
|---|---|
| Cell cycle | Chromosome segregation, mitotic cell cycle phase transition, nuclear division |
| DNA repair | Double-strand break repair, recombinational repair, DNA recombination |
| Gene regulation | Regulation of DNA metabolic process, regulation of cell cycle phase transition |
| RNA processing | RNA splicing (via transesterification reactions) |

All terms are significant at **BH-adjusted p < 0.05**.

### Visualizations

| File | Description |
|---|---|
| `go_barplot.png` | Top 20 GO terms ranked by gene count, colored by adjusted p-value |
| `go_dotplot.png` | Top 20 GO terms: x = GeneRatio, size = count, color = adjusted p-value |
| `go_cnetplot.png` | Gene-concept network linking GO terms to individual gene symbols |

---

## NRF1 Motif

```
Consensus:  GCGC[ACGT]{2}GCGC
Example:    GCGCATGCGC
Regex used: GCGC[ACGT]{2}GCGC  (case-insensitive)
```

---

## Citation

If you use this pipeline, please cite:

- Yu G, et al. (2012) clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. *OMICS*, 16(5):284–287.
- Carlson M (2019). org.Hs.eg.db: Genome wide annotation for Human. R package.

---

## License

MIT
