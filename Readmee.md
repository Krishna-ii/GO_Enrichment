**Project Overview**
- **Purpose:**: Small promoter/TSS processing pipeline used to extract promoter regions and identify NRF1-related genes and GO results.

**Repository Layout**
- **`hg38.chrom.sizes`**: Chromosome sizes used for coordinate operations.
- **`hg38.fa.fai`**: FASTA index for reference genome used by `samtools faidx`.
- **`tss_clean.bed`, `tss_final.bed`, `tss_sorted_clean.bed`**: TSS (transcription start site) intermediate files (cleaned, finalized, and sorted).
- **`promoters.bed`, `promoter_500bp_fixed.bed`**: BED files for promoter regions (raw and 500bp-fixed versions).
- **`promoter_sequences.fa`**: FASTA sequences extracted for promoters.
- **`nrf1_gene_list.txt`, `nrf1_genes.txt`, `nrf1_hits.txt`, `nrf1_final_genes.txt`**: Files listing NRF1 candidate genes and motif/hit summaries.
- **`codes/`**: Scripts used to run steps (not detailed here).
- **`result/go_results.csv`**: GO enrichment results output.

**Step-by-step Instructions**
- **Step 1 — Prepare TSS files:**: Clean and deduplicate raw TSS calls, then finalize and sort.

  Typical commands (example):

  ```bash
  # remove unwanted columns/headers and sort
  sort -k1,1 -k2,2n tss_clean.bed > tss_sorted_clean.bed
  # produce final TSS set (tool-specific; placeholder)
  cp tss_sorted_clean.bed tss_final.bed
  ```

- **Step 2 — Generate promoter regions:**: Create promoter regions (e.g., 500 bp upstream of TSS) and fix coordinates to chromosome bounds using `hg38.chrom.sizes`.

  Example using `bedtools slop`:

  ```bash
  # create 500bp upstream promoter regions (strand-aware)
  bedtools slop -i tss_final.bed -g hg38.chrom.sizes -l 500 -r 0 -s > promoter_500bp_fixed.bed
  # (Optional) merge or filter to produce a canonical promoters.bed
  sort -k1,1 -k2,2n promoter_500bp_fixed.bed > promoters.bed
  ```

- **Step 3 — Extract promoter sequences:**: Use `samtools faidx` or `bedtools getfasta` to extract sequences with the reference FASTA indexed by `hg38.fa.fai`.

  ```bash
  # extract sequences for promoters
  bedtools getfasta -fi /path/to/hg38.fa -bed promoter_500bp_fixed.bed -s -name -fo promoter_sequences.fa
  ```

- **Step 4 — Motif scan / NRF1 identification:**: Scan `promoter_sequences.fa` for NRF1 motifs (tool examples: FIMO, HOMER). Results are collected into `nrf1_hits.txt` and summarized into `nrf1_genes.txt` / `nrf1_final_genes.txt`.

  Example (FIMO):

  ```bash
  # run FIMO (MEME suite) with NRF1 PWM
  fimo --oc fimo_out nrf1_pwm.meme promoter_sequences.fa
  # aggregate hits
  awk '{print $2}' fimo_out/fimo.tsv | sort | uniq > nrf1_hits.txt
  ```

- **Step 5 — Collate gene lists and finalize:**: Map promoter coordinates back to gene IDs (using your annotation), filter and produce final NRF1 gene lists.

  Example (placeholder):

  ```bash
  # map promoter entries to genes (tool/annotation dependent)
  codes/map_promoters_to_genes.sh promoter_500bp_fixed.bed > nrf1_gene_list.txt
  # produce final curated gene list
  sort nrf1_gene_list.txt | uniq > nrf1_final_genes.txt
  ```

- **Step 6 — GO enrichment analysis:**: Run GO enrichment on `nrf1_final_genes.txt` (tools: clusterProfiler, g:Profiler, topGO), output results to `result/go_results.csv`.

  Example (R / clusterProfiler):

  ```r
  # in R
  library(clusterProfiler)
  genes <- scan('nrf1_final_genes.txt', what = "")
  ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL')
  write.csv(as.data.frame(ego), file='result/go_results.csv')
  ```

**Notes & Tips**
- **Reference paths:**: Replace `/path/to/hg38.fa` with your local reference FASTA path.
- **Tools used:**: `bedtools`, `samtools`, `fimo` (MEME), `R` with `clusterProfiler` are common; install them as needed.
- **Reproducibility:**: Keep the exact commands used in `codes/` so steps can be re-run; capture versions with `tool --version`.

**Files to inspect for verification**
- Promoter coordinates: [promoter_500bp_fixed.bed](promoter_500bp_fixed.bed)
- Promoter sequences: [promoter_sequences.fa](promoter_sequences.fa)
- NRF1 results: [nrf1_hits.txt](nrf1_hits.txt), [nrf1_final_genes.txt](nrf1_final_genes.txt)
- GO results: [result/go_results.csv](result/go_results.csv)

If you want, I can:
- create separate per-step README files (one per major step),
- or expand `codes/` with runnable scripts and example inputs.
