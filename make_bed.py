import gzip

input_file = "human_gene_annotation.tsv.gz"
output_file = "tss.bed"

with gzip.open(input_file, 'rt') as f, open(output_file, 'w') as out:
    header = f.readline().strip().split('\t')  # skip header

    written = 0
    skipped = 0

    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < 8:
            skipped += 1
            continue

        chrom  = "chr" + cols[4]   # chromosome_name → add "chr" prefix (e.g. 1 → chr1)
        start  = cols[7]           # transcription_start_site
        strand = cols[5]           # strand
        gene   = cols[6]           # external_gene_name

        # Skip if chromosome is weird (e.g. patch contigs) or gene name empty
        if not gene.strip() or not chrom.strip():
            skipped += 1
            continue

        try:
            tss = int(start)
        except ValueError:
            skipped += 1
            continue

        tss_end = tss + 1
        name = f"{chrom}@{tss}-{tss_end}|{gene}"
        out.write(f"{chrom}\t{tss}\t{tss_end}\t{name}\t.\t{strand}\n")
        written += 1

print(f"Done! Written: {written} entries, Skipped: {skipped}")
print(f"BED file saved to: {output_file}")