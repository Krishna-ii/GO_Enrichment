import re

pattern = re.compile(r'GCGC[ACGT]{2}GCGC', re.IGNORECASE)
input_fa = "promoter_sequences.fa"
output_genes = "nrf1_genes.txt"

genes = set()
current_gene = None
count = 0

with open(input_fa) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            # Extract gene name from header like:
            # >chr1@11869-11870|DDX11L2::chr1:11369-11870(+)
            try:
                current_gene = line.split("|")[1].split("::")[0]
            except IndexError:
                current_gene = None
        else:
            if current_gene and pattern.search(line):
                genes.add(current_gene)
                count += 1

with open(output_genes, "w") as out:
    for g in sorted(genes):
        out.write(g + "\n")

print(f"Found {len(genes)} unique genes with NRF1 motif")
print(f"Total hits: {count}")