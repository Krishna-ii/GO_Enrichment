library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# ── 1. Load gene list ──────────────────────────────────────────
genes <- readLines("nrf1_genes.txt")
genes <- genes[genes != "" & !is.na(genes)]
cat("Total genes loaded:", length(genes), "\n")

# ── 2. Map gene symbols → Entrez IDs ──────────────────────────
gene_ids <- bitr(genes,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)
cat("Successfully mapped:", nrow(gene_ids), "genes\n")

# ── 3. GO Enrichment (Biological Process) ─────────────────────
ego <- enrichGO(gene          = gene_ids$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

cat("Enriched GO terms found:", nrow(as.data.frame(ego)), "\n")

# ── 4. Save results table ──────────────────────────────────────
write.csv(as.data.frame(ego), "go_results.csv", row.names = FALSE)
cat("Saved: go_results.csv\n")

# ── 5. Dotplot ─────────────────────────────────────────────────
png("go_dotplot.png", width = 1400, height = 900, res = 130)
print(dotplot(ego,
              showCategory = 20,
              title = "GO Biological Process — NRF1 Target Genes"))
dev.off()
cat("Saved: go_dotplot.png\n")

# ── 6. Barplot ─────────────────────────────────────────────────
png("go_barplot.png", width = 1400, height = 900, res = 130)
print(barplot(ego,
              showCategory = 20,
              title = "GO Enrichment — NRF1 Target Genes"))
dev.off()
cat("Saved: go_barplot.png\n")

# ── 7. Gene-Concept Network ────────────────────────────────────
png("go_cnetplot.png", width = 1600, height = 1200, res = 130)
print(cnetplot(ego, showCategory = 10))
dev.off()
cat("Saved: go_cnetplot.png\n")

cat("\n✓ All done!\n")