# GSE2034_analysis_final.R
#install packges
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(
  c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db"),
  ask    = FALSE,
  update = FALSE
)

install.packages(c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer"))

# load libraries
library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)

res <- read.delim(
  "GSE2034.top.table.tsv",
  header          = TRUE,
  sep             = "\t",
  stringsAsFactors = FALSE
)

cat("Total probes loaded:", nrow(res), "\n")

#identify significant values
deg <- res[
  !is.na(res$adj.P.Val) &
  !is.na(res$logFC)     &
  res$adj.P.Val < 0.05  &
  abs(res$logFC) >= 1,
]

cat("\nSignificant DEGs (adj.P.Val < 0.05 & |logFC| >= 1):", nrow(deg), "\n")
cat("  Upregulated  (logFC > 0):", sum(deg$logFC > 0), "\n")
cat("  Downregulated (logFC < 0):", sum(deg$logFC < 0), "\n")

# Save DEG table
write.csv(deg, "GSE2034_DEGs.csv", row.names = FALSE)
cat("DEG table saved → GSE2034_DEGs.csv\n")

#Volcano plot
res$significance <- "Not significant"
res$significance[res$adj.P.Val < 0.05 & res$logFC >=  1] <- "Upregulated"
res$significance[res$adj.P.Val < 0.05 & res$logFC <= -1] <- "Downregulated"
res$significance <- factor(res$significance,
                           levels = c("Upregulated", "Downregulated", "Not significant"))

top10 <- head(deg[order(deg$adj.P.Val), ], 10)

p_volcano <- ggplot(res, aes(x = logFC,
                              y = -log10(adj.P.Val),
                              color = significance)) +
  geom_point(alpha = 0.55, size = 1.2, stroke = 0) +
  scale_color_manual(
    values = c(
      "Upregulated"     = "#E63946",
      "Downregulated"   = "#457B9D",
      "Not significant" = "#AAAAAA"
    )
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             linewidth = 0.7, colour = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             linewidth = 0.7, colour = "black") +
  geom_text_repel(
    data         = top10,
    aes(label    = Gene.symbol),
    size         = 2.8,
    max.overlaps = 25,
    colour       = "#222222"
  ) +
  labs(
    title    = "Volcano Plot \u2013 GSE2034 (Breast Cancer)",
    subtitle = "35 DEGs highlighted  |  adj.P.Val < 0.05 & |log\u2082FC| \u2265 1",
    x        = "log\u2082 Fold Change",
    y        = "-log\u2081\u2080 (adjusted p-value)",
    color    = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "#555555", size = 10),
    legend.position = "top"
  )

ggsave("volcano_plot_GSE2034.png", plot = p_volcano,
       width = 9, height = 6, dpi = 180)
cat("Volcano plot saved \u2192 volcano_plot_GSE2034.png\n")


#heatmap
deg_clean <- deg[!is.na(deg$Gene.symbol) & deg$Gene.symbol != "", ]
deg_clean$label <- paste0(
  sub("///.*", "", deg_clean$Gene.symbol),  
  " (", deg_clean$ID, ")"
)
deg_clean <- deg_clean[order(deg_clean$logFC), ]

hm_mat           <- matrix(deg_clean$logFC, ncol = 1)
rownames(hm_mat) <- deg_clean$label
colnames(hm_mat) <- "log\u2082FC"

hm_colors <- colorRampPalette(
  rev(brewer.pal(n = 9, name = "RdBu"))
)(100)

pheatmap(
  hm_mat,
  color         = hm_colors,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  show_colnames = TRUE,
  fontsize_row  = 7,
  fontsize_col  = 9,
  main          = "Heatmap \u2013 All 35 DEGs (log\u2082FC)\nGSE2034 Breast Cancer",
  filename      = "heatmap_top50_GSE2034.png",
  width         = 6,
  height        = 14
)
cat("Heatmap saved \u2192 heatmap_top50_GSE2034.png\n")

#enrichGO
gene_symbols_raw <- deg$Gene.symbol
gene_symbols_raw <- gene_symbols_raw[!is.na(gene_symbols_raw) &
                                       gene_symbols_raw != ""]
gene_symbols <- unique(trimws(
  unlist(strsplit(gene_symbols_raw, "///"))
))
gene_symbols <- gene_symbols[gene_symbols != "" & gene_symbols != "nan"]
cat("\nUnique gene symbols for enrichment:", length(gene_symbols), "\n")

gene_map <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = "org.Hs.eg.db"
)
entrez_ids <- unique(gene_map$ENTREZID)
cat("Mapped to Entrez IDs:", length(entrez_ids), "\n")

ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

cat("\nTop GO:BP terms:\n")
print(head(as.data.frame(ego)[, c("Description", "GeneRatio", "p.adjust")], 10))

p_go <- dotplot(ego, showCategory = 15) +
  ggtitle("GO Enrichment: Biological Processes \u2013 GSE2034") +
  theme(plot.title = element_text(face = "bold"))

ggsave("go_enrichment_GSE2034.png", plot = p_go,
       width = 11, height = 8, dpi = 180)
cat("GO plot saved \u2192 go_enrichment_GSE2034.png\n")

#KEGG pathway analysis
ekegg <- enrichKEGG(
  gene         = entrez_ids,
  organism     = "hsa",
  pvalueCutoff = 0.10       
)

cat("\nTop KEGG pathways:\n")
print(head(as.data.frame(ekegg)[, c("Description", "GeneRatio", "p.adjust")], 10))

p_kegg <- dotplot(ekegg, showCategory = 12) +
  ggtitle("KEGG Pathway Enrichment \u2013 GSE2034") +
  theme(plot.title = element_text(face = "bold"))

ggsave("kegg_enrichment_GSE2034.png", plot = p_kegg,
       width = 11, height = 8, dpi = 180)
cat("KEGG plot saved \u2192 kegg_enrichment_GSE2034.png\n")


cat("\n\n===== ANALYSIS SUMMARY =====\n")
cat("Total probes analysed       :", nrow(res), "\n")
cat("Probes significant (FDR<0.05):", sum(res$adj.P.Val < 0.05, na.rm=TRUE), "\n")
cat("DEGs (+ |logFC| >= 1)       :", nrow(deg), "\n")
cat("  Upregulated               :", sum(deg$logFC > 0), "\n")
cat("  Downregulated             :", sum(deg$logFC < 0), "\n")
cat("\nTop upregulated genes    :",
    paste(head(deg$Gene.symbol[deg$logFC > 0][order(-deg$logFC[deg$logFC > 0])], 5),
          collapse = ", "), "\n")
cat("Top downregulated genes  :",
    paste(head(deg$Gene.symbol[deg$logFC < 0][order(deg$logFC[deg$logFC < 0])], 5),
          collapse = ", "), "\n")
cat("\nOutput files:\n")
cat("  GSE2034_DEGs.csv\n")
cat("  volcano_plot_GSE2034.png\n")
cat("  heatmap_top50_GSE2034.png\n")
cat("  go_enrichment_GSE2034.png\n")
cat("  kegg_enrichment_GSE2034.png\n")
cat("============================\n")
