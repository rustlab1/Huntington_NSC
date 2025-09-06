# Q1: % Q30 for GSM8343539 (SRR29493486)
cat("Q1) Q30 % for GSM8343539 (SRR29493486): 95.67%\n")  # from log file

# Q2: How many genes are detected in GSM8343539 (SRR29493486)?
sample_id <- "GSM8343539"
if (!sample_id %in% colnames(tbl)) stop("Sample not found in tbl.")
detected_genes <- sum(tbl[, sample_id] > 0)
cat("Q2) Genes detected in", sample_id, ":", detected_genes, "\n")

# Q3: How many DEGs (padj < 0.05 & |log2FC| > 1.5)
deg_count <- sum(!is.na(r$padj) & r$padj < 0.05 & abs(r$log2FoldChange) > 1.5)
cat("Q3) Significant DEGs (padj<0.05 & |log2FC|>1.5):", deg_count, "\n")

# Q4: Which gene is most downregulated by log2FC?
down_gene <- rownames(r)[which.min(r$log2FoldChange)]
down_fc   <- r[down_gene, "log2FoldChange"]
down_annot <- annot[down_gene, "Symbol"]
cat(sprintf("Q4) Most downregulated gene: %s (log2FC = %.2f)\n",
            ifelse(!is.null(down_annot), down_annot, down_gene), down_fc))

# Q5: Top GO BP term enriched among upregulated HTT-Mut genes
if (exists("ego_bp") && nrow(ego_bp@result) > 0) {
  top_term <- ego_bp@result[order(ego_bp@result$p.adjust), ][1, c("Description", "p.adjust")]
  cat(sprintf("Q5) Top enriched GO BP term: %s (padj = %.3g)\n", top_term$Description, top_term$p.adjust))
} else {
  cat("Q5) No significant GO BP enrichment found.\n")
}
