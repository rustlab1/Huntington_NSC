## offline input
tbl <- as.matrix(
  read.table("counts_matrix.tsv",
             header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
)
# robust read
annot <- data.table::fread(
  "gene_annotation.tsv",
  sep = "\t",
  header = TRUE,
  quote = "",         # do not treat quotes specially
  fill = TRUE,        # pad short rows
  data.table = FALSE,
  check.names = FALSE
)

rownames(annot) <- annot$GeneID

## sample selection
gsms <- "110000XXXX"                       # 1 equals HD, 0 equals ctrl, X excluded
sml <- strsplit(gsms, "")[[1]]
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[, sel, drop = FALSE]

## groups
# map 0 to ctrl and 1 to HD explicitly
gs <- factor(sml, levels = c("0", "1"), labels = c("HD", "ctrl"))
groups <- c("ctrl", "HD")
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

## DESeq2
library(DESeq2)
keep <- rowSums(tbl >= 10) >= min(table(gs))
tbl <- tbl[keep, , drop = FALSE]

ds <- DESeqDataSetFromMatrix(countData = tbl,
                             colData = sample_info,
                             design = ~ Group)
ds <- DESeq(ds, test = "Wald", sfType = "poscount")

r <- results(ds, contrast = c("Group", groups[2], groups[1]),
             alpha = 0.05, pAdjustMethod = "fdr")

## top table with annotation
ord <- order(r$padj)
topn <- head(ord, 250)
tT <- merge(as.data.frame(r[topn, ]),
            annot, by.x = 0, by.y = "GeneID", sort = FALSE)
colnames(tT)[1] <- "GeneID"
tT <- tT[, c("GeneID","padj","pvalue","lfcSE","stat",
             "log2FoldChange","baseMean","Symbol","Description")]

#write.table(tT, file = "top250.tsv",
#            sep = "\t", quote = FALSE, row.names = FALSE)

## optional export of normalized counts for reuse
norm_counts <- counts(ds, normalized = TRUE)
#write.table(norm_counts, file = "normalized_counts.tsv",
#            sep = "\t", quote = FALSE, col.names = NA)


# volcano plot
old.pal <- palette(c("#00BFFF", "#FF3030")) # low-hi colors
par(mar=c(4,4,2,1), cex.main=1.5)
plot(r$log2FoldChange, -log10(r$padj), main=paste(groups[1], "vs", groups[2]),
     xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)
with(subset(r, padj<0.05 & abs(log2FoldChange) >= 0),
     points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)

## GO BP enrichment with clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

deg <- as.data.frame(r)
deg_ids <- rownames(deg[which(deg$padj < 0.05 & abs(deg$log2FoldChange) > 1.5), ])
bg_ids  <- rownames(deg)

ego_bp <- enrichGO(gene          = deg_ids,
                   universe      = bg_ids,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# save enrichment table
ego_tab <- as.data.frame(ego_bp)
#write.table(ego_tab, file = "enrichGO_BP_results.tsv",
 #           sep = "\t", quote = FALSE, row.names = FALSE)




barplot(ego_bp, showCategory = 20, title = "GO BP enrichment")
