suppressPackageStartupMessages({library(clusterProfiler); library(org.Mm.eg.db); library(ggplot2); library(dplyr)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript R/02_enrichment_analysis.R <gene_table.csv> <mode:go|gsea> <output_prefix>", call. = FALSE)
tab <- read.csv(args[[1]], stringsAsFactors = FALSE)
mode <- args[[2]]; prefix <- args[[3]]

if (mode == "go") {
  if (!"gene" %in% colnames(tab)) stop("GO mode requires a gene column", call. = FALSE)
  universe <- if ("in_universe" %in% colnames(tab)) tab$gene[tab$in_universe] else tab$gene
  selected <- if ("selected" %in% colnames(tab)) tab$gene[tab$selected] else tab$gene
  converted <- bitr(unique(c(selected, universe)), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  genes <- converted$ENTREZID[converted$SYMBOL %in% selected]
  background <- converted$ENTREZID[converted$SYMBOL %in% universe]
  result <- enrichGO(gene = genes, universe = background, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)
  write.csv(as.data.frame(result), paste0(prefix, "_GO_BP.csv"), row.names = FALSE)
  pdf(paste0(prefix, "_GO_BP_dotplot.pdf"), width = 6, height = 4.5); print(dotplot(result, showCategory = 10)); dev.off()
} else if (mode == "gsea") {
  if (!all(c("gene", "rank") %in% colnames(tab))) stop("GSEA mode requires gene and rank columns", call. = FALSE)
  if (!all(c("term", "term_gene") %in% colnames(tab))) stop("Append or merge term and term_gene columns for custom TERM2GENE GSEA", call. = FALSE)
  ranks <- setNames(tab$rank, tab$gene); ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  term2gene <- unique(tab[, c("term", "term_gene")]); names(term2gene) <- c("term", "gene")
  result <- GSEA(ranks, TERM2GENE = term2gene, pAdjustMethod = "BH", seed = TRUE, verbose = FALSE)
  write.csv(as.data.frame(result), paste0(prefix, "_GSEA.csv"), row.names = FALSE)
  pdf(paste0(prefix, "_GSEA_ridgeplot.pdf"), width = 7, height = 5); print(enrichplot::ridgeplot(result, showCategory = 10)); dev.off()
} else stop("mode must be go or gsea", call. = FALSE)
