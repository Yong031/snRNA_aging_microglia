suppressPackageStartupMessages({library(Seurat); library(dplyr)})
source(file.path("R", "utils.R")); source(file.path("R", "gene_sets.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript R/01_construct_state_model.R <microglia.rds> <state_gene_sets.csv> <output.rds>", call. = FALSE)
micro <- load_object(args[[1]], c("sample_id", "condition"))
sets <- read.csv(args[[2]], stringsAsFactors = FALSE)
if (!all(c("axis", "gene", "direction") %in% colnames(sets))) stop("Gene-set CSV must contain axis, gene and direction", call. = FALSE)

signed_score <- function(object, axis) {
  pos <- intersect(sets$gene[sets$axis == axis & sets$direction == "positive"], rownames(object))
  neg <- intersect(sets$gene[sets$axis == axis & sets$direction == "negative"], rownames(object))
  if (length(pos) < 2L || length(neg) < 2L) stop("Each axis requires at least two detected positive and negative genes: ", axis)
  object <- AddModuleScore(object, features = list(pos), name = paste0(axis, "_positive"), seed = 1)
  object <- AddModuleScore(object, features = list(neg), name = paste0(axis, "_negative"), seed = 1)
  object[[paste0(axis, "_score")]] <- object[[paste0(axis, "_positive1")]][,1] - object[[paste0(axis, "_negative1")]][,1]
  object
}

micro <- signed_score(micro, "priming")
micro <- signed_score(micro, "pathological")
micro$frozen_state <- classify_state(micro$pathological_score, micro$priming_score)
saveRDS(micro, args[[3]], compress = "xz")
writeLines(capture.output(sessionInfo()), paste0(args[[3]], ".sessionInfo.txt"))
