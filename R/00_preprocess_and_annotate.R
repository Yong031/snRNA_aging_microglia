suppressPackageStartupMessages({library(Seurat); library(harmony); library(ggplot2)})
source(file.path("R", "utils.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript R/00_preprocess_and_annotate.R <merged_raw_seurat.rds> <cluster_annotation.csv> <output.rds>", call. = FALSE)
}

set.seed(1)
obj <- readRDS(args[[1]])
if (!inherits(obj, "Seurat")) stop("Input must be a Seurat object", call. = FALSE)
validate_metadata(obj[[]], c("sample_id"))

# These QC defaults are starting points and must be justified for the experiment.
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj <- subset(obj, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & nCount_RNA <= 50000 & percent.mt <= 10)

# DoubletFinder should be run per library before merging. The expected rate and homotypic
# adjustment depend on recovered nuclei and sample-specific clustering, so they are not guessed here.
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "sample_id", reduction = "pca", dims.use = 1:30)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 1)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, seed.use = 1)

annotation <- read.csv(args[[2]], stringsAsFactors = FALSE)
if (!all(c("cluster", "cell_type") %in% colnames(annotation))) stop("Annotation CSV must contain cluster and cell_type columns", call. = FALSE)
mapping <- setNames(annotation$cell_type, annotation$cluster)
obj$cell_type <- unname(mapping[as.character(obj$seurat_clusters)])
if (anyNA(obj$cell_type)) stop("Some clusters are missing from the annotation CSV", call. = FALSE)

saveRDS(obj, args[[3]], compress = "xz")
writeLines(capture.output(sessionInfo()), paste0(args[[3]], ".sessionInfo.txt"))
