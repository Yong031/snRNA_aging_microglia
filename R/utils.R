classify_state <- function(pathological, priming, pathological_cutoff = 0, priming_cutoff = 1) {
  stopifnot(length(pathological) == length(priming))
  ifelse(
    pathological <= pathological_cutoff & priming < priming_cutoff, "Homeostatic-like",
    ifelse(
      pathological > pathological_cutoff & priming < priming_cutoff, "Activated/pathological",
      ifelse(pathological <= pathological_cutoff, "Primed-like", "Primed-pathological")
    )
  )
}

validate_metadata <- function(metadata, required) {
  missing <- setdiff(required, colnames(metadata))
  if (length(missing)) {
    stop("Missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

ensure_output_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

script_args <- function(description) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop("Usage: Rscript <script.R> <input_seurat.rds> <output_directory>\n", description,
         call. = FALSE)
  }
  list(input = args[[1]], output = ensure_output_dir(args[[2]]))
}

load_object <- function(path, required_metadata = character()) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  object <- readRDS(path)
  if (!inherits(object, "Seurat")) stop("Input must be a Seurat object", call. = FALSE)
  validate_metadata(object[[]], required_metadata)
  object
}

save_panel <- function(plot, output_dir, stem, width, height) {
  ggplot2::ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot,
                  width = width, height = height, units = "in", device = cairo_pdf)
  ggplot2::ggsave(file.path(output_dir, paste0(stem, ".png")), plot,
                  width = width, height = height, units = "in", dpi = 600)
}

write_source_data <- function(x, output_dir, stem) {
  utils::write.csv(x, file.path(output_dir, paste0(stem, "_source_data.csv")), row.names = FALSE)
}

mean_expression_by_sample <- function(object, genes, group_cols) {
  validate_metadata(object[[]], group_cols)
  genes <- intersect(genes, rownames(object))
  if (!length(genes)) stop("None of the requested genes are present", call. = FALSE)
  values <- Seurat::FetchData(object, vars = c(group_cols, genes))
  stats::aggregate(values[genes], values[group_cols], mean)
}

theme_manuscript <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0)
    )
}
