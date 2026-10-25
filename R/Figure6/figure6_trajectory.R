suppressPackageStartupMessages({library(monocle3); library(ggplot2)})
source(file.path("R", "utils.R"))
args <- script_args("Creates Figure 6g from a monocle3 cell_data_set containing condition and microglial_state metadata.")
cds <- readRDS(args$input)
if (!inherits(cds, "cell_data_set")) stop("Input must be a monocle3 cell_data_set", call. = FALSE)
validate_metadata(as.data.frame(colData(cds)), c("condition", "microglial_state"))
for (group in c("Sham_Aged", "SNI_Aged", "SNI_Aged_M3")) {
  cells <- rownames(subset(as.data.frame(colData(cds)), condition == group))
  panel <- plot_cells(cds[, cells], color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = TRUE) + ggtitle(group)
  save_panel(panel, args$output, paste0("Fig6g_trajectory_", group), 4.2, 4)
}
trajectory_data <- as.data.frame(colData(cds))[, c("condition", "microglial_state"), drop = FALSE]
trajectory_data$cell_id <- rownames(trajectory_data)
trajectory_data$pseudotime <- monocle3::pseudotime(cds)
write_source_data(trajectory_data, args$output, "Fig6g_pseudotime")
