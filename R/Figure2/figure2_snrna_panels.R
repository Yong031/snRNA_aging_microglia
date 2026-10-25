suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr)})
source(file.path("R", "utils.R")); source(file.path("R", "gene_sets.R"))
args <- script_args("Creates Figure 2b-e and 2i-j. Figure 2f is spatial transcriptomics and is intentionally excluded.")
obj <- load_object(args$input, c("sample_id", "condition", "cell_type"))
micro <- subset(obj, subset = cell_type == "Microglia")
validate_metadata(micro[[]], c("microglial_state", "priming_score", "pathological_score"))

# Fig. 2b: integrated microglial-state UMAP.
p2b <- DimPlot(micro, reduction = "umap", group.by = "microglial_state", label = TRUE, repel = TRUE) + NoLegend() + theme_manuscript()
save_panel(p2b, args$output, "Fig2b_integrated_microglial_states", 6, 4.8)

# Fig. 2c/d: two-dimensional state model and quadrant composition.
state_df <- micro[[]] |> mutate(frozen_state = classify_state(pathological_score, priming_score))
centroids <- state_df |> group_by(condition) |> summarise(pathological_score = mean(pathological_score), priming_score = mean(priming_score), .groups = "drop")
p2c <- ggplot(state_df, aes(pathological_score, priming_score, colour = condition)) + geom_point(alpha = 0.12, size = 0.3) +
  geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) +
  geom_point(data = centroids, size = 3) + labs(x = "Pathological status", y = "Priming status") + theme_manuscript()
save_panel(p2c, args$output, "Fig2c_state_space", 5, 4.2)
write_source_data(centroids, args$output, "Fig2c_state_space_centroids")

state_comp <- state_df |> count(sample_id, condition, frozen_state, name = "n") |> group_by(sample_id, condition) |> mutate(percent = 100*n/sum(n)) |> ungroup()
write_source_data(state_comp, args$output, "Fig2d_state_composition")
p2d <- ggplot(state_comp, aes(condition, percent, fill = frozen_state)) + stat_summary(fun = mean, geom = "bar", position = "fill") +
  labs(x = NULL, y = "Microglial cells (%)", fill = "State") + theme_manuscript() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p2d, args$output, "Fig2d_state_composition", 5, 4)

# Fig. 2e: Gpr183 localization across annotated cell types.
p2e <- DotPlot(obj, features = intersect(c("Gpr183", "Aif1", "Tmem119"), rownames(obj)), group.by = "cell_type") +
  RotatedAxis() + theme_manuscript()
save_panel(p2e, args$output, "Fig2e_Gpr183_by_cell_type", 4.8, 5)

# Fig. 2i: temporal association; one point per biological library.
if (all(c("timepoint", "Gpr183") %in% c(colnames(micro[[]]), rownames(micro)))) {
  temporal <- mean_expression_by_sample(micro, "Gpr183", c("sample_id", "timepoint"))
  write_source_data(temporal, args$output, "Fig2i_Gpr183_timecourse")
  p2i <- ggplot(temporal, aes(timepoint, Gpr183, group = sample_id)) + geom_point() + stat_summary(aes(group = 1), fun = mean, geom = "line") +
    labs(x = NULL, y = "Mean normalized Gpr183 expression") + theme_manuscript()
  save_panel(p2i, args$output, "Fig2i_Gpr183_timecourse", 4.8, 3.4)
}

# Fig. 2j: Gpr183-positive versus negative priming scores.
micro$Gpr183_status <- ifelse(FetchData(micro, "Gpr183")$Gpr183 > 0, "Gpr183_pos", "Gpr183_neg")
p2j <- ggplot(micro[[]], aes(Gpr183_status, priming_score, fill = Gpr183_status)) + geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") + labs(x = NULL, y = "Priming score") + theme_manuscript() + theme(legend.position = "none")
save_panel(p2j, args$output, "Fig2j_priming_by_Gpr183_status", 3.2, 3.5)
