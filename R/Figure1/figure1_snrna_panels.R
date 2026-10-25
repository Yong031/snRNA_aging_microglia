suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(patchwork)})
source(file.path("R", "utils.R")); source(file.path("R", "gene_sets.R"))
args <- script_args("Creates Figure 1d-k snRNA-seq panels.")
obj <- load_object(args$input, c("sample_id", "condition", "cell_type"))

# Fig. 1d: major-cell-type UMAP.
p1d <- DimPlot(obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend() + theme_manuscript()
save_panel(p1d, args$output, "Fig1d_cell_type_UMAP", 6.2, 4.8)

# Fig. 1e: microglial marker expression.
markers <- intersect(c("P2ry12", "Cx3cr1", "Tmem119", "Aif1"), rownames(obj))
p1e <- FeaturePlot(obj, features = markers, combine = TRUE, order = TRUE, min.cutoff = "q05", max.cutoff = "q95")
save_panel(p1e, args$output, "Fig1e_microglial_markers", 8, 5)

# Fig. 1f: report composition by biological sample; plot group means without cell-level inference.
composition <- obj[[]] |> count(sample_id, condition, cell_type, name = "n") |> group_by(sample_id, condition) |>
  mutate(fraction = n / sum(n)) |> ungroup() |> filter(cell_type == "Microglia")
write_source_data(composition, args$output, "Fig1f_microglial_abundance")
p1f <- ggplot(composition, aes(condition, 100 * fraction, colour = condition)) +
  geom_point(position = position_jitter(width = 0.08), size = 2) + stat_summary(fun = mean, geom = "crossbar", width = 0.45) +
  labs(x = NULL, y = "Microglia (%)") + theme_manuscript() + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p1f, args$output, "Fig1f_microglial_abundance", 4.3, 3.5)

micro <- subset(obj, subset = cell_type == "Microglia") |> add_signature_scores()

# Fig. 1g/i/k: state scores. Cell-level violins reproduce the display; sample summaries are also exported.
score_cols <- intersect(c("PrimingCore_score", "DiseaseAssociated_score", "Homeostatic_score"), colnames(micro[[]]))
long <- tidyr::pivot_longer(cbind(micro[[]][c("sample_id", "condition")], micro[[]][score_cols]), all_of(score_cols), names_to = "score", values_to = "value")
write_source_data(long, args$output, "Fig1gik_state_scores_cell_level")
sample_scores <- long |> group_by(sample_id, condition, score) |> summarise(value = mean(value), .groups = "drop")
write_source_data(sample_scores, args$output, "Fig1gik_state_scores_sample_level")
p_scores <- ggplot(long, aes(condition, value, fill = condition)) + geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") + facet_wrap(~score, scales = "free_y") +
  labs(x = NULL, y = "Module score") + theme_manuscript() + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p_scores, args$output, "Fig1g_i_k_state_scores", 7.2, 3.5)

# Fig. 1h/j: microglial program dot plots by condition.
program_genes <- unique(c(HOMEOSTATIC, PRO_INFLAMMATORY, M2_LIKE, DISEASE_ASSOCIATED))
p_dot <- DotPlot(micro, features = intersect(program_genes, rownames(micro)), group.by = "condition", dot.scale = 5) +
  RotatedAxis() + scale_colour_gradient(low = "grey90", high = "#138A9E") + theme_manuscript()
save_panel(p_dot, args$output, "Fig1h_j_microglial_programs", 13, 4.2)
