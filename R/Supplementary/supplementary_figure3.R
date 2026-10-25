suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr)})
source(file.path("R", "utils.R")); source(file.path("R", "gene_sets.R"))
args <- script_args("Creates Supplementary Figure 3a-g from the integrated microglial atlas.")
micro <- load_object(args$input, c("sample_id", "condition", "microglial_state", "priming_score", "pathological_score"))

# S3a-d: atlas by sample/state and condition-specific distributions.
save_panel(DimPlot(micro, group.by = "sample_id", raster = TRUE) + NoLegend() + theme_manuscript(), args$output, "SuppFig3a_atlas_by_sample", 6, 4.8)
save_panel(DimPlot(micro, group.by = "microglial_state", label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() + theme_manuscript(), args$output, "SuppFig3b_state_annotation", 6, 4.8)
state_markers <- c("P2ry12", "Cx3cr1", "Hexb", "Sall1", "Fa2h", "Pde4b", "Mef2c", "Tnf", "Il1b", "Nlrp3", "Apoe", "Lpl", "Ifi27l2a", "H2-Aa", "Cd74", "Top2a", "Mki67")
save_panel(DotPlot(micro, features = intersect(state_markers, rownames(micro)), group.by = "microglial_state") + RotatedAxis() + theme_manuscript(), args$output, "SuppFig3c_state_markers", 10, 4.6)
for (group in intersect(c("Sham_Young", "Sham_Aged", "SNI_Young", "SNI_Aged"), unique(micro$condition))) {
  p <- DimPlot(micro, cells.highlight = WhichCells(micro, expression = condition == group), cols = c("grey90", "#2878B5"), raster = TRUE) + ggtitle(group) + NoLegend() + theme_manuscript()
  save_panel(p, args$output, paste0("SuppFig3d_distribution_", group), 4.2, 3.8)
}

# S3e: state-program heat map source data; trajectory plotting uses the monocle3 script in Figure6.
programs <- intersect(c("Homeostatic_score", "priming_score", "aging_amplified_score", "pathological_score"), colnames(micro[[]]))
heat <- micro[[]] |> group_by(condition) |> summarise(across(all_of(programs), mean), .groups = "drop")
heat_long <- pivot_longer(heat, -condition, names_to = "program", values_to = "score") |> group_by(program) |> mutate(gene_wise_z = as.numeric(scale(score))) |> ungroup()
write_source_data(heat_long, args$output, "SuppFig3e_program_heatmap")
p_heat <- ggplot(heat_long, aes(condition, program, fill = gene_wise_z)) + geom_tile(colour = "white") + scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#E45756") + theme_manuscript() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p_heat, args$output, "SuppFig3e_program_heatmap", 5.5, 3.2)

# S3f: sample-level Ch25h/Gpr183 time course.
if ("timepoint" %in% colnames(micro[[]])) {
  temporal <- mean_expression_by_sample(micro, c("Ch25h", "Gpr183"), c("sample_id", "timepoint", intersect("sex", colnames(micro[[]]))))
  write_source_data(temporal, args$output, "SuppFig3f_Ch25h_Gpr183_timecourse")
}

# S3g: export gene-wise Spearman correlations for downstream pathway enrichment.
expr <- GetAssayData(micro, assay = DefaultAssay(micro), layer = "data")
rho <- apply(expr, 1, function(x) suppressWarnings(cor(x, micro$priming_score, method = "spearman")))
ranking <- data.frame(gene = names(rho), spearman_rho = unname(rho)) |> arrange(desc(spearman_rho))
write_source_data(ranking, args$output, "SuppFig3g_priming_axis_gene_ranking")
