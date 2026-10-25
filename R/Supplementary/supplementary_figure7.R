suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr)})
source(file.path("R", "utils.R"))
args <- script_args("Creates Supplementary Figure 7a-e/g. Trajectory panel f uses Figure6/figure6_trajectory.R.")
micro <- load_object(args$input, c("sample_id", "condition", "priming_score", "pathological_score"))

# S7a: all perturbations in the common state space.
p7a <- ggplot(micro[[]], aes(pathological_score, priming_score, colour = condition)) + geom_point(alpha = .1, size = .25) + stat_ellipse(level = .68) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + theme_manuscript()
save_panel(p7a, args$output, "SuppFig7a_all_conditions_state_space", 6.2, 4.8)

# S7b: condition-level score heat map.
programs <- intersect(c("aging_amplified_score", "gpr183_dependent_pathological_score", "homeostatic_score", "priming_score", "pathological_score"), colnames(micro[[]]))
summary <- micro[[]] |> group_by(condition) |> summarise(across(all_of(programs), mean), .groups = "drop") |> pivot_longer(-condition, names_to = "program", values_to = "score") |> group_by(program) |> mutate(z = as.numeric(scale(score))) |> ungroup()
write_source_data(summary, args$output, "SuppFig7b_program_heatmap")
p7b <- ggplot(summary, aes(condition, program, fill = z)) + geom_tile(colour = "white") + scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#E45756") + theme_manuscript() + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_panel(p7b, args$output, "SuppFig7b_program_heatmap", 7.2, 3.5)

# S7c: aged density redistribution and priming score.
aged <- subset(micro, subset = condition %in% c("Sham_Aged", "SNI_Aged", "SNI_Aged_M3"))
for (group in unique(aged$condition)) {
  p <- DimPlot(aged, cells.highlight = WhichCells(aged, expression = condition == group), cols = c("grey90", "#2A9D8F"), raster = TRUE) + ggtitle(group) + NoLegend() + theme_manuscript()
  save_panel(p, args$output, paste0("SuppFig7c_density_", group), 4, 3.7)
}

# S7d/e: export directional DEG lists for GO enrichment.
Idents(aged) <- "condition"
deg <- FindMarkers(aged, ident.1 = "SNI_Aged_M3", ident.2 = "SNI_Aged", logfc.threshold = 0, min.pct = .05) |> tibble::rownames_to_column("gene")
write_source_data(filter(deg, p_val_adj < .05, avg_log2FC < 0), args$output, "SuppFig7d_M3_downregulated_GO_input")
write_source_data(filter(deg, p_val_adj < .05, avg_log2FC > 0), args$output, "SuppFig7e_M3_upregulated_GO_input")

# S7g: Gpr183-centred correlation network.
expr <- GetAssayData(micro, layer = "data"); target <- as.numeric(expr["Gpr183", ])
rho <- apply(expr, 1, function(x) suppressWarnings(cor(x, target, method = "spearman")))
network <- data.frame(source = "Gpr183", target = names(rho), correlation = unname(rho)) |> filter(is.finite(correlation), target != "Gpr183") |> arrange(desc(abs(correlation))) |> slice_head(n = 30)
write_source_data(network, args$output, "SuppFig7g_Gpr183_regulatory_network")
