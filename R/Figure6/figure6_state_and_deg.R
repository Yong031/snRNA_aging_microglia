suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr)})
source(file.path("R", "utils.R")); source(file.path("R", "gene_sets.R"))
args <- script_args("Creates Figure 6a-d state-space, state-score, DEG and dot-plot panels.")
obj <- load_object(args$input, c("sample_id", "condition", "cell_type"))
micro <- subset(obj, subset = cell_type == "Microglia") |> add_signature_scores()
aged <- subset(micro, subset = condition %in% c("Sham_Aged", "SNI_Aged", "SNI_Aged_M3"))
validate_metadata(aged[[]], c("priming_score", "pathological_score"))

# Fig. 6a/b.
p6a <- ggplot(aged[[]], aes(pathological_score, priming_score, colour = condition)) + geom_point(alpha = .12, size = .3) +
  stat_ellipse(level = .68, linewidth = .7) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + theme_manuscript()
save_panel(p6a, args$output, "Fig6a_aged_state_space", 5.5, 4.4)
score_cols <- intersect(c("PrimingCore_score", "DiseaseAssociated_score", "Homeostatic_score"), colnames(aged[[]]))
score_df <- pivot_longer(cbind(aged[[]][c("sample_id", "condition")], aged[[]][score_cols]), all_of(score_cols), names_to = "score", values_to = "value")
write_source_data(score_df, args$output, "Fig6b_state_scores_cell_level")
p6b <- ggplot(score_df, aes(condition, value, fill = condition)) + geom_violin(trim = TRUE) + geom_boxplot(width = .1, outlier.shape = NA, fill = "white") + facet_wrap(~score, scales = "free_y") + theme_manuscript() + theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
save_panel(p6b, args$output, "Fig6b_state_scores", 8, 3.6)

# Fig. 6c: exploratory cell-level markers; use pseudobulk for biological-sample inference in a confirmatory analysis.
Idents(aged) <- "condition"
deg <- FindMarkers(aged, ident.1 = "SNI_Aged_M3", ident.2 = "SNI_Aged", logfc.threshold = 0, min.pct = .05) |>
  tibble::rownames_to_column("gene") |> mutate(neg_log10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)))
write_source_data(deg, args$output, "Fig6c_M3_vs_aCSF_DEG")
p6c <- ggplot(deg, aes(avg_log2FC, neg_log10_padj, colour = p_val_adj < .05 & abs(avg_log2FC) >= .25)) + geom_point(alpha = .55, size = .8) + scale_colour_manual(values = c("grey75", "#D95F59")) + labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~P)) + theme_manuscript() + theme(legend.position = "none")
save_panel(p6c, args$output, "Fig6c_M3_vs_aCSF_volcano", 4.8, 4)

# Fig. 6d.
p6d <- DotPlot(aged, features = intersect(unique(c(PRO_INFLAMMATORY[1:6], HOMEOSTATIC[1:6])), rownames(aged)), group.by = "condition") + RotatedAxis() + theme_manuscript()
save_panel(p6d, args$output, "Fig6d_pathological_homeostatic_programs", 7, 3.6)
