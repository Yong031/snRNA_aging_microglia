suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)})
source(file.path("R", "utils.R")); source(file.path("R", "gene_sets.R"))
args <- script_args("Creates Figure 3d/e and 3i/j microglial program panels.")
obj <- load_object(args$input, c("sample_id", "condition", "cell_type"))
micro <- subset(obj, subset = cell_type == "Microglia") |> add_signature_scores()

# Fig. 3d/e: sham, SNI and microglial Gpr183-cKO.
cko <- subset(micro, subset = condition %in% c("Sham_Young", "SNI_Young", "SNI_Young_GPR183cKO"))
p3d <- DotPlot(cko, features = intersect(unique(c(HOMEOSTATIC, PRO_INFLAMMATORY, DISEASE_ASSOCIATED)), rownames(cko)), group.by = "condition") + RotatedAxis() + theme_manuscript()
save_panel(p3d, args$output, "Fig3d_cKO_gene_programs", 11, 4)
score_cols <- intersect(c("Homeostatic_score", "ProInflammatory_score", "DiseaseAssociated_score"), colnames(cko[[]]))
score_df <- tidyr::pivot_longer(cbind(cko[[]][c("sample_id", "condition")], cko[[]][score_cols]), all_of(score_cols), names_to = "score", values_to = "value")
write_source_data(score_df, args$output, "Fig3e_cKO_state_scores")
p3e <- ggplot(score_df, aes(condition, value, fill = condition)) + geom_violin(trim = TRUE) + geom_boxplot(width = .12, outlier.shape = NA, fill = "white") + facet_wrap(~score, scales = "free_y") + theme_manuscript() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
save_panel(p3e, args$output, "Fig3e_cKO_state_scores", 7, 3.7)

# Fig. 3i/j: agonist-treated young microglia.
agonist <- subset(micro, subset = condition %in% c("Sham_Young", "Sham_Young_Agonist"))
p3i <- DotPlot(agonist, features = intersect(unique(c(HOMEOSTATIC, DISEASE_ASSOCIATED)), rownames(agonist)), group.by = "condition") + RotatedAxis() + theme_manuscript()
save_panel(p3i, args$output, "Fig3i_agonist_gene_programs", 9, 3.6)
agonist_scores <- agonist[[]][intersect(c("sample_id", "condition", "PrimingCore_score", "DiseaseAssociated_score"), colnames(agonist[[]]))]
write_source_data(agonist_scores, args$output, "Fig3j_agonist_state_scores")
p3j <- ggplot(tidyr::pivot_longer(agonist_scores, -c(sample_id, condition), names_to = "score", values_to = "value"), aes(condition, value, fill = condition)) + geom_violin(trim = TRUE) + geom_boxplot(width = .12, outlier.shape = NA, fill = "white") + facet_wrap(~score, scales = "free_y") + theme_manuscript() + theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
save_panel(p3j, args$output, "Fig3j_agonist_state_scores", 5.5, 3.5)
