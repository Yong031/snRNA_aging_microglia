suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(igraph); library(ggraph)})
source(file.path("R", "utils.R")); source(file.path("R", "gene_sets.R"))
args <- script_args("Creates Supplementary Figure 6 transcriptomic panels for Gpr183-positive versus Gpr183-negative microglia.")
obj <- load_object(args$input, c("sample_id", "condition", "cell_type", "priming_score", "pathological_score"))
micro <- subset(obj, subset = cell_type == "Microglia") |> add_signature_scores()
micro$Gpr183_status <- ifelse(FetchData(micro, "Gpr183")$Gpr183 > 0, "Gpr183_pos", "Gpr183_neg")
Idents(micro) <- "Gpr183_status"

# S6a: differential expression.
deg <- FindMarkers(micro, ident.1 = "Gpr183_pos", ident.2 = "Gpr183_neg", logfc.threshold = 0, min.pct = .05) |> tibble::rownames_to_column("gene")
write_source_data(deg, args$output, "SuppFig6a_Gpr183_positive_DEG")
p <- ggplot(deg, aes(avg_log2FC, -log10(pmax(p_val_adj, .Machine$double.xmin)), colour = p_val_adj < .05 & abs(avg_log2FC) >= .25)) + geom_point(alpha = .55, size = .8) + scale_colour_manual(values = c("grey80", "#D95F59")) + theme_manuscript() + theme(legend.position = "none")
save_panel(p, args$output, "SuppFig6a_Gpr183_positive_volcano", 5, 4)

# S6b: data-driven co-expression table; network styling can be applied in Cytoscape.
expr <- GetAssayData(micro, layer = "data")
target <- as.numeric(expr["Gpr183", ])
correlations <- apply(expr, 1, function(x) suppressWarnings(cor(x, target, method = "spearman")))
network <- data.frame(source = "Gpr183", target = names(correlations), spearman_rho = unname(correlations)) |> filter(is.finite(spearman_rho), target != "Gpr183") |> arrange(desc(abs(spearman_rho))) |> slice_head(n = 50)
write_source_data(network, args$output, "SuppFig6b_Gpr183_coexpression_network")

# S6c: functional modules; S6f-h/j-k: signature scores.
modules <- unique(c(PRO_INFLAMMATORY, DISEASE_ASSOCIATED, PRIMING_CORE, CALCIUM_MOBILIZATION))
save_panel(DotPlot(micro, features = intersect(modules, rownames(micro)), group.by = "Gpr183_status") + RotatedAxis() + theme_manuscript(), args$output, "SuppFig6c_i_functional_modules", 12, 4)
score_cols <- intersect(c("ProInflammatory_score", "DiseaseAssociated_score", "PrimingCore_score", "priming_score", "pathological_score"), colnames(micro[[]]))
score_df <- pivot_longer(cbind(micro[[]][c("sample_id", "condition", "Gpr183_status")], micro[[]][score_cols]), all_of(score_cols), names_to = "score", values_to = "value")
write_source_data(score_df, args$output, "SuppFig6f_h_j_k_signature_scores")
p_scores <- ggplot(score_df, aes(Gpr183_status, value, fill = Gpr183_status)) + geom_violin(trim = TRUE) + geom_boxplot(width = .1, outlier.shape = NA, fill = "white") + facet_wrap(~score, scales = "free_y") + theme_manuscript() + theme(legend.position = "none")
save_panel(p_scores, args$output, "SuppFig6f_h_j_k_signature_scores", 9, 5)

# S6e: frozen-state distribution for sham young versus agonist-treated young.
micro$frozen_state <- classify_state(micro$pathological_score, micro$priming_score)
frozen <- micro[[]] |> filter(condition %in% c("Sham_Young", "Sham_Young_Agonist")) |> count(sample_id, condition, frozen_state, name = "n") |> group_by(sample_id, condition) |> mutate(percent = 100*n/sum(n)) |> ungroup()
write_source_data(frozen, args$output, "SuppFig6e_frozen_state_distribution")

# S6i: normalized expression by the four core conditions.
gpr <- mean_expression_by_sample(micro, "Gpr183", c("sample_id", "condition"))
write_source_data(gpr, args$output, "SuppFig6i_Gpr183_expression")
