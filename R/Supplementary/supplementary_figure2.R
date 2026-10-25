suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr)})
source(file.path("R", "utils.R")); source(file.path("R", "gene_sets.R"))
args <- script_args("Creates Supplementary Figure 2a-i.")
obj <- load_object(args$input, c("sample_id", "condition", "cell_type"))

# S2a: post-QC distributions. A before-QC object can be supplied separately if a direct before/after panel is required.
qc <- obj[[]] |> select(sample_id, any_of(c("nFeature_RNA", "nCount_RNA", "percent.mt"))) |> pivot_longer(-sample_id, names_to = "metric", values_to = "value")
write_source_data(qc, args$output, "SuppFig2a_QC_metrics")
p <- ggplot(qc, aes("After QC", value)) + geom_violin(fill = "#2A9D8F") + facet_wrap(~metric, scales = "free_y") + labs(x = NULL, y = NULL) + theme_manuscript()
save_panel(p, args$output, "SuppFig2a_QC_metrics", 7, 3.3)

# S2b-d: clusters, annotation markers and cell-type composition.
save_panel(DimPlot(obj, group.by = "seurat_clusters", label = TRUE) + NoLegend() + theme_manuscript(), args$output, "SuppFig2b_unsupervised_clusters", 5, 4)
annotation_markers <- c("Mbp", "Plp1", "Pdgfra", "Aqp4", "Slc1a3", "Aif1", "P2ry12", "Pecam1", "Kdr", "Rgs5", "Pdgfrb", "Foxj1", "Gad1", "Slc17a6", "Chat", "Pkd2l1")
save_panel(DotPlot(obj, features = intersect(annotation_markers, rownames(obj)), group.by = "cell_type") + RotatedAxis() + theme_manuscript(), args$output, "SuppFig2c_annotation_markers", 10, 4.8)
composition <- obj[[]] |> count(sample_id, condition, cell_type, name = "n") |> group_by(sample_id, condition) |> mutate(percent = 100*n/sum(n)) |> ungroup()
write_source_data(composition, args$output, "SuppFig2d_cell_type_composition")

# S2e-i: microglial programs, age DEGs, GO input, and module scores.
micro <- subset(obj, subset = cell_type == "Microglia") |> add_signature_scores()
save_panel(DotPlot(micro, features = intersect(unique(c(HOMEOSTATIC, PRO_INFLAMMATORY, M2_LIKE, DISEASE_ASSOCIATED)), rownames(micro)), group.by = "condition") + RotatedAxis() + theme_manuscript(), args$output, "SuppFig2e_microglial_programs", 13, 4.4)
Idents(micro) <- "condition"
deg_age <- FindMarkers(micro, ident.1 = "Sham_Aged", ident.2 = "Sham_Young", logfc.threshold = 0, min.pct = .05) |> tibble::rownames_to_column("gene")
write_source_data(deg_age, args$output, "SuppFig2f_age_DEG")
deg_injury <- FindMarkers(micro, ident.1 = "SNI_Aged", ident.2 = "Sham_Aged", logfc.threshold = .1, min.pct = .05) |> tibble::rownames_to_column("gene")
write_source_data(deg_injury, args$output, "SuppFig2g_injury_DEG_for_GO")
score_cols <- grep("_score$", colnames(micro[[]]), value = TRUE)
scores <- pivot_longer(cbind(micro[[]][c("sample_id", "condition")], micro[[]][score_cols]), all_of(score_cols), names_to = "score", values_to = "value")
write_source_data(scores, args$output, "SuppFig2h_i_module_scores")
p_scores <- ggplot(scores, aes(condition, value, fill = condition)) + geom_violin(trim = TRUE) + geom_boxplot(width = .1, outlier.shape = NA, fill = "white") + facet_wrap(~score, scales = "free_y") + theme_manuscript() + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p_scores, args$output, "SuppFig2h_i_module_scores", 10, 5)
