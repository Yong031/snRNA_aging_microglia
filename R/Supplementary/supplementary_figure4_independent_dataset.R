suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)})
source(file.path("R", "utils.R"))
args <- script_args("Creates Supplementary Figure 4a-e from CRA023125/GSE286467 or another independent spinal-cord dataset.")
obj <- load_object(args$input, c("sample_id", "condition", "cell_type"))

save_panel(DimPlot(obj, group.by = "condition", raster = TRUE) + theme_manuscript(), args$output, "SuppFig4a_UMAP_condition", 5, 4)
save_panel(DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() + theme_manuscript(), args$output, "SuppFig4a_UMAP_cell_type", 5, 4)
markers <- c("Plp1", "Cnp", "Mbp", "P2ry12", "Tmem119", "Cx3cr1", "Aif1", "Gfap", "Aqp4", "Pecam1", "Rgs5", "Cd3d", "Cd79a", "Nkg7", "S100a8", "S100a9")
save_panel(DotPlot(obj, features = intersect(markers, rownames(obj)), group.by = "cell_type") + RotatedAxis() + theme_manuscript(), args$output, "SuppFig4b_annotation_markers", 11, 5)

composition <- obj[[]] |> count(sample_id, condition, cell_type, name = "n") |> group_by(sample_id, condition) |> mutate(percent = 100*n/sum(n)) |> ungroup() |> filter(cell_type == "Microglia")
write_source_data(composition, args$output, "SuppFig4c_microglial_frequency")
save_panel(FeaturePlot(obj, features = "Gpr183", order = TRUE, min.cutoff = "q05", max.cutoff = "q95") + theme_manuscript(), args$output, "SuppFig4d_Gpr183_featureplot", 5, 4)

micro <- subset(obj, subset = cell_type == "Microglia")
gpr <- mean_expression_by_sample(micro, "Gpr183", c("sample_id", "condition"))
write_source_data(gpr, args$output, "SuppFig4e_Gpr183_microglia")
p <- ggplot(gpr, aes(condition, Gpr183, fill = condition)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = .08)) + labs(x = NULL, y = "Mean normalized Gpr183 expression") + theme_manuscript() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
save_panel(p, args$output, "SuppFig4e_Gpr183_microglia", 4.5, 3.5)
