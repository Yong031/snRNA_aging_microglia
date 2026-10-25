suppressPackageStartupMessages({library(CellChat); library(Seurat); library(ggplot2)})
source(file.path("R", "utils.R"))
args <- script_args("Creates Figure 6e/f from precomputed per-condition CellChat objects stored in an RDS named list.")
cellchat_list <- readRDS(args$input)
if (!is.list(cellchat_list) || is.null(names(cellchat_list))) stop("Input must be a named list of CellChat objects", call. = FALSE)
required <- c("Sham_Aged", "SNI_Aged", "SNI_Aged_M3")
missing <- setdiff(required, names(cellchat_list)); if (length(missing)) stop("Missing CellChat groups: ", paste(missing, collapse = ", "))

# Fig. 6e: outgoing networks from Microglia. Use the same maximum edge weight for every condition.
weights <- unlist(lapply(cellchat_list[required], function(x) x@net$weight["Microglia", ]))
edge_max <- max(weights, na.rm = TRUE)
pdf(file.path(args$output, "Fig6e_microglia_outgoing_networks.pdf"), width = 12, height = 4)
par(mfrow = c(1, 3), xpd = TRUE)
for (nm in required) {
  mat <- cellchat_list[[nm]]@net$weight
  keep <- matrix(0, nrow(mat), ncol(mat), dimnames = dimnames(mat)); keep["Microglia", ] <- mat["Microglia", ]
  netVisual_circle(keep, weight.scale = TRUE, edge.weight.max = edge_max, title.name = nm)
}
dev.off()

# Fig. 6f: descriptive total outgoing probability; no error bars without biological replicates.
totals <- data.frame(condition = required, outgoing_probability = vapply(cellchat_list[required], function(x) sum(x@net$weight["Microglia", ], na.rm = TRUE), numeric(1)))
write_source_data(totals, args$output, "Fig6f_total_microglia_outgoing_probability")
p6f <- ggplot(totals, aes(condition, outgoing_probability, fill = condition)) + geom_col(width = .7) + labs(x = NULL, y = "Sum of communication probabilities") + theme_manuscript() + theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
save_panel(p6f, args$output, "Fig6f_total_microglia_outgoing_probability", 4.2, 3.5)
