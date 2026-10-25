HOMEOSTATIC <- c("P2ry12", "P2ry13", "Tmem119", "Cx3cr1", "Hexb", "Csf1r", "Sall1", "Siglech", "Fcrls", "Olfml3", "Mertk", "Tgfbr1", "Tgfbr2")
PRO_INFLAMMATORY <- c("Tnf", "Il1a", "Il1b", "Nlrp3", "Cd86", "Ccl4", "Tlr2", "Tlr4", "Nfkbia", "Relb", "Irf7")
M2_LIKE <- c("Arg1", "Mrc1", "Cd163", "Msr1", "Tgfb1", "Maf", "Mafb")
DISEASE_ASSOCIATED <- c("Trem2", "Apoe", "Lpl", "Cst7", "Clec7a", "Tyrobp", "Axl", "Itgax", "C1qb", "Ctsb", "Ctsd", "Lgals3", "Spp1", "Gpnmb")
PRIMING_CORE <- c("Cd83", "Egr1", "Nfkbiz", "Nlrp3", "Tnfaip3", "Il1b", "Tlr2", "Ccl4", "Ptgs2", "Hif1a")
CALCIUM_MOBILIZATION <- c("P2ry12", "P2rx4", "P2rx7", "Plcg2", "Orai1", "Orai2", "Calm1", "Camk1", "Camk2d", "Camk2g", "Nfatc3", "Cfl1", "Ptk2b", "Cdc42", "Msn", "Actr2")

add_signature_scores <- function(object) {
  modules <- list(HOMEOSTATIC, PRO_INFLAMMATORY, M2_LIKE, DISEASE_ASSOCIATED, PRIMING_CORE)
  names_out <- c("Homeostatic", "ProInflammatory", "M2Like", "DiseaseAssociated", "PrimingCore")
  for (i in seq_along(modules)) {
    genes <- intersect(modules[[i]], rownames(object))
    if (length(genes) >= 2) {
      object <- Seurat::AddModuleScore(object, features = list(genes), name = names_out[[i]], seed = 1)
      object[[paste0(names_out[[i]], "_score")]] <- object[[paste0(names_out[[i]], "1")]][, 1]
    }
  }
  object
}
