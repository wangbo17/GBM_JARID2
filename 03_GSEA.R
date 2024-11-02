rm(list = ls())
load(file = 'step2output.Rdata')

library(dplyr)
library(data.table)
library(fgsea)

######## GMT file ########
read_gmt <- function(file) {
  gmt_lines <- readLines(file)
  gmt_list <- list()
  for (line in gmt_lines) {
    parts <- unlist(strsplit(line, "\t"))
    gene_set_name <- parts[1]
    genes <- parts[3:length(parts)]
    genes <- genes[genes != ""]
    gmt_list[[gene_set_name]] <- genes
  }
  return(gmt_list)
}

gmt_file <- "data/TFs_ENS_1000_GTRDv19_10_gencodev27.gmt"
gmt_data <- read_gmt(gmt_file)
names(gmt_data) <- trimws(names(gmt_data), which = "right")
gmt_data <- lapply(gmt_data, trimws, which = "both")
gmt_data <- gmt_data["JARID2"]


######## GSEA ########
set.seed(17)
nes_results <- data.frame(Patient = character(), NES = numeric(), stringsAsFactors = FALSE)
for (patient in names(log2fc_list)) {
  res <- log2fc_list[[patient]]
  ranks <- res$Log2FC
  names(ranks) <- res$Gene
  ranks <- sort(ranks, decreasing = TRUE)
  fgsea_res <- fgsea(pathways = gmt_data, stats = ranks, minSize = 15, nPermSimple = 50000)
  nes <- fgsea_res$NES[1]
  nes_results <- rbind(nes_results, data.frame(Patient = patient, NES = nes))
}


save(nes_results, file = "step3output.Rdata")
