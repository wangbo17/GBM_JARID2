rm(list = ls())
load(file = 'step0output.Rdata')

library(edgeR)

######## FPKM ########
calculate_fpkm <- function(colData, countData, gene_lengths) {
  dge <- DGEList(counts = countData, group = factor(colData$Stage, levels = c("Primary", "Recurrent")))
  cpm_values <- cpm(dge, normalized.lib.sizes = FALSE)
  gene_lengths <- gene_lengths[match(rownames(countData), gene_lengths$ID), ]
  gene_lengths_kb <- gene_lengths$Length / 1000
  fpkm_values <- sweep(cpm_values, 1, gene_lengths_kb, FUN = "/")
  return(fpkm_values)
}

fpkm_countData <- calculate_fpkm(colData, countData, gene_lengths)


save(colData, countData, fpkm_countData,
     file = "step1output.Rdata")
