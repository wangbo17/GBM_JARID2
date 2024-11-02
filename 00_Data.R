rm(list = ls())

library(GenomicFeatures)
library(caret)

######## annotation ########
tx_db <- makeTxDbFromGFF("data/gencode.v46.chr_patch_hapl_scaff.annotation.gtf", format = "gtf")
exons_list_per_gene <- exonsBy(tx_db, by = "gene")

gene_lengths <- data.frame(
  ID = sub("\\..*", "", names(exons_list_per_gene)), 
  Length = sapply(exons_list_per_gene, function(exons) {
    sum(width(reduce(exons)))
  })
)


######## colData ########
meta_data <- read.csv(file = "data/meta_purity.csv", row.names = 1)
meta_data <- meta_data[
  meta_data$LibraryType == "Stranded_Total" & 
    meta_data$Local_Recurrence.. != "FALSE" & 
    as.logical(meta_data$Meets.Purity.Threshold) & 
    (meta_data$Non.Surgical.Treatment == "Radiotherapy and TMZ" | meta_data$Non.Surgical.Treatment == "Radiotherapy and TMZ +") & 
    meta_data$IDH == 0, 
]

colData <- data.frame(ID = row.names(meta_data), 
                      Stage = meta_data[, "Stage"],
                      Source = meta_data[, "Sample.Source"],
                      Patient = gsub("_.*", "", row.names(meta_data)))

######## countData ########
raw_data <- read.csv(file = "data/raw_data.csv", row.names = 1)
raw_data <- raw_data[, colnames(raw_data) %in% gene_lengths$ID]
raw_data <- raw_data[rownames(meta_data), , drop = FALSE]

countData <- na.omit(t(raw_data))


colData <- colData[colData$Source == "Stead", ]
countData <- countData[, colnames(countData) %in% colData$ID]
save(colData, countData, gene_lengths, file = "step0output.Rdata")
