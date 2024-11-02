rm(list = ls())
load(file = 'step1output.Rdata')

library(dplyr)

######## log2FC ########
log2fc_results <- data.frame()
unique_patients <- unique(colData$Patient)
for (patient in unique_patients) {
  patient_samples <- colData[colData$Patient == patient, ]
  
  primary_sample <- patient_samples$ID[patient_samples$Stage == "Primary"]
  recurrent_sample <- patient_samples$ID[patient_samples$Stage == "Recurrent"]

  primary_expr <- countData[, primary_sample]
  recurrent_expr <- countData[, recurrent_sample]
  
  log2fc <- log2((recurrent_expr + 0.01) / (primary_expr + 0.01))
  log2fc_df <- data.frame(Gene = rownames(countData), Log2FC = log2fc)
  log2fc_df$Patient <- patient
  rownames(log2fc_df) <- NULL
  log2fc_results <- rbind(log2fc_results, log2fc_df)
}


log2fc_list <- split(log2fc_results, log2fc_results$Patient)
save(log2fc_list, file = "step2output.Rdata")
