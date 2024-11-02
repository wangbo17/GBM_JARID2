rm(list = ls())
load(file = 'step1output.Rdata')
load(file = 'step2output.Rdata')
load(file = 'data/pca_gbm_idhwt_rt_tmz_local_JARID2_log2fc_all_FALSE.Rdata')

######## PC1 ########
common_genes <- intersect(rownames(pca$rotation), unique(unlist(lapply(log2fc_list, function(x) x$Gene))))

pc1_rotation <- pca$rotation[common_genes, 1]

pc1_score_df <- data.frame(Patient = character(), PC1_Score = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(log2fc_list)) {
  patient_data <- log2fc_list[[i]]
  patient_data <- patient_data[patient_data$Gene %in% common_genes, ]
  patient_data <- patient_data[match(common_genes, patient_data$Gene), ]
  pc1_score <- sum(patient_data$Log2FC * pc1_rotation)
  patient_id <- unique(patient_data$Patient)
  pc1_score_df <- rbind(pc1_score_df, data.frame(Patient = patient_id, PC1_Score = pc1_score))
}


save(pc1_score_df, file = "step4output.Rdata")
