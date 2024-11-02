rm(list = ls())
load(file = 'step3output.Rdata')
load(file = 'step4output.Rdata')

library(ggplot2)
library(ggrepel)

######## Plot ########
nes_results <- na.omit(nes_results)

merged_results <- merge(nes_results, pc1_score_df, by = "Patient")
merged_results$Responder <- ifelse(merged_results$NES > 0, "Positive", "Negative")

p <- ggplot(merged_results, aes(x = PC1_Score, y = NES, label = Patient, color = Responder)) +
  geom_point(size = 6, alpha = 0.8) + 
  geom_text_repel(size = 3, color = "black", 
                  segment.color = "grey", segment.size = 0.5) + 
  scale_color_manual(values = c("Positive" = "#1f77b4", "Negative" = "#d62728")) + 
  theme_minimal(base_size = 15) + 
  labs(x = "log2FC PC1", 
       y = "JARID2 NES") + 
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    axis.title = element_text(face = "bold", size = 12), 
    panel.grid.major = element_line(color = "grey80"), 
    panel.grid.minor = element_line(color = "grey80"), 
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  ) +
  scale_x_continuous(limits = c(-150, 150), expand = expansion(mult = c(0, 0.05))) + 
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") + 
  geom_vline(xintercept = 0, color = "black", linetype = "solid") 


ggsave("NES_vs_PC1_plot.png", plot = p, width = 10, height = 5, dpi = 600)
