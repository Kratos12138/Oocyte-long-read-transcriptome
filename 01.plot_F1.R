
# Fig1

library(ggplot2)
library(data.table)
library(tidyverse)
library(furrr)
library(ggpubr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ComplexHeatmap)


#----figure1 dominant isoform
gene_expr_quint_colors <- c( 
  "#dae8f0","#BFD8E7","#6AACD5","#2F82BE","#07529C")

p<-ggplot(dominant_dt, aes(x = dominant_fraction, fill = factor(gene_expr_quintile)) )+
  geom_histogram(alpha = 0.75,binwidth = .05) +
  scale_fill_manual(values = gene_expr_quint_colors) +
  labs(x = "Dominant isoform Expression fraction", y = "Gene count",
       title = "Dominant isoform Expression",fill = "Gene Expression \n Quintile" ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "grey40", linewidth = 0.8), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
    legend.position = "right"
  )
p
ggsave("Fig1C.dominant_bin0.05.pdf",p,width=6,height =4)


df<-fread(".../n_trans.csv")
pdf("Fig1.D.new2.pdf",width = 12,height =10)
Fig2.1 <- df %>% 
  ggplot(aes(x = gene_rank, y = n_transcripts , color=n_transcripts_category2)) + 
  geom_point(size=2) + 
  geom_label_repel(
    data=df %>% filter(oocyte==TRUE),
    aes(label = gene_name),
    force = 20, direction='both', max.overlaps = 100,
    size=2.5, show.legend = FALSE
  ) +
  scale_color_manual(values=my_cat_color3) +
  scale_y_log10(limits = c(1, NA), breaks = c(1, 10, 100)) + 
  
  scale_x_log10(limits = c(NA, 10000)) +
  theme_bw() +
  annotation_logticks() +
  labs(x="Gene rank", y="unique isoform number") +
  ggtitle("oocyte maturation genes ~ transcripts ")
Fig2.1
dev.off()


