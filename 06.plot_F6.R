
# F6C, F6HI, F6J, F6K

library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(viridis)
library(ggridges)


########### F6C ############

plot1 <- fread("../F6C.t01.SD.multi.incor.txt")

plot1 <- all %>%
  group_by(species, te) %>%
  summarize(count = n()) %>%
  ungroup()

plot1 # 1481  3

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text = element_text(size=14,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = "top",
                 #legend.justification = c(1,1),
                 legend.box.just = "top",
                 legend.margin = margin(6, 6, 6, 6), 
                 legend.text = element_text(size = 10, family="ArialMT"),
                 axis.text.x = element_text(size = 14,  family="ArialMT", angle = 0, vjust = 0.5, hjust = 0.5),
                 axis.text.y = element_text(size = 14,  family="ArialMT"))

plot1$log_num <- log2(plot1$count+1)
plot1$log_num %>% max()

p1 <- ggplot(plot1, aes(x = log_num, y = species, fill = species)) + 
  geom_density_ridges(adjust = 0.5) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  scale_fill_manual(values = c("#CE8C99","#819FC2")) +
  #geom_signif(comparisons = list(c("human", "mouse"))) +
  labs(x = "Log2(No.Isoforms+1)", y = "Density", title = "Multi-loci incorporation", fill="") +mytheme+
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1))
# geom_signif(comparisons = list(c("human", "mouse")),
#             textsize = 5, vjust = 0.5, hjust = 0.5, map_signif_level = TRUE,
#             y_position = c(6,7,8)) 

p1 

ggsave(p1, filename = "./fig/F6C.new.multi.loci.pdf",width = 6,height = 5)
fwrite(plot1, file = "./table/F6C.t01.SD.multi.incor.txt",sep = "\t")


########### F6HI ########## 

halftime_final <- fread("../F6HI.txt")

color1 <- c("#DC0000" ,"#0072B5", "#808180")
color2 <- c("#808180" ,"#0072B5", "#E18727FF")
mycolor <- c("#E3786C" ,"#A9AEAE")

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text = element_text(size=14,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = "right",
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 10, family="ArialMT"),
                 axis.text.x = element_text(size = 14,  family="ArialMT", angle = 0, vjust = 0.5, hjust = 0.5),
                 axis.text.y = element_text(size = 14,  family="ArialMT"))

halftime_final$derive <- factor(halftime_final$derive, levels = c("TE-derived","Non-TE-derived"))

derive_levels <- unique(halftime_final$derive)
ks_results1 <- data.frame(Group1 = character(), Group2 = character(),
                          p.value = numeric(), stringsAsFactors = FALSE)

for (i in 1:(length(derive_levels) - 1)) {
  for (j in (i + 1):length(derive_levels)) {
    group1 <- derive_levels[i]
    group2 <- derive_levels[j]
    ks_test <- ks.test(halftime_final$halftime[halftime_final$derive == group1],
                       halftime_final$halftime[halftime_final$derive == group2])
    ks_results1 <- rbind(ks_results1, data.frame(Group1 = group1, Group2 = group2, p.value = ks_test$p.value))
  }
}

###### box

derive_medians <- aggregate(halftime ~ derive, data = halftime_final, FUN = median)
p1 <- ggboxplot(halftime_final, x = "derive", y = "halftime", fill = "derive",
                palette = mycolor, outlier.shape = NA) +
  labs(x = "\n", y = "Pseudo-mRNA half-time", title = "wilcox.test") +
  mytheme +
  theme(strip.text = element_text(size = 18),
        axis.text.x = element_blank()) +
  stat_compare_means(comparisons = list(c("TE-derived", "Non-TE-derived")),
                     method = "wilcox.test", aes(label = ..p.adj..)) +
  geom_text(data = derive_medians, aes(x = derive, y = -Inf, label = sprintf("%.2f", halftime)),
            vjust = -1, color = "black", size = 5)

p1

###### accumulateed 
p2 <- ggplot(halftime_final, aes(x = halftime, color = derive)) +
  stat_ecdf(size = 1) +
  scale_x_continuous(limits = c(0, 15), expand = c(0, 0), breaks = seq(0, 15, by = 3)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(values = mycolor) +
  labs(x = "Pseudo-mRNA half-time (h)", y = "Cumulative fraction",title = "Mouse oocytes with ActD treatment, 17390 iso") +
  theme_minimal() +
  mytheme  +
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.position = "bottom",         
        legend.title = element_blank(),
        plot.margin = margin(1, 1, 1,1 , "cm"))  +
  geom_text(data = ks_results1, 
            aes(x = 3, y = 0.9, family="ArialMT", label = paste("p = ", sprintf("%.2e", ks_results1$p.value))),
            color = "black", size = 5, hjust = 0)

p2



########  TSS classification 
rep_classification <- unique(halftime_final$classification)
comparisons <- list()

for (i in 1:(length(rep_classification) - 1)) {
  for (j in (i + 1):length(rep_classification)) {
    comparisons <- append(comparisons, list(c(rep_classification[i] %>% as.character(),
                                              rep_classification[j] %>% as.character())))
  }
}

class_medians_s <- aggregate(halftime ~ classification, data = halftime_final, FUN = median)

p3 <- ggboxplot(halftime_final, x = "classification", y = "halftime", fill = "classification",
                palette = c("#DCC53B", "#38AF75", "#336589", "#A9AEAE"), outlier.shape = NA) +
  #facet_wrap(~coding, ncol = 10, scales = "free_y")+
  labs(x = "\n", y = "Pseudo-mRNA half-time", title = "wilcox.test") +
  mytheme +
  stat_compare_means(comparisons = list(c("exonization", "ntd"),c("ntd", "TES"),c("TSS", "ntd")),
                     method = "wilcox.test", aes(label = ..p.adj..)) +
  geom_text(data = class_medians_s, aes(x = classification, y = -Inf, label = sprintf("%.2f", halftime)),
            vjust = -1, color = "black", size = 4, angle= 270)+
  scale_y_continuous(breaks = seq(floor(min(halftime_final$halftime)),
                                  ceiling(max(halftime_final$halftime)),
                                  by = 5))+
  coord_flip()  

p3


ggsave(p1, filename="./F6H.p02.box.ACTD.pdf",width = 5,height = 6)
ggsave(p2, filename="./F6H.cucum.ACTD.pdf",width = 8,height = 8)
ggsave(p3, filename="./F6I.box.TES.ACTD.pdf",width = 15,height = 4)
fwrite(sd_final_final[,c("halftime","classification","derive")], file="./table_260128//F6HI.txt",sep = "\t")

########### F6J rbp binding ########


plot <- fread("../F6J.txt")
plot$rbpname <- factor(plot$rbpname, levels = plot$rbpname)
#plot_final <- top_10[]

plot$facet_row <- "v"
plot[1:5,]$facet_row <- 1
plot[6:10,]$facet_row <- 2

p1 <- ggplot(plot, aes(x = rbpname, y = species, size = row_mean, color = log10(as.numeric(`adj_p-value`)))) +
  geom_point() +
  scale_size_continuous(range = c(8, 15)) +
  scale_color_gradient(low="red", high="#F4D020") +
  facet_wrap(~facet_row,  nrow = 2, scales = "free")+  # 设置颜色渐变
  labs(x = "Enriched RBP", y = "", size = "Mean nomalized TPM", color = "Log10(P-adj)",
       title = "Human TES Representive 10 RBP, Bonferroni Correction, Mean nomalized TPM, 250918") +
  theme_minimal() +
  mytheme+
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(vjust = -0.5),  # 可选：调整 x 轴标题位置
    axis.title.y = element_text(vjust = 0.5),     # 可选：调整 y 轴标题位置
    legend.position = "bottom",                  # 设置图例位置
    legend.direction = "horizontal",
    strip.text = element_blank()                  # 
  )

p1

fwrite(top_10, file = "./table/F6J.t03.SD.ho.tes.RBP.txt",sep = "\t")
ggsave(p1, filename="./fig/F6J.ho.RBP.TES.padj.norm.TPM.pdf",width = 10,height = 4)




########### F6K GO analysis ########

df <- fread("../F6K.txt")
df <- df %>% select(stage, Description, LogP, Count) %>% 
  group_by(Description) %>%
  filter(LogP == min(LogP)) %>%
  ungroup() %>% 
  arrange(stage, LogP) %>% setDT()

df$Description %>% unique()

df$Description <- factor(df$Description, levels = df$Description)

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text = element_text(size=14,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = "right",
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6), 
                 legend.text = element_text(size = 10,family="ArialMT"),
                 axis.text.x = element_text(size = 14,  family="ArialMT", angle = 0, vjust = 0.5, hjust = 0.5),
                 axis.text.y= element_text(size = 12,  family="ArialMT"))

df$Description <- factor(df$Description, levels = df$Description)

p1 <- ggplot(df, aes(stage, Description)) +
  geom_point(aes(size = Count, color = LogP), shape=16) +
  scale_size_continuous(range = c(8, 15)) +
  scale_colour_gradient(low="red", high="#FFD700") +
  #facet_wrap(~type, ncol = 1) + 
  labs(size="Counts", color="log10(P-value)", x="Species", y="", title = "F6K.TE-derived host genes") + 
  mytheme +
  theme(strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.9,size = 20),
        legend.background = element_rect(fill = "transparent", color = NA) ) 

p1

fwrite(df, file = "./table/F6K.t01.SD.enrich.select.txt",sep = "\t")
ggsave(p1, filename = "./fig/F6K.enrichment_select_te.derive.pdf",width = 9, height = 6)

fwrite(df, file = "./table_260128/F6K.txt",sep = "\t")








