
#  F3B, F3C, F3D

library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(viridis)

###### F3A.category  #####

pct <- fread("../F3A.txt")

pct$Structural_Category %>% table()
pct$Structural_Category <- factor(pct$Structural_Category,
                                  levels = c("FSM", "ISM","NIC", "NNC",
                                             "Fusion", "Genic", "Antisense", "Intergenic"))
unique(pct$stage)
#    "m_GV"   "m_GVBD" "m_MI"   "m_MII" 
pct$stage <- factor(pct$stage, levels = c("m_GV"  , "m_GVBD" ,"m_MI"   ,"m_MII"))

#pct$sp <- ifelse(grepl("h_",pct$stage),"Human","Mouse")

my_colors <-  c("#D53E4F","#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD")

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=14,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = "right",
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6), 
                 legend.text = element_text(size = 14,family="ArialMT"),
                 axis.text.x= element_text(size=14,  family="ArialMT",angle = 0, vjust=0.5),
                 axis.text.y= element_text(size=14,  family="ArialMT"))

text <- pct[pct$Structural_Category %in% c("FSM","ISM","NIC","NNC"),]
text$Structural_Category <- as.character(text$Structural_Category)

p1 <- ggplot(pct, aes(x = stage, y = percentage, fill = Structural_Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  #facet_wrap(~sp)+
  #scale_y_continuous(labels = scales::percent) +
  labs(x = "Stage", y = "Percentage", title = "F3A GENCODE VM37") +
  scale_fill_manual(name = "Structural Category", values = rev(my_colors)) +
  geom_text(data = pct,  aes(x=stage, y=percentage, ymax = percentage,label = paste0(round(percentage,1), "%")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4) +  
  geom_text(aes(y = 102, label = paste0("n=", stage_all)), 
            data = pct, size = 4, vjust = -0.5, hjust = 0.5) +
  mytheme


p1

ggsave(p1, filename = "./fig/F3A.pdf",width = 6.5, height = 5.5)
#fwrite(pct, file = "./table/F3A.SD.txt",sep = "\t")
fwrite(pct[,!colnames(pct) %in% c("count")], file = "./table_260128/F3A.txt",sep = "\t")




###### F3B.stackbarplot #######

tmp <- fread("../F3B.txt")

stage_colors <- c("#C96470","#E398A2","#F9D9D4",
                  "#376BAB","#5E88C6","#8EABD4","#C9D9E8")

mytheme1 <- theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  text=element_text(size=14,  family="ArialMT"),
                  axis.title.x = element_text(vjust=-0.5, colour = "black"),
                  axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  legend.position = "right",
                  #legend.justification = c(1,1),
                  legend.box.just = "right",
                  legend.margin = margin(6, 6, 6, 6), 
                  legend.text = element_text(size = 14,family="ArialMT"),
                  axis.text.x= element_text(size=14,  family="ArialMT"),
                  axis.text.y= element_text(size=14,  family="ArialMT"))

tmp$Structural_Category <- factor(tmp$Structural_Category,
                                  levels = c("FSM", "ISM","NIC", "NNC","Fusion", "Genic", "Antisense", "Intergenic"))

tmp$ratio <- tmp$ratio_final*100

p1 <- ggplot(tmp, aes(x = Structural_Category, y = ratio, fill = stage)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Percentage of stages") +
  xlab("Structural_Category") +
  ggtitle("F3B GENCODE V48")+
  facet_wrap(~ species) +
  scale_fill_manual(values = stage_colors) +
  theme(legend.title = element_blank()) +
  #scale_y_continuous(labels = scales::percent) +
  mytheme1+
  theme(axis.text.x = element_text(size=12,  family="ArialMT", angle = 45, vjust = 1, hjust = 1))

p1

ggsave(p1, filename = "./fig/F3B.category.stack.bar.pdf",width = 8,height = 5)
#fwrite(tmp[,!c("ratio_final")], file = "./table/F3B.t01.SD.category.stack.bar.txt",sep = "\t")
fwrite(tmp[,!c("ratio_final","count")], file = "./table_260128/F3B.txt",sep = "\t")




###### F3C.stackbarplot #######
plot_samp <- fread("../F3C.txt")

plot_samp$Structural_Category %>% unique()
plot_samp <- plot_samp[!plot_samp$Structural_Category=="Genic_intron",]
plot_samp$Structural_Category <- factor(plot_samp$Structural_Category,
                                        levels = c("FSM", "ISM","NIC", "NNC","Fusion", "Genic", "Antisense", "Intergenic"))

sc_colors <- c("#D53E4F","#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD")

plot_samp <- plot_samp %>%
  group_by(stage, Structural_Category) %>%
  mutate(ratio_min = min(ratio, na.rm = TRUE),
         ratio_max = max(ratio, na.rm = TRUE),
         ratio_median = median(ratio, na.rm = TRUE) )


stage_colors <- c("#E56E6B", "#E7A36B", "#54C4D0", "#836CB0","#F3B2AF","#FF9F68","#A6DCE1")
mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=14,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = "right",
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6), 
                 legend.text = element_text(size = 14,family="ArialMT"),
                 #axis.text.x= element_text(size=14,  family="ArialMT",angle = 45，hjust=1),
                 axis.text.x = element_text(size=14,  family="ArialMT", angle = 45, vjust = 0.6, hjust = 0.5),
                 axis.text.y= element_text(size=14,  family="ArialMT"))


p3 <- ggplot(plot_samp, aes(x = stage, y = ratio_median*100, 
                            group = interaction(species, Structural_Category),
                            linetype = species, color = Structural_Category)) +
  #geom_ribbon(aes(ymin = ratio_min*100, ymax = ratio_max*100, fill = Structural_Category), alpha = 0.5) +
  geom_line(size = 1) +
  #facet_wrap(~ species) +
  scale_linetype_manual(values = c("Human Oocyte" = "solid", "Mouse Oocyte" = "dashed")) +
  labs(x = "Stage", y = "Ratio", color = "Structural Category", linetype = "Species", title = "F3C GENCODE V48") +
  theme(legend.position = "right")+
  scale_color_manual(values = rev(sc_colors)) + 
  mytheme

p3

ggsave(p3, filename = "./fig/F3C.category.dynamic.pdf",width = 8,height = 6)

sd <- plot_samp[,c("Structural_Category","stage","ratio_median")] %>% distinct(.,.keep_all = T)
sd$ratio_median <- sd$ratio_median*100
fwrite(sd, file = "./table_260128//F3C.txt",sep = "\t")


###### F3D.ISM stackbarplot #######
sd <- fread("F3D.txt")

plot_samp$subcategory %>% unique()
plot_samp$subcategory <- factor(plot_samp$subcategory,
                                levels = c("5prime_fragment","3prime_fragment", "internal_fragment", "intron_retention",  "mono-exon" ))

# colors
color_gradient <- colorRampPalette(c("#409675","white"))
colors <- color_gradient(6)
colors
subc_colors <- colors[1:5]

plot_samp <- plot_samp %>%
  group_by(stage, subcategory) %>%
  mutate(ratio_min = min(ratio, na.rm = TRUE),
         ratio_max = max(ratio, na.rm = TRUE),
         ratio_median = median(ratio, na.rm = TRUE) )

# 折线图
p4 <- ggplot(plot_samp, aes(x = stage, y = ratio_median*100, 
                            group = interaction(species, subcategory),
                            linetype = species, color = subcategory)) +
  #geom_ribbon(aes(ymin = ratio_min*100, ymax = ratio_max*100, fill = subcategory), alpha = 0.5) +
  geom_line(size = 1) +
  #facet_wrap(~ species) +
  scale_linetype_manual(values = c("Human Oocyte" = "solid", "Mouse Oocyte" = "dashed")) +
  labs(x = "Stage", y = "Ratio", color = "Structural Category", linetype = "Species", title = "F3D GENCODE V48") +
  theme(legend.position = "right")+
  scale_color_manual(values = subc_colors) +  
  mytheme

p4

# fwrite(plot_samp, file = "./table/F3D.SD.t01.ISM.dynamic.txt",sep = "\t")
# sd <- plot_samp[,c("subcategory","stage","ratio_median")] %>% distinct(.,.keep_all = T)
# sd$ratio_median <- sd$ratio_median*100
ggsave(p4, filename = "./fig/F3D.ISM.dynamic.pdf",width = 8,height = 6)
fwrite(sd, file = "./table_260128/F3D.txt",sep = "\t")


###### F3E.coding dist  #####################

plot_final <- fread("../F3E.txt")

my_colors <-  c("#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD","#D53E4F")

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=14,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = "none",
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6), 
                 legend.text = element_text(size = 14,family="ArialMT"),
                 axis.text.x= element_text(size=12,  family="ArialMT"),
                 axis.text.y= element_text(size=12,  family="ArialMT"))

Coding_prob <- plot_final[type=="Coding_prob",]
orf_length <- plot_final[type=="ORF_length",]

p1 <- ggplot(Coding_prob, aes(x = Coding_prob, y = Structural_Category, fill = Structural_Category)) + 
  geom_density_ridges() +
  facet_wrap(~species) +
  labs(x = "Coding probability (CPAT)", y = "Structural Category") + 
  mytheme +
  scale_fill_manual(values = my_colors) + 
  theme(strip.text = element_text(size = 15))

p1

orf_length <- orf_length[species=='human' & Coding_prob>0.364 | species=="mouse" & Coding_prob>=0.44,]

p2 <- ggplot(orf_length,  aes(x = ORF, y = Structural_Category, fill = Structural_Category)) +
  facet_wrap(~species)+
  geom_violin(scale = "width", adjust = 0.5, trim = TRUE) + 
  labs(x = "ORF nucleotide length (bp)", y = "Structural Category") + mytheme +
  scale_fill_manual(values = my_colors)+ theme(strip.text = element_text(size = 15) ) +
  scale_x_reverse()  

p2

pall <- ggarrange(p2,p1,ncol = 2)

pall

ggsave(pall, filename = "./fig/F3E.coding.prob.ORF.length.MO.pdf",width = 12, height = 5)

colnames(plot_final)[3] <- "value"
colnames(orf_length)[4] <- "value"

cpat <- rbind(mutate(plot_final[,c(2,3)], type="Coding_prob"),
              mutate(orf_length[,c(2,4)], type="ORF_length"))

cpat %>% dim() # 78249     3

#fwrite(cpat, file = "./table/F3E.t01.SD.cpat.MO.txt",sep = "\t")
fwrite(cpat, file = "./table_260128/F3E.txt",sep = "\t")


###### F3F.category ratio in coding   #########

pct1 <- fread("../F3F.txt")

#pct$sp <- ifelse(grepl("h_",pct$stage),"Human","Mouse")
colorRampPalette(c("black","white"))(10)

my_colors <-  c("#c1bbb8", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD")

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=14,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = "right",
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 14,family="ArialMT"),
                 axis.text.x= element_text(size=14,  family="ArialMT",angle = 45, vjust=0.5),
                 axis.text.y= element_text(size=14,  family="ArialMT"))


p5 <- ggplot(pct1, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Structural_Category)) +
  geom_rect() +
  facet_wrap(~ stage, labeller = label_parsed)  +
  scale_fill_manual(name = "Structural Category", values = rev(my_colors)) +
  coord_polar(theta = "y") +
  xlim(c(2.5, 4)) + 
  theme_void() +
  theme(legend.position = "right", text = element_text(family = "ArialMT"),
        strip.text = element_text(size = 14, face = "bold")) +
  ggtitle("F3F Ratio of isoform categories in coding transcripts\n GENCODE V48")+
  geom_text(x = 3.5, aes(y = labelPosition, label = label), size = 5,
            color = "white", family = "ArialMT", fontface = "bold") +
  geom_text(data = pct1 %>% distinct(stage, stage_all, .keep_all = TRUE),
            aes(x = 3.5, y = 25, label = paste0("n=",stage_all)), 
            size = 4, color = "black", family = "ArialMT", fontface = "bold", hjust=2)

p5

ggsave(p5, filename = "./fig/F3F.class.ratio.in.coding.pdf",width = 8, height = 6)
#fwrite(pct1, file = "./table/F3F.SD.class.ratio.in.coding.txt",sep = "\t")
fwrite(pct1[,c("stage", "Structural_Category","percentage","stage_all")],
       file = "./table_260128/F3F.txt",sep = "\t")






###### F3G.category ratio in non-coding #########

#pct$sp <- ifelse(grepl("h_",pct$stage),"Human","Mouse")
colorRampPalette(c("black","white"))(10)

my_colors <- c("#D53E4F","#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD")

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=14,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = "right",
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 14,family="ArialMT"),
                 axis.text.x= element_text(size=14,  family="ArialMT",angle = 45, vjust=0.5),
                 axis.text.y= element_text(size=14,  family="ArialMT"))

pct1 <- fread("F3G.txt")

p6 <- ggplot(pct1, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Structural_Category)) +
  facet_wrap(~ stage, ncol = 1) +
  scale_fill_manual(name = "Structural Category", values = rev(my_colors)) + 
  geom_rect() +
  coord_polar(theta = "y") +
  xlim(c(2.5, 4)) +
  theme_void() +
  theme(legend.position = "right", text = element_text(family = "ArialMT"),
        strip.text = element_text(size = 14, face = "bold")) +
  ggtitle("F3G.Ratio of isoform categories\n in non-coding transcripts\n GENCODE V48") +
  geom_text(x = 3.5, aes(y = labelPosition, label = label), size = 5,
            color = "white", family = "ArialMT", fontface = "bold") +
  geom_text(data = pct1 %>% distinct(stage, stage_all, .keep_all = TRUE),
            aes(x = 3.5, y = 25, label = paste0("n=", stage_all)), 
            size = 4, color = "black", family = "ArialMT", fontface = "bold", hjust=2)

p6


ggsave(p6, filename = "./fig/F3G.class.ratio.in.non-coding.pdf",width = 8, height = 6)
#fwrite(pct1, file = "./table/F3G.SD.class.ratio.in.non-coding.txt",sep = "\t")
fwrite(pct1[,c("stage","Structural_Category","percentage","stage_all")], 
       file = "./table_260128/F3G.txt",sep = "\t")



