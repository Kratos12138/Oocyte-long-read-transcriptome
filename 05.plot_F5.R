
# F5B, F5C, F5D, F5E, F5F

library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(viridis)

####### F5B overview #####

plot_df <- fread("../F5B.txt")
plot_df
dim(plot_df)  # 36   7

plot_df$Type <- factor(plot_df$Type, levels = c("All","Known","Novel","NIC","NNC","TE-derived"))
plot_df$classification <- factor(plot_df$classification, levels = c("TSS","TES","exonization","ntd"))

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
                 legend.text = element_text(size = 10,family="ArialMT"),
                 axis.text.x = element_text(size = 14,  family="ArialMT", angle = 45, vjust = 1, hjust = 1),
                 axis.text.y= element_text(size = 14,  family="ArialMT"))

p1 <- ggplot(plot_df, aes(x = Type, y = Percentage, fill = classification, label = Ratio)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = Percentage, label = Ratio), data = plot_df, size = 5, position = position_stack(vjust = 0.5), color = "white") +
  geom_text(aes(y = 100, label = paste0("n=",Total_count)), data = distinct(plot_df, Type,.keep_all = T),
            size = 4, vjust = -0.5, hjust = 0.5) +
  scale_fill_manual(values = c("#DCC53B", "#38AF75", "#336589", "#7E7D7E")) +
  facet_wrap(~ species, scales = "free_y", ncol = 2) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Type", y = "Percentage (%)", fill = "Classification",title = "F5B GENCODE V48 & VM37") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 15, family = "ArialMT")) +
  mytheme


p1

plot_df %>% dim() # 43  7

ggsave(p1, filename="./fig/F5B.HO.pdf",width = 12,height = 6)
fwrite(plot_df[,!colnames(plot_df) %in% c("Count","Ratio")],
       file = "./table_260128/F5B.txt",sep = "\t")


####### F5C position coding ######
sum <- fread("../F5C.txt")

sum$classification <- factor(sum$classification, levels = c("TSS", "TES", "exonization"))
sum$coding <- factor(sum$coding, levels = c("Coding", "Non-coding"))

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
                 axis.text.x= element_text(size=14,  family="ArialMT", angle = 45, vjust = 0.8, hjust = 0.8),
                 axis.text.y= element_text(size=14,  family="ArialMT"))

p2s <- ggplot(sum, aes(x = coding, y = percentage, fill = classification)) +
  facet_wrap(~species, scales="free_x")+
  geom_bar(stat = "identity", position = "stack", width = 0.8, color = "black", size = 0.3) +  
  scale_fill_manual(values = c("#DCC53B", "#38AF75", "#336589"), name = "", labels = c("TSS", "TES", "exonization")) + 
  labs(x = "", y = "Percentage", title = "Human vs Mouse oocytes") +  
  geom_text(aes(label = paste0(round(percentage * 100, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 4) +  
  mytheme +  # 自定义主题
  theme(strip.text = element_text(size = 12)) + 
  geom_text(data = sum %>% distinct(species,coding,.keep_all = T), 
            aes(x = coding, y = 1.02, label = paste0("n=", total)), 
            vjust = 0, size = 4.5, color = "black") 

print(p2s)

ggsave(p2s, filename = "./fig/F5C.te_derive.coding.pdf",width = 6, height = 6)

sum$percentage <- sum$percentage*100
fwrite(sum[,!colnames(sum) %in% c("count")],
       file = "./table_260128//F5C.txt",sep = "\t")



####### F5D ho family #####

result <- fread("../F5D.txt")
result$repFamily

### set others, set colors
specified_repFamily <- c("Alu", "MIR", "L1", "ERVL-MaLR", "L2", "ERV1", 
                         "ERVL", "TcMar-Tigger", "hAT-Charlie", "CR1", "ERVK", "Simple_repeat")

result[!repFamily %in% specified_repFamily,]$repFamily <- "Others"
result$repFamily %>% unique() %>% length() # 13

result$repFamily <- factor(result$repFamily, levels= c("Alu", "MIR", "L1", "ERVL-MaLR", "L2", "ERV1","ERVL",
                                                               "TcMar-Tigger", "hAT-Charlie", "CR1", "ERVK",
                                                               "Simple_repeat","Others"))

result[result$classification=="TES",] %>% arrange(-percentage)

# coul <- brewer.pal(11, "RdYlBu") 
# coul <- colorRampPalette(coul)(length(unique(result$repFamily)))

coul <- c("#9D1E28", "#B23137", "#C85D49", "#DF8B5C", "#F5BA70", "#E7C28E",
          "#DACDAE", "#CAD6CC", "#BADCE9", "#96B1D1", "#7286B8", "#505CA0", "#2D348C")

custom_colors <- c(
  "Alu" = "#9D1E28",
  "MIR" = "#B23137",
  "L1" = "#C85D49",
  "ERVL-MaLR" = "#DF8B5C",
  "L2" = "#F5BA70",
  "ERV1" = "#E7C28E",
  "ERVL" = "#DACDAE",
  "TcMar-Tigger" = "#CAD6CC",
  "hAT-Chailie" = "#BADCE9",
  "CR1" = "#96B1D1",
  "ERVK" = "#7286B8",
  "Simple repeat" = "#505CA0",
  "Others" = "#2D348C"
)

result$classification <- factor(result$classification, levels = c("TSS","TES","exonization"))

p3 <- ggplot(result, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = repFamily)) +
  geom_rect() +
  facet_wrap(~ classification, labeller = label_parsed) +
  geom_text(x = 3.5, aes(y = labelPosition, label = label), size = 5,
            color = "white", family = "ArialMT", fontface = "bold") +
  # Add total label in the center of each segment
  geom_text(data = result %>% distinct(classification, total, .keep_all = TRUE),
            aes(x = 3.5, y = 25, label = paste0("n=",total)), 
            size = 5, color = "black", family = "ArialMT", fontface = "bold", hjust=2) +
  #scale_fill_manual(values = coul) +
  scale_fill_manual(values = custom_colors) +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) + theme_void() +
  theme(legend.position = "right", text = element_text(family = "ArialMT"),
        strip.text = element_text(size = 14, face = "bold")) +
  ggtitle("F5D V48")


p3

ggsave(p3, filename = "./fig/F5D.family.dis.HO.pdf", colormodel="cmyk", width = 13,height = 5)
#fwrite(result, file = "./table/F5D.family.dis.HO.txt",sep = "\t")

fwrite(result[,c("repFamily","classification","percentage","total")], 
       file = "./table_260128//F5D.txt",sep = "\t")



####### F5E mo family #####
result <- fread("../F5E.txt")
result$repFamily
result[result$classification=="TES",] %>% arrange(-percentage)

### set others, set colors
specified_repFamily <- c("ERVK", "Alu", "Unknown", "ERVL-MaLR", "B2", "B4", "L1",
                         "hAT-Charlie", "ERVL", "ERV1", "L1-dep", "MIR", "LTR")

result[!repFamily %in% specified_repFamily,]$repFamily <- "Others"
result$repFamily %>% unique() %>% length() # 14

result$repFamily <- factor(result$repFamily, levels= c("ERVK", "Alu", "Unknown", "ERVL-MaLR", "B2", "B4", "L1",
                                                               "hAT-Charlie", "ERVL", "ERV1", "L1-dep", "MIR", "LTR","Others"))


coul <- brewer.pal(11, "RdYlBu") 
coul <- colorRampPalette(coul)(length(unique(result$repFamily)))

result$classification <- factor(result$classification, levels = c("TSS","TES","exonization"))

p3 <- ggplot(result, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = repFamily)) +
  geom_rect() +
  facet_wrap(~ classification, labeller = label_parsed) +
  geom_text(x = 3.5, aes(y = labelPosition, label = label), size = 5,
            color = "white", family = "ArialMT", fontface = "bold") +
  # Add total label in the center of each segment
  geom_text(data = result %>% distinct(classification, total, .keep_all = TRUE),
            aes(x = 3.5, y = 25, label = paste0("n=",total)), 
            size = 5, color = "black", family = "ArialMT", fontface = "bold", hjust=2) +
  scale_fill_manual(values = coul) +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) + theme_void() +
  theme(legend.position = "right", text = element_text(family = "ArialMT"),
        strip.text = element_text(size = 14, face = "bold")) +
  ggtitle("F5E GENCODE VM37")


p3

ggsave(p3, filename = "./fig/F5E.family.dis.MO.pdf", colormodel = "cmyk", width = 13,height = 5)
#fwrite(result, file = "./table/F5E.family.dis.MO.txt",sep = "\t")
fwrite(result[,c("repFamily","classification","percentage","total")],
       file = "./table_260128/F5E.txt",sep = "\t")


####### F5F class ##########

pct0 <- fread("../F5F.txt")

pct0

pct0$classification <- factor(pct0$classification, levels = c("TSS","TES","exonization"))

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=14,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = "bottom",
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6), 
                 legend.text = element_text(size = 14,family="ArialMT"),
                 axis.text.x= element_text(size=14,  family="ArialMT", angle = 45, vjust = 0.8, hjust = 0.8),
                 axis.text.y= element_text(size=14,  family="ArialMT"))

p4 <- ggplot(pct0, aes(x = classification, y = percentage, fill = repClass)) +
  facet_wrap(~species, scales = "free",ncol = 1)+
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(x = "", y = "Percentage",title = "Derive locus, total, repClass GENCODE V48 & VM37") +
  scale_fill_manual(name = "repClass", values = c("LTR" = "#FFFFB3", "DNA" = "#FB8072", "SINE" = "#8DD3C7",  "LINE" = "#BEBADA", 
                                                  "Others" = "#80B1D3", "Unknown" = "#FCCDE5")) +
  geom_text(data = pct0,  aes(x = classification, y = percentage, ymax = percentage, label = paste0(round(percentage,1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "black", size = 4) +  
  geom_text(aes(y = 101, label = paste0("n=", stage_all)), 
            data = pct0, size = 5, vjust = -0.5, hjust = 0.5) +
  guides(fill = guide_legend(reverse = F)) +
  theme(strip.text = element_text(size = 18)) +
  mytheme

p4

ggsave(p4, filename = "./fig/F5F.class.dis.HO.MO.pdf",width = 5,height = 12)
fwrite(pct0[,!colnames(pct0) %in% c("count")],
       file = "./table_260128/F5F.txt",sep = "\t")

