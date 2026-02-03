
# F2B, F2C, F2E

library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(viridis)
library(ggridges)


#####  F2B #######

pct <- fread("../F2B.SD.txt")

pct$Structural_Category %>% table()
pct$Structural_Category <- factor(pct$Structural_Category,
                                  levels = c("FSM", "ISM","NIC", "NNC",
                                             "Fusion", "Genic", "Antisense", "Intergenic"))
unique(pct$stage)
#    "h_GV"   "h_MI"   "h_MII" 
pct$stage <- factor(pct$stage, levels = c("h_GV", "h_MI", "h_MII"))

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
  labs(x = "Stage", y = "Percentage", title = "F2B GENCODE V48") +
  scale_fill_manual(name = "Structural Category", values = rev(my_colors)) +
  geom_text(data = pct,  aes(x=stage, y=percentage, ymax = percentage,label = paste0(round(percentage,1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "white", size = 4) + 
  geom_text(aes(y = 102, label = paste0("n=", stage_all)), 
            data = pct, size = 4, vjust = -0.5, hjust = 0.5) +
  mytheme


p1

ggsave(p1, filename = "./fig/F2B.pdf",width = 6, height = 6)
fwrite(pct[,!colnames(pct) %in% c("count")], file = "./table_260128/F2B.txt",sep = "\t")


#####  F2C #######

plot_final <- fread("../F2C.txt")
plot_final$Structural_Category <- factor(plot_final$Structural_Category,
                                         levels = rev(c("FSM", "ISM","NIC", "NNC","Fusion", "Genic", "Antisense", "Intergenic")))



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

Coding_prob <- plot_final[type="Coding_prob",]
orf_length <- plot_final[type="ORF_length",]

p2 <- ggplot(Coding_prob, aes(x = Coding_prob, y = Structural_Category, fill = Structural_Category)) + 
  geom_density_ridges() +
  #facet_wrap(~species) +
  labs(x = "Coding probability (CPAT)", y = "Structural Category",title = "F2C GENCODE V48") + 
  mytheme +
  scale_fill_manual(values = my_colors) + 
  theme(strip.text = element_text(size = 15))

p2

orf_length <- plot_final[Coding_prob>0.364,]

p2s <- ggplot(orf_length,  aes(x = ORF, y = Structural_Category, fill = Structural_Category)) +
  #facet_wrap(~species)+
  geom_violin(scale = "width", adjust = 0.5, trim = TRUE) + 
  labs(x = "ORF nucleotide length (bp)", y = "Structural Category",title = "F2C GENCODE V48") + mytheme +
  scale_fill_manual(values = my_colors)+ theme(strip.text = element_text(size = 15) ) +
  scale_x_reverse()  

p2s

pall <- ggarrange(p2s,p2,ncol = 2)

pall

coding_prob <- plot_final[!is.na(plot_final$Coding_prob),]
coding_prob$Coding_prob <- coding_prob$Coding_prob %>% as.numeric()
coding_prob %>% group_by(Structural_Category) %>% summarise(min_coding_prob=min(Coding_prob),
                                                            mean_coding_prob=mean(Coding_prob),
                                                            median_coding_prob=median(Coding_prob),
                                                            max_coding_prob=max(Coding_prob))


ggsave(pall, filename = "./fig/F2C.coding.prob.ORF.length.pdf",width = 12, height = 5)

colnames(plot_final)[2] <- "value"
colnames(orf_length)[3] <- "value"

cpat <- rbind(mutate(plot_final[,c(1:2)], type="Coding_prob"),
              mutate(orf_length[,c(1,3)], type="ORF_length"))

cpat %>% dim() # 82164     3

fwrite(cpat, file = "./table_260128/F2C.txt",sep = "\t")

#####  F2E.ISM  ########

plot_df1 <- fread("../F2E.txt")
plot_df1$percentage <- round(plot_df1$percentage,1)
str(plot_df1)

colors <- c("#489474","#62A188","#86B6A3","#ACCCBF","#D5E5DF")

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
                 axis.text.x = element_text(size=14,  family="ArialMT", angle = 0, vjust = 0.5, hjust = 0.5),
                 #axis.text.x = element_blank(),
                 axis.text.y = element_text(size=14,  family="ArialMT"))


p3 <- ggplot(plot_df1, aes(x = stage, y = percentage, fill = subcategory)) +
  #facet_wrap(~species,scales = "free")+
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(#title = "Structural Category and Type Percentage by Stage",
    x = "Stage",
    y = "Percentage",title = "F2E GENCODE V48") +
  scale_fill_manual(name = "subcategory", values = colors) +
  geom_text(aes(label = paste0(percentage, "%")),
            position = position_stack(vjust = 0.5), 
            color = "white", size = 4) +  
  mytheme+
  theme(strip.text = element_text(size = 20))+
  geom_text(aes(label = count_total, y=102),
            #position = position_stack(vjust = 0.5),  
            color = "black", size = 4, data = plot_df1 %>% distinct(.,stage,.keep_all = T))  

p3

plot_df1 %>% dim() # 15  5

ggsave(p3, filename = "./fig/F2E.ISM.h.pdf",width = 6,height = 5)
fwrite(plot_df1[,!colnames(plot_df1) %in% c("count")], file = "./table_260128/F2E.txt",sep = "\t")


