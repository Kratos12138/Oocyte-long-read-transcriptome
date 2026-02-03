#----------DTU DEXSeq-------------

# F4

rm(list=ls())
setwd("")
outputpath=(".../fig4")
dir.create(outputpath,recursive = T)
library(ggplot2)
library(data.table)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(IsoformSwitchAnalyzeR)
library(ggpubr)



count<-fread(".../normalized_counts.csv")
count<-count%>%column_to_rownames(var="iso_id")
aslist <- importRdata(isoformCountMatrix=count,
                      designMatrix=sampleTable,
                      isoformExonAnnoation="human.gtf",
                      isoformNtFasta="human.fasta",
                      fixStringTieAnnotationProblem = TRUE,
                      showProgress = FALSE)
SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aslist,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE)

dtulist <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=TRUE
)
result <- dtulist $isoformSwitchAnalysis

significant_results <- result[result$padj < 0.05 &abs(result$dIF)>0.1, ]

significant_results<-significant_results%>%
  dplyr::rename( iso_id= isoform_id)%>%
  left_join(transcript_to_gene, by = "iso_id")

significant_results$cluster <- ifelse(significant_results$dIF > 0, "MII", "GV_MI")


write.csv(result,"fig4//dtulist_h_GV_MI.vs.MII.csv",row.names = F)
write.csv(significant_results,"fig4//dtulist_h_GV_MI.vs.MII_padj<0.05.dIF>0.1.csv",row.names = F)
saveRDS(dtulist, "fig4//dtulist_h.rds")


dte <- fread("fig4//DTE_human_MII.vs.GV+MI.txt")
deg <- fread("fig4//DEG_human_MII.vs.GV+MI.tpm>5.txt") 
dtu <- fread("fig4//dtulist_h_GV_MI.vs.MII_padj<0.05.dIF>0.1.csv")
dtu$gene_name %>%unique()%>%length()

library(VennDiagram)

venn1 <- venn.diagram(
  x = list(dte[cluster=="GV_MI",]$gene_name , deg[cluster=="GV_MI",]$gene_id, dtu[cluster=="GV_MI"]$gene_name ),
  fill = c("#7FB996", "#E4F48A","gray"),
  cat.col = rep("black", 3),
  category.names = c("GV_DTE", "GV_DEG", "GV_DTU"),
  cat.pos = c(-30, 30, 150), # Updated to have three values
  alpha = 0.6, lwd = 2, lty = 'blank', scaled = TRUE, cex = 1,
  fontfamily = "ArialMT", cat.cex = 1, main.cex = 1, cat.fontfamily = "ArialMT",
  filename = NULL, 
  main = list(label = "Overlap of GV-MI_up", fontfamily = "ArialMT"),
  main.fontfamily = "ArialMT",
  print.mode = c("percent", "raw"), inverted=T
)
pdf("fig4//p01.DTE.DEG.DTU.Human.pdf",width = 12,height = 6)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
# Specify the positions for each Venn diagram
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(venn1)
popViewport()


venn2 <- venn.diagram(
  x = list(dte[cluster=="MII",]$gene_name , deg[cluster=="MII",]$gene_id, dtu[cluster=="MII"]$gene_name ),
  fill = c("#7FB996", "#E4F48A", "gray"),
  cat.col = rep("black", 3),
  category.names = c("MII_DTE", "MII_DEG", "MII_DTU"),
  cat.pos = c(-30, 30, 150), # Updated to have three values
  alpha = 0.6, lwd = 2, lty = 'blank', scaled = TRUE, cex = 1,
  fontfamily = "ArialMT", cat.cex = 1, main.cex = 1, cat.fontfamily = "ArialMT",
  filename = NULL, 
  main = list(label = "Overlap of MII_up ", fontfamily = "ArialMT"),
  main.fontfamily = "ArialMT",
  print.mode = c("percent", "raw")
)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(venn2)
popViewport()
dev.off()


#----Fig4B
header<-fread("isoform.ann.txt")
header$novelty = as.character(header$Structural_Category)
header$novelty[header$novelty %in% c("Antisense", "Genic", "Fusion","Intergenic")] = "Other"
header$novelty = factor(header$novelty,levels=c("FSM", "ISM", "NIC", "NNC", "Other"))

idx = match(SwitchList_part1$isoformFeatures$isoform_id  , header$iso_rename )
SwitchList_part1$isoformFeatures$gene_name<-header$gene_name[idx]#gene_name
SwitchList_part1$isoformFeatures$iso_biotype<-header$Structural_Category[idx]
this_df =SwitchList_part1$isoformFeatures %>% 
  left_join(header %>% dplyr::select(isoform_id=iso_rename, novelty)) %>%
  filter(isoform_switch_q_value < .05,abs(dIF) > .1) %>% group_by(novelty) %>% summarise(switches=n_distinct(isoform_id))


this_df<-this_df %>% arrange(desc(novelty)) %>%
  mutate(prop = switches / sum(switches) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

Fig_cato = ggplot(this_df, aes(x="", y=prop, fill=novelty)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y",start=0) + 
  theme_void() + theme(legend.position = 'none', plot.title = element_text(hjust=.5)) +
  scale_fill_manual(values=colorVector) + 
  geom_text(aes(y=ypos, label = paste0(novelty,"\n(N=", switches, ")")), color = "BLACK", size=5) 
Fig_cato
ggsave(Fig_cato, file="F4B.DTU.type_count>0.1.pdf",width=6,height=6)


data_to_label1 =SwitchList_part1 $isoformFeatures %>%
  filter((-log10(isoform_switch_q_value) > 9 & (abs(dIF) > .3)))

gene_to_label<-SwitchList_part1$isoformFeatures[SwitchList_part1$isoformFeatures$gene_name %in%gene_to ,]%>%
  filter(( isoform_switch_q_value < 0.05 & (abs(dIF) > .1)))

data_to_label1<-rbind(data_to_label1,gene_to_label)

plot_data <- SwitchList_part1$isoformFeatures %>% 
  left_join(header %>% dplyr::select(isoform_id=iso_rename, novelty))

plot_data <- plot_data %>%
  mutate(color = case_when(
    isoform_switch_q_value > 0.05 | abs(dIF) < 0.1 ~ 'grey',
    TRUE ~ novelty  
  ))

Fig1 = ggplot(plot_data, 
              aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes(color=color), 
    size=2
    , alpha=.5
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed',color='red') + # default cutoff
  labs(x='difference in isoform fraction (dIF)', y='-log10 ( Isoform Switch Q Value )') +
  theme_bw() + xlim(-1,1) + scale_color_manual(values=colorVector) + 
  ggrepel::geom_text_repel(data = data_to_label1 %>% filter(dIF < 0),aes(label=gene_name, segment.size	= .1),size=3,force = 10, max.overlaps = 50,nudge_y = 5, nudge_x = -.1) + 
  ggrepel::geom_text_repel(data = data_to_label1 %>% filter(dIF > 0),aes(label=gene_name, segment.size	= .1),size=3,force = 10, max.overlaps = 50,nudge_y = 5, nudge_x = .1) + 
  theme(legend.position = 'none') +
  geom_vline(xintercept = c(0.1, -0.1),lty=1,color='grey')

ggsave(Fig1, file="F4B.DTUvolcano.pdf",width=6,height=4)



#----------
### Extract AA sequences
dtulist$aaSequence = NULL
isoformFeatures_part1 = dtulist$isoformFeatures

### Add DTE/DGE to switchList
idx = match(SwitchList_part1$isoformFeatures$isoform_id, DTE_results$iso_id)
SwitchList_part1$isoformFeatures$iso_q_value = DTE_results$padj[idx]

idx = match(SwitchList_part1$isoformFeatures$isoform_id, header$iso_rename)
SwitchList_part1$isoformFeatures$gene_id = header$associated_gene[idx]
SwitchList_part1$isoformFeatures$gene_name = header$gene_name[idx]

(DGE_results$gene_id%in%header$gene_name )%>%table()
idx = match(SwitchList_part1$isoformFeatures$gene_name, DGE_results$gene_id)
SwitchList_part1$isoformFeatures$gene_q_value = DGE_results$padj[idx]



dtulist$isoformFeatures = isoformFeatures_part1 %>%
  as_tibble() %>%
  group_by(gene_id) %>%
  mutate(
    isoform_switch_q_value = if_else(any(
      # our actual filtering criteria - genes with DTU, DTE, or DGE
      (isoform_switch_q_value < 0.05 & dIF > 0.1) | iso_q_value < 0.05 | gene_q_value < 0.05
    ), 0, 1),
    dIF = 1
  ) %>%
  ungroup() %>%
  as.data.frame()


#------ORF
SwitchList_part2 <- analyzeORF(
  SwitchList_part1,
  orfMethod = "longest",  
  showProgress=FALSE
)
orf_isoforms = SwitchList_part2 $orfAnalysis %>% as_tibble() %>%
  drop_na(orfTransciptStart) %>%
  pull(isoform_id)

head(SwitchList_part2$orfAnalysis, 3)
SwitchList_part2 <- extractSequence(
  switchAnalyzeRlist = SwitchList_part2,
  pathToOutput       = "/fig4/IsoformSwitchAnalyzeR_webtools",
  extractNTseq       = TRUE,
  extractAAseq       = TRUE,
  removeShortAAseq   = TRUE,
  removeLongAAseq    = FALSE,
  onlySwitchingGenes = TRUE,
  alsoSplitFastaFile =TRUE
)


SwitchList_part2 <-  analyzeCPC2(
  switchAnalyzeRlist   =SwitchList_part2,
  pathToCPC2resultFile = "IsoformSwitchAnalyzeR_webtools/cpc2_results.txt.txt",
  codingCutoff         = 0.364, # the coding potential cutoff we suggested for human
  removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)

SwitchList_part2<- analyzeSignalP(
  switchAnalyzeRlist       = SwitchList_part2,
  pathToSignalPresultFile  = "IsoformSwitchAnalyzeR_webtools/SignalIP_results.txt"
)
#> Added signal peptide information to 42 (1.96%) transcripts

SwitchList_part2 <- analyzePFAM(
  switchAnalyzeRlist   = SwitchList_part2,
  pathToPFAMresultFile = "IsoformSwitchAnalyzeR_webtools/result_pfam_scan_human.txt",
  showProgress=FALSE)

SwitchList_part2<- analyzeIUPred2A(
  switchAnalyzeRlist        = SwitchList_part2,
  pathToIUPred2AresultFile ="IsoformSwitchAnalyzeR_webtools/IUPred2A.result",
  showProgress = FALSE
)

#------------AS
SwitchList_part2<-readRDS("SwitchList_part2.rds")

SwitchList_part2<- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchList_part2,
  quiet=TRUE
)


pdf("Fig4D.SplicingEnrichment_plot_human.dIF>0.1.pdf", width=7,height=5)
SplicingEnrichment<-extractSplicingEnrichment(
  SwitchList_part2,dIFcutoff = 0.1,onlySigIsoforms = T,
  countGenes = F,
  returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off()


pdf("Fig4G.RAB27A.pdf", width=7,height=5)
switchPlot(SwitchList_part2, gene="RAB27A",plotTopology=FALSE)
dev.off()


SwitchList_part2 <- analyzeSwitchConsequences(
  SwitchList_part2,
  onlySigIsoforms = T,
  dIFcutoff = 0,
  consequencesToAnalyze = c('tss','tts','last_exon','isoform_length','exon_number','intron_structure','ORF_length', 
                            '5_utr_seq_similarity', '5_utr_length', '3_utr_seq_similarity', '3_utr_length',
                            'coding_potential','ORF_seq_similarity','NMD_status',
                            'domains_identified','signal_peptide_identified'))


pdf(file='Fig4F.ConsequenceEnrichment_human.pdf', width=7,height=6)
switch_consequences <- extractConsequenceEnrichment(
  SwitchList_part2,dIFcutoff = 0.1,countGenes = F,minEventsForPlotting = 10,
  returnResult = T # if TRUE returns a data.frame with the summary statistics
)
dev.off()

saveRDS(SwitchList_part2,"SwitchList_part2.rds")