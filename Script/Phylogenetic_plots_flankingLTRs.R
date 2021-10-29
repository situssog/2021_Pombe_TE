#!/usr/bin/env Rscript
rm(list=ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("Biostrings"))

# install.packages("remotes") 
# remotes::install_github("YuLab-SMU/ggtree")

library(BiocManager)
library(Biostrings)

library("ggtree")
library("ggplot2")
library("tidyverse")
library("dplyr")

library(ggstance)
library(ggtree)
library(reshape2)

library(treeio)
library(ape)
library(tidytree)
library(aplot)
library(gdtools)

library(cowplot)


getwd()
# setwd("C:/Users/sertu336/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Phylogenies/LTR_completed_allSeq/IQTreebb")

tree_file<-paste0("tf_1500.treefile")
tree_ed<-treeio::read.iqtree(tree_file)
fig_tree<- ggtree(tree_ed) + 
  theme_tree2() #+
#geom_label2(aes(subset=!isTip, label=node), size=2, color="darkred", alpha=0.5) 
#
fig_tree

tree_ed@data$UFboot[is.na(tree_ed@data$UFboot)]<-0
tree_ed@data$UFboot[tree_ed@data$UFboot<90]<-0
tree_ed@data$UFboot[tree_ed@data$UFboot>=90]<-90

# tree_ed@data$UFboot[is.na(tree_ed@data$UFboot)]<-0
# tree_ed@data$UFboot[tree_ed@data$UFboot<90]<-0
# tree_ed@data$UFboot[tree_ed@data$UFboot>=90 & tree_ed@data$UFboot<95]<-0
# tree_ed@data$UFboot[tree_ed@data$UFboot>=95]<-95

ggtree(tree_ed) +
  theme_tree2() +
  geom_nodepoint(aes(alpha=as.factor(UFboot)), color="darkgreen", size=1) +
  geom_treescale() +
  scale_alpha_manual(values=c(0,0.6))+
  # theme(legend.position="right") +
  theme(legend.position="none") #+
  # ggsave("LTR_completed_tree.png", width = 7, height = 6)


alig=1500
# ltr_alig = readDNAStringSet("alig_all_tf_masked_conSeq_minLen3000_break_plusRef.fasta")
ltr_alig = readDNAStringSet(paste0("alig_all_tf_masked_conSeq_minLen", alig,"_break.fasta"))

head(names(ltr_alig))

## DO NOT RUN the look - File already produced

for (seq_n in seq(1,length(ltr_alig))){
# for (seq_n in seq(1,4)){
  print(seq_n)
  # if (!(names(ltr_alig[seq_n]) %in% c("TF2_I", "TF1_107"))){
    TF2<-(!(as.matrix(ltr_alig[grep("1_Pomberef",names(ltr_alig))][[1]])==as.matrix(ltr_alig[seq_n][[1]])))*1
    TF2[as.matrix(ltr_alig[seq_n][[1]])=="-"]<-2
    #print(sum(TF2))
    TF1<-(!(as.matrix(ltr_alig[grep("_49_JB918_EBC111_III",names(ltr_alig))][[1]])==as.matrix(ltr_alig[seq_n][[1]])))*1
    TF1[as.matrix(ltr_alig[seq_n][[1]])=="-"]<-2
    #print(sum(TF1))
    #sample<-rep(names(ltr_alig[seq_n]), width(ltr_alig[seq_n]))
    # sample<-rep(strsplit((strsplit(str_remove_all(names(ltr_alig[seq_n]), "_R_"), "_I*_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*")[[1]]), "_AB325691_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*")[[1]], width(ltr_alig[seq_n]))
    # sample<-rep(strsplit((strsplit(str_remove_all(names(ltr_alig[seq_n]), "_R_"), "_I*_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*")[[1]]), "_AB325691_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*")[[1]], width(ltr_alig[seq_n]))
    sample<-rep(str_remove_all(str_remove_all(str_remove_all(str_remove_all(names(ltr_alig[seq_n]), "_R_"), "_I*_[0-9]*_[0-9]*_[0-9]*"), "_AB325691_[0-9]*_[0-9]*_[0-9]*"), "_[0-9]*$"), width(ltr_alig[seq_n]))
    #table_seq<-rbind(table_seq, data.frame(sample,TF2, TF1))
    write.table(data.frame(sample,TF2, TF1), paste0("table_seq_all_",alig,".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
  #} 
}

table_seq_dif<-read.table(paste0("table_seq_all_",alig,".txt"), F)
names(table_seq_dif)<-c("sample", "TF2", "TF1")

head(table_seq_dif)

seq_ID_original<-tree_ed@phylo$tip.label
samples_order<-rev(seq_ID_original)
samples_order2<-str_remove_all(samples_order, "_R_")
samples_order3<-str_remove_all(samples_order2, "_I*_[0-9]*_[0-9]*_[0-9]*")
samples_order3<-str_remove_all(samples_order3, "_AB325691_[0-9]*_[0-9]*_[0-9]*")
samples_order4<-str_remove_all(samples_order3, "_[0-9]*$")
tree_ed@phylo$tip.label<-rev(samples_order4)

table_seq_dif <- table_seq_dif %>% 
  full_join(data.frame(sample=tree_ed@phylo$tip.label, 
                       seq_ID_original), by="sample") 


# samples_order2<-sapply(strsplit(samples_order2, "_AB325691_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)
# samples_order3<-sapply(strsplit(samples_order2, "_I*_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)
# tree_ed@phylo$tip.label<-rev(samples_order3)

# table_seq_dif2<-table_seq_dif %>% 
#   mutate(pos=factor(rep(seq(1,width(ltr_alig[1])), length(ltr_alig)), 
#         levels=seq(1,width(ltr_alig[1]))), 
#         sample = factor(sample, levels=samples_order), 
#          TF2=factor(TF2), 
#          TF1=factor(TF1))

head(table_seq_dif)
table_seq_dif2<-table_seq_dif %>% 
  mutate(pos=rep(seq(1,width(ltr_alig[1])), length(ltr_alig)), 
         sample = factor(sample, levels=samples_order4),
         TF2=factor(TF2), 
         TF1=factor(TF1))
head(table_seq_dif2)



# table_seq_dif2 %>% 
#   ggplot(aes(pos, sample)) +
#   geom_tile(aes(fill = TF2)) +
#   scale_fill_manual(values=c("white", "red", "gray90"), expand = c(0, 0)) +
#   theme_classic() + 
#   theme(axis.line = element_line(colour = "black"), 
#         axis.text.x=element_blank(), 
#         axis.text.y = element_blank(), 
#         legend.position="right", 
#         axis.ticks = element_blank()) +
#   ggsave("variants_relativeto_TF2.png", width = 14, height = 10)


LTR_tree<-ggtree(tree_ed) +
  theme_tree2() +
  geom_nodepoint(aes(alpha=as.factor(UFboot)), color="darkgreen", size=1.3) +
  #geom_treescale() +
  scale_alpha_manual(values=c(0,0.9))+
  theme(legend.position="right") 
  

table_seq_dif3<-table_seq_dif2 %>% 
  filter(TF2=="1") %>%
  mutate(col_LTR=ifelse(pos<590,"blue", ifelse(pos>7063, "blue", "red")))



Plot_heatmap_seq<-LTR_tree + 
  geom_facet(panel = "Seq", data = table_seq_dif3, geom = geom_point, 
           mapping=aes(x = pos, color = col_LTR), shape = '.') +
  # theme_tree2(legend.position=c(.05, .85)) +
  theme_tree2(legend.position="none") +
  scale_color_manual(values = c("red", "blue")) 

Plot_heatmap_seq
  
facet_widths(Plot_heatmap_seq, widths = c(1, 5)) #+
  # ggsave("LTR_completed_tree_TF2.png", width = 10, height = 8)




table_seq_dif3<-table_seq_dif2 %>% 
  filter(TF1=="1") %>%
  mutate(col_LTR=ifelse(pos<590,"blue", ifelse(pos>7063, "blue", "red")))


Plot_heatmap_seq<-LTR_tree + 
  geom_facet(panel = "Seq", data = table_seq_dif3, geom = geom_point, 
             mapping=aes(x = pos, color = col_LTR), shape = '.') +
  # theme_tree2(legend.position=c(.05, .85)) +
  theme_tree2(legend.position="none") +
  scale_color_manual(values = c("red", "blue")) 

facet_widths(Plot_heatmap_seq, widths = c(1, 5)) +
  ggsave("LTR_completed_tree_TF1.png", width = 10, height = 8)




### figure by TF group:

## DO NOT RUN the look - File already produced

for (seq_n in seq(1,length(ltr_alig))){
  # for (seq_n in seq(1,4)){
  print(seq_n)
  # if (!(names(ltr_alig[seq_n]) %in% c("TF2_I", "TF1_107"))){
  # TF2<-(!(as.matrix(ltr_alig[269][[1]])==as.matrix(ltr_alig[seq_n][[1]])))*1
  TF2<-(!(as.matrix(ltr_alig[grep("1_Pomberef",names(ltr_alig))][[1]])==as.matrix(ltr_alig[seq_n][[1]])))*1
  #print(sum(TF2))
  # TF1<-(!(as.matrix(ltr_alig[2][[1]])==as.matrix(ltr_alig[seq_n][[1]])))*1
  TF1<-(!(as.matrix(ltr_alig[grep("49_JB918_EBC111_III",names(ltr_alig))][[1]])==as.matrix(ltr_alig[seq_n][[1]])))*1
  #print(sum(TF1))
  gaps<-(as.matrix(ltr_alig[seq_n][[1]])=="-")*1
  # sample<-rep(names(ltr_alig[seq_n]), width(ltr_alig[seq_n]))
  # sample<-rep(strsplit((strsplit(str_remove_all(names(ltr_alig[seq_n]), "_R_"), "_I*_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*")[[1]]), "_AB325691_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*")[[1]], width(ltr_alig[seq_n]))
  sample<-rep(str_remove_all(str_remove_all(str_remove_all(str_remove_all(names(ltr_alig[seq_n]), "_R_"), "_I*_[0-9]*_[0-9]*_[0-9]*"), "_AB325691_[0-9]*_[0-9]*_[0-9]*"), "_[0-9]*$"), width(ltr_alig[seq_n]))
  #table_seq<-rbind(table_seq, data.frame(sample,TF2, TF1))
  write.table(data.frame(sample,TF2, TF1, gaps), paste0("table_seq_bygroup",alig,".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
  #} 
}

table_seq_dif_bygroup<-read.table(paste0("table_seq_bygroup",alig,".txt"), F)
names(table_seq_dif_bygroup)<-c("sample", "TF2", "TF1", "gaps")

head(table_seq_dif_bygroup)

table_seq_dif_bygroup <- table_seq_dif_bygroup %>% 
  full_join(data.frame(sample=tree_ed@phylo$tip.label, 
                       seq_ID_original), by="sample") 


table_seq_dif2_bygroup<-table_seq_dif_bygroup %>% 
  mutate(pos=rep(seq(1,width(ltr_alig[1])), length(ltr_alig)), 
         sample = factor(sample, levels = samples_order4), 
         bygroup= factor(ifelse((TF2==1 & TF1==1), "Both", 
                          ifelse((TF2==1 & TF1==0), "TF1", 
                                  ifelse((TF2==0 & TF1==1), "TF2", "NVar"))), 
                         levels=c("NVar", "TF2", "TF1", "Both")),
         TF2=factor(TF2), 
         TF1=factor(TF1)) %>% 
  filter(bygroup!="NVar") %>%
  filter(gaps!=1)


Plot_heatmap_seq<-LTR_tree +
  geom_facet(panel = "Seq", data = table_seq_dif2_bygroup, geom = geom_point, 
             mapping=aes(x = pos, color = bygroup), shape = '.') +
  # theme_tree2(legend.position=c(.05, .85)) +
  theme_tree2(legend.position="left") +
  scale_color_manual(values = c("red", "blue", "black")) +
  scale_x_continuous(breaks = c(0,590, 2000,4000,6000, 7063))

facet_widths(Plot_heatmap_seq, widths = c(1, 5)) #+
  # ggsave("LTR_completed_tree_diff_bygroup.png", width = 10, height = 8)
Plot_heatmap_seq

head(table_seq_dif_bygroup)
head(table_seq_dif2_bygroup)
head(table_seq_dif2)

#### adding distances between solo LTR pairs:
# the input for this section was produced in the solo_LTR analysis:

# getwd()
# setwd("C:/Users/sertu336/LRZ Sync+Share/TEs/all_Samples/Phylogenies/LTR_completed/IQTree_bb/")

table_distance_soloLTR<-read.table("../../Solo_LTR_completed_allSeq/IQTreebb/table_distance_soloLTR_pairs.txt", T)

names(table_distance_soloLTR)<-c("sample", "p1", "p2", "distance", "TF1_TF1","TF2_TF2","group")

tem_table_distance_soloLTR_sample<-str_remove_all(table_distance_soloLTR$sample, "_R_")
tem_table_distance_soloLTR_sample<-str_remove_all(tem_table_distance_soloLTR_sample, "_I*_[0-9]*_[0-9]*_[0-9]*")
tem_table_distance_soloLTR_sample<-str_remove_all(tem_table_distance_soloLTR_sample, "_AB325691_[0-9]*_[0-9]*_[0-9]*")
tem_table_distance_soloLTR_sample<-str_remove_all(tem_table_distance_soloLTR_sample, "_[0-9]*$")
table_distance_soloLTR$sample<-tem_table_distance_soloLTR_sample

# table_distance_soloLTR2<-table_distance_soloLTR %>% 
#   filter(distance<0.2)


LTR_tree2<-LTR_tree + 
  geom_facet(panel = "K2P", data = table_distance_soloLTR, geom = ggstance::geom_barh, 
             mapping=aes(x = (distance+0.0001), color = group, fill = group), stat = "identity") +
  scale_color_manual(values = c("blue", "orange", "red")) +
  scale_fill_manual(values = c("blue", "orange", "red")) 
  
facet_widths(LTR_tree2, widths = c(3, 1)) #+
  # ggsave("LTR_completed_tree_divergence_soloLTR_bygroup.png", width = 6, height = 6)


  
### both plots:

LTR_tree 

head(table_distance_soloLTR2)
head(table_seq_dif2_bygroup)

library(dplyr)
library(tidytree) 

d <- filter(LTR_tree, isTip) %>% select(c(label, y))


# dd1 <- left_join(table_distance_soloLTR2 %>% mutate(label=sample), d, by='label') 
dd11 <- left_join(table_distance_soloLTR %>% mutate(label=sample), d, by='label') 
dd2 <- left_join(table_seq_dif2_bygroup %>% mutate(label=sample), d, by='label')

# p1 <- ggplot(dd1, aes(y, distance+0.001)) + 
#   geom_col(aes(fill=group, colour=group)) + 
#   coord_flip() + 
#   scale_color_manual(values = c("blue", "orange", "red")) +
#   scale_fill_manual(values = c("blue", "orange", "red")) +
#   theme_tree2() +
#   theme(legend.position='none') 

p11 <- ggplot(dd11, aes(y, distance+0.001)) + 
  geom_col(aes(fill=group, colour=group)) + 
  coord_flip() + 
  scale_color_manual(values = c("blue", "orange", "red")) +
  scale_fill_manual(values = c("blue", "orange", "red")) +
  theme_tree2() +
  theme(legend.position='none') 
  

# Plot_heatmap_seq<-LTR_tree + 
#   geom_facet(panel = "Seq", data = table_seq_dif2_bygroup, geom = geom_point, 
#              mapping=aes(x = pos, color = bygroup), shape = '.') +
#   # theme_tree2(legend.position=c(.05, .85)) +
#   theme_tree2(legend.position="left") +
#   scale_color_manual(values = c("red", "blue", "black")) +
#   scale_x_continuous(breaks = c(0,485, 2000,4000,6000, 6524))


p2 <- ggplot(dd2, aes(x=pos, y=y)) + 
  geom_tile(aes(fill=bygroup, colour=bygroup)) + 
  #scale_fill_viridis_c() + 
  scale_fill_manual(values = c("red", "blue", "black")) +
  scale_colour_manual(values = c("red", "blue", "black")) +
  scale_x_continuous(breaks = c(0,590, 2000,4000,6000, 7063)) +
  theme_tree2() + 
  theme(legend.position='none')


library(ggtree)
g<- LTR_tree + theme(legend.position='none')
# p1 <- p1 + ylim2(LTR_tree) 
p11 <- p11 + ylim2(LTR_tree) 
p2 <- p2 + ylim2(LTR_tree)

library(cowplot)
# plot_grid(g, p2, p1, ncol=3, align='h', 
#           rel_widths = c(0.15, 0.5, 0.20)) +
#   ggsave("tree_divergence_seq.png", width = 8, height = 7)

plot_grid(g, p2, p11, ncol=3, align='h', 
          rel_widths = c(0.15, 0.5, 0.1)) #+
  # ggsave("tree_divergence_seq_all.png", width = 8, height = 7)

p2


### plot extracting admixed samples:

#samples_order<-rev(tree_ed@phylo$tip.label)

tree_ed2<-read.tree(tree_file)
plot(tree_ed2,show.tip.label = FALSE,  show.node.label = FALSE, align.tip.label=TRUE)
ape::nodelabels(frame = "none", bg = "none")

cladeTF2 <- tree_subset(tree_ed, node=1666, levels_back=0)
TF2_nodes<-as_tibble(cladeTF2) %>% 
  filter(group==1) %>% 
  filter(label %in% samples_order4) %>% 
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

tree_ed2 <- groupOTU(tree_ed, TF2_nodes)
ggtree(tree_ed2) +
  # geom_tiplab() +
  geom_tippoint(aes(colour=group)) #+
# ggsave("test_tree.svg", width = 49, height = 49)


# TF2_nodes<-samples_order3[1:which(samples_order3=="18_JB758")]

# tree_ed2 <- groupOTU(tree_ed, TF2_nodes)
# ggtree(tree_ed2) +
#   # geom_tiplab() +
#   geom_tippoint(aes(colour=group)) #+
#   # ggsave("test_tree.svg", width = 49, height = 49)

head(table_distance_soloLTR)

admixed_samples<-table_distance_soloLTR %>%
  filter((group=="TF1_TF2") |
  ((group=="TF1_TF1") & (sample %in% TF2_nodes)) |
  ((group=="TF2_TF2") & (!(sample %in% TF2_nodes)))) %>% 
  ungroup() %>% 
  select(sample) %>% 
  unlist() %>% 
  as.vector()



tree_ed2 <- groupOTU(tree_ed, admixed_samples)
ggtree(tree_ed2) +
  # geom_tiplab() +
  geom_tippoint(aes(colour=group, alpha=group)) +
  scale_colour_manual(values=c("white", "orange")) +
  scale_alpha_manual(values=c(0, 0.8)) +
  geom_treescale() +
  theme(legend.position='none') #+
  # ggsave("admixed_sequences.png", width = 7, height = 7)


table_seq_dif3_bygroup<-table_seq_dif_bygroup %>% 
  mutate(pos=rep(seq(1,width(ltr_alig[1])), length(ltr_alig)), 
         sample = factor(sample, levels=c(samples_order4)), 
         bygroup= factor(ifelse((TF2==1 & TF1==1), "Both", 
                                ifelse((TF2==1 & TF1==0), "TF1", 
                                       ifelse((TF2==0 & TF1==1), "TF2", "NVar"))), 
                         levels=c("NVar", "TF2", "TF1", "Both", "Gap")),
         TF2=factor(TF2), 
         TF1=factor(TF1)) %>% 
  filter(bygroup!="NVar")

table_seq_dif3_bygroup$bygroup[table_seq_dif3_bygroup$gaps==1]<-"Gap"

p1<-table_seq_dif3_bygroup %>% 
  # mutate(sample=factor(sample,levels=c(samples_order))) %>% 
  filter(sample %in% admixed_samples) %>% 
  ggplot(aes(pos,sample)) +
  geom_tile(aes(fill=bygroup, colour=bygroup)) +
  scale_colour_manual(values=c("red", "blue", "black", "gray80")) +
  scale_fill_manual(values=c("red", "blue", "black", "gray80")) +
  # scale_x_continuous(breaks = c(0,485, 2000,4000,6000, 6524)) +
  scale_x_continuous(breaks = c(0, 590, 7063)) +
  theme_classic() +
  # theme(legend.position = "none") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=6, colour="black")) 



p2<-table_distance_soloLTR %>% 
  mutate(sample=factor(sample,levels=c(samples_order4))) %>% 
  filter(sample %in% admixed_samples) %>% 
  ggplot(aes(sample,distance, colour=group)) +
  geom_point() +
  scale_colour_manual(values=c("blue", "orange", "red")) +
  xlab("") +
  ylab("K2P distance") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.y = element_blank()) +
  # theme(legend.position = "none") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size=12, colour="black"))

library(cowplot) 
plot_grid(p1, p2, ncol=2, align='h', 
          rel_widths = c(0.6, 0.2)) #+
  # ggsave("blocks_admixSeq.png", width = 14, height = 10)



table_seq_dif3_bygroup %>% 
  # mutate(sample=factor(sample,levels=c(samples_order))) %>% 
  filter(sample %in% admixed_samples) %>% 
  ungroup() %>% 
  select(pos) %>% 
  arrange(pos) %>% 
  unique() %>% 
  filter(pos>590) %>% #590
  filter(pos>7063) %>% head()#7063


  

pos_lab<-table_seq_dif3_bygroup %>% 
  # mutate(sample=factor(sample,levels=c(samples_order))) %>% 
  filter(sample %in% admixed_samples) %>% 
  ungroup() %>% 
  select(pos) %>% 
  arrange(pos) %>% 
  unique() %>% 
  mutate(pos_l=ifelse(pos %in% c(590, 7063), "LTR", "")) 


head(pos_lab)
p1<-table_seq_dif3_bygroup %>% 
  # mutate(sample=factor(sample,levels=c(samples_order))) %>% 
  filter(sample %in% admixed_samples) %>% 
  mutate(pos_label=ifelse(pos %in% c(590, 7063),pos,"")) %>% 
  ggplot(aes(factor(pos),sample)) +
  geom_tile(aes(fill=bygroup, colour=bygroup)) +
  scale_colour_manual(values=c("red", "blue", "black", "gray80")) +
  scale_fill_manual(values=c("red", "blue", "black", "gray80")) +
  #scale_x_continuous(breaks = c(0,485, 2000,4000,6000, 6524)) +
  #scale_x_discrete(breaks=pos, labels=as.character(pos_lab$pos_l)) %>% 
  scale_x_discrete(labels=c("5' LTR","3' LTR"), breaks = c(590, 7063)) +
  xlab("") +
  ylab("Sequence") +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=6, colour="black")) 

plot_grid(p1, p2, ncol=2, align='h', 
          rel_widths = c(0.6, 0.2)) #+
  # ggsave("blocks_admixSeq_noFixpos.png", width = 10, height = 10)


data.frame(admixed_samples) %>% 
  mutate(sample=str_replace_all(admixed_samples, "[0-9]*_J", "J")) %>% 
  mutate(sample=str_replace_all(sample, "[0-9]*_P", "P")) %>% 
  # dim()
  ungroup() %>% 
  select(sample) %>% 
  unlist() %>% 
  as.vector() %>% 
  unique()



# there are 137 sequences with evidence of admixture
# there are admixed sequences in 35 strains

data.frame(admixed_samples) %>% 
  mutate(sample=str_replace_all(admixed_samples, "[0-9]*_J", "J")) %>% 
  mutate(sample=str_replace_all(sample, "[0-9]*_P", "P")) %>% 
  mutate(sample=str_replace_all(sample, "_[0-9]$", "")) %>% 
  group_by(sample) %>% 
  summarise(N_seq=n())

# ancestral block distribution:
# ancestral_data<-read.table("C:/Users/sertu336/LRZ Sync+Share/TEs/all_Samples/Ancestral_blocks/ancestralhap_57ILL_LR.txt",T)
# ancestral_data<-read.table("C:/Users/Sergio/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Ancestral_blocks/ancestralhap_57ILL_LR.txt",T)

# ancestral_data<-read.table("C:/Users/sertu336/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Ancestral_blocks/ancestralhap_57ILL_LR.txt",T)
                         
ancestral_data<-read.table("C:/Users/ru43sej/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Ancestral_blocks/ancestralhap_57ILL_LR.txt",T)

ILL_samples<-as.vector(unique(ancestral_data$sample)[grep("ILL",unique(ancestral_data$sample))])

# proportion of ancestral populations per sample

table_anc_pro<-ancestral_data %>% 
  filter(!(sample %in% ILL_samples)) %>% 
  group_by(sample, Nor_PC1_pol_sim) %>% 
  summarise(number_bins = n()) %>% 
  filter(Nor_PC1_pol_sim!=0.5) %>% 
  spread(Nor_PC1_pol_sim, number_bins, fill=1) %>%
  # group_by(sample) %>% 
  mutate(total=sum(`0`+`1`), 
         proportion_sk=`1`/total) %>% 
  arrange(proportion_sk,sample) 
  
table_anc_pro$sample<-factor(as.vector(table_anc_pro$sample), levels=c(as.vector(table_anc_pro$sample)))
head(table_anc_pro)

ggplot(table_anc_pro, aes(sample, proportion_sk)) +
  geom_bar(stat = "identity", fill="steelblue") +
  theme_classic() +
  # scale_y_continuous(expand = c(0,0))+
  ylim(0,1)+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y = element_text(size=14, colour="black"))  #+
  # ggsave("../../../Ancestral_blocks/anc_proportions_LReads.png", width = 7, height = 4)


order_samplebyAnProp<-as.vector(table_anc_pro$sample)
order_samplebyAnProp<-c("Pomberef", order_samplebyAnProp)

table_n_mix_seq<-data.frame(admixed_samples) %>% 
  mutate(sample=str_replace_all(admixed_samples, "[0-9]*_J", "J")) %>% 
  mutate(sample=str_replace_all(sample, "[0-9]*_P", "P")) %>% 
  mutate(sample=str_replace_all(sample, "_[0-9]$", "")) %>% 
  group_by(sample) %>% 
  summarise(N_seq=n()) %>% 
  mutate(sample=factor(sample, levels=c(order_samplebyAnProp))) 

table_n_mix_seq<-rbind(table_n_mix_seq, 
      data.frame(sample=order_samplebyAnProp[!(order_samplebyAnProp %in% 
                                           table_n_mix_seq$sample)], 
           N_seq=0))

table_n_mix_seq %>% 
  mutate(sample=factor(sample, levels=c(order_samplebyAnProp))) %>% 
  ggplot(aes(sample, N_seq)) +
  geom_bar(stat="identity", fill="black") +
  ylab("Num. Admixed Seq") +
  theme_classic() +
  # scale_y_continuous(expand = c(0,0))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("admixedSeq_perSample.png", width = 7, height = 4)


### divergences between completed TE elements (CORE region):
sequences_coreLTR<-ape::read.dna("alig_all_CORE_tf_masked_conSeq_minLen1500_break.fasta", format = "fasta")

distance_table_coreTF<-ape::dist.dna(sequences_coreLTR, model = "K80", pairwise.deletion = TRUE, as.matrix=TRUE) 

# colnames(distance_table_coreTF)
# rownames(distance_table_coreTF)

distance_table_coreTF2<-str_remove_all(colnames(distance_table_coreTF),
                                       "_R_")
distance_table_coreTF3<-str_remove_all(distance_table_coreTF2, 
                                       "_I*_[0-9]*_[0-9]*_[0-9]*")
distance_table_coreTF3<-str_remove_all(distance_table_coreTF3, 
                               "_AB325691_[0-9]*_[0-9]*_[0-9]*")
distance_table_coreTF4<-str_remove_all(distance_table_coreTF3, "_[0-9]*$")
colnames(distance_table_coreTF)<-distance_table_coreTF4



distance_table_coreTF2<-str_remove_all(rownames(distance_table_coreTF),
                                       "_R_")
distance_table_coreTF3<-str_remove_all(distance_table_coreTF2, 
                                       "_I*_[0-9]*_[0-9]*_[0-9]*")
distance_table_coreTF3<-str_remove_all(distance_table_coreTF3, 
                                       "_AB325691_[0-9]*_[0-9]*_[0-9]*")
distance_table_coreTF4<-str_remove_all(distance_table_coreTF3, "_[0-9]*$")
rownames(distance_table_coreTF)<-distance_table_coreTF4



head(colnames(distance_table_coreTF))
# temporal<-str_remove_all(colnames(distance_table_coreTF), "_R_")
# temporal<-sapply(strsplit(temporal, "_I*_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)
# colnames(distance_table_coreTF)<-sapply(strsplit(temporal, "_AB325691_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)
# 
# temporal<-str_remove_all(rownames(distance_table_coreTF), "_R_")
# temporal<-sapply(strsplit(temporal, "_I*_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)
# rownames(distance_table_coreTF)<-sapply(strsplit(temporal, "_AB325691_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)

TF1_nodes<-samples_order4[!(samples_order4 %in% TF2_nodes)]

library(reshape2)
distance_table_coreTF2<-distance_table_coreTF[TF2_nodes, TF2_nodes]
distance_table_coreTF2<-setNames(melt(distance_table_coreTF2), c('sample1', 'sample2', 'distance'))

distance_table_coreTF2<-distance_table_coreTF2 %>% 
  filter(sample1 != sample2) %>% 
  mutate(group=factor(ifelse(((sample1 %in% admixed_samples) | (sample2 %in% admixed_samples)), "admixed_TF2", "TF2"), levels=c("TF1", "TF2", "admixed_TF1", "admixed_TF2"))) 

distance_table_coreTF1<-distance_table_coreTF[TF1_nodes, TF1_nodes]
distance_table_coreTF1<-setNames(melt(distance_table_coreTF1), c('sample1', 'sample2', 'distance'))

distance_table_coreTF1<-distance_table_coreTF1 %>% 
  filter(sample1 != sample2) %>% 
  mutate(group=factor(ifelse(((sample1 %in% admixed_samples) | (sample2 %in% admixed_samples)), "admixed_TF1", "TF1"), levels=c("TF1", "TF2", "admixed_TF1", "admixed_TF2"))) 

distance_table_coreTF1 %>% 
  filter(sample1== "72_JB1110_EBC121_1") %>%
  filter(sample2== "98_JB874_1") 


rbind(distance_table_coreTF1, distance_table_coreTF2) %>% 
  #filter(is.na(distance))
  ggplot(aes(distance*100, fill=group)) +
  # geom_density() +
  # geom_histogram(position="dodge") +
  geom_histogram() +
  xlab("K2P distance (%)") +
  #xlim(c(-1,80))+
  scale_fill_manual(values = c("blue", "red", "orange", "orange")) +
  #scale_x_continuous(breaks= seq(0,60,5))+
  facet_grid(group ~ ., scale="free") +
  theme_classic() #+
  # ggsave("divergence_coreTF.png", width = 4, height = 5)



### divergences between completed TE elements (including flanking LTR):
sequences_LTR<-ape::read.dna("alig_all_tf_masked_conSeq_minLen1500_break.fasta", format = "fasta")

distance_table_TF<-ape::dist.dna(sequences_LTR, model = "K80", pairwise.deletion = TRUE, as.matrix=TRUE) 

# colnames(distance_table_coreTF)
# rownames(distance_table_coreTF)

temporal<-str_remove_all(colnames(distance_table_TF), "_R_")
temporal<-str_remove_all(temporal, "_I*_[0-9]*_[0-9]*_[0-9]*")
temporal<-str_remove_all(temporal, "_AB325691_[0-9]*_[0-9]*_[0-9]*")
temporal<-str_remove_all(temporal, "_[0-9]*$")
colnames(distance_table_TF)<-temporal

temporal<-str_remove_all(rownames(distance_table_TF), "_R_")
temporal<-str_remove_all(temporal, "_I*_[0-9]*_[0-9]*_[0-9]*")
temporal<-str_remove_all(temporal, "_AB325691_[0-9]*_[0-9]*_[0-9]*")
temporal<-str_remove_all(temporal, "_[0-9]*$")
rownames(distance_table_TF)<-temporal

# temporal<-str_remove_all(colnames(distance_table_TF), "_R_")
# temporal<-sapply(strsplit(temporal, "_I*_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)
# colnames(distance_table_TF)<-sapply(strsplit(temporal, "_AB325691_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)
# 
# temporal<-str_remove_all(rownames(distance_table_TF), "_R_")
# temporal<-sapply(strsplit(temporal, "_I*_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)
# rownames(distance_table_TF)<-sapply(strsplit(temporal, "_AB325691_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)

TF1_nodes<-samples_order4[!(samples_order4 %in% TF2_nodes)]
TF2_nodes
library(reshape2)

distance_table_TF2<-distance_table_TF[TF2_nodes, TF2_nodes]
distance_table_TF2<-setNames(melt(distance_table_TF2), c('sample1', 'sample2', 'distance'))

distance_table_TF2<-distance_table_TF2 %>% 
  filter(sample1 != sample2) %>% 
  mutate(group=ifelse(((sample1 %in% admixed_samples) | (sample2 %in% admixed_samples)), "admixed_TF2", "TF2")) 

distance_table_TF1<-distance_table_TF[TF1_nodes, TF1_nodes]
distance_table_TF1<-setNames(melt(distance_table_TF1), c('sample1', 'sample2', 'distance'))

distance_table_TF1<-distance_table_TF1 %>% 
  filter(sample1 != sample2) %>% 
  mutate(group=ifelse(((sample1 %in% admixed_samples) | (sample2 %in% admixed_samples)), "admixed_TF1", "TF1")) 

rbind(distance_table_TF1, distance_table_TF2) %>% 
  mutate(group=factor(group, levels=c("TF1", "TF2", "admixed_TF1", "admixed_TF2"))) %>% 
  ggplot(aes(distance*100, fill=group)) +
  # geom_density() +
  # geom_histogram(position="dodge") +
  geom_histogram() +
  xlab("K2P distance (%)") +
  scale_fill_manual(values = c("blue", "red", "orange", "orange")) +
  #scale_x_continuous(breaks= seq(0,60,5))+
  facet_grid(group ~ ., scale="free") +
  theme_classic() #+
  # ggsave("divergence_wholeTF.png", width = 4, height = 5)

######

sequences<-ape::read.dna("test.fasta", format = "fasta")
sequences
distance_table<-ape::dist.dna(sequences, model = "K80", pairwise.deletion = TRUE, as.matrix=TRUE) 
distance_table



### extracting sequences per TF group:

sequences_TF<-ape::read.dna("alig_all_tf_masked_conSeq_minLen1500_break.fasta", format = "fasta")

tree_file<-paste0("tf_1500.treefile")
tree_ed_ori<-treeio::read.iqtree(tree_file)
fig_tree_ori<- ggtree(tree_ed_ori) + 
  theme_tree2() #+

samples_order_ori<-rev(tree_ed_ori@phylo$tip.label)

cladeTF2_ori <- tree_subset(tree_ed_ori, node=1666, levels_back=0)
TF2_nodes_ori<-as_tibble(cladeTF2_ori) %>% 
  filter(group==1) %>% 
  filter(label %in% samples_order_ori) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()


samples_ID_ori<-rev(tree_ed_ori@phylo$tip.label)
samples_ID<-str_remove_all(samples_ID_ori, "_R_")
samples_ID<-sapply(strsplit(samples_ID, "_AB325691_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)
samples_ID<-sapply(strsplit(samples_ID, "_I*_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)

table_names<-data.frame(samples_ID_ori, samples_ID)

non_mixed_TF2<-table_names %>% 
  filter(samples_ID_ori %in% TF2_nodes_ori) %>% 
  filter(!(samples_ID %in% admixed_samples)) %>% 
  ungroup() %>% 
  select(samples_ID_ori) %>% 
  unlist() %>% 
  as.vector()

ltr_alig = readDNAStringSet(paste0("alig_all_tf_masked_conSeq_minLen", alig,"_break.fasta"))
ltr_alig[names(ltr_alig) %in% non_mixed_TF2]

writeXStringSet(ltr_alig[names(ltr_alig) %in% non_mixed_TF2], "TF2_non_mixed.fasta", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")


non_mixed_TF1<-table_names %>% 
  filter(!(samples_ID_ori %in% TF2_nodes_ori)) %>% 
  filter(!(samples_ID %in% admixed_samples)) %>% 
  ungroup() %>% 
  select(samples_ID_ori) %>% 
  unlist() %>% 
  as.vector()

ltr_alig = readDNAStringSet(paste0("alig_all_tf_masked_conSeq_minLen", alig,"_break.fasta"))
ltr_alig[names(ltr_alig) %in% non_mixed_TF1]

writeXStringSet(ltr_alig[names(ltr_alig) %in% non_mixed_TF1], "TF1_non_mixed.fasta", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")

## Those fasta files were aligned with MAFFT to produce the files: alig_TF1_non_mixed.fasta and alig_TF2_non_mixed.fasta




### Annotate large sequences with clade information: 

# tree_ed@phylo$tip.label the same as rev(samples_order4)
# TF2_nodes
# TF1_nodes

head(table_seq_dif2)
head(table_seq_dif2_bygroup)

#
ancestralhapl<-read.table("../../../Annotation/ancestralhap_57ILL_LR.txt", T) 
ancestralhapl<-ancestralhapl[-(grep("ILL", ancestralhapl$sample)),]

pro_AncPop <- ancestralhapl %>% 
  filter(Nor_PC1_pol_sim!=0.5) %>% 
  group_by(sample) %>% 
  summarise(proportion_sk=sum(Nor_PC1_pol_sim)/n()) 

order_samples_pro_AncPop <-pro_AncPop %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB879", 0.0007,proportion_sk)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB873", 0.4533954726,prop_sk_ed)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB1206", 0.8513139695,prop_sk_ed)) %>%   arrange(prop_sk_ed) %>% 
  ungroup() %>% 
  select(sample) %>% 
  unlist() %>% 
  as.vector()

pro_AncPop<-pro_AncPop %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB879", 0.0007,proportion_sk)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB873", 0.4533954726,prop_sk_ed)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB1206", 0.8513139695,prop_sk_ed)) %>%   arrange(prop_sk_ed) %>% 
  mutate(anc_prop=1-prop_sk_ed) %>% 
  select(sample, anc_prop)


order_samples_pro_AncPop <-c("LTR_loc", "ref_loc", "Pomberef", order_samples_pro_AncPop)


LTR_clade_table<-read.table("../../Solo_AllLTR_allSeq/IQTreebb/LTR_seq_table.txt", T)


seq_changeSerie<-c("110_JB840_2", "135_JB22_2", "16_JB879_2", "188_JB1110_EBC121_2", "223_JB873_EBC095_2", "226_JB929_2", "226_JB929_3", "40_JB934_EBC115_2", "52_Pomberef_2", "79_JB874_3")

table_seq_dif4_bygroup<-table_seq_dif2_bygroup %>% 
  extract(sample, "sim_ID_seq", "(.*)_.", remove = F) %>% 
  extract(sample, "serie_ID_seq", ".*_(.)", remove = F) %>% 
  mutate(reverse_seq=(substr(seq_ID_original, 1, 3)=="_R_")*1) %>% 
  mutate(serie_ID_seq=if_else(sample %in% seq_changeSerie, as.numeric(serie_ID_seq)-1, as.numeric(serie_ID_seq))) %>% 
  select(-pos, -gaps, -TF2, -TF1, -bygroup) %>% 
  unique() 


table_seq_dif2_bygroup %>% 
  extract(sample, "sim_ID_seq", "(.*)_.", remove = F) %>% 
  head()

# table_seq_dif4_bygroup %>% 
#   tbl_df %>%
#   print(n = Inf)

head(LTR_clade_table)

d5_LTR<-c()
d3_LTR<-c()
for (line in seq(1,dim(table_seq_dif4_bygroup)[1])){
  # line=13
  seq_ID_used=as.vector(table_seq_dif4_bygroup$sim_ID_seq[line])
  series_used<-as.numeric(table_seq_dif4_bygroup$serie_ID_seq[line])
  if(table_seq_dif4_bygroup$reverse_seq[line]==0){
    d5_LTR_vector<-LTR_clade_table %>% 
      filter(sim_ID_seq==seq_ID_used) %>% 
      filter(seq_series==series_used) %>%
      select(cladeID) %>% 
      unlist() %>% 
      as.vector()
    if(length(d5_LTR_vector)==1){
      d5_LTR[line]<-d5_LTR_vector
    } else if (length(d5_LTR_vector)==0){
      d5_LTR[line]<-"NoLTR"
    } else {
      d5_LTR[line]<-"Multiple_Hits"
    }
    d3_LTR_vector<-LTR_clade_table %>% 
      filter(sim_ID_seq==seq_ID_used) %>% 
      filter(seq_series==series_used+1) %>%
      select(cladeID) %>% 
      unlist() %>% 
      as.vector()
    if(length(d3_LTR_vector)==1){
      d3_LTR[line]<-d3_LTR_vector
    } else if (length(d3_LTR_vector)==0){
      d3_LTR[line]<-"NoLTR"
    } else {
      d3_LTR[line]<-"Multiple_Hits"
    }
  } else {
    d5_LTR_vector<-LTR_clade_table %>% 
      filter(sim_ID_seq==seq_ID_used) %>% 
      filter(seq_series_ed==series_used) %>%
      select(cladeID) %>% 
      unlist() %>% 
      as.vector()
    if(length(d5_LTR_vector)==1){
      d5_LTR[line]<-d5_LTR_vector
    } else if (length(d5_LTR_vector)==0){
      d5_LTR[line]<-"NoLTR"
    } else {
      d5_LTR[line]<-"Multiple_Hits"
    }
    d3_LTR_vector<-LTR_clade_table %>% 
      filter(sim_ID_seq==seq_ID_used) %>% 
      filter(seq_series_ed==series_used+1) %>%
      select(cladeID) %>% 
      unlist() %>% 
      as.vector()
    if(length(d3_LTR_vector)==1){
      d3_LTR[line]<-d3_LTR_vector
    } else if (length(d3_LTR_vector)==0){
      d3_LTR[line]<-"NoLTR"
    } else {
      d3_LTR[line]<-"Multiple_Hits"
    }
  }
}


table_seq_dif4_bygroup %>% head()
LTR_clade_table %>% head()

table_seq_dif4_bygroup_ed<-table_seq_dif4_bygroup %>% 
  cbind(d5_LTR, d3_LTR) %>% 
  mutate(sample_serie=sample) %>%
  select(-sample) %>% 
  left_join(LTR_clade_table %>% 
              select(sim_ID_seq, sample, chr, 
                     pos, dis, len_seq, cluster, ancHaplo_mean) %>% 
              unique(), 
            by="sim_ID_seq") %>%
  mutate(core=if_else(sample_serie %in% TF2_nodes, "beta", "alpha")) %>% 
  mutate(sample=factor(sample, levels=order_samples_pro_AncPop))
  # ggplot(aes(cluster)) +
  # geom_histogram(binwidth = 1)
  # tbl_df %>%
  # print(n = Inf)



  
table_seq_dif4_bygroup_ed %>% 
  mutate(pos_group=if_else(pos>4410000, 2,1)) %>% 
  # filter(cluster==514) %>% 
  filter(cluster==608) %>% 
  filter(pos_group==2) %>% 
  mutate(d5_ancHaplo_mean=if_else(ancHaplo_mean==1, "Sk", "Sp"), 
         d3_ancHaplo_mean=if_else(ancHaplo_mean==1, "Sk", "Sp")) %>% 
  # select(sample_serie) %>% 
  # dim()
  select(sample, d5_ancHaplo_mean, d5_LTR, core, d3_LTR, d3_ancHaplo_mean) %>% 
  gather(key="region", value = "haplotype", -sample)  %>% 
  mutate(region=factor(region, levels=c("d5_ancHaplo_mean", "d5_LTR", "core", "d3_LTR", "d3_ancHaplo_mean")), 
         haplotype=factor(haplotype, levels=c("Sp", "Sk", "alpha", "subalpha", "beta"))) %>% 
  full_join(table_seq_dif4_bygroup_ed %>% 
              mutate(pos_group=if_else(pos>4410000, 2,1)) %>% 
              filter(cluster==514) %>% 
              filter(pos_group==2) %>% 
              select(sample, ancHaplo_mean) %>% 
              unique(), 
            by="sample") %>% 
  mutate(ancHaplo_mean_ed=if_else(ancHaplo_mean==0, "Sp", "Sk")) %>%
  ggplot(aes(region,sample))+
  geom_tile(aes(fill=haplotype)) +
  labs(x="Region", y="Sample") + 
  scale_fill_manual(values=c("darkred", "steelblue", "blue","orange", "red"))+
  facet_grid(ancHaplo_mean_ed ~ ., scale="free", space="free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=14, colour="black"), 
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("cluster514_annotation.png", width = 7, height = 7)

  
table_seq_dif4_bygroup_ed %>% 
  filter(d5_LTR!="NoLTR" & d3_LTR!="NoLTR" ) %>% 
  filter(d5_LTR!="gamma_iota" & d3_LTR!="gamma_iota" ) %>% 
  filter(d5_LTR!="other" & d3_LTR!="other" ) %>% 
  unite(genotype, c("d5_LTR", "core", "d3_LTR")) %>% 
  group_by(genotype) %>% 
  summarise(Count=n()) %>% 
  ggplot(aes(x = reorder(genotype, -Count), y = Count)) +
  geom_bar(stat = "identity")+
  geom_text(aes(x = reorder(genotype, -Count), 
                y = Count+20, label=Count), size=6)+
  labs(x="Genotype", y="Count") + 
  theme_classic()+
  theme(strip.text.y = element_text(angle = 0))+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("Completed_LTR_elements_genotypes_counts.png", width = 7, height = 5)


table_seq_dif4_bygroup_ed %>% 
  filter(d5_LTR!="NoLTR" & d3_LTR!="NoLTR" ) %>% 
  filter(d5_LTR!="gamma_iota" & d3_LTR!="gamma_iota" ) %>% 
  filter(d5_LTR!="other" & d3_LTR!="other" ) %>% 
  unite(genotype, c("d5_LTR", "core", "d3_LTR")) %>% 
  group_by(genotype) %>% 
  summarise(N_clusters=length(unique(cluster))) %>%
  # filter(N_clusters>1) %>%
  ggplot(aes(x = reorder(genotype, -N_clusters), y = N_clusters)) +
  geom_bar(stat = "identity")+
  geom_text(aes(x = reorder(genotype, -N_clusters), 
                y = N_clusters+20, label=N_clusters), size=6)+
  labs(x="Genotype", y="Count Num. Clusters") + 
  theme_classic()+
  theme(strip.text.y = element_text(angle = 0))+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("Completed_LTR_elements_genotypes_countsCluster.png", width = 7, height = 5)
  # ggsave("Completed_LTR_elements_genotypes_countsCluster_min2.png", width = 7, height = 5)



common_genotypes<-c("alpha_alpha_alpha",
                    "beta_beta_beta",
                    "alpha_beta_alpha",
                    "alpha_beta_beta",
                    "beta_alpha_beta",
                    "beta_beta_alpha")

common_genotypes<-c("alpha_alpha_alpha",
                    "beta_beta_beta",
                    "alpha_beta_alpha")

table_seq_dif4_bygroup_ed %>% 
  filter(d5_LTR!="NoLTR" & d3_LTR!="NoLTR" ) %>% 
  filter(d5_LTR!="gamma_iota" & d3_LTR!="gamma_iota" ) %>% 
  filter(d5_LTR!="other" & d3_LTR!="other" ) %>% 
  unite(genotype, c("d5_LTR", "core", "d3_LTR")) %>%
  filter(genotype %in% common_genotypes) %>%
  group_by(genotype, sample) %>% 
  summarise(N_clusters=length(unique(cluster))) %>%
  ungroup() %>%
  mutate(genotype=factor(genotype, levels=common_genotypes)) %>%
  ggplot(aes(x = sample, y = N_clusters)) +
  geom_bar(stat = "identity")+
  # geom_text(aes(x = reorder(genotype, -N_clusters), 
  #               y = N_clusters+20, label=N_clusters), size=6)+
  labs(x="Sample", y="Num. Clusters") + 
  facet_grid(genotype ~ ., scales = "free")+
  theme_classic()+
  theme(strip.text.y = element_text(angle = 0))+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("Completed_LTR_elements_genotypes_countsCluster_perSample2.png", width = 10, height = 6)


no_clonal_strains<-c("Pomberef", 
                     "JB879", 
                     "JB760_EBC074", 
                     "JB938", 
                     "JB869", 
                     "JB4_EBC069", 
                     "JB918_EBC111", 
                     "JB1110_EBC121", 
                     "JB873_EBC095", 
                     "JB929", 
                     "JB934_EBC115", 
                     "JB943", 
                     "JB900_EBC131", 
                     "JB854", 
                     "JB1180", 
                     "JB858_EBC087", 
                     "JB842_EBC080", 
                     "JB840", 
                     "JB872_EBC094", 
                     "JB853_EBC085", 
                     "JB939_EBC119", 
                     "JB874", 
                     "JB1205_EBC137", 
                     "JB1197_EBC135", 
                     "JB953", 
                     "JB837", 
                     "JB1206_EBC138", 
                     "JB864", 
                     "JB758", 
                     "DY34373", 
                     "DY39827")

no_clonal_strains<-c("Pomberef", 
                     "JB879", 
                     "JB760_EBC074", 
                     "JB938", 
                     "JB869", 
                     "JB4_EBC069", 
                     "JB918_EBC111", 
                     "JB1110_EBC121", 
                     "JB873_EBC095", 
                     "JB929", 
                     "JB934_EBC115", 
                     "JB943", 
                     "JB900_EBC131", 
                     "JB854", 
                     "JB1180", 
                     "JB858_EBC087", 
                     "JB842_EBC080", 
                     "JB840", 
                     "JB872_EBC094", 
                     "JB853_EBC085", 
                     "JB939_EBC119", 
                     "JB874", 
                     "JB1205_EBC137", 
                     "JB1197_EBC135", 
                     "JB953", 
                     "JB837", 
                     "JB1206_EBC138", 
                     "JB864", 
                     "JB758", 
                     "DY34373")


singleton_table_byGenotype<-table_seq_dif4_bygroup_ed %>% 
  filter(sample %in% no_clonal_strains) %>%
  filter(d5_LTR!="NoLTR" & d3_LTR!="NoLTR" ) %>% 
  filter(d5_LTR!="gamma_iota" & d3_LTR!="gamma_iota" ) %>% 
  filter(d5_LTR!="other" & d3_LTR!="other" ) %>% 
  unite(genotype, c("d5_LTR", "core", "d3_LTR")) %>%
  filter(genotype %in% common_genotypes) %>%
  group_by(genotype, cluster, sample) %>% 
  summarise(N_seq=n()) %>%
  group_by(genotype, cluster) %>% 
  summarise(N_samples_cluster=length(unique(sample)),
            singleton_sample=sample[1]) %>%
  filter(N_samples_cluster==1) %>%
  ungroup() %>%
  mutate(sample=factor(singleton_sample, 
                       levels=order_samples_pro_AncPop)) %>%
  group_by(genotype, sample) %>% 
  summarise(N_singleton_clusters=length(unique(cluster))) %>%
  ungroup()

singleton_table_byGenotype<-rbind(singleton_table_byGenotype, data.frame(genotype="alpha_alpha_alpha", 
                   sample=no_clonal_strains[!(no_clonal_strains %in% singleton_table_byGenotype$sample)], 
                   N_singleton_clusters=0))

singleton_table_byGenotype %>%
  mutate(sample=factor(sample, levels=order_samples_pro_AncPop)) %>%
  mutate(genotype=factor(genotype, levels=common_genotypes)) %>%
  ggplot(aes(x = sample, y = N_singleton_clusters, fill=genotype)) +
  geom_bar(stat = "identity")+
  # geom_text(aes(x = reorder(genotype, -N_clusters), 
  #               y = N_clusters+20, label=N_clusters), size=6)+
  labs(x="Sample", y="Num. Singleton. Clusters") + 
  facet_grid(genotype ~ ., scales = "free")+
  theme_classic()+
  theme(strip.text.y = element_text(angle = 0))+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=14, colour="black"),
        legend.position = "none") #+
  # ggsave("Completed_LTR_elements_singletons_genotypes_countsCluster_perSample.png", width = 10, height = 4)
  # ggsave("Completed_LTR_elements_singletons_genotypes_countsCluster_perSample2.png", width = 10, height = 4)


common_genotypes_sim<-c("alpha_alpha_alpha",
                    "beta_beta_beta",
                    "alpha_beta_alpha")

singleton_table_byGenotype %>% 
  filter(genotype %in% common_genotypes_sim) %>%
  mutate(sample=factor(sample, levels=order_samples_pro_AncPop)) %>%
  mutate(genotype=factor(genotype, levels=common_genotypes)) %>%
  ggplot(aes(x = sample, y = N_singleton_clusters, fill=genotype)) +
  geom_bar(stat = "identity")+
  # geom_text(aes(x = reorder(genotype, -N_clusters), 
  #               y = N_clusters+20, label=N_clusters), size=6)+
  labs(x="Sample", y="Num. Singleton. Clusters") + 
  facet_grid(genotype ~ ., scales = "free")+
  theme_classic()+
  theme(strip.text.y = element_text(angle = 0, size=14))+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=14, colour="black"),
        legend.position = "none", 
        axis.title = element_text(size=14)) #+
  # ggsave("Completed_LTR_elements_singletons_genotypes_countsCluster_perSample2.png", width = 10, height = 6)



# correlation ancestral admixute proportion Vs num singletons:
# Lineal Regression

no_clonal_strains<-c("JB22_EBC2", 
                     "JB879", 
                     "JB760_EBC074", 
                     "JB938", 
                     "JB869", 
                     "JB4_EBC069", 
                     "JB918_EBC111", 
                     "JB1110_EBC121", 
                     "JB873_EBC095", 
                     "JB929", 
                     "JB934_EBC115", 
                     "JB943", 
                     "JB900_EBC131", 
                     "JB854", 
                     "JB1180", 
                     "JB858_EBC087", 
                     "JB842_EBC080", 
                     "JB840", 
                     "JB872_EBC094", 
                     "JB853_EBC085", 
                     "JB939_EBC119", 
                     "JB874", 
                     "JB1205_EBC137", 
                     "JB1197_EBC135", 
                     "JB953", 
                     "JB837", 
                     "JB1206_EBC138", 
                     "JB864", 
                     "JB758", 
                     "DY34373", 
                     "DY39827")


no_clonal_strains<-c("JB22_EBC2", 
                     "JB879", 
                     "JB760_EBC074", 
                     "JB938", 
                     "JB869", 
                     "JB4_EBC069", 
                     "JB918_EBC111", 
                     "JB1110_EBC121", 
                     "JB873_EBC095", 
                     "JB929", 
                     "JB934_EBC115", 
                     "JB943",
                     "JB900_EBC131", 
                     "JB854", 
                     "JB1180", 
                     "JB858_EBC087", 
                     "JB842_EBC080", 
                     "JB840", 
                     "JB872_EBC094", 
                     "JB853_EBC085", 
                     "JB939_EBC119", 
                     "JB874", 
                     "JB1205_EBC137", 
                     "JB1197_EBC135", 
                     "JB953", 
                     "JB837", 
                     "JB1206_EBC138", 
                     "JB864", 
                     "JB758", 
                     "DY34373")

common_genotypes_sim<-c("alpha_alpha_alpha",
                        "beta_beta_beta")

singleton_table_byGenotype<-table_seq_dif4_bygroup_ed %>% 
  filter(sample %in% no_clonal_strains) %>%
  unite(genotype, c("d5_LTR", "core", "d3_LTR")) %>%
  filter(genotype %in% common_genotypes_sim) %>%
  group_by(genotype, cluster, sample) %>% 
  summarise(N_seq=n()) %>%
  group_by(genotype, cluster) %>% 
  summarise(N_samples_cluster=length(unique(sample)),
            singleton_sample=sample[1]) %>%
  filter(N_samples_cluster==1) %>%
  ungroup() %>%
  mutate(sample=factor(singleton_sample, 
                       levels=order_samples_pro_AncPop)) %>%
  group_by(genotype, sample) %>% 
  summarise(N_singleton_clusters=length(unique(cluster))) %>%
  ungroup()


# singleton_table_byGenotype<-rbind(singleton_table_byGenotype, 
#                                   data.frame(genotype="alpha_alpha_alpha", 
#                                              sample=no_clonal_strains[!(no_clonal_strains %in% singleton_table_byGenotype$sample)], N_singleton_clusters=0), 
#                                   data.frame(genotype="beta_beta_beta", 
#                                              sample=no_clonal_strains[!(no_clonal_strains %in% singleton_table_byGenotype$sample)], N_singleton_clusters=0))

singleton_table_byGenotype<-rbind(singleton_table_byGenotype, 
                                  data.frame(genotype="alpha_alpha_alpha", 
    sample=no_clonal_strains[!(no_clonal_strains %in% singleton_table_byGenotype$sample[singleton_table_byGenotype$genotype=="alpha_alpha_alpha"])], N_singleton_clusters=0), 
    data.frame(genotype="beta_beta_beta", 
    sample=no_clonal_strains[!(no_clonal_strains %in% singleton_table_byGenotype$sample[singleton_table_byGenotype$genotype=="beta_beta_beta"])], N_singleton_clusters=0))


statistical_test<-singleton_table_byGenotype %>% 
  ## Plot per Genotype: 
  # left_join(pro_AncPop, by="sample") %>% 
  # mutate(admixed_proportion=if_else(anc_prop<=0.5, 
  #                                   anc_prop, 1-anc_prop)) %>%
  # data.frame()
  ### total:
  group_by(sample) %>%
  summarise(total_N_singleton_clusters=sum(N_singleton_clusters)) %>%
  filter(total_N_singleton_clusters<35) %>%
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop)) %>%
  data.frame()

m1a<-lm(total_N_singleton_clusters~admixed_proportion, data=statistical_test)
result_stat<-summary(m1a)


library(ggpmisc)
singleton_table_byGenotype %>% 
  ## Plot per Genotype: 
  # left_join(pro_AncPop, by="sample") %>% 
  # mutate(admixed_proportion=if_else(anc_prop<=0.5, 
  #                                   anc_prop, 1-anc_prop)) %>%
  # ggplot(aes(admixed_proportion, N_singleton_clusters, colour=genotype))+
  # geom_point()+
  # geom_smooth(method='lm', formula= y~x, se=FALSE)+
  # stat_poly_eq(formula = y~x, 
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #              parse = TRUE)
  ### Plot with total:
  group_by(sample) %>%
  summarise(total_N_singleton_clusters=sum(N_singleton_clusters)) %>%
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop)) %>%
  ggplot(aes(admixed_proportion, total_N_singleton_clusters))+
  geom_point(alpha=0.5)+
  geom_smooth(method='lm', formula= y~x, se=FALSE)+
  annotate(geom="text", x=0.1, y=30, 
           label=paste0("Adj.R2 = ", 
                        format(round(result_stat$adj.r.squared, 2), nsmall = 2),
                        "\n",
                        "p-value = ",
                        format(round(result_stat$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="black", size=3)+
  # stat_poly_eq(formula = y~x,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE) +
  labs(x="Admixture proportion", y="Num. Singleton. Clusters") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("regressionAdmixtureProportionVsNSingletons.png", width = 6, height = 5)
  # ggsave("regressionAdmixtureProportionVsNSingletons2.png", width = 6, height = 5)



# statistical test per genotype (whole alpha or beta)

statistical_test<-singleton_table_byGenotype %>% 
  left_join(pro_AncPop, by="sample") %>%
  mutate(admixed_proportion=if_else(anc_prop<=0.5,
                                    anc_prop, 1-anc_prop)) %>%
  data.frame()

m1alpha<-lm(N_singleton_clusters~admixed_proportion, 
  data=statistical_test[statistical_test$genotype=="alpha_alpha_alpha",])
Anova(m1alpha, Type="III")
ralpha<-summary(m1alpha)

m1beta<-lm(N_singleton_clusters~admixed_proportion, 
            data=statistical_test[statistical_test$genotype=="beta_beta_beta",])
Anova(m1beta, Type="III")
rbeta<-summary(m1beta)


library(ggpmisc)
singleton_table_byGenotype %>% 
  ## Plot per Genotype: 
  left_join(pro_AncPop, by="sample") %>%
  mutate(admixed_proportion=if_else(anc_prop<=0.5,
                                    anc_prop, 1-anc_prop)) %>%
  ggplot(aes(admixed_proportion, N_singleton_clusters, colour=genotype))+
  geom_point(alpha=0.6)+
  geom_smooth(method='lm', formula= y~x, se=FALSE)+
  # stat_poly_eq(formula = y~x,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE)
  annotate(geom="text", x=0.1, y=30, 
           label=paste0("Adj.R2 = ", 
           format(round(ralpha$adj.r.squared, 2), nsmall = 2),
           "\n", "p-value = ", format(round(ralpha$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="red", size=3)+
  annotate(geom="text", x=0.2, y=30, 
           label=paste0("Adj.R2 = ", 
                        format(round(rbeta$adj.r.squared, 2), nsmall = 2),
                        "\n", "p-value = ", format(round(rbeta$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="blue", size=3)+
  scale_color_manual(values=c("red", "blue"))+
  labs(x="Admixture proportion", y="Num. Singleton. Clusters") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("regressionAdmixtureProportionVsNSingletons_perGenotype.png", width = 8, height = 5)
  # ggsave("regressionAdmixtureProportionVsNSingletons_perGenotype2.png", width = 8, height = 5)


#### Number of singletons per chromosome (density):

no_clonal_strains<-c("JB22_EBC2", 
                     "JB879", 
                     "JB760_EBC074", 
                     "JB938", 
                     "JB869", 
                     "JB4_EBC069", 
                     "JB918_EBC111", 
                     "JB1110_EBC121", 
                     "JB873_EBC095", 
                     "JB929", 
                     "JB934_EBC115", 
                     "JB943", 
                     "JB900_EBC131", 
                     "JB854", 
                     "JB1180", 
                     "JB858_EBC087", 
                     "JB842_EBC080", 
                     "JB840", 
                     "JB872_EBC094", 
                     "JB853_EBC085", 
                     "JB939_EBC119", 
                     "JB874", 
                     "JB1205_EBC137", 
                     "JB1197_EBC135", 
                     "JB953", 
                     "JB837", 
                     "JB1206_EBC138", 
                     "JB864", 
                     "JB758", 
                     "DY34373", 
                     "DY39827")


# common_genotypes_sim<-c("alpha_alpha_alpha","beta_beta_beta")

# common_genotypes_sim<-c("alpha_alpha_alpha")
# common_genotypes_sim<-c("beta_beta_beta")
# common_genotypes_sim<-c("alpha_beta_alpha")

# common_genotypes_sim<-c("alpha_alpha_alpha","beta_beta_beta","alpha_beta_alpha")

window_size=100000
chr<-c()
window<-c()
N_singletons<-c()

for (chr_used in c("I", "II", "III")) {
# chr_used<-"I"
temporal<-table_seq_dif4_bygroup_ed %>% 
  filter(sample %in% no_clonal_strains) %>%
  unite(genotype, c("d5_LTR", "core", "d3_LTR")) %>% 
  # filter(len_seq>4000) %>%
  filter(genotype %in% common_genotypes_sim) %>%
  group_by(genotype, cluster, sample, chr) %>% 
  summarise(N_seq=n(), pos_min=min(pos)) %>%
  group_by(genotype, cluster, chr) %>% 
  summarise(pos_minPerCluster=min(pos_min), 
            N_samples_cluster=length(unique(sample)),
            singleton_sample=sample[1]) %>% 
  filter(N_samples_cluster==1) %>%
  filter(chr!="AB325691")
temporal2 <- temporal %>%
  filter(chr==chr_used)
# last_point<-
windows<-seq(0,6000000,window_size)[seq(0,6000000,window_size)<max(temporal2$pos_minPerCluster)]
for (windowsID in seq(1,length(windows-1))){
# windowsID<-1
singletons<-temporal2 %>% 
  filter(pos_minPerCluster>=windows[windowsID] & 
           pos_minPerCluster<windows[windowsID+1]) %>%
  dim()
chr<-c(chr, chr_used)
window<-c(window, windowsID)
N_singletons<-c(N_singletons, singletons[1])
}
}

singletons_Window_table<-data.frame(chr,
           window,
           N_singletons) 

ggplot(singletons_Window_table, aes(chr, N_singletons)) +
  geom_boxplot() +
  geom_point(alpha=0.2) +  
  labs(x="Chromosome", y="Num. Singleton/100kb") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("SingletonsPerChrPer100kb.png", width = 8, height = 5)
  # ggsave("SingletonsPerChrPer100kb_alpha.png", width = 8, height = 5)
  # ggsave("SingletonsPerChrPer100kb_beta.png", width = 8, height = 5)
  ggsave("SingletonsPerChrPer100kb_mix.png", width = 8, height = 5)
  

m1<-lm(N_singletons~chr, 
           data=singletons_Window_table)
Anova(m1, Type="III")
summary(m1)


# singletons_Window_table_total<-data.frame(chr,
#                                     window,
#                                     N_singletons) 
# singletons_Window_table_alpha<-data.frame(chr,
#                                           window,
#                                           N_singletons) 
# singletons_Window_table_beta<-data.frame(chr,
#                                           window,
#                                           N_singletons) 
# singletons_Window_table_mix<-data.frame(chr,
#                                          window,
#                                          N_singletons) 

rbind(data.frame(singletons_Window_table_total, group="total"), 
      data.frame(singletons_Window_table_alpha, group="alpha_alpha_alpha"), 
      data.frame(singletons_Window_table_beta, group="beta_beta_beta"), 
      data.frame(singletons_Window_table_mix, group="alpha_beta_alpha")) %>%
  ggplot(aes(window/10, N_singletons, group=group, colour=group))+
  geom_line(aes(group=group))+
  geom_point()+
  scale_colour_manual(values=c("black", "blue", "red", "orange"))+
  #scale_alpha_manual(values=c(1, 0.4, 0.4, 0.4))+
  # facet_grid(. ~ chr, space="free", scale="free")+
  labs(x="Pos (Mb)", y="Num. Singleton/100kb") + 
  scale_x_continuous(breaks = seq(0,6,0.5))+
  scale_y_continuous(breaks = seq(0,14,2))+
  facet_grid(group ~ chr, space="free_x", scale="free_x") +
  theme_classic()+
  # theme_linedraw()+
  theme(panel.grid.major = element_line(colour = "gray90"),
        axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        legend.position = "none") #+
  # ggsave("SingletonsalongChromosomePer100kb.png", width = 12, height = 7)
  


###  

tree_file<-paste0("tf_1500.treefile")
tree_ed<-treeio::read.iqtree(tree_file)

fig_tree<- ggtree(tree_ed) + 
  theme_tree2() #+

test_tree<-as.phylo(as_tibble(tree_ed))




subsample_table<-table_seq_dif4_bygroup_ed %>% 
  filter(sample %in% no_clonal_strains) %>%
  unite(genotype, c("d5_LTR", "core", "d3_LTR")) %>%
  # filter(genotype=="beta_beta_beta") %>%
  # filter(genotype=="alpha_alpha_alpha") %>%
  filter(genotype=="alpha_beta_alpha") %>%
  mutate(label=seq_ID_original) %>%
  select(label, sim_ID_seq, sample_serie, sample, 
         chr, pos, len_seq, cluster, ancHaplo_mean)

subtree <- drop.tip(test_tree, test_tree$tip.label[!(test_tree$tip.label %in% subsample_table$label)])

ggtree(subtree)  %<+% subsample_table +
  # geom_tippoint(aes(colour=sample))+
  geom_tippoint(aes(colour=factor(ancHaplo_mean)))+
  geom_tiplab(aes(label=sample))+
  geom_nodelab()+
  # geom_tiplab(aes(label=cladeID))+
  xlim(c(0, 0.08))+
  # geom_label(aes(label=label))+
  theme_tree2() 

# writeXStringSet(ltr_alig[names(ltr_alig) %in% subsample_table$label], "alpha_beta_alpha.fasta", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")

#

sub_table_seq_dif2_bygroup<-table_seq_dif2_bygroup %>%
  filter(seq_ID_original %in% subsample_table$label) %>%
  mutate(sample2=factor(seq_ID_original, levels=test_tree$tip.label)) %>%
  mutate(sample=factor(sample, levels=unique(sample)))  %>%
  filter(sample!="42_JB943_1") %>%
  filter(bygroup!="Both") #%>%
  # filter(pos<7550)
  # filter(pos<7545)

labels_table<-sub_table_seq_dif2_bygroup %>%
  ungroup() %>%
  select(pos) %>%
  arrange(pos) %>%
  unique() %>%
  mutate(label_x=if_else(pos %in% c(0,583, 2001,4006,6293, 7106), as.character(pos), " ")) 

sub_table_seq_dif2_bygroup %>%
  ggplot(aes(factor(pos), sample2)) +
  geom_tile(aes(fill=bygroup))+
  scale_fill_manual(values = c("red", "blue", "black")) +
  scale_x_discrete(labels=labels_table$label_x)+
  labs(x="Position", y="Sequence") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size=9, colour="black"),
        axis.text.y = element_text(size=9, colour="black"),
        legend.position = "none") +
  # theme(axis.text.x = element_text(angle = 45, 
  #                                  hjust = 1, size=9, colour="black"), 
  #       axis.text.y = element_blank(), 
  #       legend.position = "none") +
  ggsave("alignment_beta_alpha_beta2.png", width = 18, height = 10)
  # ggsave("alignment_alpha_alpha_alpha.png", width = 18, height = 10)
  # ggsave("alignment_beta_beta_beta.png", width = 18, height = 10)
  














##### analyses of sequences in small windows

alig=1500
ltr_alig = readDNAStringSet(paste0("alig_all_tf_masked_conSeq_minLen1500_break.fasta"))


pairwise_diff<-c()
sample_1<-c()
sample_2<-c()
start<-c()
step_size<-50
start_pos<-1

for (seq_n1 in seq(1,length(ltr_alig))){
  print(seq_n1)
  #seq_n<-20
  for (seq_n2 in seq(1,length(ltr_alig))){
    seq1<-as.matrix(ltr_alig[seq_n1][[1]])[seq(start_pos,start_pos+step_size-1),]
    gaps_seq<-sum(seq1=="-")
    if (gaps_seq<(step_size*0.5)){
      if(seq_n2>seq_n1){
        seq2<-as.matrix(ltr_alig[seq_n2][[1]])[seq(start_pos,start_pos+step_size-1),]
        gaps_seq2<-sum(seq2=="-")
        if(gaps_seq2<(step_size*0.5)){
          SS<-sum(seq1!=seq2)
          pairwise_diff<-c(pairwise_diff, SS)
          sample_1<-c(sample_1, seq_n1)
          sample_2<-c(sample_2, seq_n2)
          start<-c(start, start_pos)
        }
      }
    }
  }
}



data.frame(start,sample_1,sample_2,pairwise_diff) %>% 
  ggplot(aes(pairwise_diff))+
  geom_histogram(binwidth = 1)

data.frame(start,sample_1,sample_2,pairwise_diff) %>% 
  ggplot(aes(sample_1, sample_2))+
  geom_tile(aes(fill=factor(pairwise_diff)))

data.frame(start,sample_1,sample_2,pairwise_diff) %>% 
  filter(pairwise_diff>30)





# BE CAREFUL THIS WILL DELETE EXISTING FILES!!!

alig=1500
ltr_alig = readDNAStringSet("alig_all_tf_masked_conSeq_minLen1500_break.fasta")

# for (step_size in c(30,50,80,100,130,150)){
#  for (diffAllowed in c(3,5,10,15)){
for (step_size in c(30)){
  for (diffAllowed in c(3)){
    clusters<-c()
    # step_size<-50
    # diffAllowed<-5
    # file_name<-paste0("cluster_seq_stepSize", step_size, "_diff", diffAllowed, ".txt") # excluding gaps
    file_name<-paste0("cluster_seq_stepSize", step_size, "_diff", diffAllowed, "_withGaps.txt")
    if (file.exists(file_name)) {
      #Delete file if it exists
      file.remove(file_name)
    }
    for(start_pos in seq(1, length(as.matrix(ltr_alig[1][[1]]))-step_size, step_size)){
      print(start_pos)
      for (seq_n1 in seq(1,length(ltr_alig))){
        #print(seq_n1)
        seq1<-as.matrix(ltr_alig[seq_n1][[1]])[seq(start_pos,start_pos+step_size-1),]
        gaps_seq<-sum(seq1=="-")
        if (gaps_seq>(step_size*0.5)){
          # sample<-c(sample, seq_n1)
          # sample_group<-c(sample_group, 0)
          # SS_incluster<-c(SS_incluster, gaps_seq)
          # gaps<-c(gaps, gaps_seq)
          # start<-c(start, start_pos)
          values<-data.frame(seq_n1, 0, gaps_seq, gaps_seq, start_pos)
          write.table(values, file=file_name, append = T, quote = F, 
                      row.names = F, col.names = F)
        } else {
          if(length(clusters)==0){
            clusters<-seq_n1
            # sample<-c(sample, seq_n1)
            # sample_group<-c(sample_group, seq_n1)
            # SS_incluster<-c(SS_incluster, 0)
            # gaps<-c(gaps, gaps_seq)
            # start<-c(start, start_pos)
            values<-data.frame(seq_n1, seq_n1, 0, gaps_seq, start_pos)
            write.table(values, file=file_name, append = T, quote = F, 
                        row.names = F, col.names = F)
          } else {
            sample_group_seq1<-0
            for (group_seq in (clusters)){
              if(sample_group_seq1==0){
                seq2<-as.matrix(ltr_alig[group_seq][[1]])[seq(start_pos,start_pos+step_size-1),]
                # total_variant_sites<-sum(seq1!=seq2 & seq1!="-" & seq2!="-") # excluding gaps
                total_variant_sites<-sum(seq1!=seq2)
                if(total_variant_sites<(diffAllowed+1)){
                  sample_group_seq1<-group_seq
                  # sample<-c(sample, seq_n1)
                  # sample_group<-c(sample_group, group_seq)
                  # SS_incluster<-c(SS_incluster, total_variant_sites)
                  # gaps<-c(gaps, gaps_seq)
                  # start<-c(start, start_pos)
                  values<-data.frame(seq_n1, group_seq, total_variant_sites, gaps_seq, start_pos)
                  write.table(values, file=file_name, append = T, quote = F, 
                              row.names = F, col.names = F)
                }
              } else {
                break
              }
            }
            if(sample_group_seq1==0){
              clusters<-c(clusters, seq_n1)
              # sample<-c(sample, seq_n1)
              # sample_group<-c(sample_group, seq_n1)
              # SS_incluster<-c(SS_incluster, 0)
              # gaps<-c(gaps, gaps_seq)
              # start<-c(start, start_pos)
              values<-data.frame(seq_n1, seq_n1, 0, gaps_seq, start_pos)
              write.table(values, file=file_name, append = T, quote = F, 
                          row.names = F, col.names = F)
            }
          }
        }
      }
      # clusters<-data.frame(sample, sample_group, SS_incluster, gaps, start) %>% 
      #   group_by(sample_group, start) %>%
      #   summarise(n_seq=n()) %>% 
      #   filter(n_seq>3 & sample_group!=0) %>% 
      #   ungroup() %>%
      #   select(sample_group) %>%
      #   unlist() %>% 
      #   as.vector()
      print("done")
      clusters<-read.table(file_name)
      names(clusters)<-c("sample", "sample_group", "SS_incluster", "gaps", "start")
      clusters<-clusters %>%
        group_by(sample_group, start) %>%
        summarise(n_seq=n()) %>%
        filter(n_seq>3, sample_group!=0) %>%
        #ungroup() %>%
        select(sample_group) %>%
        unique() %>% 
        unlist() %>%
        as.vector()
      print(length(clusters))
    }
  }
}


# I manually moved the files to "table_clusters" in case I run the look by mistake:
step_size<-30
diffAllowed<-3
file_name<-paste0("table_clusters/cluster_seq_stepSize", step_size, "_diff", diffAllowed, ".txt")
file_name<-paste0("table_clusters/cluster_seq_stepSize", step_size, "_diff", diffAllowed, "_withGaps.txt")

clusters_table<-read.table(file_name)
names(clusters_table)<-c("sample", "sample_group", "SS_incluster", "gaps", "start")

clusters_table %>% 
  arrange(factor(sample_group)) %>% 
  pull(sample_group) %>% 
  unique()
  head()
  
order_colours<-c(0, 2, 17, 8, 407, 463, 599, 12, 13, 77, 87, 88, 90, 92, 165, 215, 217, 219, 900, 254, 265,  335, 346, 358, 359, 367, 370, 389, 409, 410, 411, 415, 416, 424, 426, 429, 430, 431, 433, 444, 446, 458, 459, 461, 462, 1011, 1078, 277,605, 639, 644, 656, 659, 664, 667, 671, 849, 872, 874, 880, 890, 891, 895, 896, 898, 917, 931, 937, 951, 955, 958, 959, 960, 972, 977, 979, 986, 999, 1000, 1001, 1005, 1019, 1035, 1061, 1066, 1075, 1077, 1088, 1090, 1091, 1095, 1108, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160)

clusters_table %>% 
  merge(clusters_table %>% 
          group_by(sample_group, start) %>%
          summarise(n_seq=n()), 
        by=c("start", "sample_group")) %>% 
  filter(n_seq>3) %>%
  #filter(sample_group!=0) %>% 
  mutate(sample_group2=factor(sample_group, levels=c(order_colours))) %>% 
  ggplot(aes(start, sample))+
  geom_tile(aes(fill=factor(sample_group2)))+
  geom_vline(xintercept = 590, color="black")+
  geom_vline(xintercept = 7063, color="black")+
  # scale_fill_grey()+
  scale_fill_brewer(palette="Paired") +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  xlab("Position") +
  ylab("Sequence") +
  # theme_minimal()+
  theme_classic() +
  theme(legend.position = "none", 
    axis.title = element_text(size=14, colour="black", face = "bold"), 
    axis.text.x = element_text(size=12, colour="black"), 
    axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("LTR_completed_clusterSeq_stepSize", step_size, "_diff", diffAllowed, ".png"),
  #        width = 10, height = 6)



clusters_table %>% 
  merge(clusters_table %>% 
          group_by(sample_group, start) %>%
          summarise(n_seq=n()), 
        by=c("start", "sample_group")) %>% 
  #filter(n_seq>2) %>%
  #filter(sample_group!=0) %>% 
  # filter(sample==521) %>% # reference sequence
  # filter(sample==335) %>% # 49_JB918_EBC111_III
  filter(sample %in% c(521, 335)) %>% # 
  ggplot(aes(start, factor(sample, levels=c(335, 521))))+
  geom_tile(aes(fill=factor(sample_group)))+
  geom_vline(xintercept = 590, color="black")+
  geom_vline(xintercept = 7063, color="black")+
  scale_fill_brewer(palette="Paired") +
  scale_x_continuous(expand = c(0,0))+
  #scale_y_discrete(expand = c(0,0))+
  scale_y_discrete(labels=c("521" = "TF2", "335" = "TF1"), expand = c(0,0))+
  xlab("Position") +
  ylab(" ") +
  # theme_minimal()+
  theme_classic() +
  theme(legend.position = "none", 
    axis.title = element_text(size=14, colour="black", face = "bold"), 
    axis.text.x = element_text(size=12, colour="black"), 
    axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("LTR_completed_clusterSeq_Ref_stepSize",
  #               step_size, "_diff", diffAllowed, ".svg"),width = 10, height = 1)


ltr_alig[grep("1_Pomberef",names(ltr_alig))]
grep("1_Pomberef",names(ltr_alig))
grep("_49_JB918_EBC111_III",names(ltr_alig))
grep("TF1_107",names(ltr_alig))

#

clusters_table %>% 
  group_by(sample, sample_group) %>%
  summarise(n_windows=n()) %>% 
  filter(sample_group==0) %>%
  select(-sample_group)


order_colours2<-c(960, 0, 2, 17, 8, 407, 463, 599, 12, 13, 77, 87, 88, 90, 92, 165, 215, 217, 219, 900, 254, 265,  335, 346, 358, 359, 367, 370, 389, 409, 410, 411, 415, 416, 424, 426, 429, 430, 431, 433, 444, 446, 458, 459, 461, 462, 1011, 1078, 277,605, 639, 644, 656, 659, 664, 667, 671, 849, 872, 874, 880, 890, 891, 895, 896, 898, 917, 931, 937, 951, 955, 958, 959, 972, 977, 979, 986, 999, 1000, 1001, 1005, 1019, 1035, 1061, 1066, 1075, 1077, 1088, 1090, 1091, 1095, 1108, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160)

clusters_table %>% 
  left_join(clusters_table %>% 
              group_by(sample, sample_group) %>%
              summarise(n_windows=n()) %>% 
              filter(sample_group==0) %>%
              select(-sample_group), 
            by="sample") %>% 
  filter(n_windows<100) %>% 
  # pull(sample) %>%
  # unique() %>%
  # length()
  group_by(start, sample_group) %>%
  summarise(n_seq=n()) %>% 
  filter(n_seq>4) %>%
  #rbind(data.frame(start=1, sample_group=10, n_seq=1)) %>%
  # group_by(sample_group) %>%
  # summarise(num_seq=n()) %>%
  # arrange(-num_seq)
  # mutate(sample_g=factor(sample_group, levels=c(10, 6, 254,407,459,463,0))) %>%
  # mutate(sample_g=factor(sample_group, levels=c(2, 10, 254, 407, 891, 900, 1078, 1112, 1145, 0))) %>% # with windows of 30 pb and 3 dif mix. 
  # mutate(sample_g=factor(sample_group, levels=c(1078, 2, 463, 17, 407, 900, 599, 1011, 8, 217, 254, 896, 960, 0))) %>% # with windows of 30 pb and 3 dif mix but including gaps
  mutate(sample_g=factor(sample_group, levels=order_colours2)) %>% # with windows of 30 pb and 3 dif mix but including gaps
  filter(sample_group!=0) %>%
  ggplot(aes(x=start, y=n_seq/941, fill=sample_g)) +
  geom_bar(stat="identity")+
  geom_vline(xintercept = 590, color="black")+
  geom_vline(xintercept = 7063, color="black")+
  scale_fill_brewer(palette="Paired") +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0.05,0))+
  xlab("Position") +
  ylab("Proportion") +
  # theme_minimal()+
  theme_classic() +
  theme(legend.position = "none", 
        axis.title = element_text(size=14, colour="black", face = "bold"), 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("LTR_completed_clusterSeq_proportions_stepSize",
                # step_size, "_diff", diffAllowed, ".png"),width = 8, height = 1.5)


clusters_table_TF2<-clusters_table %>% 
  # filter(sample_group %in% c(6, 254,407,459,463,0)) %>% 
  # filter(sample_group %in% c(2, 10, 254, 407, 891, 900, 1078, 1112, 1145, 0)) %>%  #with windows of 30 pb and 3 dif mix. 
  filter(sample_group %in% c(1078, 2, 463, 17, 407, 900, 599, 1011, 0)) %>%  #with windows of 30 pb and 3 dif mix but with gaps
  filter(sample==521) %>%
  mutate(sample_group_TF2=sample_group) %>%
  select(start, sample_group_TF2) 

clusters_table_TF1<-clusters_table %>% 
  # filter(sample_group %in% c(6, 254,407,459,463,0)) %>% 
  # filter(sample_group %in% c(2, 10, 254, 407, 891, 900, 1078, 1112, 1145, 0)) %>%  #with windows of 30 pb and 3 dif mix. 
  filter(sample_group %in% c(1078, 2, 463, 17, 407, 900, 599, 1011, 0)) %>%  #with windows of 30 pb and 3 dif mix but with gaps
  filter(sample==335) %>%
  mutate(sample_group_TF1=sample_group) %>%
  select(start, sample_group_TF1) 


clusters_table_TF2 %>% 
  left_join(clusters_table_TF1, by="start") %>%
  filter(sample_group_TF2!=sample_group_TF1, sample_group_TF2!=0, sample_group_TF1!=0) %>%
  dim() %>%
  first()


length(unique(clusters_table$sample))

clusters_table %>% head()


# DO NOT RUN, IT WIL OVERWRITE FILE. READ EXISTING FILE BELOW.

if (file.exists("diff_haplotypes.txt")) {
  #Delete file if it exists
  file.remove("diff_haplotypes.txt")
}
for( seq_number in seq(1, length(unique(clusters_table$sample)))){
  print(seq_number)
  values<-clusters_table %>% 
    filter(sample_group %in% c(6, 254,407,459,463,0)) %>% 
    filter(sample==seq_number) %>%
    select(start, sample_group) %>%
    left_join(clusters_table_TF2, by="start") %>% 
    left_join(clusters_table_TF1, by="start") %>% 
    summarise(diff_TF2=sum(((sample_group!=sample_group_TF2) & sample_group!=0 & 
                              sample_group_TF2!=0)), 
              diff_TF1=sum(((sample_group!=sample_group_TF1) & sample_group!=0 & 
                              sample_group_TF1!=0))) %>% 
    mutate(sample_ID=seq_number)
  write.table(values, file="diff_haplotypes.txt", append = T, quote = F, 
              row.names = F, col.names = F)
}


diff_haplotype_table<-read.table("diff_haplotypes.txt")
names(diff_haplotype_table)<-c("diff_TF2", "diff_TF1", "sample")


haplotypes<-c(521, 335)
sample<-c()
haplotype_ID<-c()

# remove sample with large gaps:
completed_samples<-clusters_table %>% 
  # mutate(gap_window=(gaps>10)*1) %>%
  mutate(gap_window=(gaps>15)*1) %>% #with windows of 30 pb and 3 dif mix. 
  filter(gap_window!=1) %>% 
  group_by(sample) %>%
  summarise(N_windows=n()) %>%
  # ggplot(aes(N_windows)) +
  # geom_histogram()
  # filter(N_windows>250) %>%
  filter(N_windows>150) %>% # with windows of 30 pb and 3 dif mix. 
  pull(sample)


# The I run a loop in which for each window I take a reference haplotype and compare other sequences to the reference. If the tested sample differ in less than two windows, it is considered the same haplotype. If the tested sequence is not present in the list of haplotypes, it is added as a new haplotype to test following samples. The whole loop will repoduce a list of the common haplotypes present along the whole genome. 

clusters_table %>% head()
tail(completed_samples)

for (seq_number in completed_samples){
  # seq_number<-1
  print(seq_number)
  haplotype_found<-0
  for (haplotype_used in haplotypes) {
    # haplotype_used<-261
    # haplotype_used<-583
    if (haplotype_found==0){
      clusters_table_ref<-clusters_table %>% 
        # filter(sample_group %in% c(6, 254,407,459,463,0)) %>% 
        # mutate(sample_group_sim=ifelse(sample_group %in% c(6, 254,407,459,463,0), 
        #                                sample_group, (-1))) %>%
        # mutate(sample_group_sim=ifelse(sample_group %in% c(2, 10, 254, 407, 891, 900, 1078, 1112, 1145, 0), 
        #                                sample_group, (-1))) %>% # with windows of 30 pb and 3 dif mix. 
        mutate(sample_group_sim=ifelse(sample_group %in% c(1078, 2, 463, 17, 407, 900, 599, 1011, 0), 
                                       sample_group, (-1))) %>% # including gaps
        filter(sample==haplotype_used) %>%
        mutate(sample_group_ref=sample_group_sim) %>%
        select(start, sample_group_ref) 
      diff_sample<-clusters_table %>% 
        #filter(sample_group %in% c(6, 254,407,459,463,0)) %>% 
        # mutate(sample_group_sim=ifelse(sample_group %in% c(6, 254,407,459,463,0), 
        #                                sample_group, (-1))) %>%
        # mutate(sample_group_sim=ifelse(sample_group %in% c(2, 10, 254, 407, 891, 900, 1078, 1112, 1145, 0), 
        #                                sample_group, (-1))) %>%
        mutate(sample_group_sim=ifelse(sample_group %in% c(1078, 2, 463, 17, 407, 900, 599, 1011, 0), 
                                       sample_group, (-1))) %>% # including gaps
        filter(sample==seq_number) %>%
        select(start, sample_group_sim) %>%
        left_join(clusters_table_ref, by="start") %>%
        filter(sample_group_sim!=sample_group_ref,
               #sample_group_sim!=0, 
               sample_group_sim!=(-1), 
               #sample_group_ref!=0, 
               sample_group_ref!=(-1)) %>%
        dim() %>%
        first()
      if (diff_sample<2){
        sample<-c(sample,seq_number)
        haplotype_ID<-c(haplotype_ID, haplotype_used)
        haplotype_found<-1
      }
    }
  }
  if (haplotype_found==0){
    haplotypes<-c(haplotypes, seq_number)
    sample<-c(sample,seq_number)
    haplotype_ID<-c(haplotype_ID, seq_number)
  }
}

CommonHaplotypes_table<-data.frame(sample, haplotype_ID) 

# write.table(CommonHaplotypes_table, "table_CommonHaplotypes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
CommonHaplotypes_table<-read.table("table_CommonHaplotypes.txt", T)
CommonHaplotypes_table %>%
  head()

# summary statistics of the different haplotypes: 
sort_list_haplotypes<-CommonHaplotypes_table %>% 
  left_join(clusters_table %>% 
              group_by(sample) %>% 
              summarise(N_windows_NoGaps=n()), 
            by="sample") %>% 
  group_by(haplotype_ID) %>% 
  summarise(N_sequences_perHaplotype=n(), 
            sample_ref=sample[(N_windows_NoGaps==max(N_windows_NoGaps))][1]) %>%
  filter(N_sequences_perHaplotype>5) %>%
  arrange(-N_sequences_perHaplotype)

sort_list_haplotypes

sample_order<-CommonHaplotypes_table %>% 
  left_join(clusters_table %>% 
              group_by(sample) %>% 
              summarise(N_windows_NoGaps=n()), 
            by="sample") %>% 
  group_by(haplotype_ID) %>% 
  summarise(N_sequences_perHaplotype=n(), 
            sample_ref=sample[(N_windows_NoGaps==max(N_windows_NoGaps))][1]) %>%
  filter(N_sequences_perHaplotype>5) %>% 
  arrange(-N_sequences_perHaplotype) %>% 
  pull(sample_ref) %>% 
  rev()

sample_order<-CommonHaplotypes_table %>%
  group_by(haplotype_ID) %>%
  summarise(N_sequences_perHaplotype=n()) %>%
  filter(N_sequences_perHaplotype>5) %>%
  arrange(-N_sequences_perHaplotype) %>%
  pull(haplotype_ID) %>%
  rev()




table_plot<-clusters_table %>% 
  left_join(clusters_table %>% 
          group_by(sample_group, start) %>%
          summarise(n_seq=n()), 
        by=c("start", "sample_group")) %>% 
  left_join(clusters_table %>% 
              group_by(sample, sample_group) %>%
              summarise(n_windows=n()) %>% 
              filter(sample_group==0) %>%
              select(-sample_group), 
            by="sample") %>%
  # filter(n_windows<80) %>%
  #filter(sample_group %in% c(6, 254,407,459,463,0)) %>% 
  #filter(n_seq>2) %>%
  #filter(sample_group!=0) %>% 
  # filter(sample==261) %>% # reference sequence
  # filter(sample==1) %>% # TF2_I
  # filter(sample==583) %>% # TF1_107
  filter(sample %in% sample_order) %>% 
  # mutate(sample_or=factor(sample, levels=sample_order)) %>%
  mutate(sample_or=factor(sample))

p <- ggplot(table_plot, aes(start, sample_or))+
  geom_tile(aes(fill=factor(sample_group)))+
  # geom_vline(xintercept = 590, color="black")+
  # geom_vline(xintercept = 7063, color="black")+
  geom_vline(xintercept = 590, color="blue")+
  geom_vline(xintercept = 7063, color="blue")+
  scale_fill_grey()+
  # scale_fill_brewer(palette="Paired") +
  scale_x_continuous(expand = c(0,0))+
  #scale_y_discrete(expand = c(0,0))+
  #scale_y_discrete(labels=c("261" = "TF2", "583" = "TF1"), expand = c(0,0))+
  # scale_y_discrete(labels=c("335" = "TF2", "521" = "TF1"), expand = c(0,0))+
  scale_y_discrete(limits = rev(levels(table_plot$sample_or)), 
                                labels=c("335" = "TF1", "521" = "TF2"), expand = c(0,0))+
  xlab("Position") +
  ylab(" ") +
  # theme_minimal()+
  theme_classic() +
  theme(legend.position = "none", 
        axis.title = element_text(size=14, colour="black", face = "bold"), 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) 

p
p + ggsave(paste0("02_LTR_completed_CommonHaplotypes", 
              step_size, "_diff", diffAllowed, "_greyScale.png"), width = 8, height = 4)
p + ggsave(paste0("02_LTR_completed_CommonHaplotypes", 
                  step_size, "_diff", diffAllowed, "_greyScale.svg"), width = 8, height = 4)
# ggsave(paste0(p, "02_LTR_completed_CommonHaplotypes",
#               step_size, "_diff", diffAllowed, ".png"),width = 8, height = 4)
# ggsave(paste0(p, "02_LTR_completed_CommonHaplotypes",
              # step_size, "_diff", diffAllowed, ".svg"),width = 8, height = 4)



head(CommonHaplotypes_table)
clusters_table %>% 
  head()

table_plot_averageG<-clusters_table %>% 
  left_join(CommonHaplotypes_table, by="sample") %>%
  group_by(haplotype_ID, start, sample_group) %>% 
  summarise(N_seq_group=n()) %>% 
  group_by(haplotype_ID, start) %>% 
  summarise(mean_sample_group=sample_group[which(N_seq_group==max(N_seq_group))]) %>% 
  left_join(CommonHaplotypes_table %>%
              group_by(haplotype_ID) %>%
              summarise(N_sequences_perHaplotype=n()), by="haplotype_ID") %>% 
  filter(N_sequences_perHaplotype>5) %>% 
  mutate(haplotype_ID=factor(haplotype_ID))

head(table_plot_averageG)

library("RColorBrewer")

ggplot(table_plot_averageG, aes(start, haplotype_ID))+
  geom_tile(aes(fill=factor(mean_sample_group)))+
  geom_vline(xintercept = 590, color="black")+
  geom_vline(xintercept = 7063, color="black")+
  # scale_fill_brewer(palette="Paired") +
  scale_x_continuous(expand = c(0,0))+
  #scale_y_discrete(expand = c(0,0))+
  #scale_y_discrete(labels=c("261" = "TF2", "583" = "TF1"), expand = c(0,0))+
  # scale_y_discrete(labels=c("335" = "TF2", "521" = "TF1"), expand = c(0,0))+
  scale_y_discrete(limits = rev(levels(table_plot_averageG$haplotype_ID)), 
                   labels=c("335" = "TF1", "521" = "TF2"), expand = c(0,0))+
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#FB9A99", "#E31A1C", "#33A02C", "#FDBF6F", "#CAB2D6", "#FF7F00", "#FFFF99", "#B15928", "#6A3D9A")) +
  xlab("Position") +
  ylab(" ") +
  # theme_minimal()+
  theme_classic() +
  theme(legend.position = "none", 
        axis.title = element_text(size=14, colour="black", face = "bold"), 
        axis.text.x = element_text(size=12, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("02_LTR_completed_CommonHaplotypes_averageG",
  #             step_size, "_diff", diffAllowed, ".png"),width = 8, height = 4)
  # ggsave(paste0("02_LTR_completed_CommonHaplotypes_averageG",
              # step_size, "_diff", diffAllowed, ".svg"),width = 8, height = 4)


# Figure using TF1 and TF2 as reference:
table_plot_averageG %>%
  left_join(table_plot_averageG %>%
              ungroup() %>%
              filter(haplotype_ID == "335") %>%
              mutate(haplotype_TF1=mean_sample_group) %>%
              select(start, haplotype_TF1), by="start") %>%
  mutate(new_haplotypeCode=factor(ifelse(mean_sample_group==0, 0, 
                                  ifelse(mean_sample_group==haplotype_TF1, "TF1", mean_sample_group )))) %>%
  ggplot(aes(start, haplotype_ID))+
  geom_tile(aes(fill=factor(new_haplotypeCode)))+
  geom_vline(xintercept = 590, color="black")+
  geom_vline(xintercept = 7063, color="black")+
  # scale_fill_brewer(palette="Paired") +
  scale_x_continuous(expand = c(0,0))+
  #scale_y_discrete(expand = c(0,0))+
  #scale_y_discrete(labels=c("261" = "TF2", "583" = "TF1"), expand = c(0,0))+
  # scale_y_discrete(labels=c("335" = "TF2", "521" = "TF1"), expand = c(0,0))+
  scale_y_discrete(limits = rev(levels(table_plot_averageG$haplotype_ID)), 
                   labels=c("335" = "TF1", "521" = "TF2"),
                   expand = c(0,0))+
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#E31A1C", "#FB9A99", "#33A02C", "#FDBF6F", "#CAB2D6", "#FF7F00", "#FFFF99", "#B15928", "#6A3D9A")) +
  xlab("Position") +
  ylab(" ") +
  # theme_minimal()+
  theme_classic() +
  theme(legend.position = "none", 
      axis.title = element_text(size=14, colour="black", face = "bold"), 
      axis.text.x = element_text(size=12, colour="black"), 
      axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave(paste0("02_LTR_completed_CommonHaplotypes_averageG",
  #             step_size, "_diff", diffAllowed, "_refTF1.png"),width = 8, height = 4)
  ggsave(paste0("02_LTR_completed_CommonHaplotypes_averageG",
  step_size, "_diff", diffAllowed, "_refTF1.svg"),width = 8, height = 4)
  


# loading cluster information per sequence:
# ordered_table2<-read.table("../../../Annotation/annotation_clusters.txt", T, sep = " ")
ordered_table2<-read.table("../../../Annotation/annotation_clusters.txt", T, sep = ",")
head(ordered_table2)

alig=1500
ltr_alig = readDNAStringSet("alig_all_tf_masked_conSeq_minLen1500_break.fasta")


sampleN<-c()
seq_ID_list<-c()
for (seqAlig in seq(1,length(ltr_alig))){
#print(seqAlig)
# seqAlig<-669
  sampleN<-c(sampleN, seqAlig)
  seq_ID_found<-paste0("^", str_remove_all(str_remove_all(str_remove_all(str_remove_all(names(ltr_alig[seqAlig]), "_R_"), "_I*_[0-9]*_[0-9]*_[0-9]*"), "_AB325691_[0-9]*_[0-9]*_[0-9]*"), "_[0-9]*_[0-9]*$"), "_[J,L,I]*_")
  #print(ordered_table2$seq_ID[grep(seq_ID_found,ordered_table2$seq_ID)])
  if(length(ordered_table2$seq_ID[grep(seq_ID_found,ordered_table2$seq_ID)])==0){
    print(seqAlig)
    seq_ID_list<-c(seq_ID_list, NA)
  } else {
    seq_ID_list<-c(seq_ID_list, ordered_table2$seq_ID[grep(seq_ID_found,ordered_table2$seq_ID)])
  }
}


stats_perHaplotype<-CommonHaplotypes_table %>% 
  left_join(data.frame(sample=sampleN, 
                       seq_ID=seq_ID_list), by="sample") %>% 
  left_join(clusters_table %>% 
              group_by(sample, sample_group) %>%
              summarise(n_windows=n()) %>% 
              filter(sample_group==0) %>%
              select(-sample_group), 
            by="sample") %>%
  # filter(n_windows<80) %>%
  select(seq_ID, haplotype_ID) %>% 
  left_join(ordered_table2, by="seq_ID") %>%
  # filter(haplotype_ID==1)
  # filter(len_froSeq>3000) %>% 
  filter(haplotype_ID %in% sample_order) %>% 
  group_by(haplotype_ID) %>% 
  summarise(N_seq=n(), 
            #N_seq_min_3500kb=sum(len_froSeq>3500),
            N_samples=length(unique(sample)), 
            N_clusters=length(unique(cluster))) %>%
  arrange(-haplotype_ID)


anc_pro_order_samples<-ordered_table2 %>% 
  select(sample, anc_prop) %>% 
  filter(!(is.na(anc_prop))) %>% 
  unique() %>% 
  arrange(-anc_prop) %>% 
  pull(sample)

ordered_table2 %>% 
  select(sample, anc_prop) %>% 
  filter(!(is.na(anc_prop))) %>% 
  unique() 

CommonHaplotypes_table %>%
  left_join(data.frame(sample=sampleN, 
                       seq_ID=seq_ID_list), by="sample") %>% 
  filter(!(is.na(seq_ID))) %>%
  left_join(clusters_table %>% 
              group_by(sample, sample_group) %>%
              summarise(n_windows=n()) %>% 
              filter(sample_group==0) %>%
              select(-sample_group), 
            by="sample") %>% 
  # filter(n_windows<80) %>%
  select(seq_ID, haplotype_ID) %>% 
  left_join(ordered_table2, by="seq_ID") %>% 
  # filter(haplotype_ID==1)
  # filter(len_froSeq>3000) %>% 
  filter(haplotype_ID %in% sample_order) %>% 
  left_join(stats_perHaplotype, by="haplotype_ID") %>% head()
  filter(N_clusters>2) %>% 
  group_by(haplotype_ID, cluster, sample, anc_prop) %>% 
  summarise(N_seq=n()) %>% 
  group_by(haplotype_ID, sample, anc_prop) %>% 
  summarise(N_clusters=length(unique(cluster))) %>% 
  rbind(data.frame(haplotype_ID=c(335,335,335), 
        sample=c("JB864", 
                 "DY34373", 
                 "DY39827"), 
        anc_prop=c(0.13955071, 
                   0.07607192, 
                   0.07564192))) %>%
  mutate(haplotype_group=factor(ifelse(haplotype_ID %in% c(335, 521), as.character(haplotype_ID), "Other"), levels = c("335", "521", "Other"))) %>% 
  #mutate(sample_or=reorder(sample, -anc_prop)) %>% 
  mutate(sample_or=factor(sample, levels=anc_pro_order_samples)) %>% 
  ggplot(aes(x = sample_or, y = N_clusters, fill=factor(haplotype_ID))) +
  geom_bar(stat = "identity")+
  # geom_text(aes(x = reorder(genotype, -N_clusters), 
  #               y = N_clusters+20, label=N_clusters), size=6)+
  labs(x="Sample", y="Num. Clusters") + 
  #scale_fill_brewer(palette="Dark2") +
  scale_fill_manual(values=c("#922B21", "#BA4A00", "#1B4F72", 
                             "#9B59B6", "#E74C3C", "#1E8449", 
                             "#F1C40F", "#7E5109", "#283747"))+
  facet_grid(haplotype_group ~ ., scales = "free")+
  theme_classic()+
  theme(strip.text.y = element_text(angle = 0))+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size=12, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("03_Completed_LTR_elements_genotypes_countsCluster_perhaplotypes.svg", width = 10, height = 4)


table_seq_dif4_bygroup_ed %>%
  select(sim_ID_seq, d5_LTR, core, d3_LTR) %>%
  unique() %>%
  head()
#   
CommonHaplotypes_table %>%
  left_join(data.frame(sample=sampleN,
                       seq_ID=seq_ID_list), by="sample") %>%
  filter(!(is.na(seq_ID))) %>%
  left_join(clusters_table %>%
              group_by(sample, sample_group) %>%
              summarise(n_windows=n()) %>%
              filter(sample_group==0) %>%
              select(-sample_group),
            by="sample") %>%
  extract(seq_ID, "serial_ID_seq", "([0-9]*)_.", remove = F) %>% 
  filter(n_windows<80) %>%
  select(seq_ID, haplotype_ID, serial_ID_seq) %>%
  left_join(ordered_table2, by="seq_ID") %>%
  # filter(haplotype_ID %in% sample_order) %>%
  mutate(sim_ID_seq=paste0(serial_ID_seq, "_", sample)) %>%
  left_join(stats_perHaplotype, by="haplotype_ID") %>%
  left_join(table_seq_dif4_bygroup_ed %>%
              select(sim_ID_seq, d5_LTR, core, d3_LTR) %>%
              unique(), by="sim_ID_seq" ) %>%
  # mutate(equal_flankingLTR=(d5_LTR==d3_LTR)*1) %>%
  # filter(equal_flankingLTR==0) %>%
  # head()
  group_by(haplotype_ID, d5_LTR, core, d3_LTR) %>%
  summarise(n_clusters=length(unique(cluster))) %>%
  mutate(equal_flankingLTR=(d5_LTR==d3_LTR)*1) %>%
  # filter(haplotype_ID==521)
  filter(haplotype_ID %in% c(335,521,434,599, 896, 900,918,1002,1011,1024,1078)) %>%
  # filter(n_clusters!=1) %>%
  data.frame()
  dim()
  head()
  
  

# 
# "2_DY34373_II_4414188_90_4795_1_4795"
#   





getwd()

seq_ref_name<-c()
for (ID in stats_perHaplotype$haplotype_ID){
  seq_ref_name<-c(seq_ref_name, names(ltr_alig)[ID])
}


stats_perHaplotype<-cbind(stats_perHaplotype, seq_ref_name)

stats_perHaplotype %>% 
  left_join(table_seq_dif4_bygroup_ed %>% 
              mutate(seq_ref_name=seq_ID_original), 
            by="seq_ref_name") %>% 
  #arrange(chr, pos) %>% 
  ungroup() %>% 
  select(haplotype_ID, N_seq, N_samples, N_clusters, d5_LTR, core, d3_LTR) %>%
  arrange(haplotype_ID) %>%
  unique()


table_seq_dif4_bygroup_ed %>% 
  mutate(seq_ref_name=seq_ID_original) %>% 
  head()



seq_ref_name2<-c()
for (ID in CommonHaplotypes_table$sample){
  seq_ref_name2<-c(seq_ref_name2, names(ltr_alig)[ID])
}

library(dplyr)
library(tidytree) 
library(cowplot)

library(ggtree)

# LTR_tree_haplotypes<-ggtree(tree_ed) +
#   theme_tree2() +
#   geom_nodepoint(aes(alpha=as.factor(UFboot)), color="darkgreen", size=1) +
#   geom_treescale() +
#   scale_alpha_manual(values=c(0,0.8))+
#   theme(legend.position="right") 

LTR_tree_haplotypes<-ggtree(tree_ed, aes(color=factor(UFboot))) +
  theme_tree2() +
  geom_treescale() +
  scale_colour_manual(na.value = "white", values=c("black", "orange")) +
  #scale_alpha_manual(values=c(0,0.8))+
  theme(legend.position="right") 

g<- LTR_tree_haplotypes + theme(legend.position='none')

d <- filter(LTR_tree_haplotypes, isTip) %>% select(c(label, y))

p11 <-data.frame(seq_ID_original=seq_ref_name2, CommonHaplotypes_table) %>%
  select(-sample) %>% 
  left_join(table_seq_dif4_bygroup_ed, by="seq_ID_original") %>% 
  mutate(label=sample_serie) %>% 
  select(label, haplotype_ID, cluster) %>% 
  left_join(d, by='label') %>% 
  filter(haplotype_ID %in% c(335,434,521,
                             896,918,1002,
                             1011,1024)) %>% 
  ggplot(aes(y, factor(haplotype_ID)))+
  geom_tile(aes(fill=factor(haplotype_ID))) +
  scale_fill_manual(values=c("#922B21", "#BA4A00", "#1B4F72", 
                             "#9B59B6", "#E74C3C", "#1E8449", 
                             "#F1C40F", "#7E5109", "#283747"))+
  coord_flip() + 
  theme_tree2() +
  theme(legend.position='none', 
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, size=12, colour="black")) 

p11 <- p11 + ylim2(LTR_tree_haplotypes) 
  
plot_grid(g, p11, ncol=2, align='h', 
          rel_widths = c(0.8, 0.2)) #+
  # ggsave("04_tree_haplotypes_distribution.png", width = 8, height = 10)
  ggsave("04_tree_haplotypes_distribution.svg", width = 8, height = 10, dpi = 450)
  
  
  
# correlation ancestral admixute proportion Vs num singletons:
# Lineal Regression
# analyses done with haplotypes

no_clonal_strains<-c("JB22_EBC2", 
                     "JB879", 
                     "JB760_EBC074", 
                     "JB938", 
                     "JB869", 
                     "JB4_EBC069", 
                     "JB918_EBC111", 
                     "JB1110_EBC121", 
                     "JB873_EBC095", 
                     "JB929", 
                     "JB934_EBC115", 
                     "JB943", 
                     "JB900_EBC131", 
                     "JB854", 
                     "JB1180", 
                     "JB858_EBC087", 
                     "JB842_EBC080", 
                     "JB840", 
                     "JB872_EBC094", 
                     "JB853_EBC085", 
                     "JB939_EBC119", 
                     "JB874", 
                     "JB1205_EBC137", 
                     "JB1197_EBC135", 
                     "JB953", 
                     "JB837", 
                     "JB1206_EBC138", 
                     "JB864", 
                     "JB758", 
                     "DY34373", 
                     "DY39827")


no_clonal_strains<-c("JB22_EBC2", 
                     "JB879", 
                     "JB760_EBC074", 
                     "JB938", 
                     "JB869", 
                     "JB4_EBC069", 
                     "JB918_EBC111", 
                     "JB1110_EBC121", 
                     "JB873_EBC095", 
                     "JB929", 
                     "JB934_EBC115", 
                     "JB943",
                     "JB900_EBC131", 
                     "JB854", 
                     "JB1180", 
                     "JB858_EBC087", 
                     "JB842_EBC080", 
                     "JB840", 
                     "JB872_EBC094", 
                     "JB853_EBC085", 
                     "JB939_EBC119", 
                     "JB874", 
                     "JB1205_EBC137", 
                     "JB1197_EBC135", 
                     "JB953", 
                     "JB837", 
                     "JB1206_EBC138", 
                     "JB864", 
                     "JB758", 
                     "DY34373")


table_seq_dif4_bygroup_ed %>% 
  head()

head(seq_ref_name2)
head(CommonHaplotypes_table)
head(table_seq_dif4_bygroup_ed)
dim(table_seq_dif4_bygroup_ed)
table_seq_dif4_bygroup_ed$cluster

singleton_table_byHaplotype<-data.frame(seq_ID_original=seq_ref_name2, CommonHaplotypes_table) %>%
  select(-sample) %>% 
  left_join(table_seq_dif4_bygroup_ed, by="seq_ID_original") %>% 
  filter(sample %in% no_clonal_strains) %>%
  filter(haplotype_ID %in% c(335,521)) %>% 
  group_by(haplotype_ID, cluster, sample) %>% 
  summarise(N_seq=n()) %>%
  group_by(haplotype_ID, cluster) %>% 
  summarise(N_samples_cluster=length(unique(sample)),
            singleton_sample=sample[1]) %>%
  # filter(N_samples_cluster==1) %>%
  ungroup() %>%
  mutate(sample=factor(singleton_sample, 
                       levels=order_samples_pro_AncPop)) %>%
  group_by(haplotype_ID, sample) %>% 
  summarise(N_singleton_clusters=length(unique(cluster))) %>%
  ungroup()


head(singleton_table_byGenotype)

singleton_table_byHaplotype<-rbind(singleton_table_byHaplotype, 
                                  data.frame(haplotype_ID=335, 
                                             sample=no_clonal_strains[!(no_clonal_strains %in% singleton_table_byHaplotype$sample[singleton_table_byHaplotype$haplotype_ID==335])], N_singleton_clusters=0), 
                                  data.frame(haplotype_ID=521, 
                                             sample=no_clonal_strains[!(no_clonal_strains %in% singleton_table_byHaplotype$sample[singleton_table_byHaplotype$haplotype_ID==521])], N_singleton_clusters=0))

head(singleton_table_byHaplotype)

statistical_test<-singleton_table_byHaplotype %>% 
  group_by(sample) %>%
  summarise(total_N_singleton_clusters=sum(N_singleton_clusters)) %>%
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop)) %>%
  data.frame()

m1a<-lm(total_N_singleton_clusters~admixed_proportion, data=statistical_test)
result_stat<-summary(m1a)
result_stat

library(ggpmisc)
singleton_table_byHaplotype %>% 
  group_by(sample) %>%
  summarise(total_N_singleton_clusters=sum(N_singleton_clusters)) %>%
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop)) %>%
  ggplot(aes(admixed_proportion, total_N_singleton_clusters))+
  geom_point(alpha=0.5, size=3)+
  geom_smooth(method='lm', formula= y~x, se=FALSE, colour="black")+
  annotate(geom="text", x=0.1, y=30, 
           label=paste0("Adj.R2 = ", 
                        format(round(result_stat$adj.r.squared, 2), nsmall = 2),
                        "\n",
                        "p-value = ",
                        format(round(result_stat$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="black", size=3)+
  # stat_poly_eq(formula = y~x,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE) +
  labs(x="Admixture proportion", 
       # y="Num. Singleton. Clusters") + 
       y="Num. Clusters") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsbyHaplotype_allSeq.png", width = 5, height = 5)
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsbyHaplotype_allSeq.svg", width = 5, height = 5)




# statistical test per genotype (whole alpha or beta)

statistical_test<-singleton_table_byHaplotype %>% 
  left_join(pro_AncPop, by="sample") %>%
  mutate(admixed_proportion=if_else(anc_prop<=0.5,
                                    anc_prop, 1-anc_prop)) %>%
  data.frame()

m1_335<-lm(N_singleton_clusters~admixed_proportion, 
            data=statistical_test[statistical_test$haplotype_ID==335,])
#Anova(m1_335, Type="III")
r_335<-summary(m1_335)
r_335

m1_521<-lm(N_singleton_clusters~admixed_proportion,
           data=statistical_test[statistical_test$haplotype_ID==521,])
#Anova(m1beta, Type="III")
r_521<-summary(m1_521)
r_521

library(ggpmisc)
singleton_table_byHaplotype %>% 
  ## Plot per Genotype: 
  left_join(pro_AncPop, by="sample") %>%
  mutate(admixed_proportion=if_else(anc_prop<=0.5,
                                    anc_prop, 1-anc_prop)) %>%
  ggplot(aes(admixed_proportion, N_singleton_clusters, colour=factor(haplotype_ID)))+
  geom_point(alpha=0.7, size=3)+
  geom_smooth(method='lm', formula= y~x, se=FALSE)+
  # stat_poly_eq(formula = y~x,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE)
  annotate(geom="text", x=0.1, y=25, 
           label=paste0("Adj.R2 = ", 
                        format(round(r_335$adj.r.squared, 2), nsmall = 2),
                        "\n", "p-value = ", 
                        format(round(r_335$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="#922B21", size=4)+
  annotate(geom="text", x=0.1, y=20, 
           label=paste0("Adj.R2 = ", 
                        format(round(r_521$adj.r.squared, 2), nsmall = 2),
                        "\n", "p-value = ", 
                        format(round(r_521$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="#1B4F72", size=4)+
  scale_color_manual(values=c("#922B21", "#1B4F72"))+
  # scale_color_manual(values=c("red", "blue"))+
  labs(x="Admixture proportion", 
       y="Num. Clusters") + 
       # y="Num. Singleton. Clusters") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        legend.position = "none",
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsbyHaplotype_TF1andTF2_allSeq.png", width = 5, height = 5)
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsbyHaplotype_TF1andTF2_allSeq.svg", width = 5 , height = 5)


# The previous plot of correlations were produced by considering singletons when there is only one sequence in one sample at the most. However, this does not consider samples with solo LTR since the alignment does not contain solo LTRs. I produced a table with haplotypes ID, and move to with with the annotation table to produce the same analyses but excluding clusters with additional solo LTRs. 

data.frame(seq_ID_original=seq_ref_name2, CommonHaplotypes_table) %>%
  select(-sample) %>% 
  left_join(table_seq_dif4_bygroup_ed, by="seq_ID_original") %>% 
  left_join(pro_AncPop, by="sample") %>%
  select(sample, chr, pos, sim_ID_seq, cluster, haplotype_ID) %>% 
  write.table("table_haplotypesIDperClusterandSample.txt", quote = FALSE, 
              sep = "\t", row.names = FALSE, col.names = TRUE)


