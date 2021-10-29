# install.packages("tidyverse",  lib="/home/sergio/R/libraries_Rackham/3.4")
# module load R/3.6.1
# R_LIBS_USER=/home/sergio/R/libraries_Rackham/3.6


# rm(list=ls())

#library("ggplot2")
#library("tidyr")
#library("plyr")
#library("dplyr")
library("tidyverse")
library(wesanderson)
library(aplot)

getwd()
# setwd("C:/Users/sertu336/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Annotation")
# setwd("C:/Users/Sergio/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Annotation")
setwd("/Users/ru43sej/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Annotation")
  
#
ancestralhapl<-read.table("ancestralhap_57ILL_LR.txt", T) 
ancestralhapl<-ancestralhapl[-(grep("ILL", ancestralhapl$sample)),]

pro_AncPop <- ancestralhapl %>% 
  filter(Nor_PC1_pol_sim!=0.5) %>% 
  group_by(sample) %>% 
  summarise(proportion_sk=sum(Nor_PC1_pol_sim)/n()) 

order_samples_pro_AncPop <-pro_AncPop %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB879", 0.0007,proportion_sk)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB873", 0.4533954726,prop_sk_ed)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB1206", 0.8513139695,prop_sk_ed)) %>%   
  arrange(prop_sk_ed) %>% 
  ungroup() %>% 
  select(sample) %>% 
  unlist() %>% 
  as.vector()

pro_AncPop<-pro_AncPop %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB879", 0.0007,proportion_sk)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB873", 0.4533954726,prop_sk_ed)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB1206", 0.8513139695,prop_sk_ed)) %>%   
  arrange(prop_sk_ed) %>% 
  mutate(anc_prop=1-prop_sk_ed) %>% 
  select(sample, anc_prop)


order_samples_pro_AncPop <-c("LTR_loc", "ref_loc", "Pomberef", order_samples_pro_AncPop)


ancestralhapl_ref<-ancestralhapl %>% 
  filter(sample=="JB22") %>% 
  mutate(sample="Pomberef") 

ancestralhapl<-rbind(ancestralhapl, ancestralhapl_ref)

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
                     "JB1174_EBC132", 
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

# heatmap:


data_small_windows<-matrix(rep(c(rep(-1, 5563146/1000), rep(-1, 4537806/1000), rep(-1, 2440870/1000)), 
                               length(unique(ancestralhapl$sample))), ncol=length(unique(ancestralhapl$sample)))

data_small_windows <- as.data.frame(data_small_windows)
names(data_small_windows) <- unique(ancestralhapl$sample)

head(data_small_windows)

data_small_windows <- data_small_windows %>% 
  mutate(chromosome=c(rep("I", 5563146/1000), rep("II", 4537806/1000), rep("III", 2440870/1000)), 
         position_bin=c(seq(0,(5563146/1000)-1), seq(0, (4537806/1000)-1), seq(0, (2440870/1000)-1)))

data_small_windows$position_bin<-as.numeric(data_small_windows$position_bin)

mark_values <- function(a,b,c,d,e){ data_small_windows[data_small_windows$chromosome==a & data_small_windows$position_bin %in% c(round((b/1000)):round((c/1000))), d] <<- e}

table_tobin<-ancestralhapl

apply(table_tobin, 1 , function(x) mark_values(x[1],as.numeric(x[2]),as.numeric(x[3]),x[4],as.numeric(x[5])))

library(reshape2)
dm <- melt(data_small_windows,id.var=c("chromosome", "position_bin"))
names(dm)<-c("chromosome_name", "start_pos", "sample", "Nor_PC1_pol_sim")

matrix_window <-dm %>% 
  mutate(windowID=paste0(chromosome_name, "_", start_pos), Nor_PC1_pol_sim=as.factor(Nor_PC1_pol_sim)) %>% 
  dplyr::select(-chromosome_name, -start_pos) %>% 
  group_by(sample) %>% 
  spread(windowID, Nor_PC1_pol_sim) %>% ungroup() %>% 
  dplyr::select(paste0(dm$chromosome_name, "_", dm$start_pos)) %>% as.matrix()

samples_ID<- dm %>% mutate(windowID=paste0(chromosome_name, "_", start_pos)) %>% 
  dplyr::select(-chromosome_name, -start_pos) %>% 
  group_by(sample) %>% 
  spread(windowID, Nor_PC1_pol_sim) %>% 
  ungroup() %>% 
  dplyr::select(sample) %>% 
  unlist() %>% as.vector()

row.names(matrix_window)<- samples_ID

dm$sample <- factor(dm$sample, levels = order_samples_pro_AncPop)

dm %>% 
  filter(Nor_PC1_pol_sim>(-1)) %>%
  # ggplot(aes(start_pos*1000, sample)) + 
  ggplot(aes(start_pos/1000, sample)) + 
  geom_tile(aes(fill = Nor_PC1_pol_sim)) + 
  scale_fill_gradient(low = "darkred", high = "steelblue", expand=c(0,0)) + 
  scale_x_continuous(breaks=seq(0.5,5.5,0.5), expand=c(0,0))+
  # scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + 
  theme_classic() + 
  #ylab("\nSample") + 
  labs(x="Position (Mb)", y="Sample") + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="none", 
        axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y = element_text(colour="black")) + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") #+
  # ggsave("heatmap_longReadData.png", width = 10, height = 6, dpi = 450)





##### TF #####

samples_table<-read.table("all_TE_annotation_refCoor_Masked_conSeq.bed", F) 

names(samples_table)<-c("sample", "or_chr",	"or_start",	"or_end",	"seq_ID",	"ref_start_chr",	"ref_start",	"dis_start",	"ref_end_chr",	"ref_end",	"dis_end",	"length_ori_seq",	"direction",	"concordant_chr",	"chr_start_ed",	"pos_start_ed",	"dis_strat_ed", "chr_end_ed",	"pos_end_ed",	"dis_end_ed")

samples_table %>% 
  filter(sample=="Pomberef") %>% 
  filter(length_ori_seq>1000) %>% 
  arrange(chr_start_ed, pos_start_ed) %>% 
  mutate(len=pos_end_ed-pos_start_ed)

samples_table %>% 
  mutate(len=pos_end_ed-pos_start_ed) %>% 
  filter(len>10000)

samples_table<-samples_table %>% 
  mutate(len_tem=pos_end_ed-pos_start_ed) %>%
  # mutate(len_tem=pos_end_ed-dis_strat_ed) %>% 
  mutate(len=ifelse(len_tem>100000, length_ori_seq, len_tem)) %>% 
  mutate(chr=chr_start_ed, 
         start=pos_start_ed, 
         dis=dis_strat_ed, 
         len_froSeq=length_ori_seq, 
         addID="NA", 
         #dir=ifelse(dis_strat_ed=="-", "minus", "plus"),
         dir=direction, 
         rep="NA") %>% 
    select(seq_ID, sample, chr, start, dis, len_froSeq, addID, dir, len, rep) 

head(samples_table)


# names(samples_table)<-c("seq_ID", "sample", "chr", "start", "dis", "len_froSeq", "addID", "dir", "len", "rep")

# or_chr	or_start	or_end	seq_ID	ref_start_chr	ref_start	dis_start	ref_end_chr	ref_end	dis_end	length_ori_seq	direction	concordant_chr	chr_start_ed	pos_start_ed	dis_strat_ed

sp_tf_loc<-read.table("location_sp_tf.txt", T)

head(sp_tf_loc)

ordered_table_sp_tf<-sp_tf_loc %>% 
  mutate(seq_ID="NA", dis=0, len_froSeq=end-start, 
         addID=sample, 
         sample="ref_loc", 
         dir="+", len=end-start, rep="NA") %>% 
  select(seq_ID, sample, chr, start, dis, len_froSeq, addID, dir, len, rep) 

### Bowen 2003 LTR sequences:
head(ordered_table_sp_tf)
sp_bowenLTR<-read.csv("../../../LiLIn_results/second/Bowen_LTR_annotation.csv", T)

head(sp_bowenLTR)

ordered_table_bowen<-sp_bowenLTR %>% 
  mutate(dis=0, len_froSeq=end-start, 
         addID=clade_affiliation, 
         dir=LTR_orientation, 
         rep="NA") %>% 
  select(seq_ID, sample, chr, start, dis, len_froSeq, addID, dir, len, rep) 


head(ordered_table_bowen)


#join tables: 

head(samples_table)
ordered_table<-rbind(samples_table, ordered_table_sp_tf, ordered_table_bowen) %>% 
  # merge(ID_table, by="sample", all.x=T) %>% 
  arrange(chr, start, sample, dis) %>% 
  filter(!(is.na(chr)))


# ordered_table<-samples_table %>% 
#   merge(ID_table, by="sample", all.x=T) %>% 
#   arrange(chr, start, sample, dis, rep) %>% 
#   filter(seq_ID!="NA") %>% 
#   filter(!(is.na(chr)))

dim(ordered_table)
head(ordered_table)

cluster<-c(1)
length_cluster<-0
min_start<-0


# it worked quite ok with 150
# also with 100 bp
for (i in seq(2,dim(ordered_table)[1])){
# i=2
  print(i)
  if(length_cluster==0){
    length_cluster<-ordered_table$len[i-1]
    min_start<-ordered_table$start[i-1]
  }
  if(ordered_table$chr[i]==ordered_table$chr[i-1]){
    # if ((ordered_table$start[i]-(min_start+length_cluster))<=(-50)){
    if ((ordered_table$start[i]-(min_start+length_cluster))<=(100)){
      print(i)
      cluster[i]<-cluster[i-1]
      length_cluster<-max(length_cluster,ordered_table$len[i])
      min_start<-min(min_start, ordered_table$start[i])
    } else {
      print(i)
      cluster[i]<-cluster[i-1]+1
      length_cluster<-ordered_table$len[i]
      min_start<-ordered_table$start[i]
    }
  } else {
    print(i)
    cluster[i]<-cluster[i-1]+1
    length_cluster<-ordered_table$len[i]
    min_start<-ordered_table$start[i]
  }
}


ordered_table2<-cbind(ordered_table, cluster) %>% 
  mutate(sample = factor(sample, levels=order_samples_pro_AncPop)) 

head(ordered_table2)
ordered_table2 %>% 
  filter(cluster==517) %>%
  ungroup() %>%
  summarise(startP=min(start),
            endP=max(start))

# ordered_table2 %>% 
#   filter(cluster==514) #%>% 
#   # filter(start>4410171)
#   # filter(sample=="Pomberef") #%>% 
#   filter(sample=="JB934_EBC115")
#   ggplot(aes(start, len_froSeq)) +
#   geom_point()


ancestralhap<-ancestralhapl %>% 
  # filter(sample %in% ordered_table$sample_JB) %>%
  mutate(Nor_PC1_pol=(Nor_PC1_pol_sim==1)*1)

head(ancestralhap)

ancestralhap_GT_start<-c()
ancestralhap_GT_end<-c()

for (line in seq(1, dim(ordered_table2)[1])){
  print(line)
  start_point<-ancestralhap[ancestralhap$sample==as.vector(ordered_table2[line,"sample"]) & ancestralhap$chromosome_name==as.vector(ordered_table2[line,"chr"]) & ancestralhap$start_pos<=as.vector(ordered_table2[line,"start"]-50) & ancestralhap$end_ed>as.vector(ordered_table2[line,"start"]+50),"Nor_PC1_pol"]
  end_point<-ancestralhap[ancestralhap$sample==as.vector(ordered_table2[line,"sample"]) & ancestralhap$chromosome_name==as.vector(ordered_table2[line,"chr"]) & ancestralhap$start_pos<= as.vector(ordered_table2[line,"start"]+ordered_table2[line,"len"]-50) & ancestralhap$end_ed> as.vector(ordered_table2[line,"start"]+ordered_table2[line,"len"]+50),"Nor_PC1_pol"]
  if (length(start_point)==1){
    ancestralhap_GT_start[line]<-start_point
  } else if (length(start_point)==0){
    ancestralhap_GT_start[line]<-NA
  } else if (length(start_point)>1){
    ancestralhap_GT_start[line]<-2
  }
  if (length(end_point)==1){
    ancestralhap_GT_end[line]<-end_point
  } else if (length(end_point)==0){
    ancestralhap_GT_end[line]<-NA
  } else if (length(end_point)>1){
    ancestralhap_GT_end[line]<-2
  }
}

ordered_table2$ancHaplo<-ancestralhap_GT_start
ordered_table2$ancHaplo_end<-ancestralhap_GT_end
# ordered_table2<-ordered_table2 %>% 
#   merge(pro_AncPop, by="sample", all.x = T)

ordered_table2<-ordered_table2 %>%
  full_join(pro_AncPop, by="sample")



ordered_table2 %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  # select(sample) %>%
  # unique() %>%
  # dim()
  mutate(start_end=(ancHaplo==ancHaplo_end)*1) %>%
  filter(cluster==514) %>%
  select(-addID, -dir, -rep, -anc_prop, -chr) %>%
  tbl_df %>%
  print(n = Inf)
  ## Number of consistent flankin sequences (7954 consistent, 138 inconsistent, 734 with NA):
  # group_by(start_end) %>%
  # summarise(count=n())
  ## Number of clusters with at lest one sample inconsistent (55):
  # filter(start_end==0) %>% 
  # select(cluster) %>% 
  # unique() %>% dim()
  ##
  group_by(cluster, sample) %>%
  summarise(inconsistent_cluster=(sum((start_end==0)*1)!=0)*1) %>%
  filter(!(is.na(inconsistent_cluster))) %>%
  group_by(cluster) %>%
  summarise(total_NSamples=n(), inconsistent_Samples=sum((inconsistent_cluster==1)*1)) %>%
  # filter(inconsistent_Samples!=0) %>%
  filter(inconsistent_Samples>1) %>% 
  tbl_df %>%
  print(n = Inf)
  head()
  # dim()
  ggplot(aes(total_NSamples, inconsistent_Samples))+
  geom_point(alpha=0.5)
# tbl_df %>%
#   print(n = Inf)


#### identify clusters in sequences:

rbind(ordered_table2 %>% 
        filter(sample=="ref_loc"), 
      ordered_table2 %>% 
        filter(sample=="Pomberef"), 
      ordered_table2 %>% 
        filter(sample=="LTR_loc")) %>% 
  arrange(chr, start) %>% 
  filter(len_froSeq>1000) %>%
  select(-rep, -dir, -dis, -ancHaplo, -anc_prop) %>%
  tbl_df %>% print(n = Inf)

clean_ordered_table2<-ordered_table2 %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) 
  
tem<-str_remove_all(clean_ordered_table2$seq_ID, "J_")
tem<-str_remove_all(tem, "S_")
tem<-str_remove_all(tem, "L_")
tem2<-str_remove_all(tem, "_[0-9]*_I*_[0-9]*_[0-9]*")
tem2<-str_remove_all(tem2, "_[0-9]*_I*_I*_[0-9]*_[0-9]*")
tem2<-str_remove_all(tem2, "_I*_[0-9]*_[0-9]*")
tem2<-str_remove_all(tem2, "_[0-9]*_[0-9]*_[0-9]*")
tem2<-str_remove_all(tem2, "_AB325691_[0-9]*_[0-9]*")
tem2<-str_remove_all(tem2, "_9_MT")
sim_ID_seq<-str_remove_all(tem2, "_[0-9]*.[0-9]_I*")

# which(clean_ordered_table2$seq_ID=="111_JB934_EBC115_J_S_09.2_I_II_1834544_1834896")
# "84_JB934_EBC115_J_S_09.1_II_I_2231899_2232278" %in% clean_ordered_table2$seq_ID
# sim_ID_seq[1977]


clean_ordered_table2$sim_ID_seq<-sim_ID_seq

head(clean_ordered_table2)
clean_ordered_table2 %>% 
  filter(sim_ID_seq=="45_JB760")

clean_ordered_table3<-clean_ordered_table2 %>% 
  select(sim_ID_seq, cluster, dir, ancHaplo, ancHaplo_end, anc_prop)
head(clean_ordered_table3)
  
####

getwd()

# setwd("C:/Users/sertu336/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Phylogenies/Solo_AllLTR_allSeq/IQTreebb")
# setwd("C:/Users/Sergio/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Phylogenies/Solo_AllLTR_allSeq/IQTreebb")
setwd("/Users/ru43sej/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Phylogenies/Solo_AllLTR_allSeq/IQTreebb")


alignment_samples<-ape::read.FASTA("alig_ltr_allSolo_ltr_minLen100_plusRefseq.fasta", type = "DNA")

ref_seq<-c("TF1_107_ltr_3_Levin_etal_1990", 
           "TF2_21_Levin_etal_1990", 
           "TF2_22_Levin_etal_1990", 
           "LTR_solo")

seq_names<-names(alignment_samples) 
seq_names<-seq_names[!(seq_names %in% ref_seq)] 

tem<-str_remove_all(seq_names, "_R_")
tem2<-str_remove_all(tem, "_I*_[0-9]*_[0-9]*_[0-9]*")
tem2<-str_remove_all(tem2, "_AB325691_[0-9]*_[0-9]*_[0-9]*")
seq_ID<-str_remove_all(tem2, "_allSolo_ltr$")

tem2<-str_remove_all(tem, "_I*_[0-9]*_[0-9]*_[0-9]*_[0-9]*_allSolo_ltr$")
tem2<-str_remove_all(tem2, "_AB325691_[0-9]*_[0-9]*_[0-9]*_[0-9]*_allSolo_ltr$")
sim_ID_seq<-tem2
sample<-str_remove_all(tem2, "^[0-9]*_")

reverse_seq<-(substr(seq_names, 1, 3)=="_R_")*1

extract_len<-function(x){
  return(x[length(x)])
}
extract_dis<-function(x){
  return(x[length(x)-1])
}
extract_pos<-function(x){
  return(x[length(x)-2])
}
extract_chr<-function(x){
  return(unlist(x[length(x)-3]))
}

tem2<-str_remove_all(tem, "_[0-9]*_allSolo_ltr$")
chr<-unlist(sapply(strsplit(tem2, "_"), FUN=extract_chr))
dis<-as.numeric(unlist(sapply(strsplit(tem2, "_"), FUN=extract_dis)))
len_seq<-as.numeric(unlist(sapply(strsplit(tem2, "_"), FUN=extract_len)))
pos<-as.numeric(unlist(sapply(strsplit(tem2, "_"), FUN=extract_pos)))

table_seq_names<-data.frame(seq_ID, sample, 
                            chr, pos, dis, 
                            len_seq, sim_ID_seq, seq_names, reverse_seq) %>% 
  full_join(clean_ordered_table3, by="sim_ID_seq") 

head(table_seq_names)


### identify the LTR clade based on Bowen:
clade_ref_info<-ordered_table2 %>% 
  filter(sample=="LTR_loc") %>% 
  filter(addID!="") %>% 
  select(cluster, addID) %>% 
  unique() %>% 
  arrange(cluster) %>% 
  group_by(cluster) %>% 
  summarise(clade=paste(addID, collapse="_")) %>% 
  mutate(sample="Pomberef")


# ordered_table_seq2<-table_seq_names %>% 
#   merge(clade_ref_info, by=c("sample", "cluster"), all.x = T)

ordered_table_seq2<-table_seq_names %>% 
  full_join(clade_ref_info, by=c("sample", "cluster"))

ordered_table_seq2$anc_prop[ordered_table_seq2$sample=="Pomberef"]<-1


###

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")

#BiocManager::install("ggtree")
#BiocManager::install("treeio")

library(BiocManager)
library(Biostrings)

library("ggtree")
library("ggplot2")
library("tidyverse")
library("dplyr")

library(ggstance)
library(reshape2)

library(treeio)
library(ape)
library(tidytree)
library(aplot)

# sessionInfo()

# getwd()
tree_file<-paste0("allSoloLTR_100.treefile")
tree_ed<-treeio::read.iqtree(tree_file)
fig_tree<- ggtree(tree_ed) + 
  theme_tree2() #+
#geom_label2(aes(subset=!isTip, label=node), size=2, color="darkred", alpha=0.5) 
#

tree_ed@data$UFboot[is.na(tree_ed@data$UFboot)]<-0
tree_ed@data$UFboot[tree_ed@data$UFboot<90]<-0
tree_ed@data$UFboot[tree_ed@data$UFboot>=90]<-90

# I tested using over 50 but the result is pretty much the same:
# tree_ed@data$UFboot[is.na(tree_ed@data$UFboot)]<-0
# tree_ed@data$UFboot[tree_ed@data$UFboot<50]<-0
# tree_ed@data$UFboot[tree_ed@data$UFboot>=50]<-50


ggtree(tree_ed, aes(color=factor(UFboot))) +
  # theme_tree2() +
  # geom_nodepoint(aes(alpha=as.factor(UFboot)), color="darkgreen", size=1) +
  geom_treescale() +
  scale_colour_manual(na.value = "white", values=c("black", "orange")) #+
  #scale_alpha_manual(values=c(0,0.6))+
  # theme(legend.position="right") +
  # theme(legend.position="none") #+
  # ggsave("00_allseqLTR_tree_MLmin90.svg", width = 5, height =8)
  # ggsave("00_allseqLTR_tree_MLmin90.png", width = 5, height =8)


  


# change sequences IDs in the phylogeny:
tem<-rev(tree_ed@phylo$tip.label)
tem<-str_remove_all(tem, "_R_")
tem2<-str_remove_all(tem, "_I*_[0-9]*_[0-9]*_[0-9]*")
tem2<-str_remove_all(tem2, "_AB325691_[0-9]*_[0-9]*_[0-9]*")
seq_ID<-str_remove_all(tem2, "_allSolo_ltr$")
tree_ed@phylo$tip.label<-rev(seq_ID)


## plot per ref Clade:

library(ggplot2)
library(ggtree)

for (clade_used in unique(ordered_table_seq2$clade)[!(is.na(unique(ordered_table_seq2$clade)))]){
  print(clade_used)
  ordered_table3<-ordered_table_seq2 %>% 
    mutate(sample=seq_ID) %>% 
    filter(clade==clade_used) 
  ggtree(tree_ed) %<+% ordered_table3 +
    geom_tippoint(aes(color=clade, alpha=clade))+
    scale_alpha_manual(na.value = 0, values=c(1))+
    scale_colour_manual(na.value = "white", values=c("orange"))+
    geom_treescale(fontsize=6, offset=1)+
    theme(legend.position = "top") #+
    # ggsave(paste0("01_allseqLTR_tree_clade_",clade_used,".png"), width=7, height = 6)
  # ggtree(tree_ed, layout="daylight") %<+% ordered_table3 +
  #   geom_tippoint(aes(color=clade, alpha=clade))+
  #   scale_alpha_manual(na.value = 0, values=c(1))+
  #   scale_colour_manual(na.value = "white", values=c("orange"))+
  #   geom_treescale(fontsize=6, offset=1)+
  #   theme(legend.position = "left") +
  #   ggsave(paste0("01_2_allseqLTR_tree_clade_",clade_used,".png"), width=7, height = 6)
}
  

# distribution of large sequences. min 1500 bp length
ordered_table3<-ordered_table_seq2 %>% 
  mutate(sample=seq_ID) %>% 
  mutate(longSeq=ifelse(len_seq>2000, "1", NA)) %>% 
  filter(longSeq=="1") 

ggtree(tree_ed, aes(color=factor(UFboot))) %<+% ordered_table3 +
  # geom_tippoint(aes(colour=longSeq, alpha=longSeq))+
  geom_tippoint(aes(alpha=longSeq), colour="red")+
  # scale_colour_manual(na.value = "white", values=c("darkgreen"))+
  scale_colour_manual(na.value = "white", values=c("black", "orange")) +
  scale_alpha_manual(na.value = 0, values=c(1))+
  geom_treescale()+
  theme(legend.position = "left") #+
  # ggsave(paste0("02_allseqLTR_tree_min2000.svg"), width = 5, height =8)
  # ggsave(paste0("02_allseqLTR_tree_min2000.png"), width = 5, height =8)



ggtree(tree_ed, aes(color=factor(UFboot))) %<+% ordered_table3 +
  # geom_tippoint(aes(colour=longSeq, alpha=longSeq))+
  geom_tippoint(aes(alpha=Cluster_514), colour="red")+
  # scale_colour_manual(na.value = "white", values=c("darkgreen"))+
  scale_colour_manual(na.value = "white", values=c("black", "orange")) +
  scale_alpha_manual(na.value = 0, values=c(1))+
  geom_treescale()+
  theme(legend.position = "left") #+
  # ggsave(paste0("02_allseqLTR_tree_min2000.svg"), width = 5, height =8)
  # ggsave(paste0("02_allseqLTR_tree_example_Cluster514.png"), width = 5, height =8)






#### Identify clades for all samples:

tree_ed2<-read.tree(tree_file)
# plot(tree_ed2,show.tip.label = FALSE,  
#      show.node.label = FALSE, align.tip.label=TRUE)
# ape::nodelabels(frame = "none", bg = "none")
# 
# png("test.png",
#   width     = 15.25,
#   height    = 30.25,
#   units     = "in",
#   res       = 600,
#   pointsize = 4)
# par(mar = c(1, 1, 1, 1),  xaxs     = "i",
#   yaxs     = "i",  cex.axis = 2,cex.lab  = 2)
# plot(tree_ed2,show.tip.label = FALSE,  
#      show.node.label = FALSE, align.tip.label=TRUE)
# ape::nodelabels(frame = "none", bg = "none")
# dev.off()


#beta
clade_beta <- tree_subset(tree_ed, node=10256, levels_back=0)
beta_nodes<-as_tibble(clade_beta) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()


#alpha
clade_alpha <- tree_subset(tree_ed, node=9558, levels_back=0)
no_alpha_nodes<-as_tibble(clade_alpha) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

alpha_nodes<-unique(seq_ID[!(seq_ID %in% no_alpha_nodes)])

#delta
clade_delta <- tree_subset(tree_ed, node=12404, levels_back=0)
delta_nodes<-as_tibble(clade_delta) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

#gamma
clade_gamma <- tree_subset(tree_ed, node=13048, levels_back=0)
gamma_nodes<-as_tibble(clade_gamma) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

#iota
clade_iota <- tree_subset(tree_ed, node=13970, levels_back=0)
iota_nodes<-as_tibble(clade_iota) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

#epsilon
clade_epsilon <- tree_subset(tree_ed, node=14807, levels_back=0)
epsilon_nodes<-as_tibble(clade_epsilon) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

#zeta
clade_zeta <- tree_subset(tree_ed, node=15386, levels_back=0)
zeta_nodes<-as_tibble(clade_zeta) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()


## Mixed clades:


#beta_delta
clade_beta_delta <- tree_subset(tree_ed, node=10253, levels_back=0)

beta_delta_nodes_tem<-as_tibble(clade_beta_delta) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

beta_delta_nodes<-unique(beta_delta_nodes_tem[!(beta_delta_nodes_tem %in% beta_nodes)])


#delta_gamma
clade_delta_gamma <- tree_subset(tree_ed, node=10247, levels_back=0)

delta_gamma_nodes_tem<-as_tibble(clade_delta_gamma) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

delta_gamma_nodes<-unique(delta_gamma_nodes_tem[!(delta_gamma_nodes_tem %in% c(beta_nodes, beta_delta_nodes, delta_nodes))])


#gamma_iota
clade_gamma_iota1 <- tree_subset(tree_ed, node=13698, levels_back=0)
clade_gamma_iota2 <- tree_subset(tree_ed, node=13700, levels_back=0)

gamma_iota_nodes_tem1<-as_tibble(clade_gamma_iota1) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

gamma_iota_nodes_tem2<-as_tibble(clade_gamma_iota2) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

gamma_iota_nodes<-unique(gamma_iota_nodes_tem1, gamma_iota_nodes_tem2)


## subiota

clade_subiota <- tree_subset(tree_ed, node=13747, levels_back=0)

subiota_nodes_tem<-as_tibble(clade_subiota) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

subiota_nodes<-unique(subiota_nodes_tem[!(subiota_nodes_tem %in% c(iota_nodes))])


## zeta_alpha

clade_zeta_alpha <- tree_subset(tree_ed, node=17497, levels_back=0)

zeta_alpha_nodes<-as_tibble(clade_zeta_alpha) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()


## subalpha
clade_subalpha <- tree_subset(tree_ed, node=9559, levels_back=0)

subalpha_nodes<-as_tibble(clade_subalpha) %>% 
  filter(group==1) %>% 
  filter(label %in% tree_ed@phylo$tip.label) %>%
  ungroup() %>% 
  select(label) %>% 
  unlist() %>% 
  as.vector()

# tree_ed2 <- groupOTU(tree_ed, subalpha_nodes)
# ggtree(tree_ed2) +
#   # geom_tiplab() +
#   geom_tippoint(aes(colour=group)) #+
# # ggsave("test_tree.svg", width = 49, height = 49)






### merge tables 
seq_ID_tree<-seq_ID
general_clades_table<-rbind(data.frame(seq_ID=beta_nodes, 
                                       cladeID="beta"), 
      data.frame(seq_ID=beta_delta_nodes, cladeID="beta_delta"),
      data.frame(seq_ID=delta_nodes, cladeID="delta"),
      data.frame(seq_ID=delta_gamma_nodes, cladeID="delta_gamma"),
      data.frame(seq_ID=gamma_nodes, cladeID="gamma"), 
      data.frame(seq_ID=gamma_iota_nodes, cladeID="gamma_iota"), 
      data.frame(seq_ID=iota_nodes, cladeID="iota"), 
      data.frame(seq_ID=subiota_nodes, cladeID="subiota"), 
      data.frame(seq_ID=epsilon_nodes, cladeID="epsilon"), 
      data.frame(seq_ID=zeta_nodes, cladeID="zeta"), 
      data.frame(seq_ID=zeta_alpha_nodes, cladeID="zeta_subalpha"), 
      data.frame(seq_ID=subalpha_nodes, cladeID="subalpha"), 
      data.frame(seq_ID=alpha_nodes, cladeID="alpha")) %>%
  filter(seq_ID %in% tree_ed@phylo$tip.label) 

general_clades_table<-rbind(general_clades_table, 
                            data.frame(seq_ID=tree_ed@phylo$tip.label[!(tree_ed@phylo$tip.label %in% general_clades_table$seq_ID)], cladeID="other"))

head(general_clades_table)

ordered_table_seq2_withClades<-ordered_table_seq2 %>% 
  full_join(general_clades_table, by="seq_ID") %>%
  mutate(sample=factor(sample, levels=order_samples_pro_AncPop)) %>%
  filter(seq_ID %in% tree_ed@phylo$tip.label) %>%
  filter(!(is.na(cluster)))

setwd("/Users/ru43sej/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Phylogenies/Solo_AllLTR_allSeq/IQTreebb")

ordered_table_seq2_withClades %>%
  filter(sample %in% no_clonal_strains) %>%
  filter(!(is.na(cladeID))) %>%
  # filter(len_seq>2000) %>%
  filter(len_seq>4500) %>%
  filter(cladeID!="other") %>%
  mutate(cladeID=factor(cladeID, levels=c("beta", "beta_delta", 
                                          "delta", "delta_gamma", "gamma", 
                                          "gamma_iota", "iota", "subiota",
                                          "epsilon", "zeta", "zeta_subalpha",
                                          "subalpha", "alpha", "other"))) %>%
  # filter(len_seq<1000) %>%
  # gg plot(aes(cladeID, fill=cladeID)) +
  ggplot(aes(cladeID)) +
  geom_bar(stat="count", fill="black") +
  #scale_fill_manual(values=wes_palette(30, name = "Darjeeling1", type = "continuous"))+
  labs(x="Clade", y="Number of LTR seq") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black"), 
        legend.position = "none") #+
  # ggsave("03_LTR_elements_perClade.png", width = 4, height = 3)
  # ggsave("03_LTR_elements_perClade_min2000.png", width = 4, height = 3)
  # ggsave("03_LTR_elements_perClade_max1000.png", width = 4, height = 3)
  # ggsave("03_LTR_elements_perClade_min1500.png", width = 3, height = 2.5, dpi = 450)
  # ggsave("03_LTR_elements_perClade_min1500.svg", width = 3, height = 2.5, dpi = 450)
# ggsave("03_LTR_elements_perClade_min4500.png", width = 3, height = 2.5, dpi = 450)
# ggsave("03_LTR_elements_perClade_min4500.svg", width = 3, height = 2.5, dpi = 450)


ordered_table_seq2_withClades %>%
  filter(chr %in% c("I", "II", "III")) %>%
  filter(sample %in% no_clonal_strains) %>%
  # filter(len_seq>2000) %>%
  # filter(len_seq<1000) %>%
  ggplot(aes(chr, fill=factor(cladeID))) +
  geom_bar(stat="count", position=position_dodge()) +
  scale_fill_manual(values=wes_palette(30, name = "Darjeeling1", type = "continuous"))+
  labs(x="Chromosome", y="Number of LTR seq") + 
  theme_classic() +
  theme(axis.text.x = element_text(size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black"), 
        legend.position = "none") +
  facet_grid(. ~ cladeID, scale="free", space = "free_x") #+
  # ggsave("03_LTR_elements_perClade_perChr.png", width = 13, height = 4)
  # ggsave("03_LTR_elements_perClade_perChr_min2000.png", width = 10, height = 4)
  # ggsave("03_LTR_elements_perClade_perChr_max1000.png", width = 11, height = 4)


ordered_table_seq2_withClades %>%
  filter(sample %in% no_clonal_strains) %>%
  # filter(anc_prop<0.2) %>%
  # filter(len_seq>2000) %>%
  # filter(len_seq<1000) %>%
  filter(chr %in% c("I", "II", "III")) %>%
  ggplot(aes(pos/1000000, fill=cladeID))+
  geom_histogram(binwidth = 0.05)+
  geom_hline(yintercept=0, color = "black")+
  scale_fill_manual(values=wes_palette(30, name = "Darjeeling1", type = "continuous"))+
  scale_x_continuous(breaks=seq(0,5.5,0.5))+
  labs(x="Pos (Mb)", y="Number of LTR seq") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                   size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black"), 
        strip.text.y = element_text(size=10, angle = 0, colour="black"),
        legend.position = "none") +
  facet_grid(cladeID ~ chr, scale="free", space = "free_x") #+
  # ggsave("03_LTR_elements_perClade_pos.png", width = 8, height = 8)
  # ggsave("03_LTR_elements_perClade_pos_min2000.png", width = 8, height = 5)
  # ggsave("03_LTR_elements_perClade_pos_max1000.png", width = 8, height = 8)


ordered_table_seq2_withClades %>%
  filter(sample %in% no_clonal_strains) %>% 
  # filter(anc_prop<0.2) %>%
  # filter(len_seq>3000) %>%
  filter(chr %in% c("I", "II", "III")) %>%
  ggplot(aes(cluster, fill=cladeID))+
  geom_histogram(binwidth = 1)+
  facet_grid(cladeID ~ chr, scale="free", space = "free_x")


ordered_table_seq2_withClades %>%
  filter(sample %in% no_clonal_strains) %>% 
  # filter(anc_prop<0.2) %>%
  # filter(len_seq>3000) %>%
  filter(chr %in% c("I", "II", "III")) %>%
  filter(cluster==514) %>%
  mutate(FL=ifelse((len_seq>3000), "min3000", "lower")) %>%
  ggplot(aes(factor(ancHaplo), fill=cladeID))+
  geom_bar(stat="count", position=position_dodge()) +
  facet_grid(. ~ FL, scale="free", space = "free_x")

ordered_table_seq2_withClades %>% 
  filter(sample %in% no_clonal_strains) %>%
  filter(chr %in% c("I", "II", "III")) %>%
  ggplot(aes(sample, fill=cladeID))+
  geom_bar(stat="count", position=position_dodge())

ordered_table_seq2_withClades %>% 
  filter(sample %in% no_clonal_strains) %>%
  filter(chr %in% c("I", "II", "III")) %>%
  # filter(len_seq>2000) %>%
  filter(len_seq<1000) %>%
  ggplot(aes(sample, fill=cladeID))+
  geom_bar(stat="count")+
  scale_fill_manual(values=wes_palette(30, name = "Darjeeling1", type = "continuous"))+
  geom_hline(yintercept=0, color = "black", alpha=0.6)+
  facet_grid(cladeID ~ ., scale="free")+
  # facet_grid(cladeID ~ . , scale="free")+
  labs(x="Sample", y="Number of Seq") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black"), 
        strip.text.y = element_text(size=10, angle = 0, colour="black"),
        legend.position = "none") #+
  # ggsave("03_LTR_elements_perClade_perSample.png", width = 7, height = 9)
  # ggsave("03_LTR_elements_perClade_perSample_min2000.png", width = 7, height = 6)
  # ggsave("03_LTR_elements_perClade_perSample_max1000.png", width = 7, height = 9)



ordered_table_seq2_withClades %>%
  # filter(sample %in% no_clonal_strains) %>%
  filter(len_seq>2000) %>%
  filter(cluster==514) %>%
  select(-cluster, -chr, -anc_prop, -dis, -seq_names, -ancHaplo, -clade) %>%
  tbl_df %>%
  print(n = Inf)


ordered_table_seq2_withClades<-ordered_table_seq2_withClades %>% 
  extract(seq_ID, "seq_series", ".*(.)", remove = F) %>% 
  mutate(ori_seq_ID=str_remove_all(seq_ID, "_[0-9]$")) 

head(ordered_table_seq2_withClades)


### distribution of LTR per cluster
ordered_table_seq2_withClades2<- ordered_table_seq2_withClades %>%
  filter(!(is.na(cluster)))



ordered_table_seq2_withClades2 %>% head()
### these lines are just to check that the direction of the strain is consistent within sequences:
  # group_by(ori_seq_ID) %>%
  # summarise(total=n(), plusD=sum((dir=="+")*1)) %>%
  # mutate(diff=ifelse((total==plusD) | (plusD==0),1,0)) %>%
  # filter(diff!=0)



### homogenise direction with reference. 
seq_ID<-c()
seq_series_ed<-c()
dir_ed<-c()
reverse_seq_ed<-c()
for (seq_used in unique(ordered_table_seq2_withClades2$ori_seq_ID)) {
  #print(seq_used)
  tem_table<-ordered_table_seq2_withClades2 %>%
    filter(ori_seq_ID==seq_used)
  if(tem_table$dir[1]=="+"){
    seq_ID<-c(seq_ID, as.vector(unlist(tem_table$seq_ID)))
    seq_series_ed<-c(seq_series_ed,as.vector(unlist(tem_table$seq_series)))
    dir_ed<-rep("+", dim(tem_table)[1])
    reverse_seq_ed<-c(reverse_seq_ed,as.vector(unlist(tem_table$reverse_seq)))
  }
  if(tem_table$dir[1]=="-"){
    tem_table2<-tem_table %>%
      arrange(seq_series)
    seq_ID<-c(seq_ID, as.vector(unlist(tem_table2$seq_ID)))
    seq_series_ed<-c(seq_series_ed,rev(as.vector(unlist(tem_table2$seq_series))))
    dir_ed<-rep("+", dim(tem_table2)[1])
    reverse_seq_ed<-c(reverse_seq_ed,as.vector(1-unlist(tem_table2$reverse_seq)))
  }
}

ordered_table_seq2_withClades2 %>%
  head()

ordered_table_seq2_withClades3<-ordered_table_seq2_withClades2 %>%
  full_join(data.frame(seq_ID,
                       seq_series_ed,
                       dir_ed,
                       reverse_seq_ed), by="seq_ID") 

# write.table(ordered_table_seq2_withClades3, "LTR_seq_table.txt", quote = F, row.names = F)
# ordered_table_seq2_withClades3<-read.table("LTR_seq_table.txt", T)

cluster<-c()
sample<-c()
LRT_series<-c()
for(cluster_used in seq(1,max(ordered_table_seq2_withClades3$cluster))){
# cluster_used=565
  print(cluster_used)
  tem_table<-ordered_table_seq2_withClades3 %>% 
    filter(cluster==cluster_used) 
  tem_samples<-unique(tem_table$sample)
  for(sample_used in tem_samples){
    # print(sample_used)
    LRT_series_value<-tem_table %>%
      filter(sample==sample_used) %>%
      arrange(pos, dis, seq_series_ed) %>%
      mutate(dir_text=ifelse(reverse_seq_ed==0,"For","Rev"),
        cladeID_dir=paste(cladeID,dir_text, sep="")) %>%
      select(cladeID_dir) %>%
      unlist() %>% 
      as.vector() %>%
      paste(sep="", collapse="/") 
    cluster<-c(cluster, cluster_used)
    sample<-c(sample, sample_used)
    LRT_series<-c(LRT_series, LRT_series_value)
  }
}


# ordered_table_seq2_withClades3 %>%
#   # filter(cluster==514) %>%
#   arrange(pos, dis, seq_series_ed) %>%
#   group_by(sample) %>%
#   summarise(ancHaplo_start=ancHaplo[1]) 
ordered_table_seq2_withClades3 %>% 
  arrange(pos, dis, seq_series_ed) %>%
  group_by(sample) %>%
  head()

data.frame(cluster, 
           sample=factor(sample, levels=order_samples_pro_AncPop),
           LRT_series) %>% 
  left_join(ordered_table_seq2_withClades3 %>%
              select(sample, cluster, chr) %>%
              unique() , by=c("sample", "cluster")) %>%
  filter(chr %in% c("I", "II", "III")) %>%
  ggplot(aes(factor(cluster),sample))+
  geom_tile(aes(fill=LRT_series)) +
  scale_fill_manual(values=replicate(450, sample(wes_palette(450, name = "Darjeeling1", type = "continuous"))))+
  labs(x="Cluster", y="Sample") + 
  facet_grid(. ~ chr, scale="free", space="free")+
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=9, colour="black"), 
        strip.text.y = element_text(size=10, angle = 0, colour="black"),
        legend.position = "none") #+
  # ggsave("04_LTR_haplotype_perSample.png", width = 10, height = 7)
  


data.frame(cluster, 
           sample=factor(sample, levels=order_samples_pro_AncPop),
           LRT_series) %>% 
  filter(sample %in% no_clonal_strains) %>%
  full_join(data.frame(cluster, 
                       sample=factor(sample, levels=order_samples_pro_AncPop),
                       LRT_series) %>% 
              filter(sample %in% no_clonal_strains) %>%
              group_by(cluster, LRT_series) %>%
              summarise(N_samples=n()), by=c("cluster", "LRT_series")) %>%
  filter(N_samples==1) %>%
  ggplot(aes(sample)) +
  geom_bar(stat="count") +
  labs(x="Sample", y="Count") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=9, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave("04_LTR_singletons_haplotype_perSample.png", width = 7, height = 4)


data.frame(cluster, 
           sample=factor(sample, levels=order_samples_pro_AncPop),
           LRT_series) %>% 
  filter(sample %in% no_clonal_strains) %>%
  full_join(data.frame(cluster, 
                       sample=factor(sample, levels=order_samples_pro_AncPop),
                       LRT_series) %>% 
              filter(sample %in% no_clonal_strains) %>%
              group_by(cluster, LRT_series) %>%
              summarise(N_samples=n()), by=c("cluster", "LRT_series")) %>%
  filter(N_samples==1) %>%
  left_join(ordered_table_seq2_withClades3 %>%
              select(sample, cluster, chr) %>%
              unique() , by=c("sample", "cluster")) %>%
  filter(chr %in% c("I", "II", "III")) %>%
  ggplot(aes(cluster)) +
  geom_bar(stat="count") +
  labs(x="Cluster", y="Count") + 
  facet_grid(. ~ chr, scale="free")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=9, colour="black"), 
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave("04_LTR_singletons_alongChr.png", width = 7, height = 4)


data.frame(cluster, 
           sample=factor(sample, levels=order_samples_pro_AncPop),
           LRT_series) %>% 
  filter(sample %in% no_clonal_strains) %>%
  full_join(data.frame(cluster, 
                       sample=factor(sample, levels=order_samples_pro_AncPop),
                       LRT_series) %>% 
              filter(sample %in% no_clonal_strains) %>%
              group_by(cluster, LRT_series) %>%
              summarise(N_samples=n()), by=c("cluster", "LRT_series")) %>%
  filter(N_samples==1) %>%
  left_join(ordered_table_seq2_withClades3 %>%
              select(sample, cluster, chr) %>%
              unique() , by=c("sample", "cluster")) %>%
  filter(chr %in% c("I", "II", "III")) %>%
  group_by(chr) %>%
  summarise(N_singletons=n()) %>%
  mutate(chr_size=ifelse(chr=="I", 5.58, 
                         ifelse(chr=="III", 2.45, 4.54)), 
         singleton_perMb=N_singletons/chr_size) %>%
  ggplot(aes(chr, singleton_perMb)) +
  # ggplot(aes(chr, N_singletons)) +
  geom_bar(stat="identity") +
  # labs(x="Chromosome", y="Num. Singletons") +
  labs(x="Chromosome", y="Num. Singletons/Mb") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black")) #+
  # ggsave("04_LTR_singletons_perChr.png", width = 4, height = 3)
  # ggsave("04_LTR_singletons_perChr_PerMb.png", width = 4, height = 3)






##### analysis per cluster:

### Cluster 514 - 

ancHaplo_mean_table<-ordered_table_seq2_withClades3 %>%
  # filter(cluster==514) %>%
  group_by(cluster,sample) %>% 
  summarise(total=n()*2,
            ancHaplo_1=sum((ancHaplo==1)*1),
            ancHaplo_end_1=sum((ancHaplo_end==1)*1), 
            series_start=sum(((seq_series_ed==1) & pos==min(pos))*ancHaplo),
            series_end=sum(((seq_series_ed==max(as.numeric(seq_series_ed))) & 
                              pos==max(pos))*ancHaplo)) %>% 
  mutate(total_1=ancHaplo_1+ancHaplo_end_1, 
         ratio_1_total=total_1/total) %>%
  mutate(ancHaplo_mean=ifelse(ratio_1_total<0.5, 0, 
                              ifelse(ratio_1_total>0.5, 1, series_start))) %>%
  # filter(ratio_1_total>0 & ratio_1_total<1) %>%
  select(sample, cluster, ancHaplo_mean)
# tbl_df %>%
# print(n = Inf)

ordered_table_seq2_withClades4<-ordered_table_seq2_withClades3 %>%
  full_join(ancHaplo_mean_table, by=c("sample", "cluster")) 

head(ordered_table_seq2_withClades4)

# write.table(ordered_table_seq2_withClades4, "LTR_seq_table.txt", quote = F, row.names = F)
# ordered_table_seq2_withClades4<-read.table("LTR_seq_table.txt", T)

sample<-c()
LRT_series<-c()
cluster_used<-514
print(cluster_used)
tem_table<-ordered_table_seq2_withClades4 %>% 
  filter(cluster==cluster_used) %>%
  # filter(pos>4411000)
  filter(pos<4411000)
tem_samples<-unique(tem_table$sample)
for(sample_used in tem_samples){
  LRT_series_value<-tem_table %>%
    filter(sample==sample_used) %>%
    arrange(pos, dis, seq_series_ed) %>%
    mutate(directionName=ifelse(reverse_seq_ed==0, "For", "Rev"),
      cladeID_dir=paste(cladeID, directionName, sep="")) %>%
    select(cladeID_dir) %>%
    unlist() %>% 
    as.vector() %>%
    paste(sep="", collapse="/") 
  sample<-c(sample, sample_used)
  LRT_series<-c(LRT_series, LRT_series_value)
}

data.frame(cluster=514,
  sample=factor(sample, levels=order_samples_pro_AncPop),
           LRT_series) %>%
  left_join(ordered_table_seq2_withClades4 %>%
              filter(cluster==514) %>%
              filter(pos>4411000) %>%
              arrange(pos, dis, seq_series_ed) %>%
              group_by(sample) %>%
              summarise(ancHaplo_start=ancHaplo_mean[1]), by="sample") %>%
  mutate(sample=factor(sample, levels=order_samples_pro_AncPop)) %>%
  mutate(ancHaplo_start=factor(ancHaplo_start, levels=c(1,0))) %>%
  ggplot(aes("Haplotype",sample))+
  geom_tile(aes(fill=LRT_series)) +
  geom_text(aes(label=LRT_series))+
  scale_fill_manual(values=replicate(10, sample(wes_palette(10, name = "Darjeeling1", type = "continuous"))))+
  #facet_grid(. ~ ancHaplo_start, scale="free") +
  facet_grid(ancHaplo_start ~ ., scale="free") +
  theme_classic()+
  theme(legend.position = "none") #+
  # ggsave("04_LTR_cluster514_firstG_perSample_AncProp.png", width = 5, height = 7)
  # ggsave("04_LTR_cluster514_secondG_perSample_AncProp.png", width = 5, height = 7)



#####


ordered_table_seq2_withClades4 %>%
  filter(cluster==514) %>%
  filter(pos>4411000) %>% 
  filter(sample=="JB760_EBC074")
  head()

ordered_table_seq2_withClades %>%
  filter(seq_ID=="149_JB22_1") 
  unique()
  head()
"149_JB22_1"


#tree_ed<-treeio::read.iqtree(tree_file)
fig_tree<- ggtree(tree_ed) + 
  theme_tree2() #+

test_tree<-as.phylo(as_tibble(tree_ed))

# cluster_used<-514
for (cluster_used in unique(ordered_table_seq2_withClades4$cluster)){
  print(cluster_used)
  test<-ordered_table_seq2_withClades4 %>%
    mutate(label_seqID=str_remove_all(seq_ID, "_.*")) %>%
    mutate(label=seq_ID) %>%
    # filter(cladeID=="beta") %>%
    filter(cluster==cluster_used) %>%
    # filter(cladeID %in% c("subalpha", "alpha")) %>%
    # filter(pos>4411000) %>%
    mutate(directionName=ifelse(reverse_seq_ed==0, "For", "Rev")) %>%
    mutate(ID_sample=paste0(label_seqID, "_", sample, "_", 
                            seq_series_ed,"_", directionName)) %>%
    select(label, ID_sample, ancHaplo_mean, cladeID)
  if(dim(test)[1]>1 & 
     sum((test$ancHaplo_mean==0)*1, na.rm=T)>2 & 
     sum((test$ancHaplo_mean==1)*1, na.rm=T)>2 ){
    if(sum(!(is.na(test$ancHaplo_mean))*1)>0){
      subtree <- drop.tip(test_tree, test_tree$tip.label[!(test_tree$tip.label %in% test$label)])
      p1<-ggtree(subtree)  %<+% test +
        geom_tippoint(aes(colour=factor(ancHaplo_mean)))+
        geom_tiplab(aes(label=ID_sample))+
        # geom_tiplab(aes(label=cladeID))+
        xlim(c(0, 1.3))+
        # geom_label(aes(label=label))+
        theme_tree2() +
        theme(legend.position = "none")
        # ggsave(paste0("05_plots_perCluster/Cluster_",cluster_used,".png"), width = 6, height = 15)
      p2<-ggtree(subtree)  %<+% test +
        geom_tippoint(aes(colour=factor(ancHaplo_mean)))+
        # geom_tiplab(aes(label=ID_sample))+
        geom_tiplab(aes(label=cladeID))+
        xlim(c(0, 1.3))+
        # geom_label(aes(label=label))+
        theme_tree2() +
        theme(legend.position = "none") 
      p1 + p2 +
        ggsave(paste0("05_plots_perCluster/Cluster_",cluster_used,".png"), width = 10, height = 15)
        # ggsave(paste0("05_plots_perCluster_clades/Cluster_",cluster_used,".png"), width = 6, height = 15)
    } else {
      p1 <- ggtree(subtree)  %<+% test +
        geom_tippoint()+
        geom_tiplab(aes(label=ID_sample))+
        # geom_tiplab(aes(label=cladeID))+
        xlim(c(0, 1.3))+
        # geom_label(aes(label=label))+
        theme_tree2()+
        theme(legend.position = "none")
        # ggsave(paste0("05_plots_perCluster/NHan_Cluster_",cluster_used,".png"), width = 6, height = 15)
      p2 <- ggtree(subtree)  %<+% test +
        geom_tippoint()+
        # geom_tiplab(aes(label=ID_sample))+
        geom_tiplab(aes(label=cladeID))+
        xlim(c(0, 1.3))+
        # geom_label(aes(label=label))+
        theme_tree2()+
        theme(legend.position = "none")
      p1 + p2 +
        ggsave(paste0("05_plots_perCluster/NHan_Cluster_",cluster_used,".png"), width = 10, height = 15)
        # ggsave(paste0("05_plots_perCluster_clades/NHan_Cluster_",cluster_used,".png"), width = 6, height = 15)
    }
  }
}


# divergence between LTR elements of the same clade between Sp and Sk ancestral groups:

ordered_table_seq2_withClades4 %>% 
  arrange(cluster) %>% 
  pull(cluster) %>%
  unique() %>%
  length()
  head()

getwd()
setwd("/Users/ru43sej/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Phylogenies/Solo_AllLTR_allSeq/IQTreebb")
# setwd("C:/Users/Sergio/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Phylogenies/Solo_AllLTR_allSeq/IQTreebb")
# setwd("C:/Users/sertu336/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Phylogenies/Solo_AllLTR_allSeq/IQTreebb")


alignment_samples<-ape::read.FASTA("alig_ltr_allSolo_ltr_minLen100_plusRefseq.fasta", type = "DNA")


# cluster_used<-639
# clade_used<-"zeta"
# 
# unique(ordered_table_seq2_withClades4$cluster)
# unique(ordered_table_seq2_withClades4$cladeID)


ordered_table_seq2_withClades4 %>% head()

for (cluster_used in unique(ordered_table_seq2_withClades4$cluster)){
  for (clade_used in unique(ordered_table_seq2_withClades4$cladeID)){
    # print(cluster_used)
    # print(clade_used)
    tem<-ordered_table_seq2_withClades4 %>% 
      filter(cluster==cluster_used) %>%
      filter(cladeID==clade_used)
    tem_sp<-tem %>% 
      filter(ancHaplo_mean==0) %>% 
      select(seq_names) %>% 
      unlist() %>% 
      as.vector()
    tem_sk<-tem %>% 
      filter(ancHaplo_mean==1) %>% 
      select(seq_names) %>% 
      unlist() %>% 
      as.vector()
    data.frame(cluster=cluster_used, 
               clade=clade_used,
               N_sp=length(tem_sp), 
               N_sk=length(tem_sk)) %>% 
      write.table(paste0("shared_variants_Sp_Sk.txt"), 
                  row.names = F, col.names = F, 
                  quote = F, sep = "\t", append = T)
  }
}

table_byCluster_Clade<-read.table("shared_variants_Sp_Sk.txt", F)
names(table_byCluster_Clade)<-c("cluster", "clade", "N_sp", "N_sk")

head(table_byCluster_Clade)


tem_list_before<-ordered_table_seq2_withClades4 %>% 
  filter(!(is.na(ancHaplo_mean))) %>%
  arrange(cluster) %>% 
  pull(cluster) %>%
  unique() 

tem_list_after<-table_byCluster_Clade %>%
  # filter(N_sp + N_sk > 0) %>%
  arrange(cluster) %>% 
  pull(cluster) %>%
  unique() 

ordered_table_seq2_withClades4 %>% 
  filter(cluster==14)

ordered_table_seq2_withClades4 %>% 
  filter(!(cluster %in% tem_list_after)) %>%
  arrange(cluster) %>% 
  pull(cluster) %>%
  unique() %>%
  length()

sum(!(tem_list_after %in% tem_list_before))


table_byCluster_Clade %>% 
  filter(N_sp + N_sk >0) %>% 
  # pull(cluster) %>%
  # unique() %>%
  # length()
  # total of 662 clusters with identified flanking backgrounds
  mutate(group=ifelse((N_sp > 0 & N_sk > 0), "Both", "one_ancGroup")) %>%
  mutate(group2=ifelse(group == "Both", "Both", ifelse(N_sk > 0, "Sk_group", "Sp_group"))) %>%
  # filter(group=="Both") %>%
  # pull(cluster) %>%
  # unique() %>%
  # length()
  # 152 clusters have at least one sequence family in both background
  # filter(group!="Both") %>%
  # pull(cluster) %>%
  # unique() %>%
  # length()
  # 581 clusters contain at least one sequence family only found in one background
  group_by(cluster) %>% 
  summarise(Both=(sum(group2 == "Both")!=0)*1,
            Sk_group = (sum(group2 == "Sk_group")!=0)*1,
            Sp_group = (sum(group2 == "Sp_group")!=0)*1) %>%
  # filter(Both == 1 & (Sk_group + Sp_group > 0) ) %>%
  # dim()
  # 71 clusters contain both at least one sequence family shared between backgrounds and at least one sequence family only found in one background
  gather("group", "count", -cluster) %>%
  group_by(group) %>% 
  summarise(n_clusters=sum(count), 
            prop=sum(count)/n())



cluster_clades_bothAncBackground<-table_byCluster_Clade %>% 
  filter(N_sp + N_sk >0) %>% 
  mutate(group=ifelse((N_sp > 0 & N_sk > 0), "Both", "one_ancGroup")) %>%
  mutate(group2=ifelse(group == "Both", "Both", ifelse(N_sk > 0, "Sk_group", "Sp_group"))) %>%
  filter(group=="Both") %>% 
  mutate(subClusters_clades=paste0(cluster,"_",clade), 
         found=1) %>%
  select(subClusters_clades, found) 
  # pull(cluster) %>%
  # unique() %>%
  # length()
  # dim()
  # head()

ordered_table_seq2_withClades4 %>% 
  mutate(subClusters_clades=paste0(cluster,"_",cladeID)) %>%
  left_join(cluster_clades_bothAncBackground, by=c("subClusters_clades")) %>%
  filter(!(is.na(found))) %>%
  filter(ancHaplo_mean %in% c(0,1)) %>%
  group_by(cluster ,cladeID, ancHaplo_mean) %>%
  summarise(NSample=length(unique(sample))) %>%
  mutate(min2Samples=(NSample>1)*1) %>%
  group_by(cluster ,cladeID) %>%
  summarise(min2Samples_perBG=sum(min2Samples)) %>%
  filter(min2Samples_perBG==2) %>%
  pull(cluster) %>%
  unique() %>%
  length()


# DO NOT RUN THIS SECTION. IT WILL OVER WRITE FILE> LOAD TABLE DIRECTLY BELOW.
min_samplesSize<-4
for (cluster_used in unique(ordered_table_seq2_withClades4$cluster)){
  for (clade_used in unique(ordered_table_seq2_withClades4$cladeID)){
    # print(cluster_used)
    # print(clade_used)
  tem<-ordered_table_seq2_withClades4 %>% 
    filter(cluster==cluster_used) %>%
    filter(cladeID==clade_used)
  tem_sp<-tem %>% 
    filter(ancHaplo_mean==0) %>% 
    select(seq_names) %>% 
    unlist() %>% 
    as.vector()
  tem_sk<-tem %>% 
    filter(ancHaplo_mean==1) %>% 
    select(seq_names) %>% 
    unlist() %>% 
    as.vector()
  if(((length(c(tem_sp, tem_sk))>min_samplesSize) & 
      (length(tem_sp>1)) & (length(tem_sk>1)))){
  # if(length(c(tem_sp, tem_sk))>min_samplesSize){
  distance_table<-ape::dist.dna(alignment_samples[c(tem_sp, tem_sk)], 
                                model = "K80", pairwise.deletion = TRUE, 
                                as.matrix=TRUE)
  # print(cluster_used)
  # print(clade_used)
  distance_table[lower.tri(distance_table, diag = TRUE)] <- NA
  data.frame(sample1=rownames(distance_table)[row(distance_table)], 
             sample2=colnames(distance_table)[col(distance_table)],
             distance=c(distance_table)) %>% 
    filter(!(is.na(distance))) %>%
    mutate(group=factor(ifelse(((sample1 %in% tem_sp) & (sample2 %in% tem_sp)), "Sp_Sp", 
                               ifelse(((sample1 %in% tem_sk) & 
                                         (sample2 %in% tem_sk)), "Sk_Sk", 
                                      "Sp_Sk")), 
                        levels=c("Sp_Sp", "Sk_Sk", "Sp_Sk"))) %>% 
    group_by(group) %>%
    summarise(mean_Dis=mean(distance), 
              var_Dis=var(distance)) %>%
    mutate(cluster=cluster_used, 
           cladeID=clade_used, 
           N_Sp=length(tem_sp), 
           N_Sk=length(tem_sk)) %>% 
    write.table(paste0("pairdistance_Sp_Sk_min", min_samplesSize, ".txt"), 
                row.names = F, col.names = F, 
                quote = F, sep = "\t", append = T)
  }
  }
}

distance_table_byGroup<-read.table("pairdistance_Sp_Sk_min4.txt", F)
names(distance_table_byGroup)<-c("group", "mean_Dis", "var_Dis", "cluster", 
                                 "cladeID", "N_Sp", "N_Sk")

# distance_table_byGroup<-read.table("pairdistance_Sp_Sk_min10.txt", F)
# names(distance_table_byGroup)<-c("group", "mean_Dis", "var_Dis", "cluster", 
#                                  "cladeID", "N_Sp", "N_Sk")


distance_table_byGroup %>% 
  ggplot(aes(factor(cluster), mean_Dis, colour=group)) +
  geom_point() +
  geom_line(aes(group=group))

distance_table_byGroup %>% 
  filter(cladeID!="other") %>% 
  mutate(group=factor(group, levels=c("Sp_Sp", "Sk_Sk", "Sp_Sk"))) %>% 
  ggplot(aes(mean_Dis, fill=group)) +
  geom_histogram(aes(y=..density..), position="dodge")+
  geom_vline(data=distance_table_byGroup %>% 
               group_by(group) %>% 
               summarise(mean_distance_Group=mean(mean_Dis)) %>% 
               mutate(group=factor(group, levels=c("Sp_Sp", "Sk_Sk", "Sp_Sk"))), 
             aes(xintercept=mean_distance_Group, color=group),
             linetype="dashed") +
  scale_fill_manual(values=c("darkred", "steelblue", "orange")) +
  scale_color_manual(values=c("darkred", "steelblue", "orange")) +
  labs(x="Mean K2P distance", y="Density") + 
  theme_classic() +
  theme(axis.text.x = element_text(size=14, colour="black"), 
        axis.text.y = element_text(size=14, colour="black")) #+
  # ggsave("05_pairdistance_Sp_Sk_min4.png", width = 6, height =4)
  # ggsave("05_pairdistance_Sp_Sk_min4.svg", width = 6, height =4, dpi = 450)
  # ggsave("05_pairdistance_Sp_Sk_min10.png", width = 6, height =4)
  # ggsave("05_pairdistance_Sp_Sk_min10.svg", width = 6, height =4, dpi = 450)



# analyses using only full-length LTRs

ordered_table_seq2_withClades5 <- ordered_table_seq2_withClades4 %>%
  filter(len_seq>4500)

min_samplesSize<-3
for (cluster_used in unique(ordered_table_seq2_withClades5$cluster)){
  for (clade_used in unique(ordered_table_seq2_withClades5$cladeID)){
    # print(cluster_used)
    # print(clade_used)
    tem<-ordered_table_seq2_withClades5 %>% 
      filter(cluster==cluster_used & cladeID==clade_used) 
    tem_sp<-tem %>% 
      filter(ancHaplo_mean==0) %>% 
      select(seq_names) %>% 
      unlist() %>% 
      as.vector()
    tem_sk<-tem %>% 
      filter(ancHaplo_mean==1) %>% 
      select(seq_names) %>% 
      unlist() %>% 
      as.vector()
    if(((length(c(tem_sp, tem_sk))>min_samplesSize) & 
       (length(tem_sp>1)) & (length(tem_sk>1)))){
      distance_table<-ape::dist.dna(alignment_samples[c(tem_sp, tem_sk)], 
                                    model = "K80", pairwise.deletion = TRUE, 
                                    as.matrix=TRUE)
      distance_table[lower.tri(distance_table, diag = TRUE)] <- NA
      data.frame(sample1=rownames(distance_table)[row(distance_table)], 
                 sample2=colnames(distance_table)[col(distance_table)],
                 distance=c(distance_table)) %>% 
        filter(!(is.na(distance))) %>%
        mutate(group=factor(ifelse(((sample1 %in% tem_sp) & (sample2 %in% tem_sp)), "Sp_Sp", 
                                   ifelse(((sample1 %in% tem_sk) & 
                                             (sample2 %in% tem_sk)), "Sk_Sk", 
                                          "Sp_Sk")), 
                            levels=c("Sp_Sp", "Sk_Sk", "Sp_Sk"))) %>% 
        group_by(group) %>%
        summarise(mean_Dis=mean(distance), 
                  var_Dis=var(distance)) %>%
        mutate(cluster=cluster_used, 
               cladeID=clade_used, 
               N_Sp=length(tem_sp), 
               N_Sk=length(tem_sk)) %>% 
        write.table(paste0("pairdistance_Sp_Sk_min", min_samplesSize, 
                           "_fullLengthLTRs_alpha_beta.txt"), 
                    row.names = F, col.names = F, 
                    quote = F, sep = "\t", append = T)
    }
  }
}

distance_table_byGroup<-read.table("pairdistance_Sp_Sk_min3_fullLengthLTRs_alpha_beta.txt", F)
names(distance_table_byGroup)<-c("group", "mean_Dis", "var_Dis", "cluster", 
                                 "cladeID", "N_Sp", "N_Sk")

unique(distance_table_byGroup$cladeID)
  
distance_table_byGroup %>% 
  mutate(group=factor(group, levels=c("Sp_Sp", "Sk_Sk", "Sp_Sk"))) %>% 
  ggplot(aes(mean_Dis, fill=group)) +
  geom_histogram(aes(y=..density..), position="dodge")+
  geom_vline(data=distance_table_byGroup %>% 
               filter(cladeID %in% c("beta")) %>% 
               group_by(group) %>% 
               summarise(mean_distance_Group=mean(mean_Dis)) %>% 
               mutate(group=factor(group, levels=c("Sp_Sp", "Sk_Sk", "Sp_Sk"))), 
             aes(xintercept=mean_distance_Group, color=group),
             linetype="dashed") +
  scale_fill_manual(values=c("darkred", "steelblue", "orange")) +
  scale_color_manual(values=c("darkred", "steelblue", "orange")) +
  labs(x="Mean K2P distance", y="Density") + 
  theme_classic() +
  theme(axis.text.x = element_text(size=14, colour="black"), 
        axis.text.y = element_text(size=14, colour="black")) #+
# ggsave("05_pairdistance_Sp_Sk_min3_fullLengthLTRs_alpha_beta.png", width = 6, height =4)
# ggsave("05_pairdistance_Sp_Sk_min3_fullLengthLTRs_alpha_beta.svg", width = 6, height =4, dpi = 450)



## distance without differentiating between Clade groups:

min_samplesSize<-10
for (cluster_used in unique(ordered_table_seq2_withClades4$cluster)){
    tem<-ordered_table_seq2_withClades4 %>% 
      filter(cluster==cluster_used)
    tem_sp<-tem %>% 
      filter(ancHaplo_mean==0) %>% 
      select(seq_names) %>% 
      unlist() %>% 
      as.vector()
    tem_sk<-tem %>% 
      filter(ancHaplo_mean==1) %>% 
      select(seq_names) %>% 
      unlist() %>% 
      as.vector()
    if(length(c(tem_sp, tem_sk))>min_samplesSize){
      distance_table<-ape::dist.dna(alignment_samples[c(tem_sp, tem_sk)], 
                                    model = "K80", pairwise.deletion = TRUE, 
                                    as.matrix=TRUE)
      distance_table[lower.tri(distance_table, diag = TRUE)] <- NA
      data.frame(sample1=rownames(distance_table)[row(distance_table)], 
                 sample2=colnames(distance_table)[col(distance_table)],
                 distance=c(distance_table)) %>% 
        filter(!(is.na(distance))) %>%
        mutate(group=factor(ifelse(((sample1 %in% tem_sp) & 
                                      (sample2 %in% tem_sp)), "Sp_Sp", 
                                   ifelse(((sample1 %in% tem_sk) & 
                                             (sample2 %in% tem_sk)), "Sk_Sk", 
                                          "Sp_Sk")), 
                            levels=c("Sp_Sp", "Sk_Sk", "Sp_Sk"))) %>% 
        group_by(group) %>%
        summarise(mean_Dis=mean(distance), 
                  var_Dis=var(distance)) %>%
        mutate(cluster=cluster_used, 
               N_Sp=length(tem_sp), 
               N_Sk=length(tem_sk)) %>% 
        write.table(paste0("pairdistanceWithoutClade_Sp_Sk_min", 
                           min_samplesSize, ".txt"), 
                    row.names = F, col.names = F, 
                    quote = F, sep = "\t", append = T)
    }
}


distanceWithoutClade_table_byGroup<-read.table("pairdistanceWithoutClade_Sp_Sk_min10.txt", F)
names(distanceWithoutClade_table_byGroup)<-c("group", "mean_Dis", "var_Dis", "cluster", 
                                 "N_Sp", "N_Sk")


head(distanceWithoutClade_table_byGroup)

distanceWithoutClade_table_byGroup %>% 
  mutate(group=factor(group, levels=c("Sp_Sp", "Sk_Sk", "Sp_Sk"))) %>% 
  ggplot(aes(mean_Dis, fill=group)) +
  geom_histogram(aes(y=..density..), position="dodge")+
  geom_vline(data=distanceWithoutClade_table_byGroup %>% 
               group_by(group) %>% 
               summarise(mean_distance_Group=mean(mean_Dis)) %>% 
               mutate(group=factor(group, levels=c("Sp_Sp", "Sk_Sk", "Sp_Sk"))), 
             aes(xintercept=mean_distance_Group, color=group),
             linetype="dashed") +
  scale_fill_manual(values=c("darkred", "steelblue", "orange")) +
  scale_color_manual(values=c("darkred", "steelblue", "orange")) +
  labs(x="Mean K2P distance", y="Density") + 
  theme_classic() +
  theme(axis.text.x = element_text(size=14, colour="black"), 
        axis.text.y = element_text(size=14, colour="black"))#+
  # ggsave("05_pairdistanceWithoutClade_Sp_Sk_min10.png", width = 6, height =4)


# Distribution of sub-alpha sequences
# in the new version cluster 504 is actually 608

ordered_table_seq2_withClades4 %>% head()

# ordered_table_seq2_withClades4 %>%
#   filter(seq_names =="154_Pomberef_II_4414174_0_4881_2_allSolo_ltr")


subalpha_seq_names_cluster514<-ordered_table_seq2_withClades4  %>%
  # filter(cluster==514) %>%
  filter(cluster==608) %>%
  arrange(sample) %>%
  select(seq_names) %>% 
  unlist() %>%
  as.vector()
  
alignment_samples<-ape::read.FASTA("alig_ltr_allSolo_ltr_minLen100_plusRefseq.fasta", type = "DNA")

alignment_samples[c(subalpha_seq_names)]

ape::write.FASTA(alignment_samples[c(subalpha_seq_names_cluster514)], "alig_ltr_allSolo_ltr_minLen100_plusRefseq_subalpha_cluster514.fasta")

subalpha_seq_names_others<-ordered_table_seq2_withClades4 %>% 
  filter(cladeID=="subalpha") %>%
  filter(!(seq_names %in% subalpha_seq_names_cluster514)) %>%
  select(seq_names) %>% 
  unlist() %>%
  as.vector()

ape::write.FASTA(alignment_samples[c(subalpha_seq_names_others)], "alig_ltr_allSolo_ltr_minLen100_plusRefseq_Othersubalpha_NoCluster514.fasta")


alpha_seq_names_others<-ordered_table_seq2_withClades4 %>% 
  filter(cladeID=="alpha") %>%
  filter(!(seq_names %in% subalpha_seq_names_cluster514)) %>%
  select(seq_names) %>% 
  unlist() %>%
  as.vector()

ape::write.FASTA(alignment_samples[c(alpha_seq_names_others)], "alig_ltr_allSolo_ltr_minLen100_plusRefseq_Otheralpha_NoCluster514.fasta")


beta_seq_names_others<-ordered_table_seq2_withClades4 %>% 
  filter(cladeID=="beta") %>%
  # filter(cluster!=514) %>%
  filter(cluster!=608) %>%
  select(seq_names) %>% 
  unlist() %>%
  as.vector()

ape::write.FASTA(alignment_samples[c(beta_seq_names_others)], "alig_ltr_allSolo_ltr_minLen100_plusRefseq_Otherbeta_NoCluster514.fasta")

# reference beta sequence:

library(Biostrings)
solo_ltr_alig = readDNAStringSet("alig_ltr_allSolo_ltr_minLen100_plusRefseq.fasta")

solo_ltr_alig_beta<-solo_ltr_alig[names(solo_ltr_alig) %in% c(beta_seq_names_others)]

table_consensusHaplotype_beta<-as.matrix(solo_ltr_alig_beta) %>%
  data.frame() %>%
  mutate(seq_name=row.names(as.matrix(solo_ltr_alig_beta))) %>%
  gather(key = "pos", "genotype", -seq_name) %>%
  group_by(pos, genotype) %>%
  summarise(N_seq=n()) %>% 
  separate(pos,into = c("tem", "posN"), sep = "X") %>%
  mutate(pos=as.numeric(posN)) %>% 
  select(-tem, -posN) %>%
  # ggplot(aes(factor(pos), N_seq, fill=genotype)) +
  # geom_bar(stat = "identity")
  group_by(pos) %>%
  summarise(total=sum(N_seq), 
            max_value1=sort(N_seq, decreasing=T)[1], 
            max_value2=ifelse(max_value1==total, 0, sort(N_seq, decreasing=T)[2]), 
            fq1=max_value1/total, 
            fq2=max_value2/total, 
            genotype1=genotype[N_seq==max_value1]) 
  # ggplot(aes(fq1, fq2)) +
  # geom_point()
  # filter(fq1<0.8)


solo_ltr_alig_alpha<-solo_ltr_alig[names(solo_ltr_alig) %in% c(alpha_seq_names_others)]

table_consensusHaplotype_alpha<-as.matrix(solo_ltr_alig_alpha) %>%
  data.frame() %>%
  mutate(seq_name=row.names(as.matrix(solo_ltr_alig_alpha))) %>%
  gather(key = "pos", "genotype", -seq_name) %>%
  group_by(pos, genotype) %>%
  summarise(N_seq=n()) %>% 
  separate(pos,into = c("tem", "posN"), sep = "X") %>%
  mutate(pos=as.numeric(posN)) %>% 
  select(-tem, -posN) %>%
  # ggplot(aes(factor(pos), N_seq, fill=genotype)) +
  # geom_bar(stat = "identity")
  group_by(pos) %>%
  summarise(total=sum(N_seq), 
            max_value1=sort(N_seq, decreasing=T)[1], 
            max_value2=ifelse(max_value1==total, 0, sort(N_seq, decreasing=T)[2]), 
            fq1=max_value1/total, 
            fq2=max_value2/total, 
            genotype1=genotype[N_seq==max_value1]) 
  # ggplot(aes(fq1, fq2)) +
  # geom_point()
  # filter(fq1<0.8)

ambiguous_pos<-rbind(table_consensusHaplotype_alpha, table_consensusHaplotype_beta) %>%
  filter(fq1<0.8) %>%
  pull(pos) %>%
  unique()

ambiguous_pos  


table_consensusHaplotype_alpha %>%
  mutate(alpha_genotype=genotype1) %>%
  select(pos, alpha_genotype) %>%
  left_join(table_consensusHaplotype_beta %>%
              mutate(beta_genotype=genotype1) %>%
              select(pos, beta_genotype), by="pos") %>%
  mutate(var_sites=factor(ifelse(alpha_genotype!=beta_genotype, 1, 0))) %>%
  filter(!(alpha_genotype=="-" & beta_genotype=="-")) %>%
  arrange_(pos) %>%
  mutate(new_pos=1:n()) %>%
  #filter(var_sites==1) %>%
  select(new_pos, var_sites, alpha_genotype, beta_genotype) %>%
  gather(key = "TF", value = "Genotype", c(-new_pos, -var_sites)) %>% 
  filter(Genotype!="-") %>%
  mutate(labels=ifelse(var_sites==1, "o", "")) %>% 
  # filter(!(new_pos %in% ambiguous_pos)) %>%
  ggplot(aes(x=new_pos, y=TF))+
  geom_tile(aes(fill=Genotype))+
  geom_text(aes(label= labels))+ 
  scale_x_continuous(breaks = seq(0,400,50), limits=c(0,356))+
  xlab("Position") + 
  ylab("Haplotype") + 
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.title = element_text(colour="black", size=16), 
        axis.text.x = element_text(colour="black", size=14), 
        axis.text.y = element_text(colour="black", size=14)) #+
  # ggsave("07_consensus_alpha_beta_Haplotypes.svg", width = 18, height = 2, dpi = 450)
  # ggsave("07_consensus_alpha_beta_Haplotypes.png", width = 18, height = 2, dpi = 450)



table_consensusHaplotype_both<-table_consensusHaplotype_alpha %>%
  mutate(alpha_genotype=genotype1) %>%
  select(pos, alpha_genotype) %>%
  left_join(table_consensusHaplotype_beta %>%
              mutate(beta_genotype=genotype1) %>%
              select(pos, beta_genotype), by="pos") %>%
  mutate(var_sites=factor(ifelse(alpha_genotype!=beta_genotype, 1, 0)))

head(table_consensusHaplotype_both)

    
solo_ltr_alig_alphaplusbeta<-solo_ltr_alig[names(solo_ltr_alig) %in% c(beta_seq_names_others, alpha_seq_names_others)]

table_genotypes_solo_ltr_alig_alphaplusbeta<-as.matrix(solo_ltr_alig_alphaplusbeta) %>%
  data.frame() %>%
  mutate(seq_name=row.names(as.matrix(solo_ltr_alig_alphaplusbeta))) %>%
  gather(key = "pos", "genotype", -seq_name) %>%
  separate(pos,into = c("tem", "posN"), sep = "X") %>%
  mutate(pos=as.numeric(posN)) %>% 
  select(-tem, -posN) 


table_genotypes_solo_ltr_alig_alphaplusbeta %>%
  left_join(table_consensusHaplotype_both, by="pos") %>%
  filter(var_sites==1) %>%
  mutate(haplotypeGroup=ifelse(genotype==alpha_genotype, "alpha", ifelse(genotype==beta_genotype, "beta", "other"))) %>%
  filter(!(pos %in% ambiguous_pos)) %>%
  mutate(seq_name_order=factor(seq_name, levels=c(beta_seq_names_others, alpha_seq_names_others))) %>%
  ggplot(aes(x=factor(pos), y=seq_name_order))+
  geom_tile(aes(fill=haplotypeGroup))+
  #scale_x_continuous(breaks = seq(0,400,50), limits=c(0,356))+
  #xlab("Position") + 
  #ylab("Haplotype") + 
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.title = element_text(colour="black", size=16), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=14, colour="black"), 
        axis.text.y = element_blank()) 


solo_ltr_alig_alphaplusbetaplussubalpha<-solo_ltr_alig[names(solo_ltr_alig) %in% c(beta_seq_names_others, alpha_seq_names_others, subalpha_seq_names_cluster514, subalpha_seq_names_others)]

table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha<-as.matrix(solo_ltr_alig_alphaplusbetaplussubalpha) %>%
  data.frame() %>%
  mutate(seq_name=row.names(as.matrix(solo_ltr_alig_alphaplusbetaplussubalpha))) %>%
  gather(key = "pos", "genotype", -seq_name) %>%
  separate(pos,into = c("tem", "posN"), sep = "X") %>%
  mutate(pos=as.numeric(posN)) %>% 
  select(-tem, -posN) 


group_table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha<-table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha %>%
  left_join(table_consensusHaplotype_both, by="pos") %>%
  filter(var_sites==1) %>%
  filter(!(pos %in% ambiguous_pos)) %>%
  mutate(haplotypeGroup=ifelse(genotype==alpha_genotype, "alpha", ifelse(genotype==beta_genotype, "beta", "other"))) %>%
  group_by(seq_name) %>%
  summarise(N_alpha_var=sum(haplotypeGroup=="alpha"), 
            N_beta_var=sum(haplotypeGroup=="beta"), 
            N_other_var=sum(haplotypeGroup=="other"),
            total=n()) %>%
  mutate(groupHaplotype=ifelse((N_alpha_var>2 & N_beta_var<3), "alpha", 
                               ifelse((N_alpha_var<3 & N_beta_var>2), "beta", 
                                      ifelse((N_alpha_var>N_beta_var), "admixed_alpha", "admixed_beta")))) %>%
  select(seq_name, groupHaplotype) 


table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha %>% 
  left_join(table_consensusHaplotype_both, by="pos") %>% 
  filter(var_sites==1) %>% 
  mutate(haplotypeGroup=ifelse(genotype==alpha_genotype, "alpha", ifelse(genotype==beta_genotype, "beta", "other"))) %>%
  filter(!(pos %in% ambiguous_pos)) %>% 
  left_join(group_table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha, by="seq_name") %>% 
  mutate(seq_name_order=factor(seq_name, levels=c(beta_seq_names_others, 
                                                  subalpha_seq_names_cluster514, 
                                                  subalpha_seq_names_others, 
                                                  alpha_seq_names_others))) %>% 
  #filter(genotype!="-") %>%
  #filter(seq_name %in% c(subalpha_seq_names_cluster514, subalpha_seq_names_others)) %>%
  #filter(groupHaplotype %in% c("admixed_alpha", "admixed_beta")) %>%
  ggplot(aes(x=factor(pos), y=seq_name_order))+
  geom_tile(aes(fill=haplotypeGroup))+
  # scale_fill_manual(values=c("#7B241C", "#1A5276", "#F39C12")) +
  scale_fill_manual(values=c("#E74C3C", "#3498DB", "#F39C12")) +
  #scale_x_continuous(breaks = seq(0,1000,50)) +
  #scale_x_continuous(breaks = seq(0,400,50), limits=c(0,356))+
  xlab("Variant site") +
  ylab("Sequence") +
  facet_grid(groupHaplotype ~ ., scale="free", space="free") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.title = element_text(colour="black", size=16), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=12, colour="black"), 
        axis.text.y = element_blank()) #+
  # ggsave("07_allSubalpha_plusAlphaBetaHaplotypes.svg", width = 12, height = 8, dpi = 450)
  # ggsave("07_allSubalpha_plusAlphaBetaHaplotypes.png", width = 12, height = 8, dpi = 450)
  # ggsave("07_allSubalpha.svg", width = 9, height = 4, dpi = 450)
  # ggsave("07_allSubalpha.png", width = 9, height = 4, dpi = 450)
  # ggsave("07_allSubalpha_VS.svg", width = 12, height = 6, dpi = 450)
  # ggsave("07_allSubalpha_VS.png", width = 12, height = 6, dpi = 450)



order_subalpha_514<-table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha %>%
  left_join(table_consensusHaplotype_both, by="pos") %>%
  filter(var_sites==1) %>%
  mutate(haplotypeGroup=ifelse(genotype==alpha_genotype, "alpha", ifelse(genotype==beta_genotype, "beta", "other"))) %>%
  filter(!(pos %in% ambiguous_pos)) %>%
  left_join(group_table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha, by="seq_name") %>%
  filter(seq_name %in% c(subalpha_seq_names_cluster514)) %>%
  group_by(seq_name) %>%
  summarise(total_variants=n(), 
            N_other_var=sum(haplotypeGroup=="other")) %>%
  mutate(prop_other=N_other_var/total_variants) %>%
  arrange(prop_other) %>%
  pull(seq_name) 

library(tidyverse)
library(reshape2)
table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha %>%
  left_join(table_consensusHaplotype_both, by="pos") %>%
  filter(var_sites==1) %>%
  mutate(haplotypeGroup=ifelse(genotype==alpha_genotype, "alpha", ifelse(genotype==beta_genotype, "beta", "other"))) %>%
  filter(!(pos %in% ambiguous_pos)) %>%
  left_join(group_table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha, by="seq_name") %>%
  filter(seq_name %in% c(subalpha_seq_names_cluster514)) %>%
  mutate(seq_name_order=factor(seq_name, levels=order_subalpha_514)) %>%
  left_join(ordered_table_seq2_withClades4 %>% 
              mutate(seq_name=seq_names) %>%
              select(seq_name, sample, seq_series_ed), by="seq_name") %>% 
  #filter(genotype!="-") %>%
  ggplot(aes(x=factor(pos), y=seq_name_order))+
  geom_tile(aes(fill=haplotypeGroup))+
  # scale_fill_manual(values=c("#7B241C", "#1A5276", "#F39C12")) +
  scale_fill_manual(values=c("#E74C3C", "#3498DB", "#F39C12")) +
  #scale_x_continuous(breaks = seq(0,1000,50)) +
  #scale_x_continuous(breaks = seq(0,400,50), limits=c(0,356))+
  xlab("Variant site - TF2-11 Locus") +
  ylab("Haplotype") +
  facet_grid(sample ~ ., scale="free", switch = "y") +
  # facet_grid(sample ~ seq_series_ed, scale="free") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.title = element_text(colour="black", size=16), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_blank(), 
        strip.text.y.left = element_text(angle = 0, size=12, colour="black")) #+
  # ggsave("07_Subalpha_514.svg", width = 20, height = 11, dpi = 450)
  # ggsave("07_Subalpha_514.png", width = 20, height = 11, dpi = 450)


# Identify Subalpha/subBeta haplotypes:

subalpha_allAdmixedseq<-table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha %>%
  left_join(table_consensusHaplotype_both, by="pos") %>%
  filter(var_sites==1) %>%
  mutate(haplotypeGroup=ifelse(genotype==alpha_genotype, "alpha", ifelse(genotype==beta_genotype, "beta", "other"))) %>%
  filter(!(pos %in% ambiguous_pos)) %>%
  left_join(group_table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha, by="seq_name") %>%
  #filter(seq_name %in% c(subalpha_seq_names_cluster514, subalpha_seq_names_others)) %>%
  filter(groupHaplotype %in% c("admixed_alpha", "admixed_beta"))
  
subalpha_haplotype<-c()
list_unique_haplotypes<-c()

list_sequences<-unique(subalpha_allAdmixedseq$seq_name)
for (sequence_used in seq(2, length(list_sequences))) {
  print(sequence_used)
  if (sequence_used==1) {
    subalpha_haplotype<-c(sequence_used)
    list_unique_haplotypes<-c(sequence_used)
  } else {
    sample_name<-list_sequences[sequence_used]
    tem_new_sample_table<-subalpha_allAdmixedseq %>%
      filter(seq_name==sample_name) %>%
      select(pos, haplotypeGroup)
    min_diff<-200
    identified_haplotype<-c()
    for (ref_sequence_used in list_unique_haplotypes){
      ref_sample_name<-list_sequences[ref_sequence_used]
      tem_ref_sample_table<-subalpha_allAdmixedseq %>%
        filter(seq_name==ref_sample_name) %>%
        mutate(ref_haplotypeGroup=haplotypeGroup) %>%
        select(pos, ref_haplotypeGroup)
      N_diff<-tem_new_sample_table %>%
        left_join(tem_ref_sample_table, by="pos") %>%
        mutate(diff_genotype=(ref_haplotypeGroup!=haplotypeGroup)*1) %>%
        ungroup() %>%
        summarise(total_diff=sum(diff_genotype)) %>%
        pull(total_diff)
      if (N_diff<min_diff) {
        min_diff<-N_diff
        identified_haplotype<-ref_sequence_used
      }
    }
    if (min_diff<5) {
      subalpha_haplotype<-c(subalpha_haplotype, identified_haplotype)
    } else {
      subalpha_haplotype<-c(subalpha_haplotype, sequence_used)
      list_unique_haplotypes<-c(list_unique_haplotypes, sequence_used)
    }
  }
}

subalpha_haplotype<-c(2,subalpha_haplotype)
length(list_unique_haplotypes)
list_unique_haplotypes

subalpha_haplotype_or<-data.frame(list_sequences, subalpha_haplotype) %>%
  mutate(seq_names=list_sequences) %>%
  select(seq_names, subalpha_haplotype) %>%
  left_join(ordered_table_seq2_withClades4, by="seq_names") %>%
  select(seq_names, subalpha_haplotype, sample, cluster, chr, pos, len_seq) %>%
  group_by(subalpha_haplotype) %>%
  summarise(N_seq=n(), 
            N_samples=length(unique(sample)), 
            N_clusters=length(unique(cluster))) %>%
  arrange(-N_clusters) %>%
  pull(subalpha_haplotype)



table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha %>%
  left_join(table_consensusHaplotype_both, by="pos") %>%
  filter(var_sites==1) %>%
  mutate(haplotypeGroup=ifelse(genotype==alpha_genotype, "alpha", ifelse(genotype==beta_genotype, "beta", "other"))) %>%
  filter(!(pos %in% ambiguous_pos)) %>%
  left_join(group_table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha, by="seq_name") %>%
  #filter(seq_name %in% c(subalpha_seq_names_cluster514, subalpha_seq_names_others)) %>%
  filter(groupHaplotype %in% c("admixed_alpha", "admixed_beta")) %>%
  left_join(data.frame(list_sequences, subalpha_haplotype) %>%
              mutate(seq_name=list_sequences) %>%
              select(seq_name, subalpha_haplotype), by="seq_name") %>%
  mutate(subalpha_haplotype_order=factor(subalpha_haplotype, levels=subalpha_haplotype_or)) %>%
  ggplot(aes(x=factor(pos), y=seq_name))+
  geom_tile(aes(fill=haplotypeGroup))+
  # scale_fill_manual(values=c("#7B241C", "#1A5276", "#F39C12")) +
  scale_fill_manual(values=c("#E74C3C", "#3498DB", "#F39C12")) +
  #scale_x_continuous(breaks = seq(0,1000,50)) +
  #scale_x_continuous(breaks = seq(0,400,50), limits=c(0,356))+
  xlab("Variant site") +
  ylab("Sequence") +
  facet_grid(subalpha_haplotype_order ~ ., scale="free", space="free", switch = "y") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.title = element_text(colour="black", size=16), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_blank(), 
        strip.text.y.left = element_blank()) #+
  # ggsave("07_allSubalpha_byHaplotype_4SNPsmax.png", width = 12, height = 6, dpi = 450)
  # ggsave("07_allSubalpha_byHaplotype_4SNPsmax.svg", width = 12, height = 6, dpi = 450)
  # ggsave("07_allSubalphaBeta_byHaplotype_4SNPsmax.png", width = 12, height = 6, dpi = 450)


list_sequences[list_unique_haplotypes]

table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha %>%
  left_join(table_consensusHaplotype_both, by="pos") %>%
  filter(var_sites==1) %>%
  mutate(haplotypeGroup=ifelse(genotype==alpha_genotype, "alpha", ifelse(genotype==beta_genotype, "beta", "other"))) %>%
  filter(!(pos %in% ambiguous_pos)) %>%
  left_join(group_table_genotypes_solo_ltr_alig_alphaplusbetaplussubalpha, by="seq_name") %>%
  #filter(seq_name %in% c(subalpha_seq_names_cluster514, subalpha_seq_names_others)) %>%
  # filter(seq_name %in% c(subalpha_seq_names_cluster514)) %>%
  filter(groupHaplotype %in% c("admixed_alpha", "admixed_beta")) %>%
  left_join(data.frame(list_sequences, subalpha_haplotype) %>%
              mutate(seq_name=list_sequences) %>%
              select(seq_name, subalpha_haplotype), by="seq_name") %>%
  # pull(subalpha_haplotype) %>%
  # unique()
  mutate(subalpha_haplotype_order=factor(subalpha_haplotype, levels=subalpha_haplotype_or)) %>%
  filter(seq_name %in% list_sequences[list_unique_haplotypes]) %>%
  left_join(data.frame(list_sequences, subalpha_haplotype) %>%
              mutate(seq_names=list_sequences) %>%
              select(seq_names, subalpha_haplotype) %>%
              left_join(ordered_table_seq2_withClades4, by="seq_names") %>%
              select(seq_names, subalpha_haplotype, sample, cluster, chr, pos, len_seq) %>%
              group_by(subalpha_haplotype) %>%
              summarise(N_seq=n(), 
                        N_samples=length(unique(sample)), 
                        N_clusters=length(unique(cluster))) , by="subalpha_haplotype") %>%
  mutate(seq_label=paste0(subalpha_haplotype, "_", N_seq, "_", N_samples, "_", N_clusters)) %>%
  ggplot(aes(x=factor(pos), y=seq_label))+
  geom_tile(aes(fill=haplotypeGroup))+
  # scale_fill_manual(values=c("#7B241C", "#1A5276", "#F39C12")) +
  scale_fill_manual(values=c("#E74C3C", "#3498DB", "#F39C12")) +
  #scale_x_continuous(breaks = seq(0,1000,50)) +
  #scale_x_continuous(breaks = seq(0,400,50), limits=c(0,356))+
  xlab("Variant site") +
  ylab("Sequence") +
  facet_grid(subalpha_haplotype_order ~ ., scale="free", space="free", switch = "y") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.title = element_text(colour="black", size=16), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=10, colour="black"), 
        axis.text.y = element_text(colour="black", size=12), 
        strip.text.y.left = element_blank()) #+
  # ggsave("07_allSubalpha_commonHaplotypes_summary_4SNPsmax.png", width = 12, height = 4, dpi = 450)
  # ggsave("07_allSubalpha_commonHaplotypes_summary_4SNPsmax.svg", width = 12, height = 4, dpi = 450)
  # ggsave("07_allSubalphaBeta_commonHaplotypes_summary_4SNPsmax.png", width = 12, height = 6, dpi = 450)
  # ggsave("07_allSubalphaBeta_commonHaplotypes_summary_4SNPsmax.svg", width = 12, height = 6, dpi = 450)



data.frame(list_sequences, subalpha_haplotype) %>%
  mutate(seq_names=list_sequences) %>%
  select(seq_names, subalpha_haplotype) %>%
  left_join(ordered_table_seq2_withClades4, by="seq_names") %>% 
  # filter(subalpha_haplotype %in% c(1, 269, 309, 451, 490)) %>%
  filter(subalpha_haplotype %in% c(2, 659)) %>%
  mutate(sample_order=factor(sample, levels=order_samples_pro_AncPop)) %>%
  # filter(len_seq>1500) %>%
  mutate(full_LTR_1.5kb=factor((len_seq>1500)*1)) %>%
  ggplot(aes(sample_order, fill=full_LTR_1.5kb))+
  geom_bar(stat="count")+
  scale_fill_manual(values=c("grey80", "black"))+
  #scale_fill_manual(values=wes_palette(30, name = "Darjeeling1", type = "discrete"))+
  geom_hline(yintercept=0, color = "black", alpha=0.6)+
  facet_grid(subalpha_haplotype ~ ., scale="free")+
  labs(x="Sample", y="Number of Seq") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black"), 
        strip.text.y = element_text(size=10, angle = 0, colour="black"),
        legend.position = "right") #+
  # ggsave("07_allSubalpha_commonHaplotypes_summary_4SNPsmax_perSample.png", width = 8, height = 4, dpi = 450)
  # ggsave("07_allSubalpha_commonHaplotypes_summary_4SNPsmax_perSample.svg", width = 8, height = 4, dpi = 450)
  # ggsave("07_allSubalphaBeta_commonHaplotypes_summary_4SNPsmax_perSample.png", width = 8, height = 3, dpi = 450)
  # ggsave("07_allSubalphaBeta_commonHaplotypes_summary_4SNPsmax_perSample.svg", width = 8, height = 3, dpi = 450)




order_samples_pro_AncPop
  
### SFS considering family and direccion of sequences:


ordered_table_seq2_withClades3 %>% head()
hist(ordered_table_seq2_withClades3$len_seq)
ordered_table_seq2_withClades3 %>%
  filter(len_seq<700)

unique_samples<-c("Pomberef",
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



ordered_table_seq2_withClades3_SFS<-ordered_table_seq2_withClades3 %>% 
  filter(sample %in% unique_samples) %>%
  filter(!(is.na(ancHaplo))) %>%
  #filter(cluster==100) %>%
  mutate(reverse_seq_ed=if_else(reverse_seq_ed==0,"Fow", "Rev")) %>% 
  mutate(variantID=paste(cluster, sample, cladeID, reverse_seq_ed, sep="__")) %>%
  arrange(variantID) %>%
  group_by(variantID) %>% 
  mutate(Variant_SeriesPerFamily = row_number()) %>%
  mutate(variantID=paste(cluster, cladeID, reverse_seq_ed, Variant_SeriesPerFamily, 
                         sep="__")) %>%
  mutate(genotype=1)

added_samples<-c()
added_cluster<-c()
added_variantID<-c()
added_chr<-c()
added_ancHaplo<-c()

for (tem_variantID in unique(ordered_table_seq2_withClades3_SFS$variantID)){
  # tem_variantID<-"281__other__Fow__1"
  tem_data<-ordered_table_seq2_withClades3_SFS %>%
    filter(variantID==tem_variantID)
  tem_chr<-as.vector(tem_data$chr[1])
  tem_cluster<-as.vector(tem_data$cluster[1])
  print(tem_cluster)
  tem_start<-min(tem_data$pos)
  missing_samples<-as.vector(unique(ordered_table_seq2_withClades3_SFS$sample)[!(unique(ordered_table_seq2_withClades3_SFS$sample) %in% tem_data$sample)])
  for (tem_sample in missing_samples){
    # tem_sample<-"Pomberef"
    # ancestralhap %>% 
    #   filter(sample==tem_sample) %>% 
    #   filter(chromosome_name==tem_chr) %>% 
    #   arrange(start_pos) %>% 
    #   filter(start_pos<=(tem_start-50)) %>%
    #   filter(end_ed>(tem_start+50))
    start_point<-ancestralhap[ancestralhap$sample==tem_sample & ancestralhap$chromosome_name==tem_chr & ancestralhap$start_pos<=(tem_start-50) & ancestralhap$end_ed>(tem_start+50),"Nor_PC1_pol"]
    if (length(start_point)==1){
      added_ancHaplo<-c(added_ancHaplo, start_point)
    } else if (length(start_point)==0){
      added_ancHaplo<-c(added_ancHaplo, NA)
    } else if (length(start_point)>1){
      added_ancHaplo<-c(added_ancHaplo, 2)
    }
    added_samples<-c(added_samples, tem_sample)
    added_cluster<-c(added_cluster, tem_cluster)
    added_chr<-c(added_chr, tem_chr)
    added_variantID<-c(added_variantID, tem_variantID)
  }
}

sfs_table<-ordered_table_seq2_withClades3_SFS %>%
  ungroup() %>%
  select(sample, chr, cluster, ancHaplo, variantID, genotype) %>%
  rbind(data.frame(sample=added_samples, 
                   chr=added_chr, 
                   cluster=added_cluster, 
                   ancHaplo=added_ancHaplo, 
                   variantID=added_variantID,
                   genotype=0)) %>%
  ### Checking step: checking for inconsistent ancestral background
  # group_by(chr, cluster, sample, variantID, genotype) %>% 
  # summarise(NumSeq_perSample=n(), 
  #           sp_seq=sum(ancHaplo==0, na.rm = TRUE), 
  #           sk_seq=sum(ancHaplo==1, na.rm = TRUE), 
  #           na_seq=sum(is.na(ancHaplo), na.rm = TRUE)) %>%
  group_by(chr, cluster, variantID) %>% 
  summarise(num_samples=length(unique(sample)),
            num_seq=n(), 
            num_seq_sp=sum(ancHaplo==0),
            num_seq_sk=sum(ancHaplo==1),
            num_Gen1_sp=sum(ancHaplo==0 & genotype==1),
            num_Gen1_sk=sum(ancHaplo==1 & genotype==1),
            fq_sp=num_Gen1_sp/num_seq_sp, 
            fq_sk=num_Gen1_sk/num_seq_sk) %>% 
  select(chr, cluster, variantID, num_seq_sp, num_seq_sk, fq_sp, fq_sk)
  
head(sfs_table)


min_sample_size<-4
fq_sp<-c()
fq_sk<-c()
final_counts<-c()
for (tem_fq_sp in seq(0,1,0.05)) {
  for (tem_fq_sk in seq(0,1,0.05)) {
    counts<-sfs_table %>%
      filter(fq_sp>=tem_fq_sp & fq_sp<tem_fq_sp+0.04999) %>%
      filter(fq_sk>=tem_fq_sk & fq_sk<tem_fq_sk+0.04999) %>%
      filter(num_seq_sp>=min_sample_size & num_seq_sk>=min_sample_size) %>% 
      dim()
    fq_sp<-c(fq_sp, tem_fq_sp)
    fq_sk<-c(fq_sk, tem_fq_sk)
    final_counts<-c(final_counts, counts[1])
  }
}

sfs_table2<-data.frame(fq_sp, fq_sk, counts=final_counts) 
sfs_table2 %>%
  filter(counts!=0) %>%
  ggplot(aes(fq_sp, fq_sk)) + 
  geom_tile(aes(fill = counts)) + 
  geom_text(aes(label=ifelse(counts!=0, counts, "")), size=5) + 
  #geom_text(aes(label= (formatC(sfs_2d*100, format = "e", digits = 1)))) + 
  scale_fill_gradientn(colours = rev(terrain.colors(20)),
                       name="")+
  # #scale_fill_gradient(low = "white", high = "steelblue") +
  # #scale_fill_gradientn(colours = c("white", "blue", "red", "green"),
  # #trans = "log",
  # limits = c(0,100),
  # breaks=c(0,20,50,100)) +
  #limits = c(-7,4.6)) + 
  xlab("Allele Fq. Sp group") + 
  ylab("Allele Fq. Sk group") + 
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x=element_text(angle = 45, hjust = 1, size=17, colour="black"), 
        axis.text.y = element_text(colour="black", size=17), 
        axis.title=element_text(size=17),
        #strip.text.x = element_blank(), 
        legend.text = element_text(angle = 45, hjust = 1, colour="black", size = 14)) #+
  # ggsave("06_2dSFS_AncGroup_withLTRSeq_Family_Direction_MinSS4.png", width = 9, height = 9, dpi = 400)



min_sample_size<-4
fq_sp<-c()
fq_sk<-c()
final_counts<-c()
for (tem_fq_sp in seq(0,1,0.1)) {
  for (tem_fq_sk in seq(0,1,0.1)) {
    counts<-sfs_table %>%
      filter(fq_sp>=tem_fq_sp & fq_sp<tem_fq_sp+0.09999) %>%
      filter(fq_sk>=tem_fq_sk & fq_sk<tem_fq_sk+0.09999) %>%
      filter(num_seq_sp>=min_sample_size & num_seq_sk>=min_sample_size) %>% 
      dim()
    fq_sp<-c(fq_sp, tem_fq_sp)
    fq_sk<-c(fq_sk, tem_fq_sk)
    final_counts<-c(final_counts, counts[1])
  }
}

sfs_table2<-data.frame(fq_sp, fq_sk, counts=final_counts) 


sfs_table2 %>%
  mutate(percentage=counts*100/sum(sfs_table2$counts)) %>%
  filter(counts>2) %>%
  ggplot(aes(fq_sp, fq_sk)) + 
  geom_tile(aes(fill = percentage)) + 
  geom_text(aes(label=ifelse(counts!=0, counts, "")), size=5) +
  # geom_text(aes(label=ifelse(percentage>0.1, percentage, "")), size=5) +
  # geom_text(aes(label= (formatC(sfs_2d*100, format = "e", digits = 1)))) +
  scale_fill_gradientn(colours = rev(terrain.colors(20)),
                       name="")+
  # #scale_fill_gradient(low = "white", high = "steelblue") +
  # #scale_fill_gradientn(colours = c("white", "blue", "red", "green"),
  # #trans = "log",
  # limits = c(0,100),
  # breaks=c(0,20,50,100)) +
  #limits = c(-7,4.6)) + 
  xlab("Allele Fq. Sp group") + 
  ylab("Allele Fq. Sk group") + 
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x=element_text(angle = 45, hjust = 1, size=17, colour="black"), 
        axis.text.y = element_text(colour="black", size=17), 
        axis.title=element_text(size=17),
        #strip.text.x = element_blank(), 
        legend.text = element_text(angle = 45, hjust = 1, colour="black", size = 14)) #+
  # ggsave("06_2dSFS_AncGroup_withLTRSeq_Family_Direction_MinSS4_percentages.png", width = 9, height = 9, dpi = 400)
  # ggsave("06_2dSFS_AncGroup_withLTRSeq_Family_Direction_MinSS4_percentages.svg", width = 9, height = 9, dpi = 400)


####


# Supplementary table solo LTRs and flanking LTRs differing by family and strand direction:

ordered_table_seq2_withClades3 %>% head()
hist(ordered_table_seq2_withClades3$len_seq)
ordered_table_seq2_withClades3 %>%
  filter(len_seq<700)

unique_samples<-c("Pomberef",
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


# 
# ordered_table_seq2_withClades3_family<-
ordered_table_seq2_withClades3 %>%
  # filter(sample %in% unique_samples) %>%
  filter(!(is.na(ancHaplo))) %>%
  #filter(cluster==100) %>%
  mutate(reverse_seq_ed=if_else(reverse_seq_ed==0,"For", "Rev")) %>%
  mutate(variantID=paste(cladeID, reverse_seq_ed, sep="__")) %>%
  group_by(cluster, sample, variantID) %>%
  summarise(NumSeq=n()) %>%
  write.table("soloLTR_seq_table_FamilyDirection.txt", quote = F, row.names = F)


####













































# download.file("https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.0.5.tar.gz",
#               dest="dplyr_1.0.5.tar.gz")
# install.packages("dplyr_1.0.5.tar.gz",repos=NULL,type="source")

# remotes::install_github("YuLab-SMU/ggtree", force = TRUE)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("treeio", force = TRUE)

library("ggtree")
#tree<-read.tree("../phylogenies/phy_tf_minLen1000.treefile")
tree_ed<-treeio::read.iqtree("allSoloLTR_100.treefile")
#tree_ed<- as.phylo(tree)
fig_tree<- ggtree(tree_ed) + 
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
  #geom_tiplab() + 
  #geom_tippoint() + 
  theme_tree2() + 
  xlim(0, 0.5)
fig_tree

tree_ed@data$UFboot[is.na(tree_ed@data$UFboot)]<-0
tree_ed@data$UFboot[tree_ed@data$UFboot<90]<-0
tree_ed@data$UFboot[tree_ed@data$UFboot>=90 & tree_ed@data$UFboot<95]<-90
tree_ed@data$UFboot[tree_ed@data$UFboot>=95]<-95


# # change sequences IDs in the phylogeny:
# tem<-rev(tree_ed@phylo$tip.label)
# tem<-str_remove_all(tem, "_R_")
# tem2<-str_remove_all(tem, "_I*_[0-9]*_[0-9]*_[0-9]*")
# tem2<-str_remove_all(tem2, "_AB325691_[0-9]*_[0-9]*_[0-9]*")
# seq_ID<-str_remove_all(tem2, "_allSolo_ltr$")
# tree_ed@phylo$tip.label<-rev(seq_ID)



head(tree_ed@phylo$tip.label)
head(ordered_table_seq2)

dim(ordered_table_seq2)
length(tree_ed@phylo$tip.label)

cluster_ID<-data.frame(id=tree_ed@phylo$tip.label) %>%
  left_join(ordered_table_seq2 %>% 
          # mutate(id=seq_ID) %>%
          mutate(id=seq_names) %>%
          # filter(sample=="JB837") %>%
          select(id, cluster, chr, anc_prop, ancHaplo), by="id", all.x=T) %>% 
  mutate(anc_group=factor(ifelse(anc_prop>0.9,"Sp",ifelse(anc_prop<0.1,"Sk","Hyb")), 
                          levels=c("Sp", "Hyb", "Sk"))) %>%
  mutate(col=1) %>%
  arrange(match(id, tree_ed@phylo$tip.label))

head(cluster_ID)
row.names(cluster_ID) <- NULL
ordered_table_seq2 %>% 
  head()

head(tree_ed@phylo$tip.label)

gg <- data.frame(id=tree_ed@phylo$tip.label) %>%
  left_join(ordered_table_seq2 %>% 
              # mutate(id=seq_ID) %>%
              mutate(id=seq_names) %>%
              # filter(sample=="JB837") %>%
              select(id, cluster, chr, anc_prop, ancHaplo, len_seq), by="id", all.x=T) %>% 
  mutate(anc_group=factor(ifelse(anc_prop>0.90,"Sp",ifelse(anc_prop<0.15,"Sk","Hyb")), 
                          levels=c("Sp", "Hyb", "Sk"))) %>%
  # mutate(id_ed=factor(id, levels=c(tree_ed@phylo$tip.label))) %>%
  left_join(data.frame(id=tree_ed@phylo$tip.label, 
                       order_tree=seq(1,length(tree_ed@phylo$tip.label))), 
            by="id") %>%
  filter(len_seq<1000) %>%
  filter(chr %in% c("I", "II", "III")) %>%
  filter(anc_group %in% c("Sp", "Hyb", "Sk")) %>%
  ggplot(aes(order_tree, color=anc_group, fill=anc_group))+
  geom_histogram(aes(y=..density..), alpha=0.5, 
                 position="identity", bins=500)+
  scale_fill_manual(values=c("darkred", "orange", "steelblue")) +
  scale_colour_manual(values=c("darkred", "orange", "steelblue")) +
  facet_grid(anc_group ~ .)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        panel.grid.major = element_line(colour = "gray80"), 
        axis.text.y=element_text(size=10, colour="black"), 
        axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(), 
        # legend.position="none") 
        legend.position="top", 
        legend.title = element_blank())
  
ggsave("08_distribution_SoloLTRs_perAncGroup.png", plot = gg, width = 10, height = 4, dpi = 450)
ggsave("08_distribution_SoloLTRs_perAncGroup.svg", plot = gg, width = 10, height = 4, dpi = 450)


