# install.packages("tidyverse",  lib="/home/sergio/R/libraries_Rackham/3.4")
#rm(list=ls())
# end 3623

#library("ggplot2")
#library("tidyr")
#library("plyr")
#library("dplyr")
library("tidyverse")
library("ggpubr")
library(svglite)
library(ggplot2)


  
#
ancestralhapl<-read.table("ancestralhap_57ILL_LR.txt", T) 
ancestralhapl<-ancestralhapl[-(grep("ILL", ancestralhapl$sample)),]

pro_AncPop <- ancestralhapl %>% 
  filter(Nor_PC1_pol_sim!=0.5) %>% 
  group_by(sample) %>% 
  summarise(proportion_sk=sum(Nor_PC1_pol_sim)/n()) 

order_samples_pro_AncPop <- pro_AncPop %>%
  mutate(prop_sk_ed=ifelse(sample=="JB22", 0.0006501950,proportion_sk)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB879", 0.0007,prop_sk_ed)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB873", 0.4533954726,prop_sk_ed)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB1206", 0.8513139695,prop_sk_ed)) %>%    
  arrange(prop_sk_ed) %>% 
  ungroup() %>% 
  select(sample) %>% 
  unlist() %>% 
  as.vector()

pro_AncPop<-pro_AncPop %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB22", 0.0006501950,proportion_sk)) %>% 
  mutate(prop_sk_ed=ifelse(sample=="JB879", 0.0007,prop_sk_ed)) %>% 
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



pro_AncPop2<-pro_AncPop %>% 
  filter(sample=="JB22") %>% 
  mutate(sample="Pomberef") 

pro_AncPop2<-rbind(pro_AncPop, pro_AncPop2)

pro_AncPop2 %>% 
  mutate(sample2=factor(sample, levels=order_samples_pro_AncPop)) %>%
  ggplot(aes(sample2, 1-anc_prop))+
  geom_bar(stat="identity") +
  labs(x="Sample", y="Proportion (Sk)") + 
  scale_y_continuous(limits = c(0,1))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black")) #+
  # ggsave("anc_proportions_LReads.png", width = 7, height = 3)
  # ggsave("anc_proportions_LReads.svg", width = 7, height = 3)

##### TF #####

samples_table_or<-read.table("all_TE_annotation_refCoor_Masked_conSeq.bed", F) 

samples_table_EDTA_missingSeq<-read.table("EDTA_all_missingSeq_TE_annotation__refCoor_Masked_conSeq.bed", F, sep = "\t") 

names(samples_table_or)<-c("sample", "or_chr",	"or_start",	"or_end",	"seq_ID",	"ref_start_chr",	"ref_start",	"dis_start",	"ref_end_chr",	"ref_end",	"dis_end",	"length_ori_seq",	"direction",	"concordant_chr",	"chr_start_ed",	"pos_start_ed",	"dis_strat_ed", "chr_end_ed",	"pos_end_ed",	"dis_end_ed")


names(samples_table_EDTA_missingSeq)<-c("sample", "or_chr",	"or_start",	"or_end",	"seq_ID",	"ref_start_chr",	"ref_start",	"dis_start",	"ref_end_chr",	"ref_end",	"dis_end",	"length_ori_seq",	"direction",	"concordant_chr",	"chr_start_ed",	"pos_start_ed",	"dis_strat_ed", "chr_end_ed",	"pos_end_ed",	"dis_end_ed")

head(samples_table_EDTA_missingSeq)

samples_table_or<-rbind(samples_table_or, 
                     samples_table_EDTA_missingSeq) %>%
  arrange(sample, chr_start_ed, pos_start_ed, dis_strat_ed)

samples_table_EDTA_missingSeq %>%
  dim()

# samples_table$length_ori_seq <- as.integer(samples_table$length_ori_seq)

samples_table_or %>%
  filter(is.na(length_ori_seq))

samples_table_or$length_ori_seq

head(samples_table_or)
dim(samples_table_or)
str(samples_table_or)

samples_table_or %>%
  filter(length_ori_seq>200) %>% 
  filter(length_ori_seq<800) %>% 
  dim()
  

samples_table_or %>% 
  filter(sample=="Pomberef") %>% 
  filter(length_ori_seq>1000) %>% 
  arrange(chr_start_ed, pos_start_ed) %>% 
  mutate(len=pos_end_ed-pos_start_ed)

samples_table_or %>% 
  mutate(len=pos_end_ed-pos_start_ed) %>% 
  filter(len>100000)

samples_table<-samples_table_or %>% 
  # mutate(len_tem=pos_end_ed-dis_strat_ed) %>%
  mutate(len_tem=pos_end_ed-pos_start_ed) %>%
  # mutate(len=ifelse(len_tem>10000, length_ori_seq, len_tem)) %>% 
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

# samples_table_or %>% head()
# ordered_table2 %>%
#   filter(chr=="II") %>%
#   filter(start>4436600 & start<4450600)

# names(samples_table)<-c("seq_ID", "sample", "chr", "start", "dis", "len_froSeq", "addID", "dir", "len", "rep")

# or_chr	or_start	or_end	seq_ID	ref_start_chr	ref_start	dis_start	ref_end_chr	ref_end	dis_end	length_ori_seq	direction	concordant_chr	chr_start_ed	pos_start_ed	dis_strat_ed

# sp_tf_loc<-read.table("location_sp_tf.txt", T)
sp_tf_loc<-read.csv("../../../LiLIn_results/second/location_sp_tf_Bowen_revised5-200830-for-Tusso.csv", T)

head(sp_tf_loc)

sp_tf_loc

ordered_table_sp_tf<-sp_tf_loc %>% 
  mutate(seq_ID="NA", dis=0, len_froSeq=end-start, 
         addID=sample, 
         sample="ref_loc", 
         dir="+", len=end-start, rep="NA") %>% 
  select(seq_ID, sample, chr, start, dis, len_froSeq, addID, dir, len, rep) 

ordered_table_sp_tf %>% 
  filter(len_froSeq<0)

### Bowen 2003 LTR sequences:
head(ordered_table_sp_tf)
# sp_bowenLTR<-read.csv("../../../LiLIn_results/second/Bowen_LTR_annotation.csv", T)
sp_bowenLTR<-read.csv("../../../LiLIn_results/second/Bowen_revised5-200830-for-Tusso.csv", T)

head(sp_bowenLTR)

ordered_table_bowen<-sp_bowenLTR %>% 
  mutate(dis=0, len_froSeq=end-start, 
         addID=clade_affiliation, 
         dir=LTR_orientation, 
         rep="NA") %>% 
  select(seq_ID, sample, chr, start, dis, len_froSeq, addID, dir, len, rep) 


head(ordered_table_bowen)

ordered_table<-rbind(samples_table, ordered_table_sp_tf, ordered_table_bowen) %>% 
  # merge(ID_table, by="sample", all.x=T) %>% 
  arrange(chr, start, sample, dis) %>% 
  filter(!(is.na(chr)))


dim(ordered_table)
head(ordered_table)



# The following look identifies the cluster for each sequence based on overlap

cluster<-c(1)
length_cluster<-0
min_start<-0

# it worked quite ok with 150
# also with 100 bp
for (i in seq(2,dim(ordered_table)[1])){
# i=2
  print(i)
  #print(ordered_table[i,])
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

ancestralhap<-ancestralhapl %>% 
  # filter(sample %in% ordered_table$sample_JB) %>%
  mutate(Nor_PC1_pol=(Nor_PC1_pol_sim==1)*1)

head(ancestralhap)


# identify ancestral background in flanking sequences
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
ordered_table2<-ordered_table2 %>% 
  left_join(pro_AncPop, by="sample")

ordered_table2$anc_prop[ordered_table2$sample=="Pomberef"]<-1

ordered_table2 %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  # filter((sample %in% no_clonal_strains)) %>%
  # select(sample) %>%
  # unique() %>%
  # dim()
  mutate(start_end=(ancHaplo==ancHaplo_end)*1) %>%
  ## Number of consistent flanquin sequences (7954 consistent, 138 inconsistent, 734 with NA):
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
  ### there are 55 clusters with at least one inconsistent sample, but in most of the clusters it is only one or two samples (46). High number are only observed in 9 clusters (more than 2 samples). 
  # filter(inconsistent_Samples>0 & inconsistent_Samples<3) %>% dim()
  # filter(inconsistent_Samples>2) %>% dim()
  filter(inconsistent_Samples!=0) %>%
  ggplot(aes(total_NSamples, inconsistent_Samples))+
  geom_point(alpha=0.2) +
  scale_y_continuous(breaks=seq(1,14,1))+
  labs(x="N. Samples per Cluster", 
       y="N. Samples with Inconsistent 
       flanking Seq") + 
  theme_classic() #+
  # ggsave("04_inconsistentFlankingSeqs.png", width = 5, height = 4)
  # ggsave("04_inconsistentFlankingSeqs_NoClonalSamples.png", width = 5, height = 4)



# consensus ancestral haplotype per cluster:

ancHaplo_mean_table<-ordered_table2 %>%
  # filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  group_by(cluster,sample) %>% 
  summarise(total=n()*2,
            total_na=sum(is.na(ancHaplo)*1),
            total_na_end=sum(is.na(ancHaplo_end)*1),
            ancHaplo_1=sum((ancHaplo==1)*1, na.rm=TRUE),
            ancHaplo_end_1=sum((ancHaplo_end==1)*1, na.rm=TRUE)) %>%
  mutate(total_1=ancHaplo_1+ancHaplo_end_1, 
         ratio_1_total=total_1/total, 
         total_na_both=total_na+total_na_end) %>%
  mutate(ancHaplo_mean=ifelse(total_na_both==total, NA, 
                              ifelse(ratio_1_total<0.5, 0, 
                              ifelse(ratio_1_total>0.5, 1, NA)))) %>% 
  select(sample, cluster, ancHaplo_mean) 



ordered_table2 <- ordered_table2 %>% 
  full_join(ancHaplo_mean_table, by=c("sample", "cluster"))

# write.table(ordered_table2, file="annotation_clusters.txt", quote = F, row.names = F)

# ordered_table2 %>%
#   filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
#   filter(chr!="AB325691") %>%
#   mutate(ancHaplo_mean_final=ifelse(ancHaplo_mean==0, "Sp",
#                                     ifelse(ancHaplo_mean==1, "Sk", NA))) %>%
#   ungroup() %>%
#   select(seq_ID, sample, chr, start, dis, len_froSeq, dir, cluster, ancHaplo_mean_final) %>%
#   arrange(cluster, sample, start) %>%
#   write.table(file="annotation_clusters_supplementary_table.txt", quote = F, row.names = F)

ordered_table2 %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  filter(chr!="AB325691") %>%
  # filter(len_froSeq>1500) %>%
  filter(len_froSeq>200) %>%
  filter(len_froSeq<600) %>%
  # filter(len_froSeq>4500) %>%
  dim()

# 8616 sequences in total
# 7182 sequences longer than 200 bp
# 5537 solo-LTRs
# 924 full-length LTRs

ordered_table2 %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  filter(chr!="AB325691") %>%
  # filter(len_froSeq>200) %>%
  # filter(len_froSeq>4500) %>%
  pull(cluster) %>%
  unique() %>%
  length()

# 713 clusters in total
# 682 clusters if only sequences longer than 200 bp are considered
# 437 clusters with full-length LTRs


ordered_table2 %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  filter(chr!="AB325691") %>%
  # filter(len_froSeq>1500) %>%
  # filter(len_froSeq>200) %>%
  # filter(len_froSeq<600) %>%
  filter(len_froSeq>4500) %>%
  group_by(sample) %>% 
  summarize(NSeq=n()) %>%
  pull(NSeq) %>%
  max()

ordered_table2 %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  filter(chr!="AB325691") %>%
  filter(ancHaplo_mean %in% c(0,1)) %>%
  pull(cluster) %>%
  unique() %>%
  length()
  head()

# ordered_table2 %>% 
#   filter(chr!="AB325691") %>% 
#   filter(seq_ID!="NA") %>%
#   filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
#   mutate(ancetralHaplotype=ifelse(ancHaplo_mean==0, "Sp", 
#                 ifelse(ancHaplo_mean==1, "Sk", NA))) %>%
#   arrange(cluster, sample, start, dis) %>%
#   write.table(file="annotation_clusters.txt", quote = F, row.names = F, sep = ",")

# ordered_table2_end_start <- ordered_table2 

dim(ordered_table2_end_start)
dim(ordered_table2 )

max(ordered_table2_end_start$cluster)
max(ordered_table2$cluster )


###

ordered_table2 %>% 
  filter(seq_ID!="NA") %>% 
  # filter(cluster==200)
  group_by(cluster) %>% 
  summarize(chr=chr[1],
            start=min(start),
            total_seq=n(), 
            total_samples=length(unique(sample)), 
            min_length=min(len_froSeq), 
            max_length=max(len_froSeq)) %>%
  filter(total_samples>6)
  #filter(min_length<2500)
  #filter(chr=="III" & min_length<2000)

# ordered_table %>% 
#   filter(seq_ID!="NA") %>% 
#   merge((ordered_table %>% 
#            filter(seq_ID!="NA") %>% 
#            group_by(cluster) %>% 
#            summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
#   mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
#   write.table("tf_new_sorted_table_list_samples_masked_conSeq.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)


rbind(ordered_table2 %>% 
    filter(sample=="ref_loc"), 
  ordered_table2 %>% 
    filter(sample=="Pomberef"), 
  ordered_table2 %>% 
    filter(sample=="LTR_loc")) %>% 
  arrange(chr, start) %>% 
  filter(len_froSeq>2000) %>%
  select(-rep, -dir, -dis, -ancHaplo, -anc_prop) %>%
  tbl_df %>% print(n = Inf)
  # filter(chr=="I", start>3416193)
  filter(chr=="III")

rbind(ordered_table2 %>% 
        filter(sample=="Pomberef"), 
      ordered_table2 %>% 
        filter(sample=="LTR_loc")) %>% 
  mutate(LTRs=ifelse(len>5000, 4, 
                     ifelse(len>500,2,1))) %>%
  filter(len>100) %>%
  group_by(cluster, chr) %>%
  # summarise(count_seq=n(), 
  #           Pomberef=sum(sample=="Pomberef"), 
  #           Bowen_etal=sum(sample=="LTR_loc")) %>%
  summarise(count_seq=n(), 
            Pomberef=sum((sample=="Pomberef")*LTRs), 
            Bowen_etal=sum((sample=="LTR_loc")*LTRs)) %>%
  filter(chr!="AB325691") %>%
  arrange(chr, cluster) %>% 
  mutate(diff_num_LTRs_B_Pomberef=Bowen_etal-Pomberef) %>%
  ggplot(aes(factor(cluster), diff_num_LTRs_B_Pomberef, 
             colour=factor(Pomberef))) +
  geom_point(alpha=0.4)+
  scale_y_continuous(breaks=seq(-10,10,1))+
  labs(colour="Num. Seq Ref", x="Cluster", 
       y="Dif. Num. Seq LTRs
       (Bowen et al. - Ref)") + 
  facet_grid(. ~ chr, scale = "free", space="free") +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        legend.position = "bottom") #+
  # ggsave("06_LTR_elementsRefPombeVsBowenetal.png", width = 7, height = 4)
  # ggsave("06_LTR_elementsRefPombeVsBowenetal_min100.png", width = 7, height = 4)
  


used_sample="JB22_EBC2" #"JB22_EBC2" #"Pomberef" #"JB22"
used_sample2="LTR_loc" #"LTR_loc"
rbind(ordered_table2 %>% 
        # filter(sample=="Pomberef"), 
        filter(sample==used_sample), 
      ordered_table2 %>% 
        filter(sample==used_sample2)) %>% 
  mutate(LTRs=ifelse(len>5000, 4, 
                     ifelse(len>500,2,1))) %>%
  filter(len>150) %>%
  group_by(cluster, chr) %>% 
  # summarise(count_seq=n(),
  #           Pomberef=sum(sample=="Pomberef"),
  #           Bowen_etal=sum(sample=="LTR_loc")) %>%
  summarise(count_seq=n(),
            raw_Pomberef=sum(sample==used_sample),
            raw_Bowen_etal=sum(sample==used_sample2),
            Pomberef=sum((sample==used_sample)*LTRs),
            Bowen_etal=sum((sample==used_sample2)*LTRs)) %>%
  filter(chr!="AB325691") %>%
  # filter(raw_Pomberef1!=0) %>%
  # filter(raw_Bowen_etal!=0) %>%
  arrange(chr, cluster) %>% 
  mutate(diff_num_LTRs_B_Pomberef=Bowen_etal-Pomberef, 
         raw_diff_num_LTRs_B_Pomberef=raw_Bowen_etal-raw_Pomberef) %>%
  # filter(diff_num_LTRs_B_Pomberef<0) %>%
  # head()
  # filter(diff_num_LTRs_B_Pomberef == 0) %>%
  # dim()
  # head()
  ggplot(aes(diff_num_LTRs_B_Pomberef)) +
  # geom_histogram(aes(y=..density..))+
  geom_histogram()+
  # xlim(c(-4,2))+
  # scale_x_continuous(limits = c(-4.5,2.5),breaks=seq(-10,10,1))+
  scale_y_continuous(breaks=seq(0,160,20))+
  scale_x_continuous(breaks = seq(-10,10,1))+
  # scale_y_continuous(limits = c(0,3.5), breaks=seq(0,3.5,0.5))+
  labs(x="Dif. Num. Seq LTRs per Cluster (Bowen etal - Ref)") + 
  #facet_grid(. ~ chr, scale = "free", space="free") +
  theme_classic() #+
  # ggsave("07_LTR_elementsRefPombeVsBowenetal_his.png", width = 4, height = 3)
  # ggsave("07_LTR_elementsRefPombeVsBowenetal_his_min150.png", width = 4, height = 3)
# ggsave("07_LTR_elementsRefPombeVsBowenetal_his_min150.png", width = 4, height = 3)
# ggsave("07_LTR_elementsRefPombeVsBowenetal_his_min150.svg", width = 4, height = 3)
  # ggsave("07_LTR_elementsRefPombeVsBowenetal_den.png")
  # ggsave("07_LTR_elementsRefPombeVsBowenetal_den_min150.png")

  
  
  

### Comparison with Jeffares 2015:
jeffares_table<-read.csv("Jeffares_2015.csv", T)

unique_samples_conJeff<-c("JB879",
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
                          "JB758", "JB22_EBC2")


unique_samples_conJeff2<-data.frame(sample=unique_samples_conJeff, ed_sample=sapply(strsplit(unique_samples_conJeff, "_EBC[0-9]*"), "[", 1))
unique_samples_conJeff2

ordered_jeffares_table<-
  jeffares_table %>% 
  gather("sample", "genotype", c(-chr, -start, -stop)) %>%
  filter(genotype!=0) %>% 
  mutate(seq_ID=sample, 
         end=stop,
         dis=0, 
         len_froSeq=end-start, 
         len=end-start,
         addID="Jeffares", 
         dir="NA", 
         rep="NA") %>% 
  filter(sample %in% unique_samples_conJeff2$ed_sample) %>%
  filter(!(sample %in% c("JB1110", "JB900", "JB1174"))) %>% # I excluded these samples which are sample that did not show consisten clustering in the phylogeny. 
  select(seq_ID, sample, chr, start, dis, len_froSeq, addID, dir, len, rep) 


unique(ordered_jeffares_table$sample)

ordered_table_comJeffares<-rbind(samples_table %>% 
                                   filter(sample %in% 
                                            unique_samples_conJeff2$sample) %>%
                                   merge(unique_samples_conJeff2, by="sample") %>%
                                   filter(ed_sample %in% 
                                            ordered_jeffares_table$sample) %>% 
                                   mutate(sample=ed_sample) %>% 
                                   select(seq_ID, sample, 
                                          chr, start, dis,
                                          len_froSeq, addID, 
                                          dir, len, rep), 
                                 ordered_jeffares_table) %>% 
  # merge(ID_table, by="sample", all.x=T) %>% 
  arrange(chr, start, sample, dis) %>% 
  filter(!(is.na(chr)))


cluster<-c(1)
length_cluster<-0
min_start<-0


# it worked quite ok with 150
# also with 100 bp
for (i in seq(2,dim(ordered_table_comJeffares)[1])){
  # i=2
  print(i)
  if(length_cluster==0){
    length_cluster<-ordered_table_comJeffares$len[i-1]
    min_start<-ordered_table_comJeffares$start[i-1]
  }
  if(ordered_table_comJeffares$chr[i]==ordered_table_comJeffares$chr[i-1]){
    # if ((ordered_table_comJeffares$start[i]-(min_start+length_cluster))<=(-50)){
    if ((ordered_table_comJeffares$start[i]-(min_start+length_cluster))<=(100)){
      print(i)
      cluster[i]<-cluster[i-1]
      length_cluster<-max(length_cluster,ordered_table_comJeffares$len[i])
      min_start<-min(min_start, ordered_table_comJeffares$start[i])
    } else {
      print(i)
      cluster[i]<-cluster[i-1]+1
      length_cluster<-ordered_table_comJeffares$len[i]
      min_start<-ordered_table_comJeffares$start[i]
    }
  } else {
    print(i)
    cluster[i]<-cluster[i-1]+1
    length_cluster<-ordered_table_comJeffares$len[i]
    min_start<-ordered_table_comJeffares$start[i]
  }
}


ordered_table_comJeffares2<-cbind(ordered_table_comJeffares, cluster) 


max(ordered_table_comJeffares2$cluster)
ordered_table_comJeffares2$addID

ordered_table_comJeffares2 %>% 
  filter(chr!="AB325691") %>% 
  #filter(start>100000) %>%
  group_by(sample, cluster, chr) %>%
  summarise(pos=min(start), 
            count_seq=n(),
            raw_sample=sum(addID!="Jeffares"),
            raw_jeffares=sum(addID=="Jeffares")) %>%
  mutate(diff_seq=raw_sample-raw_jeffares, 
         presence_diff=ifelse(diff_seq<0, "Jeff", 
                              ifelse(diff_seq==0, "Equal", "LR")),
         presence_group=ifelse(raw_sample!=0 & raw_jeffares!=0, "Both", 
                               ifelse(raw_sample!=0 & raw_jeffares==0, "LR", "Jeff"))) %>% 
  # filter(pos>70000) %>%
  # ### plot with distribution of in consistent clusters
  # filter(presence_diff == "Equal")%>%
  # ggplot(aes(pos, presence_diff, colour=factor(raw_sample)))+
  # geom_point()+
  # facet_grid(sample ~ chr, scale="free", space="free")
  # plot with total number of inconsisten cluster per sample:
  group_by(sample, presence_group) %>%
  summarise(count_clusters=n()) %>%
  # group_by(sample) %>%
  # summarise(total_clusters=sum(count_clusters), 
  #           inconsistent=sum((presence_group!="Both")*count_clusters)) %>%
  # mutate(fraction=inconsistent*100/total_clusters) %>%
  # arrange(fraction) %>% 
  # data.frame()
  ggplot(aes(sample, count_clusters, fill=presence_group))+
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Sample", y="Number of clusters") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black")) #+
  # ggsave("08_LTR_elementsLRVsJeffaresetal.png", width = 6, height = 3)
  # ggsave("08_LTR_elementsLRVsJeffaresetal.svg", width = 6, height = 3)


# plot size distribution and number of clusters:
# ordered_table_plot<-ordered_table %>% 
#   #filter(seq_ID!="NA") %>% 
#   merge((ordered_table %>% 
#            filter(seq_ID!="NA") %>% 
#            group_by(cluster) %>% 
#            summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
#   mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
#   filter(len<6000)  

ordered_table_plot <- ordered_table2 %>%
  merge((ordered_table2 %>% 
           group_by(cluster) %>% 
           summarize(IDcluster=paste(cluster[1],
                                     min(start, na.rm = T), 
                                     sep="_"))), 
        by="cluster") %>% 
  mutate(IDcluster = fct_reorder(IDcluster, cluster)) #%>%
  # filter(len_froSeq<6000) 

ordered_table2 %>% 
  filter(len_froSeq<0)

for(chr_used in c("I", "II", "III")){
ordered_table_plot %>% 
  mutate(colour_point=ifelse(seq_ID=="NA", "Ref", 
                             ifelse(sample=="ref_loc","Ref_predicted", 
                                    ifelse(((cluster %% 2) != 0), "Samples1", "Sample2")))) %>% 
  mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>% 
  filter(chr==chr_used) %>%
  #filter(!(is.na(len=NA))) %>% 
  ggplot(aes(IDcluster, len_froSeq)) +
# ordered_table_plot %>%
#   filter(seq_ID!="NA") %>% 
  geom_line(aes(group=IDcluster)) +
  geom_jitter(aes(colour = colour_point, alpha=colour_point),
              size=2) +
  # geom_line(aes(group=cluster), alpha=0.4)+
  # geom_jitter(aes(colour = as.factor(cluster)), size=1, alpha=0.8) +
  # geom_jitter(data=ordered_table_plot %>%
  #               #mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>% 
  #               filter(sample=="ref"),
  #             # aes(colour = as.factor(cluster)),
  #             colour="blue",
  #             size=2, alpha=0.6) +
  # geom_point(data=ordered_table_plot %>%
  #              #mutate(cluster=factor(cluster)) %>%
  #              #mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
  #              filter(sample=="ref_loc"),
  #            # aes(colour = cluster),
  #            colour="red",
  #            size=2, alpha=0.8) +
  xlab("Cluster ID") +
  ylab("Length (bp)") +
  scale_color_manual(values=c("red", "blue", "gray", "black")) +
  scale_alpha_manual(values=c(0.9,0.9,0.4, 0.4)) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.text.y=element_text(colour="black", size=14),
        axis.title=element_text(colour="black", size=14),
        legend.position="none") +
  facet_grid(. ~ chr, scale="free", space="free") 
  ggsave(paste0("01_size_distribution_tf_", chr_used,".png"), width = 11, height = 3)
  # ggsave("01_size_distribution_tf_III.png", width = 11, height = 3)
  # ggsave("01_size_distribution_tf_III.svg", width = 11, height = 3)
  # ggsave("01_size_distribution_tf_II.png", width = 11, height = 3)
  # ggsave("01_size_distribution_tf_II.svg", width = 11, height = 3)
  # ggsave("01_size_distribution_tf_I.png", width = 11, height = 3)
  # ggsave("01_size_distribution_tf_I.svg", width = 11, height = 3)
  # ggsave("01_size_distribution_tf.png", width = 11, height = 3)
  # ggsave("01_size_distribution_tf.svg", width = 11, height = 3)
}


ordered_table2 %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  filter(len_froSeq<6000) %>%
  ggplot(aes(len_froSeq)) +
  geom_histogram() +
  xlab("Sequence Length") +
  ylab("Count") +
  theme_classic() +
  theme(axis.text.x=element_text(colour="black", size=12), 
        axis.text.y=element_text(colour="black", size=12),
        axis.title=element_text(colour="black", size=14)) 
# ggsave("01_size_ditribution_hist.png", width = 5, height = 4)
# ggsave("01_size_ditribution_hist.svg", width = 5, height = 4)


ordered_table2 %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  # filter(len_froSeq<7000) %>% 
  group_by(sample) %>% 
  summarise(short_seq=sum((len_froSeq<=1500)*1), 
            full_seq=sum((len_froSeq>1500)*1)) %>% 
  mutate(sample=factor(sample, levels=order_samples_pro_AncPop)) %>% 
  # head()
  # ggplot(aes(sample, full_seq)) +
  # geom_bar(stat="identity")+
  # labs(x="Sample", y="Number of Sequences") +
  # theme_classic() +
  # theme(axis.text.x = element_text(angle = 45,
  #                                  hjust = 1,
  #                                  size=9, colour="black"),
  #       axis.text.y = element_text(size=9, colour="black")) +
  # ggsave("01_Num_fulllengthLTR_perSample.png", width = 6, height = 3)
  gather("seq_group", "num_seq", -sample) %>% 
  # head()
  ggplot(aes(sample, num_seq, fill=seq_group)) +
  geom_bar(stat="identity")+
  scale_fill_grey()+
  labs(x="Sample", y="Number of Sequences") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size=9, colour="black"),
        axis.text.y = element_text(size=9, colour="black")) #+
  # ggsave("01_Num_seq_perSample.png", width = 7, height = 3)

#
#

ordered_table2 %>% 
  filter(cluster=="582")
    


ordered_table2 %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  # filter(len_froSeq<7000) %>% 
  group_by(sample, cluster) %>% 
  # summarise(fulllengthpresent=(sum(len_froSeq>2000)>0)*1,
  summarise(fulllengthpresent=(sum(len_froSeq>4500)>0)*1, 
            present_cluster=ifelse(fulllengthpresent==1,0,1)) %>% 
  group_by(sample) %>%
  summarise(short_clusters=sum(present_cluster),
            fulllength_min4500=sum(fulllengthpresent)) %>% 
            # fulllength_min2000=sum(fulllengthpresent)) %>% 
  gather("cluster_group", "num_cluster", -sample) %>% 
  mutate(sample=factor(sample, levels=order_samples_pro_AncPop)) %>% 
  # head()
  ggplot(aes(sample, num_cluster, fill=cluster_group)) +
  geom_bar(stat="identity", position="dodge")+
  # geom_bar(stat="identity")+
  scale_fill_grey()+
  labs(x="Sample", y="Number of Clusters") + 
  scale_y_continuous(breaks=seq(0,200,20)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black")) #+
  # ggsave("01_Num_clusters_perSample.png", width = 8, height = 3)
  # ggsave("01_Num_clusters_perSample_dodge.png", width = 9, height = 3)
  

ordered_table2 %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  pull(cluster) %>%
  unique() %>%
  length()

ordered_table2 %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  group_by(sample) %>% 
  summarise(fulllengthpresent=(sum(len_froSeq>2000)),
            solo_seq=(sum(len_froSeq<2000)),
            total=n()) %>%
  data.frame()
  ungroup() %>%
  summarise(max_full_length=max(fulllengthpresent), 
            min_full_length=min(fulllengthpresent), 
            max_solo=max(solo_seq), 
            min_solo=min(solo_seq), 
            max_total=max(total), 
            min_total=min(total))



ordered_table2 %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>% 
  filter(!(chr %in% c("AB325691", "NoLoc"))) %>%
  group_by(chr,cluster,sample) %>%
  summarise(cluster_found=1, 
            full_length=(sum(len_froSeq>2000)>0)*1) %>%
  group_by(chr, sample) %>%
  summarise(Num_clusters=sum(cluster_found), 
            Num_cluster_FL=sum(full_length)) %>%
  mutate(chr_size=ifelse(chr=="I", 5.58, 
                         ifelse(chr=="III", 2.45, 4.54)), 
         cluster_perMB=Num_clusters/chr_size, 
         cluster_perMB_FL=Num_cluster_FL/chr_size) %>%
  ggplot(aes(sample, cluster_perMB, fill=chr)) +
  # ggplot(aes(sample, cluster_perMB_FL, fill=chr)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="Sample", y="Number of Clusters/Mb") + 
  scale_y_continuous(breaks=seq(0,50,1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black")) #+
  # ggsave("01_Num_clusters_perSample_dodge_perChrMb.png", width = 9, height = 3)
  # ggsave("01_Num_clusters_perSample_dodge_perChrMb_FLmin2000.png", width = 9, height = 3)



ordered_table2 %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  group_by(sample) %>% 
  summarise(num_cluster=length(unique(cluster)), 
            total_seq=n(),
            fulllength=sum((len_froSeq>4500)), 
            solo_LTR_Fragments=sum((len_froSeq<900 & len_froSeq>200))) %>%
            # fulllength_min2000=sum((len_froSeq>2000)), 
            # solo_LTR_min800=sum((len_froSeq<800))) %>%
  gather("group_v", "value", -sample) %>% 
  mutate(sample=factor(sample, levels=order_samples_pro_AncPop)) %>% 
  mutate(group_v_ed=factor(group_v, levels=c("num_cluster", "total_seq", 
                                             "solo_LTR_Fragments", "fulllength"))) %>%
  # filter(group_v_ed=="solo_LTR_min800") %>%
  filter(group_v_ed!="total_seq") %>%
  # head()
  ggplot(aes(sample, value, fill=group_v_ed)) +
  geom_bar(stat="identity", position="dodge")+
  #geom_bar(stat="identity")+
  # scale_fill_grey()+
  scale_fill_manual(values = c("grey80", "grey50", "black"))+
  # scale_fill_manual(values = c("black"))+
  labs(x="Sample", y="Number of seq./Clusters") + 
  scale_y_continuous(breaks=seq(0,400,50)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black")) #+
# ggsave("01_Num_seq_perSample_dodge.svg", width = 9, height = 3)
# ggsave("01_Num_seq_perSample_dodge.png", width = 9, height = 3)

ordered_table2 %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  group_by(sample) %>% 
  summarise(num_cluster=length(unique(cluster)), 
            total_seq=n(),
            fulllength_min2000=sum((len_froSeq>4500)), 
            solo_LTR_min800=sum((len_froSeq<900))) %>%
  gather("group_v", "value", -sample) %>% 
  mutate(sample=factor(sample, levels=order_samples_pro_AncPop)) %>% 
  mutate(group_v_ed=factor(group_v, levels=c("num_cluster", "total_seq", 
                                             "solo_LTR_min800", "fulllength_min2000"))) %>%
  filter(group_v_ed=="total_seq") %>%
  ggplot(aes(sample, value, fill=group_v_ed)) +
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values = c("black"))+
  labs(x="Sample", y="Number of seq.") + 
  scale_y_continuous(breaks=seq(0,400,50)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black"),
        legend.position = "none") 
  # ggsave("NP_01_Num_seq_perSample.svg", width = 9, height = 3)
  # ggsave("01_Num_seq_perSample.png", width = 9, height = 3)

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



statistical_test<-ordered_table2 %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  group_by(sample) %>% 
  summarise(num_cluster=length(unique(cluster)), 
            total_seq=n(),
            fulllength_min2000=sum((len_froSeq>4500))) %>%
            # fulllength_min2000=sum((len_froSeq>2000))) %>%
  left_join(pro_AncPop2, by="sample") %>%
  filter(sample %in% unique_samples) %>% 
  mutate(sk_anc_prop=1-anc_prop) %>%
  data.frame()
  
m1a<-lm(total_seq~sk_anc_prop, data=statistical_test)
result_stat<-summary(m1a)
result_stat


library(ggpmisc)

ordered_table2 %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  group_by(sample) %>% 
  summarise(num_cluster=length(unique(cluster)), 
            total_seq=n(),
            fulllength_min2000=sum((len_froSeq>2000))) %>%
  left_join(pro_AncPop2, by="sample") %>%
  filter(sample %in% unique_samples) %>% 
  mutate(sk_anc_prop=1-anc_prop) %>%
  ggplot(aes(sk_anc_prop, total_seq))+
  geom_point(alpha=0.5, size=3)+
  geom_smooth(method='lm', formula= y~x, se=FALSE, colour="black")+
  annotate(geom="text", x=0.1, y=270,
           label=paste0("Adj.R2 = ",
                        format(round(result_stat$adj.r.squared, 2), nsmall = 2),
                        "\n",
                        "p-value = ",
                        format(round(result_stat$coefficients["sk_anc_prop","Pr(>|t|)"], 3), nsmall = 3)),
           color="black", size=4)+
  # stat_poly_eq(formula = y~x,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE) +
  labs(x="Sk Proportion", 
       y="Num. Total Seq.") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+  
# ggsave("05_SKProportionVstotalNSeq.png", width = 5, height = 5)
# ggsave("05_SKProportionVstotalNSeq.svg", width = 5, height = 5)





# ordered_table2 %>% 
#   filter(seq_ID!="NA") %>%
#   filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
#   group_by(sample) %>% 
#   summarise(num_cluster=length(unique(cluster)), 
#             total_seq=n(),
#             fulllength_min2000=sum((len_froSeq>2000))) %>%
#   left_join(pro_AncPop2, by="sample") %>%
#   filter(sample %in% unique_samples) %>%
#   ggplot(aes(1-anc_prop, total_seq)) +
#   geom_point()
#   # geom_bar(stat="identity", position="dodge")+
#   #geom_bar(stat="identity")+
#   scale_fill_grey()+
#   labs(x="Sample", y="Number of seq./Clusters") + 
#   scale_y_continuous(breaks=seq(0,400,50)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, 
#                                    hjust = 1, 
#                                    size=9, colour="black"), 
#         axis.text.y = element_text(size=9, colour="black")) #+
# # ggsave("01_Num_seq_perSample_dodge.svg", width = 9, height = 3)
# # ggsave("01_Num_seq_perSample_dodge.png", width = 9, height = 3)


ordered_table2 %>% 
  filter(chr!="AB325691") %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  filter(chr!="NoLoc") %>% 
  arrange(chr, cluster, min_start) %>% 
  group_by(chr, cluster) %>% 
  summarise(min_start=min(start)) %>% 
  ungroup() %>% 
  arrange(chr, cluster, min_start) %>% 
  group_by(chr) %>% 
  mutate(distanceCluster = min_start - lag(min_start, default = min_start[1])) %>% 
  filter(distanceCluster!=0) %>% 
  # filter(distanceCluster<130000) %>%
  ggplot(aes(distanceCluster)) +
  geom_histogram() +
  xlab("Distance to previous cluster") +
  ylab("Count") +
  facet_grid(. ~ chr, scale="free")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                               hjust = 1, 
                                               size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14),
        axis.title=element_text(colour="black", size=14)) #+
  # ggsave("01_distanceBetweenclusters_perChr.png", width = 9, height = 3)


ordered_table2 %>% 
  filter(chr!="AB325691") %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  filter(chr!="NoLoc") %>% 
  arrange(chr, cluster, min_start) %>% 
  group_by(chr, cluster) %>% 
  summarise(min_start=min(start)) %>% 
  ungroup() %>% 
  arrange(chr, cluster, min_start) %>% 
  group_by(chr) %>% 
  mutate(distanceCluster = min_start - lag(min_start, default = min_start[1])) %>% 
  filter(distanceCluster!=0) %>% 
  filter(distanceCluster<5000) %>% 
  ggplot(aes(distanceCluster)) +
  geom_histogram(binwidth=100) +
  xlab("Distance to previous cluster") +
  ylab("Count") +
  scale_x_continuous(limits = c(0,5000))+
  facet_grid(. ~ chr)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14),
        axis.title=element_text(colour="black", size=14)) #+
  # ggsave("01_distanceBetweenclusters_perChr_below5kb.png", width = 9, height = 3)



ordered_table2 %>% 
  filter(chr!="AB325691") %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  filter(chr!="NoLoc") %>% 
  # group_by(chr, cluster) %>%
  # summarize(position=min(start)/1000000) %>%
  # ggplot(aes(position)) +
  ggplot(aes(start/1000000)) +
  geom_histogram(binwidth = 0.05)+
  xlab("Position (Mb)") +
  # ylab("Num. Clusters") +
  ylab("Num. sequences") +
  scale_x_continuous(breaks = seq(0,20,0.5))+
  facet_grid(. ~ chr, scale="free", space="free") +
  theme_classic() +
  theme(axis.text.x = element_text(colour="black", size=14), 
        axis.text.y=element_text(colour="black", size=14),
        axis.title=element_text(colour="black", size=14)) #+
# ggsave("02_tf_clustersAlongGenome.png", width = 11, height = 3)
# ggsave("02_tf_clustersAlongGenome.svg", width = 11, height = 3)
ggsave("02_tf_SeqAlongGenome.png", width = 11, height = 3)



clustes_min1000<-ordered_table2 %>% 
  filter(seq_ID!="NA") %>%
  filter(len_froSeq>1000) %>% 
  ungroup() %>% 
  select(cluster) %>% 
  unique() %>% 
  unlist() %>% 
  as.vector()


ordered_table2 %>% 
  filter(chr!="AB325691") %>% 
  filter(chr %in% c("I", "II", "III")) %>% 
  filter(seq_ID!="NA") %>%
  merge((ordered_table2  %>% 
           filter(seq_ID!="NA") %>%
           group_by(cluster) %>% 
           summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
  mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
  # filter(len_froSeq<6000) %>%
  # filter(len_froSeq>1000) %>%
  # filter(cluster %in% clustes_min1000) %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  group_by(chr, IDcluster, sample) %>% 
  summarize(n_copies=length(seq_ID), 
            dif_length=(max(len_froSeq)-min(len_froSeq))) %>% 
  group_by(IDcluster) %>%
  summarize(max_nSeq=max(n_copies)) %>%
  ggplot(aes(max_nSeq)) +
  geom_histogram(binwidth = 1)+
  xlab("Max number of seq. per sample") +
  ylab("Num. Clusters") +
  scale_x_continuous(breaks = seq(1,20,2))+
  scale_y_continuous(breaks = seq(0,1000,50))+
  theme_classic() +
  theme(axis.text.x = element_text(colour="black", size=14), 
        axis.text.y=element_text(colour="black", size=14),
        axis.title=element_text(colour="black", size=14)) #+
# ggsave("02_tf_maxnCopies_persample.png", width = 5, height = 4)
# ggsave("02_tf_maxnCopies_persample.svg", width = 5, height = 4)

ordered_table2 %>% 
  filter(chr!="AB325691") %>% 
  filter(chr %in% c("I", "II", "III")) %>% 
  filter(seq_ID!="NA") %>%
  merge((ordered_table2  %>% 
           filter(seq_ID!="NA") %>%
           group_by(cluster) %>% 
           summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
  mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
  # filter(len_froSeq<6000) %>%
  # filter(len_froSeq>1000) %>%
  # filter(cluster %in% clustes_min1000) %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  filter(len_froSeq>4500) %>%
  group_by(chr, IDcluster, sample) %>% 
  summarize(n_copies=length(seq_ID), 
            dif_length=(max(len_froSeq)-min(len_froSeq))) %>% 
  group_by(IDcluster) %>%
  summarize(max_nSeq=max(n_copies)) %>%
  ggplot(aes(max_nSeq)) +
  geom_histogram(binwidth = 1)+
  xlab("Max number of seq. per sample") +
  ylab("Num. Clusters") +
  scale_x_continuous(breaks = seq(1,20,2))+
  scale_y_continuous(breaks = seq(0,1000,50))+
  theme_classic() +
  theme(axis.text.x = element_text(colour="black", size=14), 
        axis.text.y=element_text(colour="black", size=14),
        axis.title=element_text(colour="black", size=14)) #+
# ggsave("02_tf_maxnCopies_persample_min4500.png", width = 5, height = 4)
# ggsave("02_tf_maxnCopies_persample_min4500.svg", width = 5, height = 4)




# plot heat map with number of copies per cluster and sample:
ordered_table2 %>% 
  filter(chr!="AB325691") %>% 
  filter(chr %in% c("I", "II", "III")) %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  merge((ordered_table2  %>% 
           filter(seq_ID!="NA") %>%
           filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
           group_by(cluster) %>% 
           summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
  mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
  # filter(len_froSeq<6000) %>%
  # filter(len_froSeq>4500) %>%
  # filter(cluster %in% clustes_min1000) %>%
  filter(sample!="LTR_loc") %>% 
  group_by(chr, IDcluster, sample) %>% 
  summarize(n_copies=length(seq_ID), 
            dif_length=(max(len_froSeq)-min(len_froSeq))) %>% 
  mutate(text_label=ifelse(dif_length==0,"",as.character(dif_length))) %>% 
  mutate(sample=factor(sample, levels=order_samples_pro_AncPop)) %>%
  # filter(n_copies<9) %>%
  mutate(n_copies_category=factor(ifelse(n_copies<2, "1", 
                                  ifelse(n_copies<7, "2-6", ">6")),levels = c("1", "2-6", ">6"))) %>%
  ggplot(aes(IDcluster, sample)) +
  # geom_tile(aes(fill=factor(n_copies))) +
  geom_tile(aes(fill=n_copies_category)) +
  # scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#4E84C4", "#293352")) +
  scale_fill_manual(values=c("gray80", "#00AFBB", "red", "#E7B800", "#FC4E07", "#52854C", "#4E84C4", "#293352")) +
  xlab("Cluster ID") +
  ylab("Sample") +
  labs(fill="Num. Copies")+
  theme_classic() +
  theme(axis.text.x=element_blank(), 
        axis.text.y=element_text(colour="black", size=14),
        axis.title=element_text(colour="black", size=14),
        legend.position="top") +
  facet_grid(. ~ chr, scale="free", space="free") #+
# ggsave("02_heatmap_tf_nCopies_persampleR1.png", width = 14, height = 8)
# ggsave("02_heatmap_tf_nCopies_persampleR1.svg", width = 14, height = 8)
# ggsave("02_heatmap_tf_nCopies_persample.png", width = 14, height = 8)
# ggsave("02_heatmap_tf_nCopies_persample.svg", width = 14, height = 8)





ordered_table2 %>% 
  filter(chr!="AB325691") %>% 
  filter(chr %in% c("I", "II", "III")) %>% 
  filter(seq_ID!="NA") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  merge((ordered_table2  %>% 
           filter(seq_ID!="NA") %>%
           filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
           group_by(cluster) %>% 
           summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
  mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
  # filter(len_froSeq<6000) %>%
  filter(len_froSeq>4500) %>%
  # filter(cluster %in% clustes_min1000) %>%
  group_by(chr, IDcluster, sample) %>% 
  summarize(n_copies=length(seq_ID), 
            dif_length=(max(len_froSeq)-min(len_froSeq))) %>% 
  mutate(text_label=ifelse(dif_length==0,"",as.character(dif_length))) %>% 
  mutate(sample=factor(sample, levels=order_samples_pro_AncPop)) %>%
  # filter(n_copies<9) %>%
  mutate(n_copies_category=factor(ifelse(n_copies<2, "1", 
                                         ifelse(n_copies<5, "2-4", ">4")),levels = c("1", "2-4", ">4"))) %>%
  ggplot(aes(IDcluster, sample)) +
  # geom_tile(aes(fill=factor(n_copies))) +
  geom_tile(aes(fill=n_copies_category)) +
  # scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#52854C", "#4E84C4", "#293352")) +
  scale_fill_manual(values=c("gray70", "#00AFBB", "red", "#E7B800", "#FC4E07", "#52854C", "#4E84C4", "#293352")) +
  xlab("Cluster ID") +
  ylab("Sample") +
  labs(fill="Num. Copies")+
  theme_classic() +
  theme(axis.text.x=element_blank(), 
        axis.text.y=element_text(colour="black", size=14),
        axis.title=element_text(colour="black", size=14),
        legend.position="top") +
  facet_grid(. ~ chr, scale="free", space="free") #+
# ggsave("02_heatmap_tf_nCopies_persample_min4500.png", width = 14, height = 8)
# ggsave("02_heatmap_tf_nCopies_persample_min4500.svg", width = 14, height = 8)



ordered_table2 %>% 
  filter(chr!="AB325691") %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  # filter(len_froSeq>1000) %>%
  merge((ordered_table2 %>% 
           group_by(cluster) %>% 
           summarize(min_start=min(start))), by="cluster") %>% 
  group_by(cluster, chr, min_start, sample) %>% 
  filter(sample!="ref_loc") %>% 
  summarize(n_copies=n(), max_len=max(len_froSeq)) %>% 
  group_by(cluster, chr, min_start) %>% 
  summarize(n_samples=n(), range_len=paste0(min(max_len), "_", max(max_len))) %>% 
  mutate(min_start=min_start/1000000) %>%
  # filter(n_samples>10)
  ggplot(aes(min_start, n_samples)) +
  geom_point(size=2, alpha=0.5) +
  scale_x_continuous(breaks=seq(0,5.5,0.5), labels=seq(0,5.5,0.5))+
  scale_y_continuous(breaks=seq(0,40,2), labels=seq(0,40,2))+
  xlab("Pos.") +
  ylab("Number Samples") +
  theme_classic() +
  theme(axis.text.x=element_text(size=14, angle = 45, hjust = 1, colour="black"), 
        axis.text.y=element_text(colour="black", size=14)) +
  facet_grid(. ~ chr, scale="free", space="free") #+
  # ggsave("03_SamplesperCluster_tf_min1000.png", width = 12, height = 3)
  # ggsave("03_SamplesperCluster_tf.png", width = 12, height = 3)

unique(ordered_table2$sample)

ordered_table2 %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  # filter(len_froSeq>1500) %>%
  filter(len_froSeq>4500) %>%
  left_join((ordered_table2 %>% 
           group_by(cluster) %>% 
           summarize(min_start=min(start))), by="cluster") %>% 
  # filter(cluster=="515") %>%
  # filter((sample %in% no_clonal_strains)) %>%
  # arrange(sample) %>%
  # dim()
  # # print()
  # data.frame()
  # dim()
  # filter(len_froSeq>1000) 
  # filter(sample=="JB900_EBC131")
  # filter(len_froSeq<1500) %>% 
  # filter(sample!="ref_loc") %>% 
  group_by(cluster, sample) %>% 
  summarize(n_repeats=n()) %>% 
  filter((sample %in% no_clonal_strains)) %>%
  group_by(cluster) %>% 
  summarize(n_samples=n()) %>% 
  # filter(n_samples>27) %>%
  # filter(n_samples>9) 
  # pull(n_samples) %>%
  # max()
  # head()
  # # head()
  # dim()
  ggplot(aes(n_samples)) +
  geom_histogram(binwidth=1) +
  # geom_histogram(aes(y=..density..), binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,35,5), 38), 
                     labels=c(seq(0,35,5), 38)) +
  # ylim(c(0,300))+
  ylab("Count Clusters") +
  # ylab("Density") +
  xlab("Number Samples") +
  theme_classic() +
  theme(axis.text.x=element_text(size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14), 
        axis.title=element_text(colour="black", size=14)) #+
  # ggsave("04_fq_homologous_tf_counts.png", width = 7, height = 5)
  # ggsave("04_fq_homologous_tf_counts_min1500.png", width = 7, height = 5)
  # ggsave("04_fq_homologous_tf_density.png", width = 7, height = 5)


SFS_cluster_all_ordered_table2<-ordered_table2 %>% 
  # filter(len_froSeq>1500) %>%
  left_join((ordered_table2 %>% 
           group_by(cluster) %>% 
           summarize(min_start=min(start))), by="cluster", all.x = T) %>% 
  group_by(cluster, sample) %>% 
  summarize(n_repeats=n()) %>% 
  filter((sample %in% no_clonal_strains)) %>% 
  group_by(cluster) %>% 
  summarize(n_samples=n()) %>%
  mutate(group="all_LTR_seq")

SFS_cluster_fulllengthLTR_ordered_table2<-ordered_table2 %>% 
  filter(len_froSeq>4500) %>%
  # filter(len_froSeq>1500) %>%
  left_join((ordered_table2 %>% 
           group_by(cluster) %>% 
           summarize(min_start=min(start))), by="cluster", all.x = T) %>% 
  group_by(cluster, sample) %>% 
  summarize(n_repeats=n()) %>% 
  filter((sample %in% no_clonal_strains)) %>% 
  group_by(cluster) %>% 
  summarize(n_samples=n()) %>%
  mutate(group="full_length_LTR")


rbind(SFS_cluster_fulllengthLTR_ordered_table2, 
      SFS_cluster_all_ordered_table2) %>%
  mutate(group=factor(group, levels=c("full_length_LTR", "all_LTR_seq"))) %>%
  # # filter(group=="all_LTR_seq") %>%
  # filter(group=="full_length_LTR") %>%
  # # filter(n_samples < 2) %>%
  # filter(n_samples > 28) %>%
  # dim()
  # head()
  ggplot(aes(n_samples, fill=group)) +
  geom_histogram(binwidth=1, position="dodge") +
  scale_fill_grey()+
  scale_x_continuous(breaks=c(seq(0,35,5), 38), 
                     labels=c(seq(0,35,5), 38)) +
  scale_y_continuous(breaks=c(seq(0,300,50)))+
  #ylim(c(0,300))+
  ylab("Count Clusters") +
  xlab("Number Samples") +
  theme_classic() +
  theme(axis.text.x=element_text(size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14), 
        axis.title=element_text(colour="black", size=14),
        legend.position = "bottom") #+
  # ggsave("NP_04_fq_homologous_tf_counts_bothAll_FullLengthLTRs.png", width = 6, height = 4, dpi = 450)
  # ggsave("NP_04_fq_homologous_tf_counts_bothAll_FullLengthLTRs.svg", width = 6, height = 4, dpi = 450)
  
  

#library(gridExtra)
#library(ggExtra)
library(cowplot)

plot_table<-ordered_table2 %>% 
  merge(ordered_table2 %>% 
          group_by(cluster) %>% 
          summarize(longTE=(sum(len_froSeq>4500)>0)*1), by="cluster", all.x = T) %>% 
          # summarize(longTE=(sum(len_froSeq>1500)>0)*1), by="cluster", all.x = T) %>% 
  merge((ordered_table2 %>% 
           group_by(cluster) %>% 
           summarize(min_start=min(start))), by="cluster", all.x = T) %>%
  filter(longTE==1) %>%
  filter(sample!="ref_loc") %>% 
  group_by(cluster, sample) %>% 
  summarize(n_repeats=n(), n_repeats_1500=sum(len_froSeq>4500)) %>%
  # summarize(n_repeats=n(), n_repeats_1500=sum(len_froSeq>1500)) %>%
  #filter(n_repeats==1) %>%
  filter((sample %in% no_clonal_strains)) %>% 
  group_by(cluster) %>% 
  summarize(n_samples=n(), n_samples_1500=sum(n_repeats_1500>0))

scatter <- ggplot(plot_table, aes(n_samples, n_samples_1500)) +
  geom_hex(bins = 32) +
  # geom_bin2d(bins = 32) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()+
  scale_x_continuous(breaks=c(seq(0,35,5), 35), 
                     labels=c(seq(0,35,5), 35)) +
  scale_y_continuous(breaks=c(seq(0,35,5), 35), 
                     labels=c(seq(0,35,5), 35)) +
  ylab("Number Samples with seq > 1500 bp") +
  xlab("Number Samples") +
  theme_classic() +
  theme(axis.text.x=element_text(size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14), 
        axis.title=element_text(colour="black", size=14), 
        legend.position = "top", 
        axis.title.x=element_blank()) +
  border() +
  rremove("legend")

hist_bottom <- ggplot(plot_table, aes(n_samples)) +
  geom_histogram(binwidth=1) +
  # geom_histogram(aes(y=..density..), binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,35,5), 38), 
                     labels=c(seq(0,35,5), 38)) +
  #ylim(c(0,300))+
  ylab("Count Clusters") +
  # ylab("Density") +
  xlab("Number Samples") +
  theme_classic() +
  theme(axis.text.x=element_text(size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14), 
        axis.title=element_text(colour="black", size=14)) 

hist_right <-ggplot(plot_table, aes(n_samples_1500)) +
  geom_histogram(binwidth=1) +
  # geom_histogram(aes(y=..density..), binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,35,5), 38), 
                     labels=c(seq(0,35,5), 38)) +
  #ylim(c(0,300))+
  ylab("Count Clusters") +
  # ylab("Density") +
  xlab("Number Samples > 1500 bp") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust = 1, size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14), 
        axis.title=element_text(colour="black", size=14), 
        axis.title.y=element_blank()) +
  rotate() 

plot_grid(scatter, hist_right, hist_bottom, NULL, ncol = 2, align = "hv", 
          rel_widths = c(3, 1), rel_heights = c(3, 1)) #+
  # ggsave("./plots/04_fq_plusDensity_homologous_tf_counts_clusterswith1_min1500.png",
  # width = 10, height = 10)
  
ggplot(plot_table, aes(n_samples, n_samples_1500)) +
  geom_hex(bins = 32) +
  # geom_bin2d(bins = 32) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()+
  scale_x_continuous(breaks=c(seq(0,35,5), 35), 
                     labels=c(seq(0,35,5), 35)) +
  scale_y_continuous(breaks=c(seq(0,35,5), 35), 
                     labels=c(seq(0,35,5), 35)) +
  ylab("Number Samples with seq > 1500 bp") +
  xlab("Number Samples") +
  theme_classic() +
  theme(axis.text.x=element_text(size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14), 
        axis.title=element_text(colour="black", size=14), 
        legend.position = "top", 
        axis.title.x=element_blank()) #+
  # ggsave("./plots/04_Density_homologous_tf_counts_clusterswith1_min1500_legend.png",
  # width = 6, height = 6)


plot_table2<-ordered_table2 %>% 
  merge(ordered_table2 %>% 
          group_by(cluster) %>% 
          summarize(longTE=(sum(len_froSeq>4500)>0)*1), by="cluster", all.x = T) %>% 
  merge((ordered_table2 %>% 
           group_by(cluster) %>% 
           summarize(min_start=min(start))), by="cluster", all.x = T) %>%
  filter(longTE==1) %>%
  filter(sample!="ref_loc") %>% 
  group_by(cluster, sample) %>% 
  summarize(n_repeats=n(), n_repeats_4500=sum(len_froSeq>4500)) %>%
  #filter(n_repeats==1) %>%
  filter((sample %in% no_clonal_strains)) %>% 
  group_by(cluster) %>% 
  summarize(n_samples=n(), n_samples_4500=sum(n_repeats_4500>0))

# scatter <- 
ggplot(plot_table2, aes(n_samples, n_samples_4500)) +
  geom_hex(bins = 32) +
  # geom_bin2d(bins = 32) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()+
  scale_x_continuous(breaks=c(seq(0,35,5), 35), 
                     labels=c(seq(0,35,5), 35)) +
  scale_y_continuous(breaks=c(seq(0,35,5), 35), 
                     labels=c(seq(0,35,5), 35)) +
  ylab("Number Samples with seq > 4500 bp") +
  xlab("Number Samples") +
  theme_classic() +
  theme(axis.text.x=element_text(size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14), 
        axis.title=element_text(colour="black", size=14), 
        legend.position = "top", 
        axis.title.x=element_blank()) +
  border() 

ordered_table2 %>% 
  merge(ordered_table2 %>% 
          group_by(cluster) %>% 
          summarize(longTE=(sum(len_froSeq>1500)>0)*1), by="cluster", all.x = T) %>% 
  merge((ordered_table2 %>% 
           group_by(cluster) %>% 
           summarize(min_start=min(start))), by="cluster", all.x = T) %>%
  filter(longTE==1) %>%
  filter(sample!="ref_loc") %>% 
  group_by(cluster, sample) %>% 
  summarize(n_repeats=n(), n_repeats_1500=sum(len_froSeq>1500)) %>%
  filter(n_repeats==1) %>%
  filter((sample %in% no_clonal_strains)) %>% 
  group_by(cluster) %>% 
  summarize(n_samples=n(), n_samples_1500=sum(n_repeats_1500>0))  %>% 
  ggplot(aes(n_samples, n_samples_1500)) +
  geom_hex(bins = 32) +
  # geom_bin2d(bins = 32) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()+
  scale_x_continuous(breaks=c(seq(0,35,5), 35), 
                     labels=c(seq(0,35,5), 35)) +
  scale_y_continuous(breaks=c(seq(0,35,5), 35), 
                     labels=c(seq(0,35,5), 35)) +
  ylab("Number Samples with seq > 1500 bp") +
  xlab("Number Samples") +
  ylim(c(-1,32)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14), 
        axis.title=element_text(colour="black", size=14), 
        legend.position = "top", 
        axis.title.x=element_blank()) #+
  # ggsave("./plots/04_homologous_tf_counts_clusterswith1_min1500_legend_max1SeqperCluster.png",
  # width = 4, height = 4)
  

# empty <- ggplot()+geom_point(aes(1,1), colour="white")+
#   theme(axis.ticks=element_blank(), 
#         panel.background=element_blank(), 
#         axis.text.x=element_blank(), axis.text.y=element_blank(),           
#         axis.title.x=element_blank(), axis.title.y=element_blank())
#   
# grid.arrange(scatter, hist_right, hist_bottom, empty, ncol=2, nrow=2, 
#              widths=c(4, 1), heights=c(4, 1))



ordered_table2 %>% 
  filter(cluster==551) %>%
  # filter(start<256000) %>%
  select(sample, seq_ID, len_froSeq, chr, start, dis) %>%
  # ggplot(aes(start,len_froSeq))+
  # geom_point()
  # filter(len_froSeq>4000) %>%
  tbl_df %>% print(n=Inf)


### Comparison between Clones:

library(VennDiagram)
library(svglite)

ordered_table2 %>% head()
sample1="Pomberef"; sample2="JB22_EBC2"
# sample1="Pomberef"; sample2="JB22_EBC2"
# sample1="JB22_EBC2"; sample2="JB22"

# sample1="JB760"; sample2="JB760_EBC074"

# sample1="JB4"; sample2="JB4_EBC069"

# sample1="JB873"; sample2="JB873_EBC095"

# sample1="JB1206"; sample2="JB1206_EBC138"

# sample1="JB1174_EBC132"; sample2="JB900_EBC131"

for(run in c(1)){
comparison_venn_table<-ordered_table2 %>%
  filter(sample %in% c(sample1, sample2)) %>%
  group_by(cluster) %>%
  summarise(N_seq=n(), 
            Nseq_sample1=sum(sample==sample1), 
            Nseq_sample2=sum(sample==sample2), 
            common=min(Nseq_sample2,Nseq_sample1), 
            N_fullLenghtLTR_sample1=sum(((sample==sample1) & (len_froSeq>1500))), 
            N_fullLenghtLTR_sample2=sum(((sample==sample2) & (len_froSeq>1500)))) %>%
  mutate(diff=Nseq_sample2-Nseq_sample1,
         diff_fullLenght=N_fullLenghtLTR_sample2-N_fullLenghtLTR_sample1,
         group=if_else(diff==0, "equal", 
                 if_else(diff>0, "higher_sample2", "higher_sample1"))) %>% 
  # filter(diff!=0)
  # ungroup() %>%
  # summarise(N_total_seq_sample1=sum(Nseq_sample1), 
  #           N_total_seq_sample2=sum(Nseq_sample2)) 
  group_by(group) %>% 
  summarise(total_common=sum(common),
            sum_sample1=sum(Nseq_sample1), 
            sum_sample2=sum(Nseq_sample2), 
            extra_seq=sum(abs(diff))) 
common<-sum(comparison_venn_table$total_common)
extra_sample1<-as.vector(unlist(comparison_venn_table[comparison_venn_table$group=="higher_sample1", "extra_seq"]))

extra_sample2<-as.vector(unlist(comparison_venn_table[comparison_venn_table$group=="higher_sample2", "extra_seq"]))

svg(filename=paste0("02_vennD_",sample1,"_vs_",sample2,".svg"), 
    width=3, 
    height=3, 
    pointsize=12)
grid.newpage()
draw.pairwise.venn(extra_sample1+common, 
                   extra_sample2+common,
                   common, 
                   category = c("JB900_EBC132", "JB900_EBC131"), 
                   #category = c(sample1, sample2), 
                   lty = rep("blank", 2), 
                   fill = c("light blue", "pink"), 
                   alpha = rep(0.5, 2), 
                   cat.pos = c(-30,30), 
                   cat.dist = rep(0.025, 2)) 
dev.off()
png(filename=paste0("02_vennD_",sample1,"_vs_",sample2,".png"), 
    units="in", width=3, height=3, pointsize=12, res=400)
grid.newpage()
draw.pairwise.venn(extra_sample1+common, 
                   extra_sample2+common,
                   common, 
                   category = c(sample1, sample2), 
                   lty = rep("blank", 2), 
                   fill = c("light blue", "pink"), 
                   alpha = rep(0.5, 2), 
                   cat.pos = c(-30,30), 
                   cat.dist = rep(0.025, 2)) 
dev.off()
}


  
###

library("ggtree")
library(treeio)
#tree<-read.tree("../phylogenies/phy_tf_minLen1000.treefile")
tree_ed<-treeio::read.iqtree("../Phylogenies/LTR_completed/IQTree_bb/tf_3000.treefile")
#tree_ed<- as.phylo(tree)
fig_tree<- ggtree(tree_ed) + 
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
  #geom_tiplab() + 
  #geom_tippoint() + 
  theme_tree2()  
  #xlim(0, 0.5)
# fig_tree

tree_ed@data$UFboot[is.na(tree_ed@data$UFboot)]<-0
tree_ed@data$UFboot[tree_ed@data$UFboot<90]<-0
tree_ed@data$UFboot[tree_ed@data$UFboot>=90 & tree_ed@data$UFboot<95]<-90
tree_ed@data$UFboot[tree_ed@data$UFboot>=95]<-95


samples_order<-rev(tree_ed@phylo$tip.label)
samples_order2<-str_remove_all(samples_order, "_R_")
samples_order2<-sapply(strsplit(samples_order2, "_AB325691_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)
samples_order3<-sapply(strsplit(samples_order2, "_I*_[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*"), "[", 1)
samples_order3

tree_ed@phylo$tip.label<-rev(samples_order3)

head(ordered_table2) 

head(tree_ed@phylo$tip.label)
ordered_table2 %>% 
  filter(len_froSeq>3000) %>% 
  head()

ordered_table3<-ordered_table2 %>% 
  mutate(tem_id=seq_ID) %>% 
  separate(tem_id, into=c("seqNum"), sep="_") %>% 
  mutate(id=paste0(seqNum, "_", sample)) 

cluster_ID<-data.frame(id=tree_ed@phylo$tip.label) %>% 
  merge(ordered_table3 %>% 
          select(id, cluster, chr, anc_prop, ancHaplo), by="id", all.x=T) %>% 
  mutate(anc_group=factor(ifelse(anc_prop>0.85,"Sp",ifelse(anc_prop<0.15,"Sk","Hyb")), 
                          levels=c("Sp", "Hyb", "Sk"))) %>%
  arrange(match(id, tree_ed@phylo$tip.label))

row.names(cluster_ID) <- NULL

head(cluster_ID)

lables_tf<-ordered_table3 %>% 
  select(chr, cluster) %>%
  group_by(chr) %>% 
  dplyr::summarise(min_clus=min(cluster), 
                   max_clus=max(cluster), 
                   #multiple30=ifelse(cluster %% 30 == 0, cluster, 0)) %>% 
                   mid_1=((max_clus-min_clus)/4)+min_clus, 
                   mid_3=((max_clus-min_clus)*3/4)+min_clus, 
                   middle=round((max_clus-min_clus)/2)+min_clus) %>%
  ungroup() %>% 
  select(-chr, -min_clus, -max_clus) %>%  
  #filter(multiple30!=0) %>% 
  #select(-chr) %>%  
  unlist() %>% 
  as.vector() %>% 
  round() %>% 
  sort() %>% 
  as.numeric()

p <-ggtree(tree_ed) +
  #aes(color=as.factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)))), 
  #
  theme_tree2() +
  geom_nodepoint(aes(alpha=as.factor(UFboot)), color="darkgreen", size=2) +
  #scale_alpha_manual(values=c(0,0.5,1)) + 
  # geom_nodepoint(aes(color=factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)), levels=c(90,80,0)), 
  #                    alpha=factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0))), levels=c(90,80,0)), 
  #                               size=3) +
  geom_treescale() +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
  #geom_tiplab() + 
  #geom_tippoint(color="black", shape=4, size=1.5, alpha=0.2) +
  #geom_tree() +
  # geom_nodepoint(aes(color=as.factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)))), 
  #                alpha=0.5, 
  #                size=3) +
  #geom_tree(data=df, color='firebrick')
  #scale_color_continuous(low='darkgreen', high='red') +
  #scale_color_manual(values=c("gray", "black", "darkred")) +
  #scale_color_manual(values=c("white", "white", "darkred")) +
theme(legend.position="right")
#geom_tiplab(align=TRUE, linesize=1, col="gray90") 
# xlim(0, 0.4) 

# pp <- p %<+% cluster_ID + 
#   geom_treescale()+
#   geom_tiplab(align = TRUE, linesize = 0, alpha=0) +
#   geom_tippoint(aes(colour=anc_group,
#                     alpha=anc_group), 
#                     size=2) +
#   scale_color_manual(values=c("darkred", "orange", "steelblue")) +
#   scale_alpha_manual(values=c(0.8,0.0,0.8)) +
#   theme(legend.position = "None")

p
pp<-p
print(pp)
node_used<-identify(pp)
pp3<-pp+geom_nodelab(node_used)
#hilight(identify(pp), fill = "gray")
print(pp3)
# cluster_ID<-data.frame(id=tree_ed$tip.label) %>% 
#   merge(ordered_table2 %>% 
#           mutate(id=seq_ID) %>%
#           select(id, cluster, chr), by="id", all.x=T) %>% 
#   arrange(match(id, tree_ed$tip.label))
identify(pp)



#### colours by ancestral blocks:

head(cluster_ID)

cluster_ID %>% 
  filter(anc_group=="Sk")

p2 <- facet_plot(pp, panel="", data=cluster_ID %>% 
                   filter(chr=="III"), 
                 geom=geom_point, aes(x=cluster, 
                                      #size=factor(anc_group),
                                      colour=factor(anc_group))) + 
  # alpha=ancHaplo, 
  # size=anc_group)) + 
  #color='darkred', 
  #alpha=0.8) +
  # scale_color_manual(values=c("darkred", "orange", "steelblue")) +
  # scale_alpha_manual(values=c(0,0.3,1,0.6,0.8,0.8)) + 
  #scale_size_manual(values=c(2,1.5,2)) +
  # scale_size_manual(values=c(2,1,2)) +
  theme_classic() +
  # scale_x_continuous(breaks = c(1,2,seq(3,max(cluster_ID$cluster),3)), 
  # labels = c(1,2,seq(3,max(cluster_ID$cluster),3))) +
  # scale_y_continuous(breaks = seq(1,800,30)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        panel.grid.major = element_line(colour = "gray80"), 
        axis.text.x=element_text(size=10, angle = 45, hjust = 1, colour="black"), 
        axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(), 
        # legend.position="none") 
        legend.position="top", 
        legend.title = element_blank())


p2



# generate some random values for each tip label in the data
# Make a second plot with the original, naming the new plot "dot", 
# using the data you just created, with a point geom.

p2 <- facet_plot(pp, panel="I", data=cluster_ID %>% filter(chr=="I"), 
                 geom=geom_point, aes(x=cluster, colour=anc_group, 
                                      alpha=anc_group, 
                                      size=anc_group))+ 
  #color='darkred', 
  #alpha=0.8) +
  scale_color_manual(values=c("darkred", "orange", "steelblue")) +
  scale_alpha_manual(values=c(0,0.3,1,0.5, 0.9, 0.9)) + 
  scale_size_manual(values=c(2,1.5,2)) +
  theme_classic() +
  scale_x_continuous(breaks = lables_tf, labels = lables_tf) +
  scale_y_continuous(breaks = seq(1,1800,30)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        panel.grid.major = element_line(colour = "gray80"), 
        axis.text.x=element_text(size=12, angle = 45, hjust = 1, colour="black"), 
        axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.position="none")
#panel.grid.minor = element_line(colour = "gray80"))
p2
clade_pop<-identify(p2)
clade_pop
p3 <- facet_plot(p2, panel="II", data=cluster_ID %>% filter(chr=="II"), 
                 geom=geom_point, aes(x=cluster, colour=anc_group, 
                                      alpha=anc_group, 
                                      size=anc_group))
#color='darkred', 
# alpha=0.8) 
p4 <- facet_plot(p3, panel="III", data=cluster_ID %>% filter(chr=="III"), 
                 geom=geom_point, aes(x=cluster, colour=anc_group, 
                                      alpha=anc_group, 
                                      size=anc_group))
#color='darkred', 
# alpha=0.8) 
p4 



library(ggtree)
library(reshape2)

# load the packaged
library(grid)
library(gtable)

head(cluster_ID)

panel_size<-cluster_ID %>% 
  filter(!(is.na(chr))) %>%
  group_by(chr) %>%
  summarise(max_cluster=max(cluster), 
            min_cluster=min(cluster)) %>% 
  mutate(size_ch=max_cluster-min_cluster) %>% 
  mutate(total=cluster_ID %>% 
           filter(!(is.na(chr))) %>% 
           select(cluster) %>% 
           unlist() %>% 
           as.vector() %>% 
           max()) %>% 
  mutate(size_block=0.6*size_ch/total)


gt = ggplot_gtable(ggplot_build(p4))
#gtable_show_layout(gt) # will show you the layout - very handy function
#gt # see plot layout in table format
#gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
#gt$layout$l[grep('panel-3', gt$layout$name)] # you want to find the column specific to panel-3
#gt$layout$l[grep('panel-4', gt$layout$name)] # you want to find the column specific to panel-4
gt$widths[7] = (unlist(panel_size[panel_size$chr=="I","size_block"]))*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
gt$widths[9] = (unlist(panel_size[panel_size$chr=="II","size_block"]))*gt$widths[9]
gt$widths[11] = (unlist(panel_size[panel_size$chr=="III","size_block"]))*gt$widths[11]

grid.draw(gt) # plot with grid draw

svg("05_tree_distribution_TF.svg",height = 180/20.4, width = 180/20.4, pointsize = 7)
grid.draw(gt) # plot with grid draw
dev.off()













## flanking Sp and Sk sequence:

clusters_noNA<-ordered_table2 %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>% 
  filter(!(is.na(ancHaplo))) %>%
  ungroup() %>%
  select(cluster) %>%
  unique() %>%
  unlist() %>%
  as.vector()

ordered_table2 %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>% 
  filter(cluster %in% clusters_noNA) %>%
  # dim()
  group_by(chr, cluster, sample) %>% 
  summarise(NumSeq_perSample=n(), 
            sp_seq=sum(ancHaplo==0, na.rm = TRUE), 
            sk_seq=sum(ancHaplo==1, na.rm = TRUE), 
            na_seq=sum(is.na(ancHaplo), na.rm = TRUE)) %>% 
  mutate(fq_sp_persample=sp_seq/NumSeq_perSample, 
         fq_sk_persample=sk_seq/NumSeq_perSample) %>%
  # plot distribution of anc haplotypes within clusers and samples:
  # ggplot(aes(fq_sp_persample-fq_sk_persample))+
  # geom_histogram()
  filter((fq_sp_persample-fq_sk_persample)!=0) %>%
  mutate(con_ancGroup=ifelse((fq_sp_persample-fq_sk_persample)<0, "Sk", "Sp")) %>% 
  # # heat map with ancestral group per cluster and sample:
  # ggplot(aes(factor(cluster), sample)) +
  # geom_tile(aes(fill=factor(con_ancGroup)))+
  # labs(x="Cluster",y="Sample", fill="Anc. Group")+
  # scale_fill_manual(values = c("steelblue", "darkred"))+
  # facet_grid(. ~ chr, space="free", scale="free") +
  # theme_classic()+
  # theme(axis.text.x = element_blank())+
  # ggsave("09_anc_group_flankingSeq.png", height = 6, width = 16)
  group_by(chr, cluster) %>% 
  summarise(Total_samples=n(), 
            Num_sp_seq=sum(con_ancGroup=="Sp", na.rm = TRUE), 
            Num_sk_seq=sum(con_ancGroup=="Sk", na.rm = TRUE), 
            Total_samples_group=ifelse(Total_samples==Num_sp_seq, "Sp", 
                                       ifelse(Total_samples==Num_sk_seq, 
                                              "Sk", "Mixed"))) %>% 
  gather("Group", "Num_samples", c(-chr,-cluster, -Total_samples_group)) %>%
  filter(Group!="Num_sk_seq") %>% 
  ggplot(aes(factor(cluster), Num_samples, group=cluster)) +
  geom_line()+
  geom_point(aes(factor(cluster), Num_samples, alpha=Group, colour=Total_samples_group)) +
  scale_alpha_manual(values = c(0,0.8)) +
  scale_color_manual(values = c("black", "steelblue", "darkred"))+
  scale_y_continuous(breaks = seq(0,40,4))+
  labs(x="Cluster",y="Num Samples", colour="Anc. Group")+
  facet_grid(. ~ chr, space="free", scale="free")+
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        legend.position = "bottom")#+
  # ggsave("10_distributionAncGroupPercluster.png", height = 4, width = 10)
  # tbl_df %>% print(n=30)
  # head()
  # filter(cluster==290) %>%
  # # filter(is.na(ancHaplo)) %>%
  # # head()
  # tbl_df %>% print(n=Inf)
  # ungroup() %>%
  # select(sample) %>%
  # unique() %>%
  # unlist() %>%
  # as.vector()



# 1dSFS and 2dSFS for SNP data:

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



unique_samples<-c("JB22_EBC2",
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




library("gdsfmt")
library("SNPRelate")
library("ggplot2")
library("VariantAnnotation")
library("ggrepel")
library("tidyr")
library("plyr")
library("dplyr")
library("GenomicRanges")
library("IRanges")


ancestralhap<-read.table("../Ancestral_blocks/ancestralhap_57ILL_LR.txt", T) 
head(ancestralhap)


vcf_file_name<-"variant_SNPs_all.LRSamles.noncoding.vcf.recode.vcf"
vcf_file_name<-"variant_SNPs_all.LRSamles.allSNPs.vcf.recode.vcf"

temporal_number<-runif(1,0,100)
vcf <- readVcf(vcf_file_name)
snpgdsVCF2GDS(vcf_file_name, paste0("temporal_", temporal_number, ".gds"), method="biallelic.only")
snpgdsSummary(paste0("temporal_", temporal_number, ".gds"))
genofile <- snpgdsOpen(paste0("temporal_", temporal_number, ".gds"))
genotypes <- t(snpgdsGetGeno(genofile))

colnames(genotypes)<-read.gdsn(index.gdsn(genofile, "sample.id"))

index.gdsn(genofile, "snp.chromosome")
unique(read.gdsn(index.gdsn(genofile, "snp.chromosome")))
as.data.frame(genotypes) %>%
  head()

genotypes_table<-as.data.frame(genotypes) %>%
  mutate(ID=paste0(read.gdsn(index.gdsn(genofile, "snp.chromosome")), 
                   "_",read.gdsn(index.gdsn(genofile, "snp.position"))), 
         chr=read.gdsn(index.gdsn(genofile, "snp.chromosome")), 
         pos=read.gdsn(index.gdsn(genofile, "snp.position"))) %>% 
  dplyr::select(-JB1174_EBC132) %>%
  rowwise() %>%
  dplyr::mutate(total = sum(c_across("DY34373":"JB953"))) %>%
  ungroup() %>%
  filter(total<31) %>%
  dplyr::select(-total) %>%
  gather("sample", "GT", 1:31) 

head(genotypes_table)

ancestralhap_LRsamples<-ancestralhap %>% 
  filter(sample %in% genotypes_table$sample) %>% 
  mutate(Nor_PC1_pol=Nor_PC1_pol_sim, 
         windowID=as.factor(paste0(chromosome_name,"_", start_pos,"_", end_ed))) 

SNP_list_genotypes<-genotypes_table %>% 
  dplyr::select(chr, pos) %>% 
  mutate(ID=paste0(chr, "_", pos)) %>% 
  unique() 


windowID_list<-ancestralhap_LRsamples %>% 
  dplyr::select(chromosome_name,start_pos,end_ed,windowID) %>% 
  unique() 

head(windowID_list)

windowID<-c()
for (i in seq(1,dim(SNP_list_genotypes)[1])){
  window_found<-windowID_list %>%
    filter(chromosome_name==SNP_list_genotypes$chr[i]) %>%
    filter(start_pos<=SNP_list_genotypes$pos[i]) %>%
    filter(end_ed>=SNP_list_genotypes$pos[i]) %>%
    pull(windowID) %>%
    as.vector()
  if (length(window_found)==1){
    windowID<-c(windowID, window_found)
  } else if (length(window_found)==0) {
    windowID<-c(windowID, "No_found") 
  } else if (length(window_found)>1) {
    windowID<-c(windowID, "Multiple_found") 
  }
}

data.frame(windowID) %>%
  dplyr::filter(windowID %in% c("No_found", "Multiple_found")) %>%
  dplyr::group_by(windowID) %>%
  dplyr::summarise(total=n())

# chr_used<-"II"
# isnps <- with(SNP_list_genotypes[SNP_list_genotypes$chr==chr_used,], IRanges(pos, width=1, names=ID))
# igenes <- with(windowID_list[windowID_list$chromosome_name==chr_used,], IRanges(start_pos, end_ed, names=windowID))
# olaps <- findOverlaps(isnps, igenes)

# genotypes_andGT<-cbind(SNP_list_genotypes[queryHits(olaps),], 
#                        windowID_list[subjectHits(olaps),]) %>% 
#   select(ID, windowID) %>% 
#   left_join(genotypes_table, by="ID")

# genotypes_andGT<-genotypes_andGT %>% 
#   merge(ancestralhap_LRsamples %>% 
#           select(windowID, sample, Nor_PC1_pol), 
#         by=c("windowID", "sample"))

head(genotypes_table)
genotypes_andGT<-genotypes_table %>%
  left_join(cbind(SNP_list_genotypes, windowID) %>%
              select(ID, windowID), by="ID") %>%
  filter(!(windowID %in% c("No_found", "Multiple_found"))) %>%
  left_join(ancestralhap_LRsamples %>% 
              select(windowID, sample, Nor_PC1_pol), 
            by=c("windowID", "sample")) 

snpgdsClose(genofile)

genotypes_all_andGT<-genotypes_andGT

# write.table(genotypes_all_andGT, "genotypes_all_andGT.noncoding.txt", quote = F, row.names = F, col.names = T)
genotypes_all_andGT<-read.table("genotypes_all_andGT.noncoding.txt",T)

head(genotypes_all_andGT)

# write.table(genotypes_all_andGT, "genotypes_all_andGT.allSNPs.txt", quote = F, row.names = F, col.names = T)
# genotypes_all_andGT<-read.table("genotypes_all_andGT.allSNPs.txt",T)

genotypes_all_andGT$GT1_in0<-(genotypes_all_andGT$GT==1 & genotypes_all_andGT$Nor_PC1_pol==0)*1
genotypes_all_andGT$GT1_in1<-(genotypes_all_andGT$GT==1 & genotypes_all_andGT$Nor_PC1_pol==1)*1
genotypes_all_andGT$GT0_in0<-(genotypes_all_andGT$GT==0 & genotypes_all_andGT$Nor_PC1_pol==0)*1
genotypes_all_andGT$GT0_in1<-(genotypes_all_andGT$GT==0 & genotypes_all_andGT$Nor_PC1_pol==1)*1

# genotypes_all_andGT %>%
#   filter(sample=="JB22_EBC2") %>%
#   head()
#   pull(sample) %>%
#   unique()
#   head()

genotypes_all_andGT_pop<-genotypes_all_andGT %>% 
  select(-windowID) %>% 
  group_by(chr, pos, ID) %>% 
  summarise(total = n(), NAnc0=sum(Nor_PC1_pol==0), NAnc1=sum(Nor_PC1_pol==1), 
            NAnc0_in0=sum(GT0_in0), NAnc0_in1=sum(GT0_in1), 
            NAnc1_in0=sum(GT1_in0), NAnc1_in1=sum(GT1_in1)) %>% 
  mutate(fq_NAnc0_in0=NAnc0_in0/NAnc0, fq_NAnc0_in1=NAnc0_in1/NAnc1, 
         fq_NAnc1_in0=NAnc1_in0/NAnc0, fq_NAnc1_in1=NAnc1_in1/NAnc1) 

genotypes_all_andGT_pop[is.na(genotypes_all_andGT_pop)] <- 0

genotypes_all_andGT_pop<-genotypes_all_andGT_pop %>% 
  mutate(folded_fq_pop0=(((fq_NAnc1_in0+fq_NAnc1_in1)>=1)*1*fq_NAnc0_in0)+(((fq_NAnc1_in0+fq_NAnc1_in1)<1)*1*(fq_NAnc1_in0)), 
         folded_fq_pop1=(((fq_NAnc1_in0+fq_NAnc1_in1)>=1)*1*fq_NAnc0_in1)+(((fq_NAnc1_in0+fq_NAnc1_in1)<1)*1*(fq_NAnc1_in1))) 

genotypes_all_andGT_pop_min5<-genotypes_all_andGT_pop %>% filter(NAnc0>=5 & NAnc1>=5)
genotypes_all_andGT_pop_min5$folded_fq_pop0_group<-0
genotypes_all_andGT_pop_min5$folded_fq_pop1_group<-0

for (fq in seq(0,1,0.1)) {
  print(fq)
  genotypes_all_andGT_pop_min5$folded_fq_pop0_group[genotypes_all_andGT_pop_min5$folded_fq_pop0>=fq & genotypes_all_andGT_pop_min5$folded_fq_pop0<fq+0.0999]<-fq
  genotypes_all_andGT_pop_min5$folded_fq_pop1_group[genotypes_all_andGT_pop_min5$folded_fq_pop1>=fq & genotypes_all_andGT_pop_min5$folded_fq_pop1<fq+0.0999]<-fq
}


folded_2dSFS_data<- genotypes_all_andGT_pop_min5 %>% 
  mutate(total_Var=dim(genotypes_all_andGT_pop_min5)[1]) %>% 
  select(folded_fq_pop0_group, folded_fq_pop1_group, total_Var) %>% 
  group_by(folded_fq_pop0_group, folded_fq_pop1_group) %>%
  summarise(sfs_2d = n()/total_Var[1]) 


folded_2dSFS_data$labels<-folded_2dSFS_data$sfs_2d
# folded_2dSFS_data$labels[folded_2dSFS_data$labels<0.001]<-0

folded_2dSFS_data<-folded_2dSFS_data %>%
  filter(labels>0.001)


ggplot(folded_2dSFS_data, aes(folded_fq_pop0_group, folded_fq_pop1_group)) + 
  geom_tile(aes(fill = sfs_2d*100)) + 
  geom_text(aes(label=format(round(labels*100,1), nsmall = 1)), size=5) + 
  #geom_text(aes(label= (formatC(sfs_2d*100, format = "e", digits = 1)))) + 
  scale_fill_gradientn(colours = rev(terrain.colors(10)), 
                       name="", 
                       #scale_fill_gradient(low = "white", high = "steelblue") + 
                       #scale_fill_gradientn(colours = c("white", "blue", "red", "green"),
                       #trans = "log", 
                       limits = c(0,100), 
                       breaks=c(0,20,50,100)) + 
  #limits = c(-7,4.6)) + 
  theme_classic() + 
  xlab("Allele Fq. Sp") + 
  ylab("Allele Fq. Sk") + 
  scale_x_continuous(breaks=seq(0,1,by=0.1)) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x=element_text(angle = 45, hjust = 1, size=15, colour="black"), 
        axis.text.y = element_text(colour="black", size=15), 
        axis.title=element_text(size=15),
        strip.text.x = element_blank(), 
        legend.text = element_text(colour="black", size = 14)) #+ 
  # ggsave(paste0("2dSFS_folded_SNPs_min5Samples.allSNPs.svg"),  width = 6, height = 6, dpi = 450)
  # ggsave(paste0("2dSFS_folded_SNPs_min5Samples.allSNPs.png"),  width = 6, height = 6, dpi = 450)
  # ggsave(paste0("2dSFS_folded_SNPs_min5Samples.noncoding.svg"),  width = 6, height = 6, dpi = 450)
  # ggsave(paste0("2dSFS_folded_SNPs_min5Samples.noncoding.png"),  width = 6, height = 6, dpi = 450)


getwd()


### 2dSFS for TEs:

#removing repeated strains there are 32 samples:
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


inconsistent_clusters<-ordered_table2 %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  mutate(start_end=(ancHaplo==ancHaplo_end)*1) %>% 
  dplyr::group_by(cluster, sample) %>%
  summarise(inconsistent_cluster=(sum((start_end==0)*1)!=0)*1) %>% 
  filter(!(is.na(inconsistent_cluster))) %>%
  group_by(cluster) %>%
  summarise(total_NSamples=n(), inconsistent_Samples=sum((inconsistent_cluster==1)*1)) %>%
  filter(inconsistent_Samples>1) %>%
  ungroup() %>%
  select(cluster) %>%
  unlist() %>% as.vector()


# ordered_table3<-ordered_table2 %>% 
#   filter(sample %in% unique_samples) %>% 
#   filter(!(is.na(ancHaplo))) %>%
#   filter(!(cluster %in% inconsistent_clusters)) %>%
#   mutate(genotype=1)

ordered_table3<-ordered_table2 %>% 
  filter(sample %in% unique_samples) %>%
  filter(!(is.na(ancHaplo_mean))) %>%
  # filter(!(cluster %in% inconsistent_clusters)) %>%
  mutate(genotype=1)

head(ordered_table3)

added_samples<-c()
added_cluster<-c()
added_chr<-c()
added_ancHaplo<-c()

for (tem_cluster in unique(ordered_table3$cluster)){
  # tem_cluster<-284
  tem_data<-ordered_table3 %>%
    filter(cluster==tem_cluster)
  tem_chr<-as.vector(tem_data$chr[1])
  tem_start<-min(tem_data$start)-50
  missing_samples<-as.vector(unique(ordered_table3$sample)[!(unique(ordered_table3$sample) %in% tem_data$sample)])
  # missing_samples
  # head(tem_data)
  # tem_sample<-"Pomberef"
  for (tem_sample in missing_samples){
    start_point<-ancestralhap[ancestralhap$sample==tem_sample & ancestralhap$chromosome_name==tem_chr & ancestralhap$start_pos<=(tem_start-50) & ancestralhap$end_ed>(tem_start+50),"Nor_PC1_pol_sim"]
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
  }
}

#added_ancHaplo[added_ancHaplo == 0.5 & !(is.na(added_ancHaplo))] <- NA


sfs_table<-ordered_table3 %>%
  ungroup() %>%
  select(sample, chr, cluster, ancHaplo_mean, genotype) %>%
  rbind(data.frame(sample=added_samples, 
                   chr=added_chr, 
                   cluster=added_cluster, 
                   ancHaplo_mean=added_ancHaplo, 
                   genotype=0)) %>%
  group_by(chr, cluster, sample, genotype) %>% 
  dplyr::summarise(NumSeq_perSample=n(), 
            sp_seq=sum(ancHaplo_mean==0, na.rm = TRUE), 
            sk_seq=sum(ancHaplo_mean==1, na.rm = TRUE), 
            na_seq=sum(is.na(ancHaplo_mean), na.rm = TRUE), 
            num_genotypes_perSample=length(unique(genotype))) %>% 
  # filter(num_genotypes_perSample==1) %>%
  mutate(fq_sp_persample=sp_seq/NumSeq_perSample, 
         fq_sk_persample=sk_seq/NumSeq_perSample) %>%
  # # plot distribution of anc haplotypes within clusers and samples:
  # ggplot(aes(fq_sp_persample-fq_sk_persample))+
  # geom_histogram()
  filter((fq_sp_persample-fq_sk_persample)!=0) %>%
  mutate(con_ancGroup=ifelse((fq_sp_persample-fq_sk_persample)<0, 
                             "Sk", "Sp")) %>% 
  select(-NumSeq_perSample, -sp_seq, 
         -sk_seq, -na_seq, 
         -num_genotypes_perSample, 
         -fq_sp_persample, -fq_sk_persample) %>%
  group_by(chr, cluster, con_ancGroup) %>% 
  summarise(num_samples=n(), 
            fq=sum(genotype)/num_samples) #%>%
  # # 1d SFS
  # # filter(num_samples>2)%>%
  # filter(fq>0)%>%
  # ggplot(aes(fq, fill=con_ancGroup))+
  # geom_histogram()+
  # # geom_histogram(aes(y=..density..))+
  # scale_fill_manual(values = c("steelblue", "darkred"))+
  # # scale_y_continuous(breaks = seq(0,40,4))+
  # labs(x="Frequency",y="Count")+
  # facet_grid(con_ancGroup ~ .)+
  # theme_classic()+
  # theme(legend.position = "none")#+
  # ggsave("11_1dSFS_AncGroupPercluster.png", height = 4, width = 5)
  # tbl_df %>% print(n=30)

head(sfs_table)

# sfs_table<-ordered_table3 %>%
#   # filter(len_froSeq>1500) %>%
#   ungroup() %>%
#   select(sample, chr, cluster, ancHaplo, genotype) %>%
#   rbind(data.frame(sample=added_samples, 
#                    chr=added_chr, 
#                    cluster=added_cluster, 
#                    ancHaplo=added_ancHaplo, 
#                    genotype=0)) %>%
#   # filter(!(is.na(ancHaplo))) %>%
#   group_by(chr, cluster, sample, genotype) %>% 
#   summarise(NumSeq_perSample=n(), 
#             sp_seq=sum(ancHaplo==0, na.rm = TRUE), 
#             sk_seq=sum(ancHaplo==1, na.rm = TRUE), 
#             na_seq=sum(is.na(ancHaplo), na.rm = TRUE), 
#             num_genotypes_perSample=length(unique(genotype))) %>% 
#   # filter(num_genotypes_perSample==1) %>%
#   mutate(fq_sp_persample=sp_seq/NumSeq_perSample, 
#          fq_sk_persample=sk_seq/NumSeq_perSample) %>%
#   # # plot distribution of anc haplotypes within clusers and samples:
#   # ggplot(aes(fq_sp_persample-fq_sk_persample))+
#   # geom_histogram()
#   filter((fq_sp_persample-fq_sk_persample)!=0) %>%
#   mutate(con_ancGroup=ifelse((fq_sp_persample-fq_sk_persample)<0, 
#                              "Sk", "Sp")) %>% 
#   select(-NumSeq_perSample, -sp_seq, 
#          -sk_seq, -na_seq, 
#          -num_genotypes_perSample, 
#          -fq_sp_persample, -fq_sk_persample) %>%
#   group_by(chr, cluster, con_ancGroup) %>% 
#   summarise(num_samples=n(), 
#             fq=sum(genotype)/num_samples) #%>%
#   # # 1d SFS
#   # # filter(num_samples>2)%>%
#   # filter(fq>0)%>%
#   # ggplot(aes(fq, fill=con_ancGroup))+
#   # geom_histogram()+
#   # # geom_histogram(aes(y=..density..))+
#   # scale_fill_manual(values = c("steelblue", "darkred"))+
#   # # scale_y_continuous(breaks = seq(0,40,4))+
#   # labs(x="Frequency",y="Count")+
#   # facet_grid(con_ancGroup ~ .)+
#   # theme_classic()+
#   # theme(legend.position = "none")+
#   # ggsave("11_1dSFS_AncGroupPercluster.png", height = 4, width = 5)
#   # tbl_df %>% print(n=30) 

head(sfs_table)

fq_sp<-c()
fq_sk<-c()
final_counts<-c()
for (tem_fq_sp in seq(0,1,0.1)) {
  for (tem_fq_sk in seq(0,1,0.1)) {
    counts<-sfs_table %>%
      select(chr, cluster, con_ancGroup, fq) %>%
      spread(con_ancGroup, fq, fill=0) %>%
      filter(Sk>=tem_fq_sk & Sk<tem_fq_sk+0.0999) %>%
      filter(Sp>=tem_fq_sp & Sp<tem_fq_sp+0.0999) %>%
      dim()
    fq_sp<-c(fq_sp, tem_fq_sp)
    fq_sk<-c(fq_sk, tem_fq_sk)
    final_counts<-c(final_counts, counts[1])
  }
}

sum(final_counts)
sfs_table2<-data.frame(fq_sp, fq_sk, counts=final_counts, total=sum(final_counts)) 

sfs_table2 %>%
  mutate(proportion_counts=counts*100/total) %>%
  filter(proportion_counts>0.001) %>% 
  ggplot(aes(fq_sp, fq_sk)) + 
  # geom_tile(aes(fill = proportion_counts)) +
  # geom_text(aes(label=ifelse(counts!=0, counts, "")), size=4) +
  geom_tile(aes(fill = proportion_counts)) +
  geom_text(aes(label=ifelse(proportion_counts!=0, format(round(proportion_counts,1), nsmall = 1),
                             "")), size=5) +
  scale_fill_gradientn(colours = rev(terrain.colors(10)),
                       name="", 
                       # #scale_fill_gradient(low = "white", high = "steelblue") +
                       # #scale_fill_gradientn(colours = c("white", "blue", "red", "green"),
                       # #trans = "log",
                       limits = c(0,100),
                       breaks=c(0,20,50,100)) +
  xlab("Allele Fq. Sp group") + 
  ylab("Allele Fq. Sk group") + 
  scale_x_continuous(breaks=seq(0,1,by=0.1)) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x=element_text(angle = 45, hjust = 1, size=15, colour="black"), 
        axis.text.y = element_text(colour="black", size=15), 
        axis.title=element_text(size=15),
        strip.text.x = element_blank(), 
        legend.text = element_text(colour="black", size = 14)) #+
  # ggsave("11_2dSFS_AncGroup.png", width = 6, height = 6, dpi = 450)
  # ggsave("11_2dSFS_AncGroup_propotions.png", width = 6, height = 6, dpi = 450)
  # ggsave("11_2dSFS_AncGroup_propotions.svg", width = 6, height = 6, dpi = 450)


### I filtered the original table removing clusters in which at least one of the two groups (Sp or Sk) have a sample size lower than 4
# this leaves 668 clusters (out of 688):

length(unique(sfs_table$cluster))

head(sfs_table)
min3_samplesPerGroup_clusters<-sfs_table %>% 
  select(chr, cluster, con_ancGroup, num_samples) %>%
  spread(con_ancGroup, num_samples, fill=0) %>%
  filter(Sp>3) %>%
  filter(Sk>3) %>% 
  ungroup() %>%
  select(cluster) %>%
  unlist() %>%
  as.vector()

head(min3_samplesPerGroup_clusters)
length(min3_samplesPerGroup_clusters)

sfs_table_min4<-sfs_table %>% 
  filter(cluster %in% min3_samplesPerGroup_clusters)

fq_sp<-c()
fq_sk<-c()
final_counts<-c()
for (tem_fq_sp in seq(0,1,0.1)) {
  for (tem_fq_sk in seq(0,1,0.1)) {
    counts<-sfs_table %>%
      select(chr, cluster, con_ancGroup, fq) %>%
      spread(con_ancGroup, fq, fill=0) %>%
      filter(Sk>=tem_fq_sk & Sk<tem_fq_sk+0.0999) %>%
      filter(Sp>=tem_fq_sp & Sp<tem_fq_sp+0.0999) %>%
      dim()
    fq_sp<-c(fq_sp, tem_fq_sp)
    fq_sk<-c(fq_sk, tem_fq_sk)
    final_counts<-c(final_counts, counts[1])
  }
}


sfs_table_min4_2<-data.frame(fq_sp, fq_sk, counts=final_counts, total=sum(final_counts)) 

sfs_table_min4_2 %>% 
  mutate(proportion_counts=counts*100/total) %>%
  filter(proportion_counts>0.001) %>% 
  ggplot(aes(fq_sp, fq_sk)) + 
  # geom_tile(aes(fill = proportion_counts)) +
  # geom_text(aes(label=ifelse(counts!=0, counts, "")), size=4) +
  geom_tile(aes(fill = proportion_counts)) +
  geom_text(aes(label=ifelse(proportion_counts!=0, format(round(proportion_counts,1), nsmall = 1),
                             "")), size=5) +
  scale_fill_gradientn(colours = rev(terrain.colors(10)),
                       name="", 
  # #scale_fill_gradient(low = "white", high = "steelblue") +
  # #scale_fill_gradientn(colours = c("white", "blue", "red", "green"),
  # #trans = "log",
  limits = c(0,100),
  breaks=c(0,20,50,100)) +
  xlab("Allele Fq. Sp group") + 
  ylab("Allele Fq. Sk group") + 
  scale_x_continuous(breaks=seq(0,1,by=0.1)) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x=element_text(angle = 45, hjust = 1, size=15, colour="black"), 
        axis.text.y = element_text(colour="black", size=15), 
        axis.title=element_text(size=15),
        strip.text.x = element_blank(), 
        legend.text = element_text(colour="black", size = 14)) #+
  # ggsave("11_2dSFS_AncGroup_min4.png", width = 6, height = 6, dpi = 450)
  # ggsave("11_2dSFS_AncGroup_propotions_min4.png", width = 6, height = 6, dpi = 450)
  # ggsave("11_2dSFS_AncGroup_propotions_min4.svg", width = 6, height = 6, dpi = 450)


###

## I did a third plot including only clusters where there is a least one sequence of more than 1500 bp. I wanted to know the frequency of large sequences relative to short sequence.
# After the review comments I changed it to 4500


clusters_withmin1500<-ordered_table3 %>%
  filter(len_froSeq>4500) %>% 
  # filter(len_froSeq>1500) %>% 
  select(cluster) %>%
  unique() %>%
  unlist() %>% as.vector()


sfs_table<-ordered_table3 %>%
  ungroup() %>%
  select(sample, chr, cluster, ancHaplo_mean, genotype) %>%
  rbind(data.frame(sample=added_samples, 
                   chr=added_chr, 
                   cluster=added_cluster, 
                   ancHaplo_mean=added_ancHaplo, 
                   genotype=0)) %>%
  filter(cluster %in% clusters_withmin1500) %>%
  group_by(chr, cluster, sample, genotype) %>% 
  summarise(NumSeq_perSample=n(), 
            sp_seq=sum(ancHaplo_mean==0, na.rm = TRUE), 
            sk_seq=sum(ancHaplo_mean==1, na.rm = TRUE), 
            na_seq=sum(is.na(ancHaplo_mean), na.rm = TRUE), 
            num_genotypes_perSample=length(unique(genotype))) %>% 
  mutate(fq_sp_persample=sp_seq/NumSeq_perSample, 
         fq_sk_persample=sk_seq/NumSeq_perSample) %>%
  # # plot distribution of anc haplotypes within clusers and samples:
  # ggplot(aes(fq_sp_persample-fq_sk_persample))+
  # geom_histogram()
  filter((fq_sp_persample-fq_sk_persample)!=0) %>%
  mutate(con_ancGroup=ifelse((fq_sp_persample-fq_sk_persample)<0, 
                             "Sk", "Sp")) %>% 
  select(-NumSeq_perSample, -sp_seq, 
         -sk_seq, -na_seq, 
         -num_genotypes_perSample, 
         -fq_sp_persample, -fq_sk_persample) %>%
  group_by(chr, cluster, con_ancGroup) %>% 
  summarise(num_samples=n(), 
            fq=sum(genotype)/num_samples)


fq_sp<-c()
fq_sk<-c()
final_counts<-c()

for (tem_fq_sp in seq(0,1,0.05)) {
  for (tem_fq_sk in seq(0,1,0.05)) {
    counts<-sfs_table %>%
      select(chr, cluster, con_ancGroup, fq) %>%
      spread(con_ancGroup, fq, fill=0) %>%
      filter(Sk>=tem_fq_sk & Sk<tem_fq_sk+0.04999) %>%
      filter(Sp>=tem_fq_sp & Sp<tem_fq_sp+0.04999) %>%
      dim()
    fq_sp<-c(fq_sp, tem_fq_sp)
    fq_sk<-c(fq_sk, tem_fq_sk)
    final_counts<-c(final_counts, counts[1])
  }
}

sfs_table2<-data.frame(fq_sp, fq_sk, counts=final_counts) 
sfs_table2 %>%
  ggplot(aes(fq_sp, fq_sk)) + 
  geom_tile(aes(fill = counts)) + 
  geom_text(aes(label=ifelse(counts!=0, counts, "")), size=4) + 
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
        axis.text.x=element_text(angle = 45, hjust = 1, size=15, colour="black"), 
        axis.text.y = element_text(colour="black", size=15), 
        axis.title=element_text(size=15),
        strip.text.x = element_blank(), 
        legend.text = element_text(colour="black", size = 14)) #+
# ggsave("11_2dSFS_AncGroup_clustersWithmin1_4500bp.svg", width = 6, height = 6, dpi = 400)
# ggsave("11_2dSFS_AncGroup_clustersWithmin1_4500bp.png", width = 6, height = 6, dpi = 400)
  # ggsave("11_2dSFS_AncGroup_clustersWithmin1_1500bp.svg", width = 6, height = 6, dpi = 400)
  # ggsave("11_2dSFS_AncGroup_clustersWithmin1_1500bp.png", width = 6, height = 6, dpi = 400)







### Full-length LTRs 2dSFS:

ordered_table3<-ordered_table2 %>% 
  filter(sample %in% unique_samples) %>%
  filter(!(is.na(ancHaplo_mean))) %>%
  # filter(len_froSeq>1500) %>%
  filter(len_froSeq>4500) %>%
  # filter(!(cluster %in% inconsistent_clusters)) %>%
  mutate(genotype=1)


added_samples<-c()
added_cluster<-c()
added_chr<-c()
added_ancHaplo<-c()

for (tem_cluster in unique(ordered_table3$cluster)){
  # tem_cluster<-284
  tem_data<-ordered_table3 %>%
    filter(cluster==tem_cluster)
  tem_chr<-as.vector(tem_data$chr[1])
  tem_start<-min(tem_data$start)-50
  missing_samples<-as.vector(unique(ordered_table3$sample)[!(unique(ordered_table3$sample) %in% tem_data$sample)])
  # missing_samples
  # head(tem_data)
  # tem_sample<-"Pomberef"
  for (tem_sample in missing_samples){
    start_point<-ancestralhap[ancestralhap$sample==tem_sample & ancestralhap$chromosome_name==tem_chr & ancestralhap$start_pos<=(tem_start-50) & ancestralhap$end_ed>(tem_start+50),"Nor_PC1_pol_sim"]
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
  }
}

sfs_table<-ordered_table3 %>% 
  ungroup() %>%
  select(sample, chr, cluster, ancHaplo_mean, genotype) %>%
  rbind(data.frame(sample=added_samples, 
                   chr=added_chr, 
                   cluster=added_cluster, 
                   ancHaplo_mean=added_ancHaplo, 
                   genotype=0)) %>%
  group_by(chr, cluster, sample, genotype) %>% 
  summarise(NumSeq_perSample=n(), 
            sp_seq=sum(ancHaplo_mean==0, na.rm = TRUE), 
            sk_seq=sum(ancHaplo_mean==1, na.rm = TRUE), 
            na_seq=sum(is.na(ancHaplo_mean), na.rm = TRUE), 
            num_genotypes_perSample=length(unique(genotype))) %>% 
  # filter(num_genotypes_perSample==1) %>%
  mutate(fq_sp_persample=sp_seq/NumSeq_perSample, 
         fq_sk_persample=sk_seq/NumSeq_perSample) %>%
  # # plot distribution of anc haplotypes within clusers and samples:
  # ggplot(aes(fq_sp_persample-fq_sk_persample))+
  # geom_histogram()
  filter((fq_sp_persample-fq_sk_persample)!=0) %>%
  mutate(con_ancGroup=ifelse((fq_sp_persample-fq_sk_persample)<0, 
                             "Sk", "Sp")) %>% 
  select(-NumSeq_perSample, -sp_seq, 
         -sk_seq, -na_seq, 
         -num_genotypes_perSample, 
         -fq_sp_persample, -fq_sk_persample) %>%
  group_by(chr, cluster, con_ancGroup) %>% 
  summarise(num_samples=n(), 
            fq=sum(genotype)/num_samples) 

fq_sp<-c()
fq_sk<-c()
final_counts<-c()
for (tem_fq_sp in seq(0,1,0.1)) {
  for (tem_fq_sk in seq(0,1,0.1)) {
    counts<-sfs_table %>%
      select(chr, cluster, con_ancGroup, fq) %>%
      spread(con_ancGroup, fq, fill=0) %>%
      filter(Sk>=tem_fq_sk & Sk<tem_fq_sk+0.0999) %>%
      filter(Sp>=tem_fq_sp & Sp<tem_fq_sp+0.0999) %>%
      dim()
    fq_sp<-c(fq_sp, tem_fq_sp)
    fq_sk<-c(fq_sk, tem_fq_sk)
    final_counts<-c(final_counts, counts[1])
  }
}


sfs_table_min4_2<-data.frame(fq_sp, fq_sk, counts=final_counts, total=sum(final_counts)) 

sfs_table_min4_2 %>% 
  mutate(proportion_counts=counts*100/total) %>%
  filter(proportion_counts>0.001) %>% 
  ggplot(aes(fq_sp, fq_sk)) + 
  # geom_tile(aes(fill = proportion_counts)) +
  # geom_text(aes(label=ifelse(counts!=0, counts, "")), size=4) +
  geom_tile(aes(fill = proportion_counts)) +
  geom_text(aes(label=ifelse(proportion_counts!=0, format(round(proportion_counts,1), nsmall = 1),
                             "")), size=5) +
  scale_fill_gradientn(colours = rev(terrain.colors(10)),
                       name="", 
                       # #scale_fill_gradient(low = "white", high = "steelblue") +
                       # #scale_fill_gradientn(colours = c("white", "blue", "red", "green"),
                       # #trans = "log",
                       limits = c(0,100),
                       breaks=c(0,20,50,100)) +
  xlab("Allele Fq. Sp group") + 
  ylab("Allele Fq. Sk group") + 
  scale_x_continuous(breaks=seq(0,1,by=0.1)) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x=element_text(angle = 45, hjust = 1, size=15, colour="black"), 
        axis.text.y = element_text(colour="black", size=15), 
        axis.title=element_text(size=15),
        strip.text.x = element_blank(), 
        legend.text = element_text(colour="black", size = 14)) #+
# ggsave("11_2dSFS_AncGroup_propotions_fulllengthLTR_4500bp.png", width = 6, height = 6, dpi = 450)
# ggsave("11_2dSFS_AncGroup_propotions_fulllengthLTR_4500bp.svg", width = 6, height = 6, dpi = 450)
  # ggsave("11_2dSFS_AncGroup_fulllengthLTR_1500bp.png", width = 6, height = 6, dpi = 450)
  # ggsave("11_2dSFS_AncGroup_propotions_fulllengthLTR_1500bp.png", width = 6, height = 6, dpi = 450)
  #ggsave("11_2dSFS_AncGroup_propotions_fulllengthLTR_1500bp.svg", width = 6, height = 6, dpi = 450)



# Folded SFS for TE:

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


ordered_table3<-ordered_table2 %>% 
  filter(sample %in% unique_samples) %>%
  filter(!(is.na(ancHaplo_mean))) %>%
  mutate(genotype=1)


added_samples<-c()
added_cluster<-c()
added_chr<-c()
added_ancHaplo<-c()

for (tem_cluster in unique(ordered_table3$cluster)){
  # tem_cluster<-284
  tem_data<-ordered_table3 %>%
    filter(cluster==tem_cluster)
  tem_chr<-as.vector(tem_data$chr[1])
  tem_start<-min(tem_data$start)-50
  missing_samples<-as.vector(unique(ordered_table3$sample)[!(unique(ordered_table3$sample) %in% tem_data$sample)])
  # missing_samples
  # head(tem_data)
  # tem_sample<-"Pomberef"
  for (tem_sample in missing_samples){
    start_point<-ancestralhap[ancestralhap$sample==tem_sample & ancestralhap$chromosome_name==tem_chr & ancestralhap$start_pos<=(tem_start-50) & ancestralhap$end_ed>(tem_start+50),"Nor_PC1_pol_sim"]
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
  }
}


# ordered_table2 %>% 
#   filter(sample %in% unique_samples) %>%
#   filter(!(is.na(ancHaplo_mean))) %>%
#   group_by(cluster, sample) %>%
#   summarise(seq_n=n()) %>% 
#   pull(seq_n) %>% 
#   unique()
#   head()

sfs_table<-ordered_table3 %>% 
  ungroup() %>%
  select(sample, chr, cluster, ancHaplo_mean, genotype) %>%
  rbind(data.frame(sample=added_samples, 
                   chr=added_chr, 
                   cluster=added_cluster, 
                   ancHaplo_mean=added_ancHaplo, 
                   genotype=0)) %>%
  group_by(chr, cluster, sample, genotype) %>% 
  summarise(NumSeq_perSample=n(), 
            sp_seq=sum(ancHaplo_mean==0, na.rm = TRUE), 
            sk_seq=sum(ancHaplo_mean==1, na.rm = TRUE), 
            na_seq=sum(is.na(ancHaplo_mean), na.rm = TRUE), 
            num_genotypes_perSample=length(unique(genotype))) %>% 
  # filter(num_genotypes_perSample==1) %>%
  mutate(fq_sp_persample=sp_seq/NumSeq_perSample, 
         fq_sk_persample=sk_seq/NumSeq_perSample) %>%
  # # plot distribution of anc haplotypes within clusers and samples:
  # ggplot(aes(fq_sp_persample-fq_sk_persample))+
  # geom_histogram()
  filter((fq_sp_persample-fq_sk_persample)!=0) %>%
  mutate(con_ancGroup=ifelse((fq_sp_persample-fq_sk_persample)<0, 
                             "Sk", "Sp")) %>% 
  select(-NumSeq_perSample, -sp_seq, 
         -sk_seq, -na_seq, 
         -num_genotypes_perSample, 
         -fq_sp_persample, -fq_sk_persample) %>%
  group_by(chr, cluster, con_ancGroup) %>% 
  summarise(num_samples=n(), 
            fq=sum(genotype)/num_samples) 

head(sfs_table)
ordered_table3 %>% 
  filter(cluster==9)
  head()

folded_sfs_table<-sfs_table %>%
  group_by(chr, cluster) %>%
  summarise(num_samples_sp=sum(((con_ancGroup=="Sp")*num_samples)), 
            num_samples_sk=sum(((con_ancGroup=="Sk")*num_samples)), 
            sum_fq=sum(fq), 
            folded_fq_sp=(((sum_fq)>=1)*1*(1-sum(((con_ancGroup=="Sp")*fq))))+(((sum_fq)<1)*1*(sum(((con_ancGroup=="Sp")*fq)))), 
            folded_fq_sk=(((sum_fq)>=1)*1*(1-sum(((con_ancGroup=="Sk")*fq))))+(((sum_fq)<1)*1*(sum(((con_ancGroup=="Sk")*fq))))) %>%
  #filter(num_samples_sp>0 & num_samples_sk>0) %>%
  #filter(sum_fq!=2) %>%
  filter(num_samples_sp>3) %>%
  filter(num_samples_sk>3) 


fq_sp<-c()
fq_sk<-c()
final_counts<-c()
for (tem_fq_sp in seq(0,1,0.1)) {
  for (tem_fq_sk in seq(0,1,0.1)) {
    counts<-folded_sfs_table %>%
      select(chr, cluster, folded_fq_sp, folded_fq_sk) %>%
      filter(folded_fq_sk>=tem_fq_sk & folded_fq_sk<tem_fq_sk+0.0999) %>%
      filter(folded_fq_sp>=tem_fq_sp & folded_fq_sp<tem_fq_sp+0.0999) %>%
      dim()
    fq_sp<-c(fq_sp, tem_fq_sp)
    fq_sk<-c(fq_sk, tem_fq_sk)
    final_counts<-c(final_counts, counts[1])
  }
}

sum(final_counts)
folded_sfs_table2<-data.frame(fq_sp, fq_sk, counts=final_counts, total=sum(final_counts)) 

folded_sfs_table2 %>%
  mutate(proportion_counts=counts*100/total) %>%
  filter(proportion_counts>0.001) %>% 
  ggplot(aes(fq_sp, fq_sk)) + 
  # geom_tile(aes(fill = proportion_counts)) +
  # geom_text(aes(label=ifelse(counts!=0, counts, "")), size=4) +
  geom_tile(aes(fill = proportion_counts)) +
  geom_text(aes(label=ifelse(proportion_counts!=0, format(round(proportion_counts,1), nsmall = 1),
  "")), size=5) +
  scale_fill_gradientn(colours = rev(terrain.colors(10)),
                       name="", 
                       # #scale_fill_gradient(low = "white", high = "steelblue") +
                       # #scale_fill_gradientn(colours = c("white", "blue", "red", "green"),
                       # #trans = "log",
                       limits = c(0,100),
                       breaks=c(0,20,50,100)) +
  xlab("Allele Fq. Sp group") + 
  ylab("Allele Fq. Sk group") + 
  scale_x_continuous(breaks=seq(0,1,by=0.1)) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x=element_text(angle = 45, hjust = 1, size=15, colour="black"), 
        axis.text.y = element_text(colour="black", size=15), 
        axis.title=element_text(size=15),
        strip.text.x = element_blank(), 
        legend.text = element_text(colour="black", size = 14)) #+
  # ggsave("11_folded_2dSFS_AncGroup_min4Samples.png", width = 6, height = 6, dpi = 450)
  # ggsave("11_folded_2dSFS_AncGroup_min4Samples_propotions.png", width = 6, height = 6, dpi = 450)
  # ggsave("11_folded_2dSFS_AncGroup_min4Samples_propotions.svg", width = 6, height = 6, dpi = 450)





# Folded SFS for full-length LTRs - min length 4500 kb:

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


ordered_table3<-ordered_table2 %>%
  filter(sample %in% unique_samples) %>%
  filter(!(is.na(ancHaplo_mean))) %>%
  filter(len_froSeq>4500) %>%
  # filter(len_froSeq>1500) %>%
  mutate(genotype=1)


added_samples<-c()
added_cluster<-c()
added_chr<-c()
added_ancHaplo<-c()

for (tem_cluster in unique(ordered_table3$cluster)){
  # tem_cluster<-284
  tem_data<-ordered_table3 %>%
    filter(cluster==tem_cluster)
  tem_chr<-as.vector(tem_data$chr[1])
  tem_start<-min(tem_data$start)-50
  missing_samples<-as.vector(unique(ordered_table3$sample)[!(unique(ordered_table3$sample) %in% tem_data$sample)])
  # missing_samples
  # head(tem_data)
  # tem_sample<-"Pomberef"
  for (tem_sample in missing_samples){
    start_point<-ancestralhap[ancestralhap$sample==tem_sample & ancestralhap$chromosome_name==tem_chr & ancestralhap$start_pos<=(tem_start-50) & ancestralhap$end_ed>(tem_start+50),"Nor_PC1_pol_sim"]
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
  }
}

# ordered_table2 %>% 
#   filter(sample %in% unique_samples) %>%
#   filter(!(is.na(ancHaplo_mean))) %>%
#   group_by(cluster, sample) %>%
#   summarise(seq_n=n()) %>% 
#   pull(seq_n) %>% 
#   unique()
#   head()

sfs_table<-ordered_table3 %>% 
  ungroup() %>%
  select(sample, chr, cluster, ancHaplo_mean, genotype) %>%
  rbind(data.frame(sample=added_samples, 
                   chr=added_chr, 
                   cluster=added_cluster, 
                   ancHaplo_mean=added_ancHaplo, 
                   genotype=0)) %>%
  group_by(chr, cluster, sample, genotype) %>% 
  summarise(NumSeq_perSample=n(), 
            sp_seq=sum(ancHaplo_mean==0, na.rm = TRUE), 
            sk_seq=sum(ancHaplo_mean==1, na.rm = TRUE), 
            na_seq=sum(is.na(ancHaplo_mean), na.rm = TRUE), 
            num_genotypes_perSample=length(unique(genotype))) %>% 
  # filter(num_genotypes_perSample==1) %>%
  mutate(fq_sp_persample=sp_seq/NumSeq_perSample, 
         fq_sk_persample=sk_seq/NumSeq_perSample) %>%
  # # plot distribution of anc haplotypes within clusers and samples:
  # ggplot(aes(fq_sp_persample-fq_sk_persample))+
  # geom_histogram()
  filter((fq_sp_persample-fq_sk_persample)!=0) %>%
  mutate(con_ancGroup=ifelse((fq_sp_persample-fq_sk_persample)<0, 
                             "Sk", "Sp")) %>% 
  select(-NumSeq_perSample, -sp_seq, 
         -sk_seq, -na_seq, 
         -num_genotypes_perSample, 
         -fq_sp_persample, -fq_sk_persample) %>%
  group_by(chr, cluster, con_ancGroup) %>% 
  summarise(num_samples=n(), 
            fq=sum(genotype)/num_samples) 

head(sfs_table)

folded_sfs_table<-sfs_table %>%
  group_by(chr, cluster) %>%
  summarise(num_samples_sp=sum(((con_ancGroup=="Sp")*num_samples)), 
            num_samples_sk=sum(((con_ancGroup=="Sk")*num_samples)), 
            sum_fq=sum(fq), 
            folded_fq_sp=(((sum_fq)>=1)*1*(1-sum(((con_ancGroup=="Sp")*fq))))+(((sum_fq)<1)*1*(sum(((con_ancGroup=="Sp")*fq)))), 
            folded_fq_sk=(((sum_fq)>=1)*1*(1-sum(((con_ancGroup=="Sk")*fq))))+(((sum_fq)<1)*1*(sum(((con_ancGroup=="Sk")*fq))))) %>%
  #filter(sum_fq!=2) %>%
  filter(num_samples_sp>3) %>%
  filter(num_samples_sk>3) 


fq_sp<-c()
fq_sk<-c()
final_counts<-c()
for (tem_fq_sp in seq(0,1,0.1)) {
  for (tem_fq_sk in seq(0,1,0.1)) {
    counts<-folded_sfs_table %>%
      select(chr, cluster, folded_fq_sp, folded_fq_sk) %>%
      filter(folded_fq_sk>=tem_fq_sk & folded_fq_sk<tem_fq_sk+0.0999) %>%
      filter(folded_fq_sp>=tem_fq_sp & folded_fq_sp<tem_fq_sp+0.0999) %>%
      dim()
    fq_sp<-c(fq_sp, tem_fq_sp)
    fq_sk<-c(fq_sk, tem_fq_sk)
    final_counts<-c(final_counts, counts[1])
  }
}

sum(final_counts)
folded_sfs_table2<-data.frame(fq_sp, fq_sk, counts=final_counts, total=sum(final_counts)) 

folded_sfs_table2 %>%
  mutate(proportion_counts=counts*100/total) %>%
  filter(proportion_counts>0.001) %>% 
  ggplot(aes(fq_sp, fq_sk)) + 
  # geom_tile(aes(fill = proportion_counts)) +
  # geom_text(aes(label=ifelse(counts!=0, counts, "")), size=4) +
  geom_tile(aes(fill = proportion_counts)) +
  geom_text(aes(label=ifelse(proportion_counts!=0, format(round(proportion_counts,1), nsmall = 1),
                             "")), size=5) +
  scale_fill_gradientn(colours = rev(terrain.colors(10)),
                       name="", 
                       # #scale_fill_gradient(low = "white", high = "steelblue") +
                       # #scale_fill_gradientn(colours = c("white", "blue", "red", "green"),
                       # #trans = "log",
                       limits = c(0,100),
                       breaks=c(0,20,50,100)) +
  xlab("Allele Fq. Sp group") + 
  ylab("Allele Fq. Sk group") + 
  scale_x_continuous(breaks=seq(0,1,by=0.1), limits = c(-0.05,1)) +
  scale_y_continuous(breaks=seq(0,1,by=0.1), limits = c(-0.05,1)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position="bottom", 
        axis.text.x=element_text(angle = 45, hjust = 1, size=15, colour="black"), 
        axis.text.y = element_text(colour="black", size=15), 
        axis.title=element_text(size=15),
        strip.text.x = element_blank(), 
        legend.text = element_text(colour="black", size = 14)) #+
# ggsave("11_folded_2dSFS_AncGroup_min4Samples_propotions_fulllengthLTR_4500bp.png", width = 6, height = 6, dpi = 450)
# ggsave("11_folded_2dSFS_AncGroup_min4Samples_propotions_fulllengthLTR_4500bp.svg", width = 6, height = 6, dpi = 450)
  # ggsave("11_folded_2dSFS_AncGroup_min4Samples_fulllengthLTR_1500bp.png", width = 6, height = 6, dpi = 450)
  # ggsave("11_folded_2dSFS_AncGroup_min4Samples_propotions_fulllengthLTR_1500bp.png", width = 6, height = 6, dpi = 450)
  # ggsave("11_folded_2dSFS_AncGroup_min4Samples_propotions_fulllengthLTR_1500bp.svg", width = 6, height = 6, dpi = 450)


####

# library(phytools)
# 
# 
# tree_ed@phylo$tip.label
# # test<-data.frame(id=tree_ed$tip.label) %>% 
# test<-data.frame(id=tree_ed@phylo$tip.label) %>% 
#   merge(ordered_table2 %>% 
#           mutate(id=seq_ID) %>%
#           select(id, cluster), by="id", all.x=T) %>% 
#   arrange(match(id, tree_ed@phylo$tip.label)) %>% 
#   mutate(cluster=ifelse(is.na(cluster),0,cluster)) 
# #filter(!(is.na(cluster)))
# test_matrix<-matrix(0,dim(test)[1],max(test$cluster))
# row.names(test_matrix)<-test$id
# colnames(test_matrix)<-seq(1,max(test$cluster))
# for(line in seq(1,dim(test)[1])){
#   test_matrix[as.vector(test$id[line]),test$cluster[line]]<-1
# }
# 
# X<-test_matrix
# 
# tree<-read.tree("../phylogenies/phy_tf_minLen1000_ed.treefile")
# 
# tree<-reorder(tree,"cladewise")
# X<-X[tree_ed@phylo$tip.label,]
# #plotTree(tree,plot=FALSE)
# #plotTree(tree)
# #svg("test.svg",height = 100/25.4, width = 180/25.4, pointsize = 7)
# 
# plotTree(tree,plot=FALSE, xlim=c(0,0.5))
# obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
# 
# plotTree(tree,lwd=1,ylim=c(0,obj$y.lim[2]*1.05),xlim=c(0,obj$x.lim[2]*1.2),
#          ftype="off")
# 
# obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
# h<-max(obj$xx)
# fsize<-0.6
# for(i in 1:Ntip(tree)){
#   lines(c(obj$xx[i],h),rep(obj$yy[i],2),lty="dotted", col="grey")
#   #text(h,obj$yy[i],tree$tip.label[i],cex=fsize,pos=4,font=3,offset=0.1)
# }
# 
# # s<-max(fsize*strwidth(tree$tip.label))
# # start.x<-1.05*h+s
# 
# start.x<-1.05*h
# cols<-setNames(c("white","blue"),0:1)
# step<-((par()$usr[2]-par()$usr[1])-start.x)/ncol(X)/2
# for(i in 1:ncol(X)){
#   # for(i in 1:10){
#   text(start.x,max(obj$yy)+1,paste("",colnames(X)[i]),pos=4,srt=90, cex=0.7,offset=0)
#   for(j in 1:nrow(X)){
#     xy<-c(start.x,obj$yy[j])
#     y<-c(xy[2]-0.5,xy[2]+0.5,xy[2]+0.5,xy[2]-0.5)
#     x<-c(xy[1]-step,xy[1]-step,xy[1]+step,xy[1]+step)
#     polygon(x,y,col=cols[as.character(X[j,i])], border = NA)
#     # asp<-(par()$usr[2]-par()$usr[1])/(par()$usr[4]-par()$usr[3])*
#     #   par()$pin[2]/par()$pin[1]
#     #x<-c(xy[1]-0.5*asp,xy[1]-0.5*asp,xy[1]+0.5*asp,xy[1]+0.5*asp)
#     #polygon(x,y,col=cols[as.character(X[j,i])])
#   }
#   #start.x<-start.x+asp
#   start.x<-start.x+(step*2)
# }
# dev.off()


## 2PK distance:
install.packages("remotes")
remotes::install_github("gjuggler/ggphylo")
library(ggphylo)
library(ggplot2)
library(plyr)
library(arm)
library(reshape2)
library(dplyr)

pp
sub <- tree.extract.clade(pp, clade_pop)
library('tidytree')
# Group clades
clades <- groupClade(pp, clade_pop)
head(clades)

clades$label[clades$group==1 & clades$isTip]

alignment_samples<-ape::read.FASTA("../Alignments/alig_all_tf_masked_conSeq_minLen1000_break.fasta", type = "DNA")

library("ape")

distance_ltr<-ape::dist.dna(alignment_samples, model = "K80", variance = FALSE,
                            gamma = FALSE, pairwise.deletion = FALSE,
                            base.freq = NULL, as.matrix = T)

head(distance_ltr,1)
head(distance_ltr)
distance_ltr[lower.tri(distance_ltr, diag = T)] <- NA
distance <- melt(distance_ltr) %>% 
  filter(!(is.na(value)))

sample_pop<-clades$label[clades$group==1 & clades$isTip]
# group_id<-
distance %>% head()
#filter(Dist>=0 & Dist<0.04) %>% 
filter(Var1==seq | Var2==seq) %>% 
  mutate(id_sample=ifelse(Var1==seq, as.vector(Var2), as.vector(Var1)),
         group_sample=ifelse(value>0.1,"sk","sp")) %>% 
  dplyr::select(id_sample, group_sample) 

head(group_id)

distance_groups<-distance %>% 
  merge(group_id %>% 
          mutate(Var1=id_sample) %>% 
          dplyr::select(-id_sample), by="Var1") %>% 
  merge(group_id %>% 
          mutate(Var2=id_sample) %>% 
          dplyr::select(-id_sample), by="Var2") %>% 
  mutate(comparison=factor(ifelse(group_sample.x=="sp" & group_sample.y=="sp", "sp_sp", 
                                  ifelse(group_sample.x=="sk" & group_sample.y=="sk", 
                                         "sk_sk", "sp_sk")), 
                           levels=c("sp_sp", "sk_sk", "sp_sk"))) 



test<-ape::read.FASTA("../Alignments/alig_all_tf_masked_conSeq_minLen1000_break.fasta", type = "DNA")
test
distance_ltr<-ape::dist.dna(test, model = "K80", variance = FALSE,
                            gamma = FALSE, pairwise.deletion = FALSE,
                            base.freq = NULL, as.matrix = T)

distance_ltr

distance_ltr[lower.tri(distance_ltr, diag = T)] <- NA
distance <- melt(distance_ltr) %>% 
  filter(!(is.na(value)))

seq="01_ref_I_1465320_10_4919_4915_1_complete_ltr"
group_id<-distance %>% 
  #filter(Dist>=0 & Dist<0.04) %>% 
  filter(Var1==seq | Var2==seq) %>% 
  mutate(id_sample=ifelse(Var1==seq, as.vector(Var2), as.vector(Var1)),
         group_sample=ifelse(value>0.1,"sk","sp")) %>% 
  dplyr::select(id_sample, group_sample) 

head(group_id)

distance_groups<-distance %>% 
  merge(group_id %>% 
          mutate(Var1=id_sample) %>% 
          dplyr::select(-id_sample), by="Var1") %>% 
  merge(group_id %>% 
          mutate(Var2=id_sample) %>% 
          dplyr::select(-id_sample), by="Var2") %>% 
  mutate(comparison=factor(ifelse(group_sample.x=="sp" & group_sample.y=="sp", "sp_sp", 
                                  ifelse(group_sample.x=="sk" & group_sample.y=="sk", 
                                         "sk_sk", "sp_sk")), 
                           levels=c("sp_sp", "sk_sk", "sp_sk"))) 


distance_groups %>% 
  ggplot(aes(value, fill=comparison), alpha=0.7) +
  # geom_histogram() +
  # geom_histogram(binwidth=1) +
  geom_histogram(aes(y=..density..), position="dodge", binwidth = 0.01) +
  scale_fill_manual(values=c("darkred", "steelblue", "orange")) +
  # facet_grid(comparison ~ .)
  # scale_x_continuous(breaks=seq(0,18,1), labels=seq(0,18,1)) +
  # ylab("Count") +
  ylab("Density") +
  xlab("Distance") +
  theme_classic() +
  theme(axis.text.x=element_text(size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14), 
        axis.title=element_text(colour="black", size=14)) #+
ggsave("06_2PK_distance_groups_ltr_completed.png")


distance_groups %>% 
  ggplot(aes(value)) +
  # geom_histogram() +
  # geom_histogram(binwidth=1) +
  geom_histogram(aes(y=..density..), binwidth = 0.01, colour="black") +
  #scale_fill_manual(values=c("darkred", "steelblue", "orange")) +
  # facet_grid(comparison ~ .)
  # scale_x_continuous(breaks=seq(0,18,1), labels=seq(0,18,1)) +
  # ylab("Count") +
  ylab("Density") +
  xlab("Distance") +
  theme_classic() +
  theme(axis.text.x=element_text(size=14, colour="black"), 
        axis.text.y=element_text(colour="black", size=14), 
        axis.title=element_text(colour="black", size=14)) +
  ggsave("06_2PK_distance_ltr_completed.png")










### Correlation between ancestral admixture and number of clusters with full-length LTRs:

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

getwd()
table_haplotypeID<-read.table("/Users/ru43sej/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Phylogenies/LTR_completed_allSeq/IQTreebb/table_haplotypesIDperClusterandSample.txt", T) %>% 
  mutate(start=pos) %>% 
  select(-pos)


head(table_haplotypeID)

sim_ID_seq<-c()
for (line in seq(1, dim(ordered_table2)[1])){
  sim_ID_seq<-c(sim_ID_seq, paste0(strsplit(ordered_table2$seq_ID[line], "_")[[1]][1], 
                                   "_", 
                                   ordered_table2$sample[line]))
}


# Including only Tf1 and Tf2 haplotypes:

table_byHaplotype_full_lengthLTR <- ordered_table2 %>% 
  cbind(sim_ID_seq) %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  left_join(table_haplotypeID %>% 
              select(sample, sim_ID_seq, haplotype_ID), 
            by=c("sample", "sim_ID_seq")) %>% 
  # left_join(table_haplotypeID, by=c("sample", "chr", "cluster", "start", "sim_ID_seq")) %>% 
  filter(sample %in% no_clonal_strains) %>%
  filter(len_froSeq>4.5) %>%
  filter(haplotype_ID %in% c("335", "521")) %>% 
  group_by(cluster, sample) %>% 
  summarise(N_seqperSample=n()) %>%
  group_by(sample) %>% 
  summarise(N_clusters=n()) 

head(table_byHaplotype_full_lengthLTR)

table_byHaplotype_full_lengthLTR<-rbind(table_byHaplotype_full_lengthLTR, 
                                   data.frame(sample=no_clonal_strains[!(no_clonal_strains %in% table_byHaplotype_full_lengthLTR$sample)],
                                              N_clusters=0))


statistical_test<-table_byHaplotype_full_lengthLTR %>% 
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop))  %>% 
  data.frame()

m1a<-lm(N_clusters~admixed_proportion, data=statistical_test)
result_stat<-summary(m1a)
result_stat

# library(olsrr)
# library(stargazer)
# stargazer(m1a, type = "text")
# ols_plot_resid_qq(m1a)
# ols_test_normality(m1a)
# plot(m1a)


library(ggpmisc)
table_byHaplotype_full_lengthLTR %>% 
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop)) %>%
  ggplot(aes(admixed_proportion, N_clusters))+
  geom_point(alpha=0.5, size=3)+
  geom_smooth(method='lm', formula= y~x, se=FALSE, colour="black")+
  annotate(geom="text", x=0.1, y=50, 
           label=paste0("Adj.R2 = ", 
                        format(round(result_stat$adj.r.squared, 2), nsmall = 2),
                        "\n",
                        "p-value = ",
                        format(round(result_stat$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="black", size=3)+
  # stat_poly_eq(formula = y~x,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE) +
  labs(x="Admixture proportion", y="Num. Clusters") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+  
ggsave("05_regressionAdmixtureProportionVsNclusters_335_&_521.png", width = 5, height = 5)
ggsave("05_regressionAdmixtureProportionVsNclusters_335_&_521.svg", width = 5, height = 5)






# Correlations by haplotype:


table_byHaplotype_full_lengthLTR <- ordered_table2 %>% 
  cbind(sim_ID_seq) %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  left_join(table_haplotypeID %>% 
              select(sample, sim_ID_seq, haplotype_ID), 
            by=c("sample", "sim_ID_seq")) %>% 
  # left_join(table_haplotypeID, by=c("sample", "chr", "cluster", "start", "sim_ID_seq")) %>% 
  filter(sample %in% no_clonal_strains) %>%
  filter(len_froSeq>4.5) %>%
  filter(haplotype_ID %in% c("335", "521")) %>% 
  group_by(cluster, sample, haplotype_ID) %>% 
  summarise(N_seqperSample=n()) %>%
  group_by(sample, haplotype_ID) %>% 
  summarise(N_clusters=n()) 

table_byHaplotype_full_lengthLTR 

head(singleton_table_byHaplotype)

table_byHaplotype_full_lengthLTR <-rbind(table_byHaplotype_full_lengthLTR , 
                                   data.frame(sample=no_clonal_strains[!(no_clonal_strains %in% table_byHaplotype_full_lengthLTR$sample[table_byHaplotype_full_lengthLTR$haplotype_ID==335])],
                                              haplotype_ID=335,
                                              N_clusters=0),
                                   data.frame(sample=no_clonal_strains[!(no_clonal_strains %in% table_byHaplotype_full_lengthLTR$sample[table_byHaplotype_full_lengthLTR$haplotype_ID==521])],
                                              haplotype_ID=521,
                                              N_clusters=0))

# statistical test per genotype (whole alpha or beta)

statistical_test<-table_byHaplotype_full_lengthLTR %>% 
  left_join(pro_AncPop, by="sample") %>%
  mutate(admixed_proportion=if_else(anc_prop<=0.5,
                                    anc_prop, 1-anc_prop)) %>%
  data.frame()

m1_335<-lm(N_clusters~admixed_proportion, 
           data=statistical_test[statistical_test$haplotype_ID==335,])
#Anova(m1_335, Type="III")
r_335<-summary(m1_335)
r_335

m1_521<-lm(N_clusters~admixed_proportion, 
           data=statistical_test[statistical_test$haplotype_ID==521,])
#Anova(m1beta, Type="III")
r_521<-summary(m1_521)
r_521


library(ggpmisc)
table_byHaplotype_full_lengthLTR %>% 
  ## Plot per Genotype: 
  left_join(pro_AncPop, by="sample") %>%
  mutate(admixed_proportion=if_else(anc_prop<=0.5,
                                    anc_prop, 1-anc_prop)) %>%
  ggplot(aes(admixed_proportion, N_clusters, colour=factor(haplotype_ID)))+
  geom_point(alpha=0.5, size=3)+
  geom_smooth(method='lm', formula= y~x, se=FALSE)+
  # stat_poly_eq(formula = y~x,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE)
  annotate(geom="text", x=0.05, y=35, 
           label=paste0("Adj.R2 = ", 
                        format(round(r_335$adj.r.squared, 2), nsmall = 2),
                        "\n", "p-value = ", 
                        format(round(r_335$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="#922B21", size=3)+
  annotate(geom="text", x=0.05, y=30, 
           label=paste0("Adj.R2 = ", 
                        format(round(r_521$adj.r.squared, 2), nsmall = 2),
                        "\n", "p-value = ", 
                        format(round(r_521$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="#1B4F72", size=3)+
  scale_color_manual(values=c("#922B21", "#1B4F72"))+
  # scale_color_manual(values=c("red", "blue"))+
  labs(x="Admixture proportion", 
       y="Num. Clusters") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        legend.position = "none") #+
# ggsave("05_regressionAdmixtureProportionVsNclustersbyHaplotype_TF1orTF2.png", width = 5, height = 5)
# ggsave("05_regressionAdmixtureProportionVsNclustersbyHaplotype_TF1orTF2.svg", width = 5, height = 5)




### Correlation between ancestral admixture and singletons:

ordered_table2 %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  head()

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

getwd()
table_haplotypeID<-read.table("/Users/ru43sej/Dropbox/uppsala/Repeats_pombe_TE_Wtf/Analyses/all_Samples/Phylogenies/LTR_completed_allSeq/IQTreebb/table_haplotypesIDperClusterandSample.txt", T) %>% 
  mutate(start=pos) %>% 
  select(-pos)


head(table_haplotypeID)

sim_ID_seq<-c()
for (line in seq(1, dim(ordered_table2)[1])){
  sim_ID_seq<-c(sim_ID_seq, paste0(strsplit(ordered_table2$seq_ID[line], "_")[[1]][1], 
                                   "_", 
                                   ordered_table2$sample[line]))
}




# Using all identified haplotypes:


table_haplotypeID %>% 
  group_by(haplotype_ID) %>%
  summarise(N_clusters_haplotype=n_distinct(cluster)) %>% 
  head()
    
  
ordered_table2 %>% 
  cbind(sim_ID_seq) %>% 
  head()

singleton_table_byHaplotype <-ordered_table2 %>% 
  cbind(sim_ID_seq) %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  # left_join(table_haplotypeID, by=c("sample", "chr", "cluster", "start", "sim_ID_seq")) %>% 
  left_join(table_haplotypeID %>% 
              select(sample, sim_ID_seq, haplotype_ID), 
            by=c("sample", "sim_ID_seq")) %>% 
  # left_join(table_haplotypeID %>% 
  #             group_by(haplotype_ID) %>%
  #             summarise(N_clusters_haplotype=n_distinct(cluster)) , 
  #           by="haplotype_ID") %>% 
  left_join(ordered_table2 %>% 
              cbind(sim_ID_seq) %>% 
              filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
              # left_join(table_haplotypeID, by=c("sample", "chr", "cluster", "start", "sim_ID_seq")) %>% 
              left_join(table_haplotypeID %>% 
                          select(sample, sim_ID_seq, haplotype_ID), 
                        by=c("sample", "sim_ID_seq")) %>% 
              group_by(haplotype_ID) %>%
              summarise(N_clusters_haplotype=n_distinct(cluster)), by="haplotype_ID") %>%
  # group_by(sample, haplotype_ID) %>%
  # summarise(N_seq=n()) %>%
  # filter(haplotype_ID!="NA")%>%
  # ggplot(aes(sample, N_seq, fill=factor(haplotype_ID)))+
  # geom_bar(stat = "identity")+
  # theme_classic() +
  # theme(panel.border = element_blank(), 
  #       panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank(), 
  #       axis.line = element_line(colour = "black"), 
  #       #legend.position="bottom", 
  #       axis.text.x=element_text(angle = 45, hjust = 1, size=15, colour="black"), 
  #       axis.text.y = element_text(colour="black", size=15)) 
  filter(sample %in% no_clonal_strains) %>%
  group_by(haplotype_ID, N_clusters_haplotype, cluster, sample) %>%
  summarise(N_seqperSample=n()) %>%
  #filter(N_seqperSample==1) %>%
  # filter(sample=="JB758")
  group_by(cluster) %>% 
  summarise(N_Samples=n_distinct(sample), 
            N_haplotypes=n_distinct(haplotype_ID), 
            Samples=paste0(sample, collapse = "", sep="/"), 
            Haplotypes=paste0(haplotype_ID, collapse = "", sep="/"), 
            N_clusters=paste0(N_clusters_haplotype, collapse = "", sep="/")) %>%
  filter(N_Samples==1) %>% 
  filter(N_haplotypes==1) %>% 
  mutate(sample=str_remove_all(Samples, "/"), 
         haplotype_ID=str_remove_all(Haplotypes, "/"), 
         n_clusters=str_remove_all(N_clusters, "/")) %>% 
  select(-Samples, -Haplotypes, -N_clusters) %>% 
  filter(haplotype_ID!="NA") %>%
  filter(as.numeric(n_clusters) > 2) %>% 
  group_by(sample) %>% 
  summarise(N_singletonCluster=n())

singleton_table_byHaplotype<-rbind(singleton_table_byHaplotype, 
             data.frame(sample=no_clonal_strains[!(no_clonal_strains %in% singleton_table_byHaplotype$sample)],
                                              N_singletonCluster=0))

statistical_test<-singleton_table_byHaplotype %>% 
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop))  %>% 
  data.frame()

m1a<-lm(N_singletonCluster~admixed_proportion, data=statistical_test)
result_stat<-summary(m1a)
result_stat


library(ggpmisc)
singleton_table_byHaplotype %>% 
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop)) %>%
  ggplot(aes(admixed_proportion, N_singletonCluster))+
  geom_point(alpha=0.5, size=3)+
  geom_smooth(method='lm', formula= y~x, se=FALSE, colour="black")+
  annotate(geom="text", x=0.1, y=7, 
           label=paste0("Adj.R2 = ", 
                        format(round(result_stat$adj.r.squared, 2), nsmall = 2),
                        "\n",
                        "p-value = ",
                        format(round(result_stat$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="black", size=4)+
  # stat_poly_eq(formula = y~x,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE) +
  labs(x="Admixture proportion", 
       y="Num. Singleton. Clusters") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+  
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsbyHaplotype.png", width = 5, height = 5)
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsbyHaplotype.svg", width = 5, height = 5)



# including only 335 and 521 haplotypes:

singleton_table_byHaplotype <- ordered_table2 %>% 
  cbind(sim_ID_seq) %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  left_join(table_haplotypeID %>% 
              select(sample, sim_ID_seq, haplotype_ID), 
            by=c("sample", "sim_ID_seq")) %>% 
  # left_join(table_haplotypeID, by=c("sample", "chr", "cluster", "start", "sim_ID_seq")) %>% 
  filter(sample %in% no_clonal_strains) %>%
  group_by(haplotype_ID, cluster, sample) %>% 
  summarise(N_seqperSample=n()) %>%
  filter(N_seqperSample==1) %>% 
  group_by(cluster) %>% 
  summarise(N_Samples=length(unique(sample)), 
            N_haplotypes=length(unique(haplotype_ID)), 
            Samples=paste0(sample, collapse = "", sep="/"), 
            Haplotypes=paste0(haplotype_ID, collapse = "", sep="/")) %>%
  filter(N_Samples==1) %>% 
  filter(N_haplotypes==1) %>% 
  #group_rows() %>%
  mutate(sample=str_remove_all(Samples, "/"), 
         haplotype_ID=str_remove_all(Haplotypes, "/")) %>% 
  select(-Samples, -Haplotypes) %>% 
  filter(haplotype_ID %in% c("335", "521")) %>% 
  group_by(sample) %>% 
  summarise(N_singletonCluster=n())

head(singleton_table_byHaplotype)

singleton_table_byHaplotype<-rbind(singleton_table_byHaplotype, 
     data.frame(sample=no_clonal_strains[!(no_clonal_strains %in% singleton_table_byHaplotype$sample)],
                N_singletonCluster=0))


statistical_test<-singleton_table_byHaplotype %>% 
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop))  %>% 
  data.frame()

m1a<-lm(N_singletonCluster~admixed_proportion, data=statistical_test)
result_stat<-summary(m1a)
result_stat

library(olsrr)
library(stargazer)
stargazer(m1a, type = "text")
ols_plot_resid_qq(m1a)
ols_test_normality(m1a)
plot(m1a)


library(ggpmisc)
singleton_table_byHaplotype %>% 
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop)) %>%
  ggplot(aes(admixed_proportion, N_singletonCluster))+
  geom_point(alpha=0.5, size=3)+
  geom_smooth(method='lm', formula= y~x, se=FALSE, colour="black")+
  annotate(geom="text", x=0.1, y=6, 
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
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsbyHaplotype_335_&_521.png", width = 5, height = 5)
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsbyHaplotype_335_&_521.svg", width = 5, height = 5)


singleton_table_byHaplotype_TEs<-singleton_table_byHaplotype
head(singleton_table_byHaplotype_TEs)

# Correlations by haplotype:

singleton_table_byHaplotype <-ordered_table2 %>% 
  cbind(sim_ID_seq) %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  left_join(table_haplotypeID %>% 
              select(sample, sim_ID_seq, haplotype_ID), 
            by=c("sample", "sim_ID_seq")) %>% 
  # left_join(table_haplotypeID, by=c("sample", "chr", "cluster", "start", "sim_ID_seq")) %>% 
  filter(sample %in% no_clonal_strains) %>%
  group_by(haplotype_ID, cluster, sample) %>% 
  summarise(N_seqperSample=n()) %>%
  filter(N_seqperSample==1) %>% 
  group_by(cluster) %>% 
  summarise(N_Samples=length(unique(sample)), 
            N_haplotypes=length(unique(haplotype_ID)), 
            Samples=paste0(sample, collapse = "", sep="/"), 
            Haplotypes=paste0(haplotype_ID, collapse = "", sep="/")) %>%
  filter(N_Samples==1) %>% 
  filter(N_haplotypes==1) %>% 
  #group_rows() %>%
  mutate(sample=str_remove_all(Samples, "/"), 
         haplotype_ID=str_remove_all(Haplotypes, "/")) %>% 
  select(-Samples, -Haplotypes) %>% 
  filter(haplotype_ID %in% c("335", "521")) %>% 
  group_by(sample, haplotype_ID) %>% 
  summarise(N_singletonCluster=n()) 

head(singleton_table_byHaplotype)

singleton_table_byHaplotype<-rbind(singleton_table_byHaplotype, 
                                   data.frame(sample=no_clonal_strains[!(no_clonal_strains %in% singleton_table_byHaplotype$sample[singleton_table_byHaplotype$haplotype_ID=="335"])],
                                              haplotype_ID="335",
                                              N_singletonCluster=0),
                                   data.frame(sample=no_clonal_strains[!(no_clonal_strains %in% singleton_table_byHaplotype$sample[singleton_table_byHaplotype$haplotype_ID=="521"])],
                                              haplotype_ID="521",
                                              N_singletonCluster=0))

# statistical test per genotype (whole alpha or beta)

statistical_test<-singleton_table_byHaplotype %>% 
  left_join(pro_AncPop, by="sample") %>%
  mutate(admixed_proportion=if_else(anc_prop<=0.5,
                                    anc_prop, 1-anc_prop)) %>%
  data.frame()

m1_335<-lm(N_singletonCluster~admixed_proportion, 
           data=statistical_test[statistical_test$haplotype_ID==335,])
#Anova(m1_335, Type="III")
r_335<-summary(m1_335)
r_335

m1_521<-lm(N_singletonCluster~admixed_proportion, 
           data=statistical_test[statistical_test$haplotype_ID==521,])
#Anova(m1beta, Type="III")
r_521<-summary(m1_521)
r_521

# m1_521<-lm(N_singletonCluster~admixed_proportion, 
#            data=statistical_test[statistical_test$haplotype_ID==521 &
#                                    statistical_test$N_singletonCluster!=0,])
# #Anova(m1beta, Type="III")
# r_521<-summary(m1_521)
# r_521


library(ggpmisc)
singleton_table_byHaplotype %>% 
  ## Plot per Genotype: 
  left_join(pro_AncPop, by="sample") %>%
  mutate(admixed_proportion=if_else(anc_prop<=0.5,
                                    anc_prop, 1-anc_prop)) %>%
  ggplot(aes(admixed_proportion, N_singletonCluster, colour=factor(haplotype_ID)))+
  geom_point(alpha=0.5, size=3)+
  geom_smooth(method='lm', formula= y~x, se=FALSE)+
  # stat_poly_eq(formula = y~x,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE)
  annotate(geom="text", x=0.05, y=15, 
           label=paste0("Adj.R2 = ", 
                        format(round(r_335$adj.r.squared, 2), nsmall = 2),
                        "\n", "p-value = ", 
                        format(round(r_335$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="#922B21", size=3)+
  annotate(geom="text", x=0.05, y=12, 
           label=paste0("Adj.R2 = ", 
                        format(round(r_521$adj.r.squared, 2), nsmall = 2),
                        "\n", "p-value = ", 
                        format(round(r_521$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="#1B4F72", size=3)+
  scale_color_manual(values=c("#922B21", "#1B4F72"))+
  # scale_color_manual(values=c("red", "blue"))+
  labs(x="Admixture proportion", 
       y="Num. Singleton. Clusters") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"), 
        legend.position = "none") #+
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsbyHaplotype_TF1orTF2.png", width = 5, height = 5)
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsbyHaplotype_TF1orTF2.svg", width = 5, height = 5)



### Correlation between singleton SNPs and Ancestral admixture

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


length(no_clonal_strains)

num_samples_SNPs_noncoding<-genotypes_all_andGT %>% 
  #select(sample,ID, chr, pos,GT) %>% 
  filter(sample %in% no_clonal_strains) %>%
  group_by(ID, chr, pos) %>%
  summarise(n_lines=n(),
            n_sample=length(unique(sample)), 
            n_genotypes=length(unique(GT)), 
            total_genotype=n_sample-sum(GT)) %>%
  # filter(n_lines!=n_sample) 
  mutate(folded_total_genotypes=ifelse(total_genotype>(n_sample/2), 
                                       n_sample-total_genotype, 
                                       total_genotype), 
         folded_tf_genotypes=ifelse(total_genotype>(n_sample/2), 
                                       1, 
                                       0)) %>%
  ungroup() %>%
  select(ID, folded_total_genotypes, folded_tf_genotypes) 
  # filter(folded_total_genotypes==1) %>%
  # pull(ID)

table_singleton_SNPs<-genotypes_all_andGT %>% 
  filter(sample %in% no_clonal_strains) %>%
  left_join(num_samples_SNPs_noncoding, by="ID") %>%
  filter(folded_total_genotypes==1) %>%
  mutate(folded_genotype=ifelse(folded_tf_genotypes==0, GT, abs(1-GT))) %>%
  filter(folded_genotype==1) %>%
  group_by(sample) %>%
  summarise(num_singletons=n()) %>%
  # filter(num_singletons>7200) %>%
  left_join(pro_AncPop, by="sample") %>%
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop))

m1a<-lm(num_singletons~admixed_proportion, data=table_singleton_SNPs)
result_stat<-summary(m1a)
result_stat

library(ggpmisc)
table_singleton_SNPs %>%
  # filter(num_singletons>7200) %>%
  ggplot(aes(admixed_proportion, num_singletons))+
  geom_point(alpha=0.5, size=3)+
  geom_smooth(method='lm', formula= y~x, se=FALSE, colour="black")+
  annotate(geom="text", x=0.4, y=7400, 
           label=paste0("Adj.R2 = ", 
                        format(round(result_stat$adj.r.squared, 2), nsmall = 2),
                        "\n",
                        "p-value = ",
                        format(round(result_stat$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="black", size=4)+
  labs(x="Admixture proportion", 
       y="Num. Singleton. SNPs") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+  
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsNonCodingSNPs.png", width = 5, height = 5)
  # ggsave("05_regressionAdmixtureProportionVsNSingletonsNonCodingSNPs.svg", width = 5, height = 5)




#### Correlation between TE and SNS singletons 

table_singleton_SNPs_TEs <- table_singleton_SNPs %>%
  filter(num_singletons>7200) %>%
  left_join(singleton_table_byHaplotype_TEs, by="sample") 


m1a<-lm(num_singletons~N_singletonCluster, data=table_singleton_SNPs_TEs)
result_stat<-summary(m1a)
result_stat

result_stat$coefficients
format(round(result_stat$adj.r.squared, 2), nsmall = 2)
format(round(result_stat$coefficients["N_singletonCluster","Pr(>|t|)"], 3), nsmall = 3)

table_singleton_SNPs_TEs %>%
  filter(num_singletons>7200) %>%
  ggplot(aes(num_singletons, N_singletonCluster))+
  geom_point(alpha=0.5, size=3)+
  geom_smooth(method='lm', formula= y~x, se=FALSE, colour="black")+
  annotate(geom="text", x=6800, y=16,
           label=paste0("Adj.R2 = ",
                        format(round(result_stat$adj.r.squared, 2), nsmall = 2),
                        "\n",
                        "p-value = ",
                        format(round(result_stat$coefficients["N_singletonCluster","Pr(>|t|)"], 3), nsmall = 3)),
           color="black", size=4)+
  labs(y="Num. Singleton. Clusters", 
       x="Num. Singleton. SNPs") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+  
# ggsave("05_regressionNSingletonsNonCodingSNPsVsTEClusters.png", width = 5, height = 5)
# ggsave("05_regressionNSingletonsNonCodingSNPsVsTEClusters.svg", width = 5, height = 5)



## Correlations between total TF1 and TF2 sequences Vs Anc admixture

# Including only Tf1 and Tf2 haplotypes:

table_byHaplotype_full_lengthLTR <-ordered_table2 %>% 
  cbind(sim_ID_seq) %>%
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  left_join(table_haplotypeID %>% 
              select(sample, sim_ID_seq, haplotype_ID), 
            by=c("sample", "sim_ID_seq")) %>% 
  # left_join(table_haplotypeID, by=c("sample", "chr", "cluster", "start", "sim_ID_seq")) %>% 
  filter(sample %in% no_clonal_strains) %>%
  filter(len_froSeq>4.5) %>%
  filter(haplotype_ID %in% c("335", "521")) %>% 
  group_by(sample) %>% 
  summarise(N_seqperSample=n()) 
  
head(table_byHaplotype_full_lengthLTR)

table_byHaplotype_full_lengthLTR<-rbind(table_byHaplotype_full_lengthLTR, 
                                        data.frame(sample=no_clonal_strains[!(no_clonal_strains %in% table_byHaplotype_full_lengthLTR$sample)],
                                                   N_seqperSample=0))

statistical_test<-table_byHaplotype_full_lengthLTR %>% 
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop))  %>% 
  data.frame()

m1a<-lm(N_seqperSample~admixed_proportion, data=statistical_test)
result_stat<-summary(m1a)
result_stat

# library(olsrr)
# library(stargazer)
# stargazer(m1a, type = "text")
# ols_plot_resid_qq(m1a)
# ols_test_normality(m1a)
# plot(m1a)


library(ggpmisc)
table_byHaplotype_full_lengthLTR %>% 
  left_join(pro_AncPop, by="sample") %>% 
  mutate(admixed_proportion=if_else(anc_prop<=0.5, 
                                    anc_prop, 1-anc_prop)) %>%
  ggplot(aes(admixed_proportion, N_seqperSample))+
  geom_point(alpha=0.5, size=3)+
  geom_smooth(method='lm', formula= y~x, se=FALSE, colour="black")+
  annotate(geom="text", x=0.1, y=50, 
           label=paste0("Adj.R2 = ", 
                        format(round(result_stat$adj.r.squared, 2), nsmall = 2),
                        "\n",
                        "p-value = ",
                        format(round(result_stat$coefficients["admixed_proportion","Pr(>|t|)"], 3), nsmall = 3)),
           color="black", size=3)+
  # stat_poly_eq(formula = y~x,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE) +
  labs(x="Admixture proportion", y="Num. Full-length LTRs") + 
  theme_classic()+
  theme(axis.text.x = element_text(size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black")) #+  
# ggsave("05_regressionAdmixtureProportionVsNSeq_335_&_521.png", width = 5, height = 5)
# ggsave("05_regressionAdmixtureProportionVsNSeq_335_&_521.svg", width = 5, height = 5)





###



## Comparison between previous analyses and EDTA:

order_samples_pro_AncPop_ed <-order_samples_pro_AncPop
order_samples_pro_AncPop_ed[order_samples_pro_AncPop_ed=="JB1174_EBC132"]<-"JB900_EBC132"

summary_table_ori <- read.table("all_summary_LTRSeq_ori.txt", F) 
names(summary_table_ori) <- c("sample", "Nseq", "TotalSeqLen")
summary_table_ori <- summary_table_ori %>%
  mutate(dataset = "Extended")

summary_table_EDTA <- read.table("summary_LTRSeq_EDTA.txt", F) 
names(summary_table_EDTA) <- c("sample", "Nseq", "TotalSeqLen")
summary_table_EDTA <- summary_table_EDTA %>%
  mutate(dataset = "EDTA")

summary_table_SeqIdentifaction <- rbind(summary_table_ori, summary_table_EDTA)

summary_table_SeqIdentifaction$sample <- factor(summary_table_SeqIdentifaction$sample, levels = order_samples_pro_AncPop_ed)


summary_table_SeqIdentifaction %>%
  head()

ggplot(summary_table_SeqIdentifaction, aes(sample, Nseq, fill=dataset)) +
  geom_bar(stat = "identity", position=position_dodge())+
  scale_fill_grey()+
  scale_y_continuous(breaks = seq(0,500,20))+
  labs(x="Sample", y="Number of Seq.") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black")) #+

# ggsave("12_Comparison_EDTA_NSeq.png", width = 12, height = 5)


ggplot(summary_table_SeqIdentifaction, aes(sample, TotalSeqLen, fill=dataset)) +
  geom_bar(stat = "identity", position=position_dodge())+
  scale_fill_grey()+
  labs(x="Sample", y="Total Seq. Length") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black")) 

# ggsave("12_Comparison_EDTA_TotalSeqLength.png", width = 12, height = 5)

#


summary_table_SeqIdentifaction %>%
  select(-TotalSeqLen) %>%
  spread(dataset, Nseq) %>% 
  mutate(diif_Nseq_Extended_EDTA=Extended-EDTA) %>%
  ggplot(aes(sample, diif_Nseq_Extended_EDTA)) +
  geom_bar(stat = "identity")+
  scale_fill_grey()+
  #scale_y_continuous(breaks = seq(0,500,20))+
  labs(x="Sample", y="Diff. Num. Seq. Extended-EDTA") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black"))



summary_table_SeqIdentifaction %>%
  select(-Nseq) %>%
  spread(dataset, TotalSeqLen) %>% 
  mutate(diif_Nseq_Extended_EDTA=Extended-EDTA) %>%
  ggplot(aes(sample, diif_Nseq_Extended_EDTA)) +
  geom_bar(stat = "identity")+
  scale_fill_grey()+
  #scale_y_continuous(breaks = seq(0,500,20))+
  labs(x="Sample", y="Diff. TotalSeqLen. Extended-EDTA") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black"))


##

# sequence overlapping:

sum_table_Allcomparison_OriAnalysisSeq_EDTA <- read.table("sum_table_Allcomparison_OriAnalysisSeq_EDTA.txt", F) 
names(sum_table_Allcomparison_OriAnalysisSeq_EDTA) <- c("sample", "intersect", "NSeq", "dataset")

sum_table_Allcomparison_OriAnalysisSeq_EDTA$sample <- factor(sum_table_Allcomparison_OriAnalysisSeq_EDTA$sample, levels = order_samples_pro_AncPop_ed)


sum_table_Allcomparison_OriAnalysisSeq_EDTA %>%
  ggplot(aes(intersect, NSeq,fill=dataset))+
  geom_bar(stat = "identity", position=position_dodge())+
  scale_y_continuous(breaks = seq(0,500,40))+
  scale_fill_grey()+
  facet_wrap(~sample, scale="free")+
  # facet_grid(. ~ sample, scale="free")+
  theme_classic() +
  theme(axis.text.x = element_text(size=9, colour="black"), 
        axis.text.y = element_text(size=9, colour="black"),
        legend.position = "bottom")

# ggsave("12_Comparison_EDTA_OriAnalysisSeq_intersect.png", width = 12, height = 12)





### Length of clusters in original de-novo genomes Vs Reference genome

samples_table_ori<-read.table("all_TE_annotation_refCoor_Masked_conSeq.bed", F) 

samples_table_EDTA_missingSeq<-read.table("EDTA_all_missingSeq_TE_annotation__refCoor_Masked_conSeq.bed", F, sep = "\t") 

names(samples_table_ori)<-c("sample", "or_chr",	"or_start",	"or_end",	"seq_ID",	"ref_start_chr",	"ref_start",	"dis_start",	"ref_end_chr",	"ref_end",	"dis_end",	"length_ori_seq",	"direction",	"concordant_chr",	"chr_start_ed",	"pos_start_ed",	"dis_strat_ed", "chr_end_ed",	"pos_end_ed",	"dis_end_ed")


names(samples_table_EDTA_missingSeq)<-c("sample", "or_chr",	"or_start",	"or_end",	"seq_ID",	"ref_start_chr",	"ref_start",	"dis_start",	"ref_end_chr",	"ref_end",	"dis_end",	"length_ori_seq",	"direction",	"concordant_chr",	"chr_start_ed",	"pos_start_ed",	"dis_strat_ed", "chr_end_ed",	"pos_end_ed",	"dis_end_ed")

samples_table_ori<-rbind(samples_table_ori, 
                     samples_table_EDTA_missingSeq) %>%
  arrange(sample, chr_start_ed, pos_start_ed, dis_strat_ed)



samples_table_ori %>% 
  left_join(ordered_table2 %>% 
              filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
              select(seq_ID, cluster), by="seq_ID") %>%
  # filter(cluster==9)
  group_by(cluster, sample, or_chr) %>% 
  summarise(min_or_start=min(or_start), 
            max_or_start=max(or_start), 
            min_or_end=min(or_end), 
            max_or_end=max(or_end), 
            total_len=sum(length_ori_seq), 
            N_seq=n()) %>%
  mutate(min_start=ifelse(min_or_start<=min_or_end, min_or_start, min_or_end), 
         max_end=ifelse(max_or_start>=max_or_end, max_or_start, max_or_end)) %>%
  mutate(total_dis=max_end-min_start) %>%
  filter(N_seq!=1) %>%
  head()
  filter(cluster==11)
  head()

ordered_table2 %>% 
  filter(!(sample %in% c("LTR_loc", "ref_loc"))) %>%
  select(seq_ID, cluster) %>%
  head()






####### TF without LTRs #####

samples_table<-read.table("table_list_samples_tfwithoutltr_masked_conSeq.txt", F)
head(samples_table)
names(samples_table)<-c("seq_ID", "sample", "chr", "start", "dis", "len_froSeq", "addID", "dir", "len", "rep")

ordered_table2 -samples_table %>% 
  merge(ID_table, by="sample", all.x=T) %>% 
  arrange(chr, start, sample, dis, rep) %>% 
  filter(seq_ID!="NA") %>% 
  filter(!(is.na(chr)))


     cluster<-c(1)
     length_cluster<-0
     min_start<-0
     for (i in seq(2,dim(ordered_table2 [1])){
       if(length_cluster==0){
         length_cluster<-ordered_table2 len_froSeq[i-1]
         min_start<-ordered_table2 start[i-1]
       }
       if(ordered_table2 chr[i]==ordered_table2 chr[i-1]){
         if ((ordered_table2 start[i]-(min_start+length_cluster))<=(-50)){
           print(i)
           cluster[i]<-cluster[i-1]
           length_cluster<-max(length_cluster,ordered_table2 len_froSeq[i])
           min_start<-min(min_start, ordered_table2 start[i])
         } else {
           cluster[i]<-cluster[i-1]+1
           length_cluster<-ordered_table2 len_froSeq[i]
           min_start<-ordered_table2 start[i]
         }
       } else {
         cluster[i]<-cluster[i-1]+1
         length_cluster<-ordered_table2 len_froSeq[i]
         min_start<-ordered_table2 start[i]
       }
     }
     
     
     ordered_table2 -cbind(ordered_table2  cluster) %>% 
       mutate(sample = fct_reorder(sample, anc_prop))
     
     head(ordered_table2 
          
          ordered_table2 %>% 
            filter(seq_ID!="NA") %>% 
            group_by(cluster) %>% 
            summarize(chr=chr[1],
                      start=min(start),
                      total_seq=n(), 
                      total_samples=length(unique(sample)), 
                      min_length=min(len), 
                      max_length=max(len)) 
          
          
          ordered_table2 %>% 
            filter(seq_ID!="NA") %>% 
            filter(cluster==200)
          
          ordered_table2 %>% 
            filter(seq_ID!="NA") %>% 
            merge((ordered_table2 %>% 
                   filter(seq_ID!="NA") %>% 
                   group_by(cluster) %>% 
                   summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
            mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
            write.table("tfWithoutLTR_new_sorted_table_list_samples_masked_conSeq.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
          
          library("ggtree")
          #tree<-read.tree("../phylogenies/phy_tfWithoutLTR_minLen1000")
          tree_ed<-read.iqtree("../phylogenies/phy_tfWithoutLTR_minLen1000.treefile")
          #tree_ed<- as.phylo(tree)
          fig_tree<- ggtree(tree_ed) + 
            #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
            #geom_tiplab() + 
            #geom_tippoint() + 
            theme_tree2() +
            geom_treescale() 
          #xlim(0, 0.5)
          fig_tree
          
          # clusters:
          cluster_summary_table<-ordered_table2 %>% 
            filter(seq_ID!="NA") %>% 
            group_by(cluster) %>% 
            summarize(chr=chr[1],
                      start=min(start),
                      total_seq=n(), 
                      total_samples=length(unique(sample)), 
                      min_length=min(len), 
                      max_length=max(len)) %>% 
            filter(total_samples>=4)
          
          library(svglite)
          for(i in seq(1,dim(cluster_summary_table)[1])){
            cluster_number<-as.vector(unlist(cluster_summary_table[i,"cluster"]))
            cluster_chr<-as.vector(unlist(cluster_summary_table[i,"chr"]))
            cluster_start<-as.vector(unlist(cluster_summary_table[i,"start"]))
            cluster_total_samples<-as.vector(unlist(cluster_summary_table[i,"total_samples"]))
            cluster_total_seq<-as.vector(unlist(cluster_summary_table[i,"total_seq"]))
            cluster_min_length<-as.vector(unlist(cluster_summary_table[i,"min_length"]))
            cluster_max_length<-as.vector(unlist(cluster_summary_table[i,"max_length"]))
            cluster_table<-ordered_table2 %>% 
              mutate(cluster_set=factor(ifelse(seq_ID %in% c(ordered_table2 %>% 
                                                             filter(seq_ID!="NA") %>% 
                                                             filter(cluster==cluster_number) %>% 
                                                             select(seq_ID) %>% 
                                                             unlist() %>% 
                                                             as.vector()), 1,0), levels=c(0,1))) %>% 
              mutate(anc_group=factor(ifelse(anc_prop>0.9,"Sp",ifelse(anc_prop<0.1,"Sk","Hyb")), 
                                      levels=c("Sp", "Hyb", "Sk"))) %>% 
              select(seq_ID, cluster_set, anc_group) 
            
            fig_tree %<+% cluster_table + 
              geom_tippoint(aes(colour=anc_group, alpha=factor(cluster_set)), size=4) +
              scale_color_manual(values=c("darkred", "orange", "steelblue")) +
              scale_alpha_manual(values=c(0,0.6,0))+
              labs(title=paste0("Cluster:",cluster_number," Chr:", cluster_chr, " Pos:", cluster_start, 
                                " NSamples:", cluster_total_samples, 
                                " NSeq:", cluster_total_seq, 
                                " MinLen:", cluster_min_length, 
                                " MaxLen:", cluster_max_length))+
              theme(legend.position = "none") +
              ggsave(paste0("./tf_clusters/cluster",cluster_number,".svg"))
            fig_tree %<+% cluster_table + 
              geom_tippoint(aes(colour=anc_group, alpha=factor(cluster_set)), size=4) +
              scale_color_manual(values=c("darkred", "orange", "steelblue")) +
              scale_alpha_manual(values=c(0,0.6,0))+
              labs(title=paste0("Cluster:",cluster_number," Chr:", cluster_chr, " Pos:", cluster_start, 
                                " NSamples:", cluster_total_samples, 
                                " NSeq:", cluster_total_seq))+
              theme(legend.position = "none") +
              ggsave(paste0("./tf_clusters/cluster",cluster_number,".png"))
          }
          
          
          #####
          
          tree_ed@data$UFboot[is.na(tree_ed@data$UFboot)]<-0
          tree_ed@data$UFboot[tree_ed@data$UFboot<90]<-0
          tree_ed@data$UFboot[tree_ed@data$UFboot>=90 & tree_ed@data$UFboot<95]<-90
          tree_ed@data$UFboot[tree_ed@data$UFboot>=95]<-95
          
          cluster_ID<-data.frame(id=tree_ed@phylo$tip.label) %>% 
            merge(ordered_table2 %>% 
                  mutate(id=seq_ID) %>%
                  select(id, cluster, chr, anc_prop), by="id", all.x=T) %>% 
            mutate(anc_group=factor(ifelse(anc_prop>0.9,"Sp",ifelse(anc_prop<0.1,"Sk","Hyb")), 
                                    levels=c("Sp", "Hyb", "Sk"))) %>%
            arrange(match(id, tree_ed@phylo$tip.label))
          
          row.names(cluster_ID) <- NULL
          
          head(cluster_ID)
          lables_tf<-ordered_table2 %>% 
            select(chr, cluster) %>%
            group_by(chr) %>% 
            dplyr::summarise(min_clus=min(cluster), 
                             max_clus=max(cluster), 
                             #multiple30=ifelse(cluster %% 30 == 0, cluster, 0)) %>% 
                             mid_1=((max_clus-min_clus)/4)+min_clus, 
                             mid_3=((max_clus-min_clus)*3/4)+min_clus, 
                             middle=round((max_clus-min_clus)/2)+min_clus) %>%
            ungroup() %>% 
            select(-chr, -min_clus, -max_clus) %>%  
            #filter(multiple30!=0) %>% 
            #select(-chr) %>%  
            unlist() %>% 
            as.vector() %>% 
            round() %>% 
            sort() %>% 
            as.numeric()
          
          p <-ggtree(tree_ed) +
            #aes(color=as.factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)))), 
            #
            theme_tree2() +
            geom_nodepoint(aes(alpha=as.factor(UFboot)), color="darkgreen", size=2) +
            #scale_alpha_manual(values=c(0,0.5,1)) + 
            # geom_nodepoint(aes(color=factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)), levels=c(90,80,0)), 
            #                    alpha=factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0))), levels=c(90,80,0)), 
            #                               size=3) +
            geom_treescale() +
            #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
            #geom_tiplab() + 
            #geom_tippoint(color="black", shape=4, size=1.5, alpha=0.2) +
            #geom_tree() +
            # geom_nodepoint(aes(color=as.factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)))), 
            #                alpha=0.5, 
            #                size=3) +
            #geom_tree(data=df, color='firebrick')
            #scale_color_continuous(low='darkgreen', high='red') +
            #scale_color_manual(values=c("gray", "black", "darkred")) +
            #scale_color_manual(values=c("white", "white", "darkred")) +
          theme(legend.position="right")
          #geom_tiplab(align=TRUE, linesize=1, col="gray90") 
          # xlim(0, 0.4) 
          
          # pp <- p %<+% cluster_ID + 
          #   geom_treescale()+
          #   geom_tiplab(align = TRUE, linesize = 0, alpha=0) +
          #   geom_tippoint(aes(colour=anc_group,
          #                     alpha=anc_group), 
          #                     size=2) +
          #   scale_color_manual(values=c("darkred", "orange", "steelblue")) +
          #   scale_alpha_manual(values=c(0.8,0.0,0.8)) +
          #   theme(legend.position = "None")
          
          pp<-p
          print(pp)
          #pp3<-pp+geom_nodelab(identify(pp))
          #hilight(identify(pp), fill = "gray")
          #print(pp3)
          # cluster_ID<-data.frame(id=tree_ed$tip.label) %>% 
          #   merge(ordered_table2 %>% 
          #           mutate(id=seq_ID) %>%
          #           select(id, cluster, chr), by="id", all.x=T) %>% 
          #   arrange(match(id, tree_ed$tip.label))
          identify(pp)
          
          # generate some random values for each tip label in the data
          # Make a second plot with the original, naming the new plot "dot", 
          # using the data you just created, with a point geom.
          
          p2 <- facet_plot(pp, panel="I", data=cluster_ID %>% filter(chr=="I"), 
                           geom=geom_point, aes(x=cluster, colour=anc_group, 
                                                alpha=anc_group, 
                                                size=anc_group))+ 
            #color='darkred', 
            #alpha=0.8) +
            scale_color_manual(values=c("darkred", "orange", "steelblue")) +
            scale_alpha_manual(values=c(0,0.3,1,0.5, 0.9, 0.9)) + 
            scale_size_manual(values=c(2,1.5,2)) +
            theme_classic() +
            scale_x_continuous(breaks = lables_tf, labels = lables_tf) +
            scale_y_continuous(breaks = seq(1,800,30)) +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
                  panel.grid.major = element_line(colour = "gray80"), 
                  axis.text.x=element_text(size=12, angle = 45, hjust = 1, colour="black"), 
                  axis.text.y=element_blank(), 
                  axis.ticks.y = element_blank(), 
                  legend.position="none")
          #panel.grid.minor = element_line(colour = "gray80"))
          p2
          clade_pop<-identify(p2)
          clade_pop
          p3 <- facet_plot(p2, panel="II", data=cluster_ID %>% filter(chr=="II"), 
                           geom=geom_point, aes(x=cluster, colour=anc_group, 
                                                alpha=anc_group, 
                                                size=anc_group))
          #color='darkred', 
          # alpha=0.8) 
          p4 <- facet_plot(p3, panel="III", data=cluster_ID %>% filter(chr=="III"), 
                           geom=geom_point, aes(x=cluster, colour=anc_group, 
                                                alpha=anc_group, 
                                                size=anc_group))
          #color='darkred', 
          # alpha=0.8) 
          p4 
          
          
          
          library(ggtree)
          library(reshape2)
          
          # load the packaged
          library(grid)
          library(gtable)
          
          head(cluster_ID)
          
          panel_size<-cluster_ID %>% 
            filter(!(is.na(chr))) %>%
            group_by(chr) %>%
            summarise(max_cluster=max(cluster), 
                      min_cluster=min(cluster)) %>% 
            mutate(size_ch=max_cluster-min_cluster) %>% 
            mutate(total=cluster_ID %>% 
                   filter(!(is.na(chr))) %>% 
                   select(cluster) %>% 
                   unlist() %>% 
                   as.vector() %>% 
                   max()) %>% 
            mutate(size_block=0.4*size_ch/total)
          
          
          gt = ggplot_gtable(ggplot_build(p4))
          #gtable_show_layout(gt) # will show you the layout - very handy function
          #gt # see plot layout in table format
          #gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
          #gt$layout$l[grep('panel-3', gt$layout$name)] # you want to find the column specific to panel-3
          #gt$layout$l[grep('panel-4', gt$layout$name)] # you want to find the column specific to panel-4
          gt$widths[7] = (unlist(panel_size[panel_size$chr=="I","size_block"]))*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
          gt$widths[9] = (unlist(panel_size[panel_size$chr=="II","size_block"]))*gt$widths[9]
          gt$widths[11] = (unlist(panel_size[panel_size$chr=="III","size_block"]))*gt$widths[11]
          
          grid.draw(gt) # plot with grid draw
          
          
          svg("05_tree_distribution_TFWithoutLTRs.svg",height = 250/25.4, width = 180/25.4, pointsize = 7)
          grid.draw(gt) # plot with grid draw
          dev.off()
          
          
          #
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          # LTR:
          # install.packages("tidyverse",  lib="/home/sergio/R/libraries_Rackham/3.4")
          
          #library("ggplot2")
          #library("tidyr")
          #library("plyr")
          #library("dplyr")
          library("tidyverse")
          
          
          #
          pro_AncPop<-read.table("proportions_AncPop_161Samples_200_100.txt", T) %>% 
            select(sample_JB, an_population, number_bins_fq_total) %>% 
            spread(an_population, number_bins_fq_total) %>% 
            select(sample_JB, JB22_pop1) 
          names(pro_AncPop)<-c("sample_JB", "anc_prop")
          
          ID_table<-read.table("sample_ID.txt", T) %>% 
            merge(pro_AncPop, by="sample_JB", all.x=T)
          # 
          
          head(ID_table)
          
          
          samples_table<-read.table("table_list_samples_ltr_completed_masked_conSeq.txt", F)
          head(samples_table)
          names(samples_table)<-c("seq_ID", "sample", "chr", "start", "dis", "len_froSeq","len")
          
          ordered_table2 -samples_table %>% 
            merge(ID_table, by="sample", all.x=T) %>% 
            arrange(chr, start, sample, dis) %>% 
            filter(seq_ID!="NA") %>% 
            filter(!(is.na(chr))) %>% 
            mutate(len=350, 
                   len_froSeq=350)
          
          head(ordered_table2 
               dim(ordered_table2 [1]
                   
                   cluster<-c(1)
                   length_cluster<-0
                   min_start<-0
                   for (i in seq(2,dim(ordered_table2 [1])){
                     if(length_cluster==0){
                       length_cluster<-ordered_table2 len_froSeq[i-1]
                       min_start<-ordered_table2 start[i-1]
                     }
                     if(ordered_table2 chr[i]==ordered_table2 chr[i-1]){
                       if ((ordered_table2 start[i]-(min_start+length_cluster))<=(-50)){
                         print(i)
                         cluster[i]<-cluster[i-1]
                         length_cluster<-max(length_cluster,ordered_table2 len_froSeq[i])
                         min_start<-min(min_start, ordered_table2 start[i])
                       } else {
                         cluster[i]<-cluster[i-1]+1
                         length_cluster<-ordered_table2 len_froSeq[i]
                         min_start<-ordered_table2 start[i]
                       }
                     } else {
                       cluster[i]<-cluster[i-1]+1
                       length_cluster<-ordered_table2 len_froSeq[i]
                       min_start<-ordered_table2 start[i]
                     }
                   }
                   
                   ordered_table2 -cbind(ordered_table2  cluster)
                   
                   ordered_table2 %>% 
                     filter(seq_ID!="NA") %>% 
                     group_by(cluster) %>% 
                     summarize(chr=chr[1],
                               start=min(start),
                               total_seq=n(), 
                               total_samples=length(unique(sample)), 
                               min_length=min(len), 
                               max_length=max(len))
                   
                   ordered_table2 %>% 
                     filter(seq_ID!="NA") %>% 
                     merge((ordered_table2 %>% 
                            filter(seq_ID!="NA") %>% 
                            group_by(cluster) %>% 
                            summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
                     mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
                     write.table("ltr_completed_new_sorted_table_list_samples_masked_conSeq.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
                   
                   
                   # plot size distribution and number of clusters:
                   ordered_table2 %>% 
                     filter(seq_ID!="NA") %>% 
                     merge((ordered_table2 %>% 
                            filter(seq_ID!="NA") %>% 
                            group_by(cluster) %>% 
                            summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
                     mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>% 
                     #filter(len<450) %>% 
                     ggplot(aes((IDcluster), len)) +
                     geom_line(aes(group=cluster), alpha=0.4)+
                     geom_jitter(aes(colour = as.factor(cluster)), size=1, alpha=0.8) +
                     xlab("Cluster ID") +
                     ylab("Length (bp)") +
                     scale_color_manual(values=c(rep(c("black", "gray"),max(ordered_table2 cluster)-1), "black")) +
                     theme_classic() +
                     theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                           legend.position="none") +
                     facet_grid(. ~ chr, scale="free", space="free") #+
                   #ggsave("01_size_distribution_ltr_completed.png", width = 12, height = 4)
                   
                   
                   # plot heat map with number of copies per cluster and sample:
                   ordered_table2 %>% 
                     filter(seq_ID!="NA") %>% 
                     merge((ordered_table2 %>% 
                            filter(seq_ID!="NA") %>% 
                            group_by(cluster) %>% 
                            summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
                     mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
                     filter(len<450) %>% 
                     group_by(chr, IDcluster, sample) %>% 
                     summarize(n_copies=length(seq_ID), 
                               dif_length=(max(len)-min(len))) %>% 
                     mutate(text_label=ifelse(dif_length==0,"",as.character(dif_length))) %>% 
                     ggplot(aes(IDcluster, sample)) +
                     geom_tile(aes(fill=as.factor(n_copies))) +
                     scale_fill_manual(values=c("#FFDB6D", "#C4961A", "#F4EDCA", 
                                                "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")) +
                     #geom_text(aes(label=text_label)) +
                     xlab("Cluster ID") +
                     ylab("Sample") +
                     theme_classic() +
                     theme(axis.text.x=element_blank(), 
                           legend.position="top") +
                     facet_grid(. ~ chr, scale="free", space="free") #+
                   #ggsave("02_heatmap_ltr_completed_nCopies_persample.png", width = 12, height = 8)
                   
                   
                   ordered_table2 %>% 
                     merge((ordered_table2 %>% 
                            group_by(cluster) %>% 
                            summarize(min_start=min(start))), by="cluster") %>% 
                     group_by(cluster, chr, min_start, sample) %>% 
                     summarize(n_copies=n()) %>% 
                     group_by(cluster, chr, min_start) %>% 
                     summarize(n_samples=n()) %>% 
                     mutate(min_start=min_start/1000000) %>% 
                     ggplot(aes(min_start, n_samples)) +
                     geom_point(size=2, alpha=0.5) +
                     scale_x_continuous(breaks=seq(0,5.5,0.5), labels=seq(0,5.5,0.5))+
                     scale_y_continuous(breaks=seq(0,20,2), labels=seq(0,20,2))+
                     xlab("Pos.") +
                     ylab("Number Samples") +
                     theme_classic() +
                     theme(axis.text.x=element_text(size=14, angle = 45, hjust = 1, colour="black"), 
                           axis.text.y=element_text(colour="black", size=14)) +
                     facet_grid(. ~ chr, scale="free", space="free") #+
                   #ggsave("03_SamplesperCluster_ltr_completed.png", width = 12, height = 3)
                   
                   ordered_table2 %>% 
                     merge((ordered_table2 %>% 
                            group_by(cluster) %>% 
                            summarize(min_start=min(start))), by="cluster") %>% 
                     group_by(cluster, sample) %>% 
                     summarize(n_repeats=n()) %>% 
                     group_by(cluster) %>% 
                     summarize(n_samples=n()) %>% 
                     ggplot(aes(n_samples)) +
                     #geom_histogram(binwidth=1) +
                     geom_histogram(aes(y=..density..), binwidth=1) +
                     scale_x_continuous(breaks=seq(0,18,1), labels=seq(0,18,1)) +
                     #ylab("Count") +
                     ylab("Density") +
                     xlab("Number Samples") +
                     theme_classic() +
                     theme(axis.text.x=element_text(size=14, colour="black"), 
                           axis.text.y=element_text(colour="black", size=14), 
                           axis.title=element_text(colour="black", size=14)) #+
                   #ggsave("04_fq_homologous_ltf.png", width = 7, height = 5)
                   #ggsave("04_fq_homologous_ltr_completed_density.png", width = 7, height = 5)
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   library("ggtree")
                   tree_ed<-read.iqtree("../phylogenies/phy_ltr_complete_minLen200.treefile")
                   #tree_ed<- as.phylo(tree)
                   fig_tree<- ggtree(tree_ed) + 
                     #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
                     #geom_tiplab() + 
                     #geom_tippoint() + 
                     theme_tree2() +
                     geom_treescale() 
                   #xlim(0, 0.5)
                   fig_tree
                   
                   
                   # clusters:
                   cluster_summary_table<-ordered_table2 %>% 
                     filter(seq_ID!="NA") %>% 
                     group_by(cluster) %>% 
                     summarize(chr=chr[1],
                               start=min(start),
                               total_seq=n(), 
                               total_samples=length(unique(sample)), 
                               min_length=min(len), 
                               max_length=max(len)) %>% 
                     filter(total_samples>=4)
                   
                   library(svglite)
                   for(i in seq(1,dim(cluster_summary_table)[1])){
                     cluster_number<-as.vector(unlist(cluster_summary_table[i,"cluster"]))
                     cluster_chr<-as.vector(unlist(cluster_summary_table[i,"chr"]))
                     cluster_start<-as.vector(unlist(cluster_summary_table[i,"start"]))
                     cluster_total_samples<-as.vector(unlist(cluster_summary_table[i,"total_samples"]))
                     cluster_total_seq<-as.vector(unlist(cluster_summary_table[i,"total_seq"]))
                     cluster_min_length<-as.vector(unlist(cluster_summary_table[i,"min_length"]))
                     cluster_max_length<-as.vector(unlist(cluster_summary_table[i,"max_length"]))
                     cluster_table<-ordered_table2 %>% 
                       mutate(cluster_set=factor(ifelse(seq_ID %in% c(ordered_table2 %>% 
                                                                      filter(seq_ID!="NA") %>% 
                                                                      filter(cluster==cluster_number) %>% 
                                                                      select(seq_ID) %>% 
                                                                      unlist() %>% 
                                                                      as.vector()), 1,0), levels=c(0,1))) %>% 
                       mutate(anc_group=factor(ifelse(anc_prop>0.9,"Sp",ifelse(anc_prop<0.1,"Sk","Hyb")), 
                                               levels=c("Sp", "Hyb", "Sk"))) %>% 
                       select(seq_ID, cluster_set, anc_group) 
                     
                     fig_tree %<+% cluster_table + 
                       geom_tippoint(aes(colour=anc_group, alpha=factor(cluster_set)), size=4) +
                       scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                       scale_alpha_manual(values=c(0,0.6,0))+
                       labs(title=paste0("Cluster:",cluster_number," Chr:", cluster_chr, " Pos:", cluster_start, 
                                         " NSamples:", cluster_total_samples, 
                                         " NSeq:", cluster_total_seq, 
                                         " MinLen:", cluster_min_length, 
                                         " MaxLen:", cluster_max_length))+
                       theme(legend.position = "none") +
                       ggsave(paste0("./tf_clusters/cluster",cluster_number,".svg"))
                     fig_tree %<+% cluster_table + 
                       geom_tippoint(aes(colour=anc_group, alpha=factor(cluster_set)), size=4) +
                       scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                       scale_alpha_manual(values=c(0,0.6,0))+
                       labs(title=paste0("Cluster:",cluster_number," Chr:", cluster_chr, " Pos:", cluster_start, 
                                         " NSamples:", cluster_total_samples, 
                                         " NSeq:", cluster_total_seq))+
                       theme(legend.position = "none") +
                       ggsave(paste0("./ltr_completed_clusters/cluster",cluster_number,".png"))
                   }
                   
                   
                   
                   tree_ed@data$UFboot[is.na(tree_ed@data$UFboot)]<-0
                   tree_ed@data$UFboot[tree_ed@data$UFboot<90]<-0
                   tree_ed@data$UFboot[tree_ed@data$UFboot>=90 & tree_ed@data$UFboot<95]<-90
                   tree_ed@data$UFboot[tree_ed@data$UFboot>=95]<-95
                   
                   cluster_ID<-data.frame(id=tree_ed@phylo$tip.label) %>% 
                     merge(ordered_table2 %>% 
                           mutate(id=seq_ID) %>%
                           select(id, cluster, chr, anc_prop), by="id", all.x=T) %>% 
                     mutate(anc_group=factor(ifelse(anc_prop>0.9,"Sp",ifelse(anc_prop<0.1,"Sk","Hyb")), 
                                             levels=c("Sp", "Hyb", "Sk"))) %>%
                     arrange(match(id, tree_ed@phylo$tip.label))
                   
                   row.names(cluster_ID) <- NULL
                   
                   lables_tf<-ordered_table2 %>% 
                     select(chr, cluster) %>%
                     group_by(chr) %>% 
                     dplyr::summarise(min_clus=min(cluster), 
                                      max_clus=max(cluster), 
                                      #multiple30=ifelse(cluster %% 30 == 0, cluster, 0)) %>% 
                                      mid_1=((max_clus-min_clus)/4)+min_clus, 
                                      mid_3=((max_clus-min_clus)*3/4)+min_clus, 
                                      middle=round((max_clus-min_clus)/2)+min_clus) %>%
                     ungroup() %>% 
                     select(-chr, -min_clus, -max_clus) %>%  
                     #filter(multiple30!=0) %>% 
                     #select(-chr) %>%  
                     unlist() %>% 
                     as.vector() %>% 
                     round() %>% 
                     sort() %>% 
                     as.numeric()
                   
                   p <-ggtree(tree_ed) +
                     #aes(color=as.factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)))), 
                     #
                     theme_tree2() +
                     geom_nodepoint(aes(alpha=as.factor(UFboot)), color="darkgreen", size=2) +
                     #scale_alpha_manual(values=c(0,0.5,1)) + 
                     # geom_nodepoint(aes(color=factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)), levels=c(90,80,0)), 
                     #                    alpha=factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0))), levels=c(90,80,0)), 
                     #                               size=3) +
                     geom_treescale() +
                     #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
                     #geom_tiplab() + 
                     #geom_tippoint(color="black", shape=4, size=1.5, alpha=0.2) +
                     #geom_tree() +
                     # geom_nodepoint(aes(color=as.factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)))), 
                     #                alpha=0.5, 
                     #                size=3) +
                     #geom_tree(data=df, color='firebrick')
                     #scale_color_continuous(low='darkgreen', high='red') +
                     #scale_color_manual(values=c("gray", "black", "darkred")) +
                     #scale_color_manual(values=c("white", "white", "darkred")) +
                   theme(legend.position="right")
                   #geom_tiplab(align=TRUE, linesize=1, col="gray90") 
                   # xlim(0, 0.4) 
                   
                   # pp <- p %<+% cluster_ID + 
                   #   geom_treescale()+
                   #   geom_tiplab(align = TRUE, linesize = 0, alpha=0) +
                   #   geom_tippoint(aes(colour=anc_group,
                   #                     alpha=anc_group), 
                   #                     size=2) +
                   #   scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                   #   scale_alpha_manual(values=c(0.8,0.0,0.8)) +
                   #   theme(legend.position = "None")
                   
                   pp<-p
                   
                   pp
                   
                   # cluster_ID<-data.frame(id=tree_ed$tip.label) %>% 
                   #   merge(ordered_table2 %>% 
                   #           mutate(id=seq_ID) %>%
                   #           select(id, cluster, chr), by="id", all.x=T) %>% 
                   #   arrange(match(id, tree_ed$tip.label))
                   
                   
                   # generate some random values for each tip label in the data
                   # Make a second plot with the original, naming the new plot "dot", 
                   # using the data you just created, with a point geom.
                   
                   p2 <- facet_plot(pp, panel="I", data=cluster_ID %>% filter(chr=="I"), 
                                    geom=geom_point, aes(x=cluster, colour=anc_group, 
                                                         alpha=anc_group, 
                                                         size=anc_group))+ 
                     #color='darkred', 
                     #alpha=0.8) +
                     scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                     scale_alpha_manual(values=c(0,0.3,1,0.5, 0.9, 0.9)) + 
                     scale_size_manual(values=c(2,1.5,2)) +
                     theme_classic() +
                     scale_x_continuous(breaks = lables_tf, labels = lables_tf) +
                     scale_y_continuous(breaks = seq(1,800,30)) +
                     theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
                           panel.grid.major = element_line(colour = "gray80"), 
                           axis.text.x=element_text(size=12, angle = 45, hjust = 1, colour="black"), 
                           axis.text.y=element_blank(), 
                           axis.ticks.y = element_blank(), 
                           legend.position="none")
                   #panel.grid.minor = element_line(colour = "gray80"))
                   p2
                   p3 <- facet_plot(p2, panel="II", data=cluster_ID %>% filter(chr=="II"), 
                                    geom=geom_point, aes(x=cluster, colour=anc_group, 
                                                         alpha=anc_group, 
                                                         size=anc_group))
                   #color='darkred', 
                   # alpha=0.8) 
                   p4 <- facet_plot(p3, panel="III", data=cluster_ID %>% filter(chr=="III"), 
                                    geom=geom_point, aes(x=cluster, colour=anc_group, 
                                                         alpha=anc_group, 
                                                         size=anc_group))
                   #color='darkred', 
                   # alpha=0.8) 
                   #p4 
                   
                   
                   
                   library(ggtree)
                   library(reshape2)
                   
                   # load the packaged
                   library(grid)
                   library(gtable)
                   
                   head(cluster_ID)
                   
                   panel_size<-cluster_ID %>% 
                     filter(!(is.na(chr))) %>%
                     group_by(chr) %>%
                     summarise(max_cluster=max(cluster), 
                               min_cluster=min(cluster)) %>% 
                     mutate(size_ch=max_cluster-min_cluster) %>% 
                     mutate(total=cluster_ID %>% 
                            filter(!(is.na(chr))) %>% 
                            select(cluster) %>% 
                            unlist() %>% 
                            as.vector() %>% 
                            max()) %>% 
                     mutate(size_block=0.4*size_ch/total)
                   
                   
                   gt = ggplot_gtable(ggplot_build(p4))
                   #gtable_show_layout(gt) # will show you the layout - very handy function
                   #gt # see plot layout in table format
                   #gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
                   #gt$layout$l[grep('panel-3', gt$layout$name)] # you want to find the column specific to panel-3
                   #gt$layout$l[grep('panel-4', gt$layout$name)] # you want to find the column specific to panel-4
                   gt$widths[7] = (unlist(panel_size[panel_size$chr=="I","size_block"]))*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
                   gt$widths[9] = (unlist(panel_size[panel_size$chr=="II","size_block"]))*gt$widths[9]
                   gt$widths[11] = (unlist(panel_size[panel_size$chr=="III","size_block"]))*gt$widths[11]
                   
                   #grid.draw(gt) # plot with grid draw
                   
                   svg("05_tree_distribution_LTR_completed.svg",height = 250/25.4, width = 180/25.4, pointsize = 7)
                   grid.draw(gt) # plot with grid draw
                   dev.off()
                   
                   # 2PK distance:
                   
                   # distance<-read.csv("../2PK_distance/Distance_Data_ltr_completed.csv", T) 
                   # head(distance)
                   # 
                   # seq="01 ref I 1465320 10 4919 4915 1 complete ltr"
                   # group_id<-distance %>% 
                   #   #filter(Dist>=0 & Dist<0.04) %>% 
                   #   filter(Species.1==seq | Species.2==seq) %>% 
                   #   mutate(id_sample=ifelse(Species.1==seq, as.vector(Species.2), as.vector(Species.1)),
                   #          group_sample=ifelse(Dist>0.1,"sk","sp")) %>% 
                   #   select(id_sample, group_sample) 
                   # 
                   # head(group_id)
                   # 
                   # distance_groups<-distance %>% 
                   #   merge(group_id %>% 
                   #           mutate(Species.1=id_sample) %>% 
                   #           select(-id_sample), by="Species.1") %>% 
                   #   merge(group_id %>% 
                   #           mutate(Species.2=id_sample) %>% 
                   #           select(-id_sample), by="Species.2") %>% 
                   #   mutate(comparison=factor(ifelse(group_sample.x=="sp" & group_sample.y=="sp", "sp_sp", 
                   #                                   ifelse(group_sample.x=="sk" & group_sample.y=="sk", 
                   #                                          "sk_sk", "sp_sk")), 
                   #                            levels=c("sp_sp", "sk_sk", "sp_sk"))) 
                   # 
                   # 
                   # distance_groups %>% 
                   #   ggplot(aes(Dist, fill=comparison), alpha=0.7) +
                   #   # geom_histogram() +
                   #   # geom_histogram(binwidth=1) +
                   #   geom_histogram(aes(y=..density..), position="dodge", binwidth = 0.01) +
                   #   scale_fill_manual(values=c("darkred", "steelblue", "orange")) +
                   #   # facet_grid(comparison ~ .)
                   #   # scale_x_continuous(breaks=seq(0,18,1), labels=seq(0,18,1)) +
                   #   # ylab("Count") +
                   #   ylab("Density") +
                   #   xlab("Distance") +
                   #   theme_classic() +
                   #   theme(axis.text.x=element_text(size=14, colour="black"), 
                   #         axis.text.y=element_text(colour="black", size=14), 
                   #         axis.title=element_text(colour="black", size=14)) +
                   #   ggsave("06_2PK_distance_groups_ltr_completed.png")
                   #   
                   # 
                   # distance_groups %>% 
                   #   ggplot(aes(Dist)) +
                   #   # geom_histogram() +
                   #   # geom_histogram(binwidth=1) +
                   #   geom_histogram(aes(y=..density..), binwidth = 0.01, colour="black") +
                   #   #scale_fill_manual(values=c("darkred", "steelblue", "orange")) +
                   #   # facet_grid(comparison ~ .)
                   #   # scale_x_continuous(breaks=seq(0,18,1), labels=seq(0,18,1)) +
                   #   # ylab("Count") +
                   #   ylab("Density") +
                   #   xlab("Distance") +
                   #   theme_classic() +
                   #   theme(axis.text.x=element_text(size=14, colour="black"), 
                   #         axis.text.y=element_text(colour="black", size=14), 
                   #         axis.title=element_text(colour="black", size=14)) +
                   #   ggsave("06_2PK_distance_ltr_completed.png")
                   # 
                   
                   
                   library(ggplot2)
                   library(plyr)
                   library(arm)
                   library(reshape2)
                   library(dplyr)
                   
                   
                   test<-ape::read.FASTA("../Alignments/alig_ltr_complete_ltr_minLen200.fasta", type = "DNA")
                   test
                   distance_ltr<-ape::dist.dna(test, model = "K80", variance = FALSE,
                                               gamma = FALSE, pairwise.deletion = FALSE,
                                               base.freq = NULL, as.matrix = T)
                   distance_ltr[lower.tri(distance_ltr, diag = T)] <- NA
                   distance <- melt(distance_ltr) %>% 
                     filter(!(is.na(value)))
                   
                   seq="01_ref_I_1465320_10_4919_4915_1_complete_ltr"
                   group_id<-distance %>% 
                     #filter(Dist>=0 & Dist<0.04) %>% 
                     filter(Var1==seq | Var2==seq) %>% 
                     mutate(id_sample=ifelse(Var1==seq, as.vector(Var2), as.vector(Var1)),
                            group_sample=ifelse(value>0.1,"sk","sp")) %>% 
                     dplyr::select(id_sample, group_sample) 
                   
                   head(group_id)
                   
                   distance_groups<-distance %>% 
                     merge(group_id %>% 
                           mutate(Var1=id_sample) %>% 
                           dplyr::select(-id_sample), by="Var1") %>% 
                     merge(group_id %>% 
                           mutate(Var2=id_sample) %>% 
                           dplyr::select(-id_sample), by="Var2") %>% 
                     mutate(comparison=factor(ifelse(group_sample.x=="sp" & group_sample.y=="sp", "sp_sp", 
                                                     ifelse(group_sample.x=="sk" & group_sample.y=="sk", 
                                                            "sk_sk", "sp_sk")), 
                                              levels=c("sp_sp", "sk_sk", "sp_sk"))) 
                   
                   
                   distance_groups %>% 
                     ggplot(aes(value, fill=comparison), alpha=0.7) +
                     # geom_histogram() +
                     # geom_histogram(binwidth=1) +
                     geom_histogram(aes(y=..density..), position="dodge", binwidth = 0.01) +
                     scale_fill_manual(values=c("darkred", "steelblue", "orange")) +
                     # facet_grid(comparison ~ .)
                     # scale_x_continuous(breaks=seq(0,18,1), labels=seq(0,18,1)) +
                     # ylab("Count") +
                     ylab("Density") +
                     xlab("Distance") +
                     theme_classic() +
                     theme(axis.text.x=element_text(size=14, colour="black"), 
                           axis.text.y=element_text(colour="black", size=14), 
                           axis.title=element_text(colour="black", size=14)) #+
                   ggsave("06_2PK_distance_groups_ltr_completed.png")
                   
                   distance_groups %>% 
                     ggplot(aes(value)) +
                     # geom_histogram() +
                     # geom_histogram(binwidth=1) +
                     geom_histogram(aes(y=..density..), binwidth = 0.01, colour="black") +
                     #scale_fill_manual(values=c("darkred", "steelblue", "orange")) +
                     # facet_grid(comparison ~ .)
                     # scale_x_continuous(breaks=seq(0,18,1), labels=seq(0,18,1)) +
                     # ylab("Count") +
                     ylab("Density") +
                     xlab("Distance") +
                     theme_classic() +
                     theme(axis.text.x=element_text(size=14, colour="black"), 
                           axis.text.y=element_text(colour="black", size=14), 
                           axis.title=element_text(colour="black", size=14)) +
                     ggsave("06_2PK_distance_ltr_completed.png")
                   
                   
                   
                   
                   ##
                   ##
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   ##### soLO ltrS:
                   
                   
                   #
                   pro_AncPop<-read.table("proportions_AncPop_161Samples_200_100.txt", T) %>% 
                     select(sample_JB, an_population, number_bins_fq_total) %>% 
                     spread(an_population, number_bins_fq_total) %>% 
                     select(sample_JB, JB22_pop1) 
                   names(pro_AncPop)<-c("sample_JB", "anc_prop")
                   
                   ID_table<-read.table("sample_ID.txt", T) %>% 
                     merge(pro_AncPop, by="sample_JB", all.x=T)
                   # 
                   
                   head(ID_table)
                   
                   
                   samples_table<-read.table("table_list_samples_ltr_completePlusSolo_masked_conSeq.txt", T)
                   head(samples_table)
                   names(samples_table)<-c("seq_ID", "sample", "chr", "start", "dis", "len_froSeq","len")
                   completed_ltr<-read.table("table_list_samples_ltr_completePlusSolo_masked_conSeq_completedSeq.txt", T) %>%
                     select(idseq) %>% 
                     unlist() %>% 
                     as.vector()
                   
                   
                   # samples_table<-samples_table 
                   #   mutate(len_froSeq=ifelse(len_froSeq>600,"NA",len_froSeq),
                   #          len=ifelse(len>600,"NA",len))
                   
                   ordered_table2 -samples_table %>% 
                     merge(ID_table, by="sample", all.x=T) %>% 
                     arrange(chr, start, sample, dis) %>% 
                     filter(seq_ID!="NA") %>% 
                     filter(!(is.na(chr))) %>% 
                     mutate(completed_ltr=(seq_ID %in% completed_ltr)*1) #%>% 
                   # filter(completed_ltr==0)
                   
                   head(ordered_table2 
                        
                        cluster<-c(1)
                        length_cluster<-0
                        min_start<-0
                        for (i in seq(2,dim(ordered_table2 [1])){
                          if(length_cluster==0){
                            length_cluster<-ordered_table2 len_froSeq[i-1]
                            min_start<-ordered_table2 start[i-1]
                          }
                          if(ordered_table2 chr[i]==ordered_table2 chr[i-1]){
                            if ((ordered_table2 start[i]-(min_start+length_cluster))<=(-50)){
                              print(i)
                              cluster[i]<-cluster[i-1]
                              length_cluster<-max(length_cluster,ordered_table2 len_froSeq[i])
                              min_start<-min(min_start, ordered_table2 start[i])
                            } else {
                              cluster[i]<-cluster[i-1]+1
                              length_cluster<-ordered_table2 len_froSeq[i]
                              min_start<-ordered_table2 start[i]
                            }
                          } else {
                            cluster[i]<-cluster[i-1]+1
                            length_cluster<-ordered_table2 len_froSeq[i]
                            min_start<-ordered_table2 start[i]
                          }
                        }
                        
                        ordered_table2 -cbind(ordered_table2  cluster) %>% 
                          mutate(sample = fct_reorder(sample, anc_prop))
                        
                        max(ordered_table2 cluster)
                        
                        
                        ordered_table2 %>% 
                          filter(seq_ID!="NA") %>% 
                          group_by(cluster) %>% 
                          summarize(chr=chr[1],
                                    start=min(start),
                                    total_seq=n(), 
                                    total_samples=length(unique(sample)), 
                                    min_length=min(len), 
                                    max_length=max(len))
                        
                        ordered_table2 %>% 
                          filter(seq_ID!="NA") %>% 
                          merge((ordered_table2 %>% 
                                 filter(seq_ID!="NA") %>% 
                                 group_by(cluster) %>% 
                                 summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
                          mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
                          write.table("ltr_completedPlusSolo_new_sorted_table_list_samples_masked_conSeq.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
                        
                        
                        # plot size distribution and number of clusters:
                        ordered_table2 %>% 
                          filter(seq_ID!="NA") %>% 
                          filter(len<500) %>% 
                          merge((ordered_table2 %>% 
                                 filter(seq_ID!="NA") %>% 
                                 group_by(cluster) %>% 
                                 summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
                          mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>% 
                          #filter(len<450) %>% 
                          ggplot(aes((IDcluster), len)) +
                          geom_line(aes(group=cluster), alpha=0.4)+
                          geom_jitter(aes(colour = as.factor(cluster)), size=1, alpha=0.8) +
                          xlab("Cluster ID") +
                          ylab("Length (bp)") +
                          scale_color_manual(values=c(rep(c("black", "gray"),max(ordered_table2 cluster)-1), "black")) +
                          theme_classic() +
                          theme(axis.text.x=element_blank(), 
                                axis.text.y=element_text(colour="black", size=14), 
                                legend.position="none") +
                          facet_grid(. ~ chr, scale="free", space="free") #+
                        ggsave("01_size_distribution_ltr_SoloPlusComplete.png", width = 12, height = 4)
                        
                        
                        # plot heat map with number of copies per cluster and sample:
                        ordered_table2 %>% 
                          filter(seq_ID!="NA") %>% 
                          filter(len<500) %>% 
                          merge((ordered_table2 %>% 
                                 filter(seq_ID!="NA") %>% 
                                 group_by(cluster) %>% 
                                 summarize(IDcluster=paste(cluster[1],min(start),sep="_"))), by="cluster") %>% 
                          mutate(IDcluster = fct_reorder(IDcluster, cluster)) %>%
                          filter(len<450) %>% 
                          group_by(chr, IDcluster, sample) %>% 
                          summarize(n_copies=length(seq_ID), 
                                    dif_length=(max(len)-min(len))) %>% 
                          mutate(text_label=ifelse(dif_length==0,"",as.character(dif_length))) %>% 
                          ggplot(aes(IDcluster, sample)) +
                          geom_tile(aes(fill=as.factor(n_copies))) +
                          #scale_fill_manual(values=c("#FFDB6D", "#C4961A", "#F4EDCA", 
                          #"#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")) +
                          #geom_text(aes(label=text_label)) +
                          scale_fill_brewer(palette="Paired") +
                          xlab("Cluster ID") +
                          ylab("Sample") +
                          theme_classic() +
                          theme(axis.text.x=element_blank(), 
                                axis.text.y=element_text(colour="black", size=12), 
                                legend.position="top") +
                          facet_grid(. ~ chr, scale="free", space="free") #+
                        ggsave("02_heatmap_ltr_SoloPlusComplete_nCopies_persample.png", width = 12, height = 8)
                        
                        
                        ordered_table2 %>% 
                          filter(len<500) %>% 
                          merge((ordered_table2 %>% 
                                 group_by(cluster) %>% 
                                 summarize(min_start=min(start))), by="cluster") %>% 
                          group_by(cluster, chr, min_start, sample) %>% 
                          summarize(n_copies=n()) %>% 
                          group_by(cluster, chr, min_start) %>% 
                          summarize(n_samples=n()) %>% 
                          mutate(min_start=min_start/1000000) %>% 
                          ggplot(aes(min_start, n_samples)) +
                          geom_point(size=2, alpha=0.5) +
                          scale_x_continuous(breaks=seq(0,5.5,0.5), labels=seq(0,5.5,0.5))+
                          scale_y_continuous(breaks=seq(0,20,2), labels=seq(0,20,2))+
                          xlab("Pos.") +
                          ylab("Number Samples") +
                          theme_classic() +
                          theme(axis.text.x=element_text(size=12, angle = 45, hjust = 1, colour="black"), 
                                axis.text.y=element_text(colour="black", size=14)) +
                          facet_grid(. ~ chr, scale="free", space="free") #+
                        ggsave("03_SamplesperCluster_ltr_SoloPlusComplete.png", width = 12, height = 3)
                        
                        ordered_table2 %>% 
                          filter(completed_ltr==0) %>% 
                          #filter(len<500) %>% 
                          merge((ordered_table2 %>% 
                                 group_by(cluster) %>% 
                                 summarize(min_start=min(start))), by="cluster") %>% 
                          group_by(cluster, sample) %>% 
                          summarize(n_repeats=n()) %>% 
                          group_by(cluster) %>% 
                          summarize(n_samples=n()) %>% 
                          # filter(n_samples==1) %>% 
                          ggplot(aes(n_samples)) +
                          geom_histogram(binwidth=1) +
                          # geom_histogram(aes(y=..density..), binwidth=1) +
                          scale_x_continuous(breaks=seq(0,18,1), labels=seq(0,18,1)) +
                          ylab("Count") +
                          # ylab("Density") +
                          xlab("Number Samples") +
                          theme_classic() +
                          theme(axis.text.x=element_text(size=14, colour="black"), 
                                axis.text.y=element_text(colour="black", size=14), 
                                axis.title=element_text(colour="black", size=14)) #+
                        ggsave("04_fq_homologous_ltf_Solo_counts.png", width = 7, height = 5)
                        ggsave("04_fq_homologous_ltr_SoloPlusComplete_density.png", width = 7, height = 5)
                        
                        
                        
                        
                        library("ggtree")
                        tree_ed<-read.iqtree("../phylogenies/phy_ltr_completePlusSolo_minLen200.treefile")
                        #tree_ed<- as.phylo(tree)
                        fig_tree<- ggtree(tree_ed) + 
                          #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
                          #geom_tiplab() + 
                          #geom_tippoint() + 
                          theme_tree2() +
                          geom_treescale() 
                        #xlim(0, 0.5)
                        fig_tree
                        
                        ordered_table2 %>% head()
                        
                        # clusters:
                        cluster_summary_table<-ordered_table2 %>% 
                          filter(seq_ID!="NA") %>% 
                          group_by(cluster) %>% 
                          summarize(chr=chr[1],
                                    start=min(start),
                                    total_seq=n(), 
                                    total_samples=length(unique(sample)), 
                                    min_length=min(len), 
                                    max_length=max(len)) %>% 
                          filter(total_samples>=4)
                        
                        head(cluster_summary_table)
                        
                        library(svglite)
                        for(i in seq(1,dim(cluster_summary_table)[1])){
                          # i=5
                          cluster_number<-as.vector(unlist(cluster_summary_table[i,"cluster"]))
                          cluster_chr<-as.vector(unlist(cluster_summary_table[i,"chr"]))
                          cluster_start<-as.vector(unlist(cluster_summary_table[i,"start"]))
                          cluster_total_samples<-as.vector(unlist(cluster_summary_table[i,"total_samples"]))
                          cluster_total_seq<-as.vector(unlist(cluster_summary_table[i,"total_seq"]))
                          cluster_min_length<-as.vector(unlist(cluster_summary_table[i,"min_length"]))
                          cluster_max_length<-as.vector(unlist(cluster_summary_table[i,"max_length"]))
                          cluster_table<-ordered_table2 %>% 
                            mutate(cluster_set=factor(ifelse(seq_ID %in% c(ordered_table2 %>% 
                                                                           filter(seq_ID!="NA") %>% 
                                                                           filter(cluster==cluster_number) %>% 
                                                                           select(seq_ID) %>% 
                                                                           unlist() %>% 
                                                                           as.vector()), 1,0), levels=c(0,1))) %>% 
                            mutate(anc_group=factor(ifelse(anc_prop>0.9,"Sp",ifelse(anc_prop<0.1,"Sk","Hyb")), 
                                                    levels=c("Sp", "Hyb", "Sk"))) %>% 
                            select(seq_ID, cluster_set, anc_group) 
                          
                          fig_tree %<+% cluster_table + 
                            geom_tippoint(aes(colour=anc_group, alpha=factor(cluster_set)), size=4) +
                            scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                            scale_alpha_manual(values=c(0,0.6,0))+
                            labs(title=paste0("Cluster:",cluster_number," Chr:", cluster_chr, " Pos:", cluster_start, 
                                              " NSamples:", cluster_total_samples, 
                                              " NSeq:", cluster_total_seq, 
                                              " MinLen:", cluster_min_length, 
                                              " MaxLen:", cluster_max_length))+
                            theme(legend.position = "none") +
                            ggsave(paste0("./ltr_completePlusSolo_clusters/cluster",cluster_number,".svg"))
                          fig_tree %<+% cluster_table + 
                            geom_tippoint(aes(colour=anc_group, alpha=factor(cluster_set)), size=4) +
                            scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                            scale_alpha_manual(values=c(0,0.6,0))+
                            labs(title=paste0("Cluster:",cluster_number," Chr:", cluster_chr, " Pos:", cluster_start, 
                                              " NSamples:", cluster_total_samples, 
                                              " NSeq:", cluster_total_seq))+
                            theme(legend.position = "none") +
                            ggsave(paste0("./ltr_completePlusSolo_clusters/cluster",cluster_number,".png"))
                        }
                        
                        ## plots per sample:
                        
                        
                        cluster_summary_per_sample_table<-ordered_table2 %>%
                          #filter(sample=="ref") %>% 
                          mutate(sample_used=factor(ifelse(sample=="ref", 1,0)),
                                 anc_group=factor(ifelse(anc_prop>0.9,"Sp",ifelse(anc_prop<0.1,"Sk","Hyb")), 
                                                  levels=c("Sp", "Hyb", "Sk")), 
                                 Solo_LTR=factor(ifelse(len_froSeq>1500, 0, 1))) %>% 
                          select(seq_ID, sample, Solo_LTR) %>%
                          spread(sample, value = Solo_LTR) 
                        
                        row.names(cluster_summary_per_sample_table) <- cluster_summary_per_sample_table$seq_ID
                        head(cluster_summary_per_sample_table)
                        
                        cluster_summary_per_sample_table <- cluster_summary_per_sample_table %>% 
                          select(-seq_ID)
                        
                        # library(RColorBrewer)
                        # 
                        # brewer.pal(n = 8, name = 'Dark2')
                        # 
                        # fig_tree %<+% cluster_summary_per_sample_table + 
                        #   geom_tippoint(aes(colour=sample, alpha=Sp_group, size=4)) +
                        #   # scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                        #   scale_color_manual(values=rep(c(brewer.pal(n = 8, name = 'Dark2')),3)) +
                        #   # scale_color_brewer(palette = "Dark2")+
                        #   scale_alpha_manual(values=c(0,0.5))+
                        #   # scale_size_manual(values=c(1,0.5))+
                        #   theme(legend.position = "top")
                        
                        
                        genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
                        genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)
                        colnames(genotype) <- sub("\\.$", "", colnames(genotype))
                        p <- ggtree(beast_tree, mrsd="2013-01-01") + geom_treescale(x=2008, y=1, offset=2)
                        p <- p + geom_tiplab(size=2)
                        gheatmap(p, genotype, offset=5, width=0.5, font.size=3, 
                                 colnames_angle=-45, hjust=0) +
                          scale_fill_manual(breaks=c("HuH3N2", "pdm", "trig"), 
                                            values=c("steelblue", "firebrick", "darkgreen"), name="genotype")
                        
                        gheatmap(fig_tree, cluster_summary_per_sample_table, offset=5, width=0.5, font.size=3, 
                                 colnames_angle=-45, hjust=0) +
                          scale_fill_manual(breaks=c("0", "1", "NA"), 
                                            values=c("Orange", "Blue", "white"), name="LTR")
                        
                        fig_tree %<+% cluster_summary_per_sample_table +
                          geom_tippoint(aes(x=sample, colour=Solo_LTR, alpha=Solo_LTR), size=4) +
                          scale_alpha_manual(values=c(0.8,0.2))
                        #facet_grid(. ~ sample)
                        scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                          
                          labs(title=paste0("Cluster:",cluster_number," Chr:", cluster_chr, " Pos:", cluster_start, 
                                            " NSamples:", cluster_total_samples, 
                                            " NSeq:", cluster_total_seq, 
                                            " MinLen:", cluster_min_length, 
                                            " MaxLen:", cluster_max_length))+
                          theme(legend.position = "none") +
                          ggsave(paste0("./ltr_completePlusSolo_perSample/sample",sample_used,".svg"))
                        
                        for(i in seq(1,dim(cluster_summary_table)[1])){
                          # i=5
                          cluster_number<-as.vector(unlist(cluster_summary_table[i,"cluster"]))
                          cluster_chr<-as.vector(unlist(cluster_summary_table[i,"chr"]))
                          cluster_start<-as.vector(unlist(cluster_summary_table[i,"start"]))
                          cluster_total_samples<-as.vector(unlist(cluster_summary_table[i,"total_samples"]))
                          cluster_total_seq<-as.vector(unlist(cluster_summary_table[i,"total_seq"]))
                          cluster_min_length<-as.vector(unlist(cluster_summary_table[i,"min_length"]))
                          cluster_max_length<-as.vector(unlist(cluster_summary_table[i,"max_length"]))
                          cluster_table<-ordered_table2 %>% 
                            mutate(cluster_set=factor(ifelse(seq_ID %in% c(ordered_table2 %>% 
                                                                           filter(seq_ID!="NA") %>% 
                                                                           filter(cluster==cluster_number) %>% 
                                                                           select(seq_ID) %>% 
                                                                           unlist() %>% 
                                                                           as.vector()), 1,0), levels=c(0,1))) %>% 
                            mutate(anc_group=factor(ifelse(anc_prop>0.9,"Sp",ifelse(anc_prop<0.1,"Sk","Hyb")), 
                                                    levels=c("Sp", "Hyb", "Sk"))) %>% 
                            select(seq_ID, cluster_set, anc_group) 
                          
                          fig_tree %<+% cluster_table + 
                            geom_tippoint(aes(colour=anc_group, alpha=factor(cluster_set)), size=4) +
                            scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                            scale_alpha_manual(values=c(0,0.6,0))+
                            labs(title=paste0("Cluster:",cluster_number," Chr:", cluster_chr, " Pos:", cluster_start, 
                                              " NSamples:", cluster_total_samples, 
                                              " NSeq:", cluster_total_seq, 
                                              " MinLen:", cluster_min_length, 
                                              " MaxLen:", cluster_max_length))+
                            theme(legend.position = "none") +
                            ggsave(paste0("./ltr_completePlusSolo_perSample/sample",sample_used,".svg"))
                          fig_tree %<+% cluster_table + 
                            geom_tippoint(aes(colour=anc_group, alpha=factor(cluster_set)), size=4) +
                            scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                            scale_alpha_manual(values=c(0,0.6,0))+
                            labs(title=paste0("Cluster:",cluster_number," Chr:", cluster_chr, " Pos:", cluster_start, 
                                              " NSamples:", cluster_total_samples, 
                                              " NSeq:", cluster_total_seq))+
                            theme(legend.position = "none") +
                            ggsave(paste0("./ltr_completePlusSolo_perSample/sample",sample_used,".png"))
                        }
                        
                        
                        
                        
                        tree_ed@data$UFboot[is.na(tree_ed@data$UFboot)]<-0
                        tree_ed@data$UFboot[tree_ed@data$UFboot<90]<-0
                        tree_ed@data$UFboot[tree_ed@data$UFboot>=90 & tree_ed@data$UFboot<95]<-90
                        tree_ed@data$UFboot[tree_ed@data$UFboot>=95]<-95
                        
                        cluster_ID<-data.frame(id=tree_ed@phylo$tip.label) %>% 
                          merge(ordered_table2 %>% 
                                mutate(id=seq_ID) %>%
                                select(id, cluster, chr, anc_prop), by="id", all.x=T) %>% 
                          mutate(anc_group=factor(ifelse(anc_prop>0.9,"Sp",ifelse(anc_prop<0.1,"Sk","Hyb")), 
                                                  levels=c("Sp", "Hyb", "Sk"))) %>%
                          arrange(match(id, tree_ed@phylo$tip.label))
                        
                        row.names(cluster_ID) <- NULL
                        
                        lables_tf<-ordered_table2 %>% 
                          select(chr, cluster) %>%
                          group_by(chr) %>% 
                          dplyr::summarise(min_clus=min(cluster), 
                                           max_clus=max(cluster), 
                                           #multiple30=ifelse(cluster %% 30 == 0, cluster, 0)) %>% 
                                           mid_1=((max_clus-min_clus)/4)+min_clus, 
                                           mid_3=((max_clus-min_clus)*3/4)+min_clus, 
                                           middle=round((max_clus-min_clus)/2)+min_clus) %>%
                          ungroup() %>% 
                          select(-chr, -min_clus, -max_clus) %>%  
                          #filter(multiple30!=0) %>% 
                          #select(-chr) %>%  
                          unlist() %>% 
                          as.vector() %>% 
                          round() %>% 
                          sort() %>% 
                          as.numeric()
                        
                        p <-ggtree(tree_ed) +
                          #aes(color=as.factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)))), 
                          #
                          theme_tree2() +
                          geom_nodepoint(aes(alpha=as.factor(UFboot)), color="darkgreen", size=2) +
                          #scale_alpha_manual(values=c(0,0.5,1)) + 
                          # geom_nodepoint(aes(color=factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)), levels=c(90,80,0)), 
                          #                    alpha=factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0))), levels=c(90,80,0)), 
                          #                               size=3) +
                          geom_treescale() +
                          #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
                          #geom_tiplab() + 
                          #geom_tippoint(color="black", shape=4, size=1.5, alpha=0.2) +
                          #geom_tree() +
                          # geom_nodepoint(aes(color=as.factor(ifelse(UFboot>90,90,ifelse(UFboot>80,80,0)))), 
                          #                alpha=0.5, 
                          #                size=3) +
                          #geom_tree(data=df, color='firebrick')
                          #scale_color_continuous(low='darkgreen', high='red') +
                          #scale_color_manual(values=c("gray", "black", "darkred")) +
                          #scale_color_manual(values=c("white", "white", "darkred")) +
                        theme(legend.position="right")
                        #geom_tiplab(align=TRUE, linesize=1, col="gray90") 
                        # xlim(0, 0.4) 
                        
                        # pp <- p %<+% cluster_ID + 
                        #   geom_treescale()+
                        #   geom_tiplab(align = TRUE, linesize = 0, alpha=0) +
                        #   geom_tippoint(aes(colour=anc_group,
                        #                     alpha=anc_group), 
                        #                     size=2) +
                        #   scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                        #   scale_alpha_manual(values=c(0.8,0.0,0.8)) +
                        #   theme(legend.position = "None")
                        
                        pp<-p
                        
                        pp
                        
                        # cluster_ID<-data.frame(id=tree_ed$tip.label) %>% 
                        #   merge(ordered_table2 %>% 
                        #           mutate(id=seq_ID) %>%
                        #           select(id, cluster, chr), by="id", all.x=T) %>% 
                        #   arrange(match(id, tree_ed$tip.label))
                        
                        
                        # generate some random values for each tip label in the data
                        # Make a second plot with the original, naming the new plot "dot", 
                        # using the data you just created, with a point geom.
                        
                        p2 <- facet_plot(pp, panel="I", data=cluster_ID %>% filter(chr=="I"), 
                                         geom=geom_point, aes(x=cluster, colour=anc_group, 
                                                              alpha=anc_group, 
                                                              size=anc_group))+ 
                          #color='darkred', 
                          #alpha=0.8) +
                          scale_color_manual(values=c("darkred", "orange", "steelblue")) +
                          scale_alpha_manual(values=c(0,0.3,1,0.5, 0.9, 0.9)) + 
                          scale_size_manual(values=c(2,1.5,2)) +
                          theme_classic() +
                          scale_x_continuous(breaks = lables_tf, labels = lables_tf) +
                          scale_y_continuous(breaks = seq(1,800,30)) +
                          theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
                                panel.grid.major = element_line(colour = "gray80"), 
                                axis.text.x=element_text(size=12, angle = 45, hjust = 1, colour="black"), 
                                axis.text.y=element_blank(), 
                                axis.ticks.y = element_blank(), 
                                legend.position="none")
                        #panel.grid.minor = element_line(colour = "gray80"))
                        p2
                        p3 <- facet_plot(p2, panel="II", data=cluster_ID %>% filter(chr=="II"), 
                                         geom=geom_point, aes(x=cluster, colour=anc_group, 
                                                              alpha=anc_group, 
                                                              size=anc_group))
                        #color='darkred', 
                        # alpha=0.8) 
                        p4 <- facet_plot(p3, panel="III", data=cluster_ID %>% filter(chr=="III"), 
                                         geom=geom_point, aes(x=cluster, colour=anc_group, 
                                                              alpha=anc_group, 
                                                              size=anc_group))
                        #color='darkred', 
                        # alpha=0.8) 
                        #p4 
                        
                        
                        
                        library(ggtree)
                        library(reshape2)
                        
                        # load the packaged
                        library(grid)
                        library(gtable)
                        
                        head(cluster_ID)
                        
                        panel_size<-cluster_ID %>% 
                          filter(!(is.na(chr))) %>%
                          group_by(chr) %>%
                          summarise(max_cluster=max(cluster), 
                                    min_cluster=min(cluster)) %>% 
                          mutate(size_ch=max_cluster-min_cluster) %>% 
                          mutate(total=cluster_ID %>% 
                                 filter(!(is.na(chr))) %>% 
                                 select(cluster) %>% 
                                 unlist() %>% 
                                 as.vector() %>% 
                                 max()) %>% 
                          mutate(size_block=0.4*size_ch/total)
                        
                        
                        gt = ggplot_gtable(ggplot_build(p4))
                        #gtable_show_layout(gt) # will show you the layout - very handy function
                        #gt # see plot layout in table format
                        #gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
                        #gt$layout$l[grep('panel-3', gt$layout$name)] # you want to find the column specific to panel-3
                        #gt$layout$l[grep('panel-4', gt$layout$name)] # you want to find the column specific to panel-4
                        gt$widths[7] = (unlist(panel_size[panel_size$chr=="I","size_block"]))*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
                        gt$widths[9] = (unlist(panel_size[panel_size$chr=="II","size_block"]))*gt$widths[9]
                        gt$widths[11] = (unlist(panel_size[panel_size$chr=="III","size_block"]))*gt$widths[11]
                        
                        #grid.draw(gt) # plot with grid draw
                        
                        svg("05_tree_distribution_LTR_completed.svg",height = 250/25.4, width = 180/25.4, pointsize = 7)
                        grid.draw(gt) # plot with grid draw
                        dev.off()
                        
                        
                  
                        
                        