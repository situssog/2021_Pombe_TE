#!/usr/bin/env Rscript
rm(list=ls())


library("ggplot2")
library("tidyr")
#library("plyr")
library("dplyr")

samples_table<-read.table("all_annotation_refCoor_Masked_conSeq.txt", F)
names(samples_table)<-c("or_chr", "or_start", "or_end", "seq_ID", "ref_start_chr", "ref_start", "dis_start", "ref_end_chr", "ref_end", "dis_end")

samples_table<-samples_table %>% 
  mutate(length_ori_seq=abs(or_end-or_start))

direction<-c()
concordant_chr<-c()
chr_start_ed<-c()
pos_start_ed<-c()
dis_strat_ed<-c()
for (line in seq(1,dim(samples_table)[1])){
  if(as.character(samples_table[line,"ref_start_chr"])=="NoLoc"){
    direction<-c(direction, "-")
    concordant_chr<-c(concordant_chr, "no_start")
    chr_start_ed[line]<-as.vector(samples_table[line,"ref_end_chr"])
    pos_start_ed[line]<-samples_table[line,"ref_end"]
    dis_strat_ed[line]<-samples_table[line,"dis_end"]
  } else if(as.character(samples_table[line,"ref_end_chr"])=="NoLoc"){
    direction<-c(direction, "+")
    concordant_chr<-c(concordant_chr, "no_end")
    chr_start_ed[line]<-as.vector(samples_table[line,"ref_start_chr"])
    pos_start_ed[line]<-samples_table[line,"ref_start"]
    dis_strat_ed[line]<-samples_table[line,"dis_start"]
  } else if(as.character(samples_table[line,"ref_start_chr"])==as.character(samples_table[line,"ref_end_chr"])){
    if(samples_table[line,"ref_start"]<=samples_table[line,"ref_end"]){
      direction<-c(direction, "+")
      concordant_chr<-c(concordant_chr, 1)
      chr_start_ed[line]<-as.vector(samples_table[line,"ref_start_chr"])
      pos_start_ed[line]<-samples_table[line,"ref_start"]
      dis_strat_ed[line]<-samples_table[line,"dis_start"]
    } else {
      direction<-c(direction, "-")
      concordant_chr<-c(concordant_chr, 1)
      chr_start_ed[line]<-as.vector(samples_table[line,"ref_end_chr"])
      pos_start_ed[line]<-samples_table[line,"ref_end"]
      dis_strat_ed[line]<-samples_table[line,"dis_end"]
    }
  } else {
    if(samples_table[line,"dis_start"]<=samples_table[line,"dis_end"]){
      direction<-c(direction, "+")
      concordant_chr<-c(concordant_chr, 0)
      chr_start_ed[line]<-as.vector(samples_table[line,"ref_start_chr"])
      pos_start_ed[line]<-samples_table[line,"ref_start"]
      dis_strat_ed[line]<-samples_table[line,"dis_start"]
    } else {
      direction<-c(direction, "-")
      concordant_chr<-c(concordant_chr, 0)
      chr_start_ed[line]<-as.vector(samples_table[line,"ref_end_chr"])
      pos_start_ed[line]<-samples_table[line,"ref_end"]
      dis_strat_ed[line]<-samples_table[line,"dis_end"]
    }
  }
}

cbind(samples_table, direction, concordant_chr, chr_start_ed, pos_start_ed, dis_strat_ed) %>% 
  arrange(chr_start_ed, pos_start_ed, dis_strat_ed) %>% 
  write.table("all_annotation_refCoor_Masked_conSeq_ed.bed", quote = FALSE, row.names = FALSE, col.names = TRUE)

cbind(samples_table, direction, concordant_chr, chr_start_ed, pos_start_ed, dis_strat_ed) %>% 
  arrange(chr_start_ed, pos_start_ed, dis_strat_ed) %>% 
  mutate(location=paste(chr_start_ed, pos_start_ed, dis_strat_ed, length_ori_seq, sep="_")) %>% 
  select(seq_ID, location) %>% 
  write.table("all_annotation_refCoor_Masked_conSeq_edSim.bed", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

