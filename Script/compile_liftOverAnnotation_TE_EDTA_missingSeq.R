#!/usr/bin/env Rscript
rm(list=ls())

args <- commandArgs(TRUE)
sample=args[1]

#library("ggplot2")
library("tidyr")
#library("plyr")
library("dplyr")


samples_table<-read.table(paste0(sample, "_EDTA_missingSeq_TE_annotation__refCoor_Masked_conSeq.bed"), F)
names(samples_table)<-c("or_chr", "or_start", "or_end", "seq_ID", "ref_start_chr", "ref_start", "dis_start", "ref_end_chr", "ref_end", "dis_end", "direction")

samples_table<-samples_table %>%
  mutate(length_ori_seq=abs(or_end-or_start)) %>%
  select(or_chr, or_start, or_end, seq_ID, ref_start_chr, ref_start, dis_start, ref_end_chr, ref_end, dis_end, length_ori_seq, direction)

#direction<-c()
concordant_chr<-c()
chr_start_ed<-c()
pos_start_ed<-c()
dis_start_ed<-c()

chr_end_ed<-c()
pos_end_ed<-c()
dis_end_ed<-c()

for (line in seq(1,dim(samples_table)[1])){
  if(as.character(samples_table[line,"ref_start_chr"])=="NoLoc"){
    #direction<-c(direction, "-")
    concordant_chr<-c(concordant_chr, "no_start")
    chr_start_ed[line]<-as.vector(samples_table[line,"ref_end_chr"])
    pos_start_ed[line]<-samples_table[line,"ref_end"]
    dis_start_ed[line]<-samples_table[line,"dis_end"]
    chr_end_ed[line]<-"NA"
    pos_end_ed[line]<-pos_start_ed[line]+(max(samples_table[line,"or_start"], samples_table[line,"or_end"])-min(samples_table[line,"or_start"], samples_table[line,"or_end"]))
    dis_end_ed[line]<-0
  } else if(as.character(samples_table[line,"ref_end_chr"])=="NoLoc"){
    #direction<-c(direction, "+")
    concordant_chr<-c(concordant_chr, "no_end")
    chr_start_ed[line]<-as.vector(samples_table[line,"ref_start_chr"])
    pos_start_ed[line]<-samples_table[line,"ref_start"]
    dis_start_ed[line]<-samples_table[line,"dis_start"]
    chr_end_ed[line]<-"NA"
    pos_end_ed[line]<-pos_start_ed[line]+(max(samples_table[line,"or_start"], samples_table[line,"or_end"])-min(samples_table[line,"or_start"], samples_table[line,"or_end"]))
    dis_end_ed[line]<-0
  } else if(as.character(samples_table[line,"ref_start_chr"])==as.character(samples_table[line,"ref_end_chr"])){
    if(samples_table[line,"ref_start"]<=samples_table[line,"ref_end"]){
      #direction<-c(direction, "+")
      concordant_chr<-c(concordant_chr, 1)
      chr_start_ed[line]<-as.vector(samples_table[line,"ref_start_chr"])
      pos_start_ed[line]<-samples_table[line,"ref_start"]
      dis_start_ed[line]<-samples_table[line,"dis_start"]
      chr_end_ed[line]<-as.vector(samples_table[line,"ref_end_chr"])
      pos_end_ed[line]<-samples_table[line,"ref_end"]
      dis_end_ed[line]<-samples_table[line,"dis_end"]
    } else {
      #direction<-c(direction, "-")
      concordant_chr<-c(concordant_chr, 1)
      chr_start_ed[line]<-as.vector(samples_table[line,"ref_end_chr"])
      pos_start_ed[line]<-samples_table[line,"ref_end"]
      dis_start_ed[line]<-samples_table[line,"dis_end"]
      chr_end_ed[line]<-as.vector(samples_table[line,"ref_start_chr"])
      pos_end_ed[line]<-samples_table[line,"ref_start"]
      dis_end_ed[line]<-samples_table[line,"dis_start"]
    }
  } else {
    if(samples_table[line,"dis_start"]<=samples_table[line,"dis_end"]){
      #direction<-c(direction, "+")
      concordant_chr<-c(concordant_chr, 0)
      chr_start_ed[line]<-as.vector(samples_table[line,"ref_start_chr"])
      pos_start_ed[line]<-samples_table[line,"ref_start"]
      dis_start_ed[line]<-samples_table[line,"dis_start"]
      chr_end_ed[line]<-"NA"
      pos_end_ed[line]<-pos_start_ed[line]+(max(samples_table[line,"or_start"], samples_table[line,"or_end"])-min(samples_table[line,"or_start"], samples_table[line,"or_end"]))
      dis_end_ed[line]<-0
    } else {
      #direction<-c(direction, "-")
      concordant_chr<-c(concordant_chr, 0)
      chr_start_ed[line]<-as.vector(samples_table[line,"ref_end_chr"])
      pos_start_ed[line]<-samples_table[line,"ref_end"]
      dis_start_ed[line]<-samples_table[line,"dis_end"]
      chr_end_ed[line]<-"NA"
      pos_end_ed[line]<-pos_start_ed[line]+(max(samples_table[line,"or_start"], samples_table[line,"or_end"])-min(samples_table[line,"or_start"], samples_table[line,"or_end"]))
      dis_end_ed[line]<-0
    }
  }
}

cbind(samples_table, concordant_chr, chr_start_ed, pos_start_ed, dis_start_ed, chr_end_ed, pos_end_ed, dis_end_ed) %>%
  #arrange(chr_start_ed, pos_start_ed, dis_start_ed) %>%
  write.table(paste0(sample, "_EDTA_missingSeq_TE_annotation__refCoor_Masked_conSeq_ed.bed"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cbind(samples_table, concordant_chr, chr_start_ed, pos_start_ed, dis_start_ed) %>%
  #arrange(chr_start_ed, pos_start_ed, dis_start_ed) %>%
  mutate(location=paste(chr_start_ed, pos_start_ed, dis_start_ed, length_ori_seq, sep="_")) %>%
  select(seq_ID, location) %>%
  write.table(paste0(sample, "_EDTA_missingSeq_TE_annotation__refCoor_Masked_conSeq_edSim.bed"), sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


