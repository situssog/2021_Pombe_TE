#!/usr/bin/env Rscript
rm(list=ls())


library("tidyr")
library("dplyr")
library("ggplot2")


args <- commandArgs(TRUE)


# file_blast_table<-"blast_ltrInTFseq_1500.txt"
# file_wtf_seq<-"wtfSeq_ID_1500.txt"
# file_SPNCRNASeq<-"SPNCRNASeq_ID_1500.txt"
# file_tf_seq_ID<-"tf_seqID_1500.txt"
# file_tf_removed_seq_ID<-"tf_seqID_excluded_1500.txt"
# output1<-"annotation_tf_breakPoints_1500.txt"
# output2<-"extraTF_Seq_1500.txt"


file_blast_table<-args[1] # ex: blast_ltrInTFseq.txt
file_wtf_seq<-args[2] # ex: wtfSeq_ID.txt
file_SPNCRNASeq<-args[3] # ex: SPNCRNASeq_ID.txt
file_tf_seq_ID<-args[4] # ex: tf_seqID.txt
file_tf_removed_seq_ID<-args[5] # ex: tf_seqID_excluded.txt
output1<-args[6] # ex: annotation_tf_breakPoints.txt
output2<-args[7] # ex: extraTF_Seq.txt



table_blast<-read.table(file_blast_table, F)
names(table_blast)<-c("qseqid", "sseqid", "pident", "length", "mismatch", "gaps", "qstart", "qend", "sstart", "send", "slen", "qlen", "evalue", "bitscore")

wtf_seq<-scan(file_wtf_seq, what="", sep="\n")
SPNCRNA_seq<-scan(file_SPNCRNASeq, what="", sep="\n")
removed_seq_ID<-scan(file_tf_removed_seq_ID, what="", sep="\n")
exclused_seq<-unique(c(wtf_seq, SPNCRNA_seq, removed_seq_ID))

tfseqID<-scan(file_tf_seq_ID, what="", sep="\n")

extra_seq<-tfseqID[!(tfseqID %in% table_blast$qseqid) & !(tfseqID %in% exclused_seq)]

table_blast<-table_blast %>% 
  filter(!(qseqid %in% exclused_seq))

# table_blast %>%
#   mutate(dif=qend-qstart) %>%
#   ggplot(aes(dif)) +
#   geom_histogram()


whole_seq<-c()
new_ranges_scafold<-c()
new_ranges_start<-c()
new_ranges_end<-c()
new_ranges_size_block<-c()
new_ranges_n<-c()

for (scaffold in unique(table_blast$qseqid)){
  # print(paste0("run: ", scaffold))
  #print(n)
  #pos<-1
  sub_table<-table_blast %>% 
  filter(qseqid==scaffold) %>% 
  arrange(qstart)
  #print(dim(sub_table)[1])
##  if (dim(sub_table)[1]<3){
##    whole_seq<-c(whole_seq, scaffold)
##  } else {
    #print("more rows")
    distances<-c()
    n=1
    for (line in seq(1,dim(sub_table)[1])){
      if (line==1){
        distances<-sub_table[line,"qstart"]
        if (distances>300){
          new_ranges_scafold<-c(new_ranges_scafold,scaffold)
          new_ranges_start<-c(new_ranges_start,1)
          new_ranges_end<-c(new_ranges_end,sub_table[line,"qend"])
          new_ranges_n<-c(new_ranges_n, n)
          new_ranges_size_block<-c(new_ranges_size_block, sub_table[line,"qend"])
          n=n+1
        } else {
          #print("")
        }
      } else {
        #print(scaffold)
        distances<-sub_table[line,"qstart"]-sub_table[line-1,"qend"]
        if (distances>300){
		      if (sub_table[line-1,"qstart"]>150){
			     temp_new_ranges_start<-sub_table[line-1,"qstart"]
		      } else { 
			     temp_new_ranges_start<-1
		      }
		      if (sub_table$qlen[1]-sub_table[line,"qend"]>150){
			     temp_new_ranges_end<-sub_table[line,"qend"]
		      } else {
			     temp_new_ranges_end<-sub_table$qlen[1]
		      }
		      new_ranges_scafold<-c(new_ranges_scafold,scaffold)
        	     #new_ranges_start<-c(new_ranges_start,sub_table[line-1,"qstart"])
		      new_ranges_start<-c(new_ranges_start,temp_new_ranges_start)
        	     new_ranges_end<-c(new_ranges_end,temp_new_ranges_end)
        	     new_ranges_n<-c(new_ranges_n,n)
        	     new_ranges_size_block<-c(new_ranges_size_block, (temp_new_ranges_end-temp_new_ranges_start+1))
        	     n=n+1
        } else {
        	#print("")
        }
      }
    }
    if (sub_table$qlen[1]-sub_table[dim(sub_table)[1],"qend"]>300){
      #print("end rows")
      new_ranges_scafold<-c(new_ranges_scafold,scaffold)
      new_ranges_start<-c(new_ranges_start,sub_table[dim(sub_table)[1],"qstart"])
      new_ranges_end<-c(new_ranges_end,sub_table$qlen[1])
      new_ranges_n<-c(new_ranges_n,n)
      new_ranges_size_block<-c(new_ranges_size_block, sub_table$qlen[1]-sub_table[dim(sub_table)[1],"qstart"])
    }
}




table_blast %>% 
  filter(qseqid %in% whole_seq) %>% 
  ungroup() %>% 
  mutate(new_ranges_scafold=qseqid, 
          new_ranges_scafold_n=paste0(qseqid,"_1"), 
          new_ranges_start=1, 
          new_ranges_end=qlen, 
          new_ranges_size_block=qlen) %>% 
  select(new_ranges_scafold, 
          new_ranges_scafold_n, 
          new_ranges_start, 
          new_ranges_end, 
          new_ranges_size_block) %>% 
  unique() %>% 
  rbind(data.frame(new_ranges_scafold, new_ranges_start, new_ranges_end, new_ranges_size_block, new_ranges_n) %>% 
  mutate(new_ranges_scafold_n=paste0(new_ranges_scafold, "_", new_ranges_n)) %>% 
  select(new_ranges_scafold, 
          new_ranges_scafold_n, 
          new_ranges_start, 
          new_ranges_end, 
          new_ranges_size_block)) %>% 
  write.table(output1, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


data.frame(extra_seq) %>% 
  write.table(output2, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)







