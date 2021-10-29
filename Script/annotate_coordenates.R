#!/usr/bin/env Rscript
rm(list=ls())

args <- commandArgs(TRUE)

library("ggplot2")
library("tidyr")
library("plyr")
library("dplyr")

# sample="DY34373"
# translated_coor<-read.table("L_17_II__chain_table_DY34373_refCoor.bed",  F) 
# annotation_table<-read.table("DY34373_repeats_annotation_file.gff3", F) 

print("sample:")
print(args[1])
print("annotation file:")
print(args[2])
print("repeat annotations")
print(args[3])

#sample="EBC121"
#file="EBC121_tig00000007_I_94_638829__chain_table_EBC121_refCoor.bed"
#translated_coor<-read.table(file, F) 

sample=args[1]
translated_coor<-read.table(args[2], F) 

names(translated_coor)<-c("ref_chr", "ref_start", "ref_end", "sample_ID")

translated_coor<- translated_coor %>% 
  separate("sample_ID", c("sample_chr", "sample_pos"), "__") 


annotation_table<-read.table(args[3], F) 
names(annotation_table)<-c("sample_chr_ann", "EDTA", "type", "start_ann", "end_ann", "len", "dir", "point", "ID")

scalfold<-as.vector(translated_coor[1,"sample_chr"])

annotation_table<-annotation_table %>% 
  filter(sample_chr_ann==scalfold)

start_new_coor<-c()
end_new_coor<-c()
start_new_coor_chr<-c()
end_new_coor_chr<-c()
distance_start<-c()
distance_end<-c()
for (line in seq(1:dim(annotation_table)[1])){
  print(line)
  distance1<-0
  distance2<-0
  if (annotation_table[line,"start_ann"]<1){
      annotation_table[line,"start_ann"]<-1
    }
  max_pos<-translated_coor %>% 
          filter(sample_chr==annotation_table[line,"sample_chr_ann"]) %>% 
          ungroup() %>% 
          select(sample_pos) %>% 
          unlist() %>% 
          as.numeric() %>%
          as.vector() %>% 
          max()
  #if (annotation_table[line,"end_ann"]>max_pos){
    #annotation_table[line,"end_ann"]<-max_pos
  #}
  subLine_start<-translated_coor %>% 
    filter(sample_chr==annotation_table[line,"sample_chr_ann"] & 
      sample_pos==annotation_table[line,"start_ann"]-distance1)
  subLine_end<-translated_coor %>% 
    filter(sample_chr==annotation_table[line,"sample_chr_ann"] & 
      sample_pos==annotation_table[line,"end_ann"]+distance2)
  while (dim(subLine_start)[1]==0){
#    print(distance1)
    if (annotation_table[line,"start_ann"]-distance1<0){
      print("smaller coor")
      subLine_start<-data.frame(ref_chr=c("NoLoc"), ref_start=c("NoLoc"))
      break
    }
    distance1<-distance1+10
    subLine_start<-translated_coor %>% filter(sample_chr==annotation_table[line,"sample_chr_ann"] &  sample_pos==annotation_table[line,"start_ann"]-distance1)
  }
  if (dim(subLine_start)[1]>1){
    print("more lines start!")
  }
  start_new_coor_chr[line]<-as.vector(subLine_start[1,"ref_chr"])
  start_new_coor[line]<-subLine_start[1,"ref_start"]
  distance_start[line]<-distance1
  while (dim(subLine_end)[1]==0){
#    print(distance2)
    if ((annotation_table[line,"end_ann"]+distance2)>max_pos){
      print("larger coor")
      subLine_end<-data.frame(ref_chr=c("NoLoc"), ref_end=c("NoLoc"))
      break
    }
    distance2<-distance2+10
    subLine_end<-translated_coor %>% filter(sample_chr==annotation_table[line,"sample_chr_ann"] &  sample_pos==annotation_table[line,"end_ann"]+distance2)
  }
  if (dim(subLine_end)[1]>1){
    print("more lines end!")
  }
  end_new_coor_chr[line]<-as.vector(subLine_end[1,"ref_chr"])
  end_new_coor[line]<-subLine_end[1,"ref_end"]
  distance_end[line]<-distance2
  write.table(data.frame(annotation_table[line,], 
    start_new_coor_chr[line], start_new_coor[line], distance_start[line], 
    end_new_coor_chr[line], end_new_coor[line], distance_end[line]), paste0("temp_final_",sample,"_refCoor.txt"), 
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t", append=TRUE)
}

annotation_table %>% 
  cbind(start_new_coor_chr, start_new_coor, distance_start, end_new_coor_chr, end_new_coor, distance_end) %>% 
  write.table(paste0("final_",sample,"_refCoor.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t", append=TRUE)



