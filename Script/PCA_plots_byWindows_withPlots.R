#!/usr/bin/env Rscript
rm(list=ls())

args <- commandArgs(TRUE)


# module load R/3.4.3
# export R_LIBS_USER=/home/sergio/R/libraries/3.4
# source("https://bioconductor.org/biocLite.R")
# biocLite("gdsfmt",  lib="/home/sergio/R/libraries/3.4")
# biocLite("SNPRelate",  lib="/home/sergio/R/libraries/3.4")
# install.packages('ggplot2', lib="/home/sergio/R/libraries/3.4")
# biocLite("VariantAnnotation",  lib="/home/sergio/R/libraries/3.4")
# install.packages('ggrepel', lib="/home/sergio/R/libraries/3.4")
# install.packages('tidyr', lib="/home/sergio/R/libraries/3.4")
# install.packages('plyr', lib="/home/sergio/R/libraries/3.4")
# install.packages('dplyr', lib="/home/sergio/R/libraries/3.4")

## in rackham: 
#module load R/3.4.3 MariaDB/10.2.11
#R_LIBS_USER=/home/sergio/R/libraries_Rackham/3.4
#source("https://bioconductor.org/biocLite.R")
#biocLite("gdsfmt",  lib="/home/sergio/R/libraries_Rackham/3.4")
#biocLite("SNPRelate",  lib="/home/sergio/R/libraries_Rackham/3.4")
#install.packages('ggplot2', lib="/home/sergio/R/libraries_Rackham/3.4")
#biocLite("VariantAnnotation",  lib="/home/sergio/R/libraries_Rackham/3.4")
#install.packages('ggrepel', lib="/home/sergio/R/libraries_Rackham/3.4")
#install.packages('tidyr', lib="/home/sergio/R/libraries_Rackham/3.4")
#install.packages('plyr', lib="/home/sergio/R/libraries_Rackham/3.4")
#install.packages('dplyr', lib="/home/sergio/R/libraries_Rackham/3.4")
#install.packages('vcfR', lib="/home/sergio/R/libraries_Rackham/3.4")

library("gdsfmt")
library("SNPRelate")
library("ggplot2")
library("VariantAnnotation")
library("ggrepel")
library("tidyr")
library("plyr")
library("dplyr")

#install.packages('VariantAnnotation', lib="/home/sergio/R/x86_64-redhat-linux-gnu-library/3.4")
#biocLite("VariantAnnotation", lib="/home/sergio/R/x86_64-redhat-linux-gnu-library/3.4/test")
#install.packages('rJava', lib="/home/sergio/R/x86_64-redhat-linux-gnu-library/3.4")
#install.packages("rlang", lib="/home/sergio/R/x86_64-redhat-linux-gnu-library/3.4") 
#update.packages(lib.loc = "/home/sergio/R/x86_64-redhat-linux-gnu-library/3.4")

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sergio/R/x86_64-redhat-linux-gnu-library/3.4/mysql-5.7.21-linux-glibc2.12-x86_64/lib:/home/sergio/R/x86_64-redhat-linux-gnu-library/3.4/mysql-5.7.21-linux-glibc2.12-x86_64/bin:/home/sergio/R/x86_64-redhat-linux-gnu-library/3.4/mysql-5.7.21-linux-glibc2.12-x86_64/share


chromosome<-"I"
number_snps<-10
overlapping_windows<-5
vcf_file_name<-"chromosome_I_all57_clean_withGAPS.recode.vcf"
temporal_number<-runif(1,0,100)
#vcf <- readVcf(vcf_file_name)
start<-31
end <-40
window<-7

start<-1100
end <-1300
window<-11


chromosome<-"III"
number_snps<-200
overlapping_windows<-100
vcf_file_name<-"chromosome_III_all57_clean_withGAPS.recode.vcf"
temporal_number<-runif(1,0,100)
#vcf <- readVcf(vcf_file_name)
start<-10600
end <-10800
window<-106

start<-1100
end <-1300
window<-11

# for 57 sample
# no_mixed_group<-c("JB22", "JB938", "JB870", "JB760", "JB879", "JB869")
# for 161 samples
# no_mixed_group<-c("JB1183", "JB1195", "JB1181", "JB897", "JB759", "JB1184", "JB1182", "JB869", "JB938", "JB889", "JB885", "JB886", "JB888", "JB891", "JB760", "JB887", "JB940", "JB22", "JB945", "JB906", "JB941", "JB1170", "JB374", "JB50", "JB936", "JB761", "JB937", "JB892", "JB868", "JB1204", "JB1168", "JB1179", "JB879", "JB898", "JB1111", "JB870")

chromosome<-args[1]
number_snps<-as.numeric(args[2])
overlapping_windows<-as.numeric(args[3])
vcf_file_name<-args[4]

temporal_number<-runif(1,0,100)
vcf <- readVcf(vcf_file_name)
start<-1
end <-number_snps


varprop<-data.frame(chromosome_name=chromosome, window_number=0, 
                    start_pos=0, end_pos=0, 
                    pc1_varcom=0, 
                    pc2_varcom=0, 
                    pc3_varcom=0, 
                    pc4_varcom=0, 
                    pc5_varcom=0, 
                    pop1N=0, 
                    pop2N=0, 
                    dxy_adj=0, 
                    dxy_nom=0, 
                    dxy_fq=0, 
                    theta_pw_pop1=0, 
                    theta_pw_pop2=0, 
                    theta_pw_pop1_unco_pol=0, 
                    theta_pw_pop2_unco_pol=0, 
                    theta_w_pop1=0, 
                    theta_w_pop2=0, 
                    Taj_D_pop1=0, 
                    Taj_D_pop2=0, 
                    JB22_pop1=TRUE)
#                    pop_JB22_nomixedN=0, 
#                    pop_JB22_mixedN=0, 
#                    dxy_nom_hap_JB22=0,
#                    theta_pw_pop_JB22_nomixed=0,
#                    theta_pw_pop_JB22_mixed=0,
#                    theta_w_pop_JB22_nomixed=0,
#                    theta_w_pop_JB22_mixed=0, 
#                    JB869_pop_JB22_nomixed=TRUE, 
#                    pop_JB22_nomixed_samples="JB")


PC1_distribution<-data.frame(chromosome_name=chromosome, window_number=0, 
                    start_pos=0, end_pos=0, 
                    sample="JB000", 
                    PC1=0.0001, 
                    PC2=0.0001, 
                    max_PC1=0.0001, 
                    min_PC1=0.0001, 
                    Nor_PC1=0.0001,
                    Nor_PC1_pol=0.0001, 
                    max_PC2=0.0001, 
                    min_PC2=0.0001, 
                    Nor_PC2=0.0001)


#start<-1+(26500*overlapping_windows)
#end <-number_snps+(26500*overlapping_windows)
#
#for (window in seq(26501,floor(length(vcf)/overlapping_windows)-4)){
for (window in seq(1,floor(length(vcf)/overlapping_windows)-4)){
  print(window)
  new_vcf<-vcf[start:end]
  file_name<-paste0(chromosome, "_", window, "_", 
                    start(new_vcf[1]), "_", end(new_vcf[length(new_vcf)]), "_", 
                    end(new_vcf[length(new_vcf)])-start(new_vcf[1]))
  if(!file.exists(paste0("./vcf_files/",file_name, ".vcf"))){
    writeVcf(new_vcf, paste0("./vcf_files/",file_name, ".vcf"))
  }
  snpgdsVCF2GDS(paste0("./vcf_files/",file_name, ".vcf"), paste0("temporal_", temporal_number, ".gds"), method="biallelic.only")
  snpgdsSummary(paste0("temporal_", temporal_number, ".gds"))
  genofile <- snpgdsOpen(paste0("temporal_", temporal_number, ".gds"))
  if (length(read.gdsn(index.gdsn(genofile, "snp.rs.id")))!=0){
  ccm_pca<-snpgdsPCA(genofile,  autosome.only=FALSE)
  table_plot<-data.frame(PC1=ccm_pca$eigenvect[,1], PC2=ccm_pca$eigenvect[,2], sample=ccm_pca$sample)
  ggplot(table_plot, aes(x=PC1, y=PC2, colour=sample, label=sample)) +
    geom_point() + 
    geom_text_repel(aes(label = sample), box.padding = unit(0.1, "lines")) +
    theme(legend.position="none") +
    labs(x = paste0("PC1 ",round(ccm_pca$varprop[1], digits=3)), 
         y = paste0("PC2 ",round(ccm_pca$varprop[2], digits=3)), 
         title = file_name ) +
    ggsave(paste0("./PCA_byWindow_", number_snps, "_", overlapping_windows, "/", file_name, "_PC1_PC2", ".png"), width = 15, height = 10)
  PC1_component<-ccm_pca$eigenvect[,1]
  Nor_PC1<-(PC1_component-min(PC1_component))/(max(PC1_component)-min(PC1_component))
  pop1<-table_plot %>% mutate(Nor_PC1=Nor_PC1) %>% filter(Nor_PC1>0.75) %>% dplyr::select(sample) %>% unlist() %>% as.vector()
  pop2<-table_plot %>% mutate(Nor_PC1=Nor_PC1) %>% filter(Nor_PC1<0.25) %>% dplyr::select(sample) %>% unlist() %>% as.vector()
#  if ("JB22" %in% pop1){
#    pop_JB22_nomixed<-pop1[pop1 %in% no_mixed_group]
#    pop_JB22_mixed<-pop1[!(pop1 %in% no_mixed_group)]
#  } else if ("JB22" %in% pop2){
#    pop_JB22_nomixed<-pop2[pop2 %in% no_mixed_group]
#    pop_JB22_mixed<-pop2[!(pop2 %in% no_mixed_group)]
#  } else {
#    pop_JB22_nomixed<-c()
#    pop_JB22_mixed<-c()
#  }
  gpop1 <- snpgdsGetGeno(genofile, sample.id=pop1)
  gpop2 <- snpgdsGetGeno(genofile, sample.id=pop2)
#  if (length(pop_JB22_nomixed)>0){
#    gpop_JB22_nomixed <- snpgdsGetGeno(genofile, sample.id=pop_JB22_nomixed)
#  }
#  if (length(pop_JB22_mixed)>0){
#    gpop_JB22_mixed<- snpgdsGetGeno(genofile, sample.id=pop_JB22_mixed)
#  }
  #read.gdsn(index.gdsn(genofile, "sample.id"))
  #read.gdsn(index.gdsn(genofile, "snp.rs.id"))
  #read.gdsn(index.gdsn(genofile, "genotype"))
  #read.gdsn(index.gdsn(genofile, "snp.position"))
  #read.gdsn(index.gdsn(genofile, "snp.allele"))
  start_pos<-read.gdsn(index.gdsn(genofile, "snp.position"))[1]
  end_pos<-read.gdsn(index.gdsn(genofile, "snp.position"))
  end_pos<-end_pos[length(end_pos)]
  # dxy
  dif_pairs_adj<-c()
  for (i in seq(1,length(pop1))){
    for (j in seq(1,length(pop2))){
      distanceij<-sum(!(gpop1[i,]==gpop2[j,])*1)/(end_pos-start_pos)
      distanceij<-(-3/4)*(log(1-((3/4)*distanceij)))
      dif_pairs_adj<-c(dif_pairs_adj, (distanceij*(1/(length(pop1)))*(1/(length(pop2)))))
    }
  }
  dif_pairs_nom<-c()
  for (i in seq(1,length(pop1))){
    for (j in seq(1,length(pop2))){
      distanceij<-sum(!(gpop1[i,]==gpop2[j,])*1)/(end_pos-start_pos)
      dif_pairs_nom<-c(dif_pairs_nom, (distanceij*(1/(length(pop1)))*(1/(length(pop2)))))
    }
  }
  dxy_adj<-((length(pop1)+length(pop2))/(length(pop1)+length(pop2)-1))*(sum(dif_pairs_adj))
  dxy_nom<-((length(pop1)+length(pop2))/(length(pop1)+length(pop2)-1))*(sum(dif_pairs_nom))
  dif_pairs_fq<-c()
  for (i in seq(1,length(gpop1[1,]))){
    a_pop1<-sum(gpop1[,i])/length(pop1)
    b_pop1<-(1-a_pop1)
    a_pop2<-sum(gpop2[,i])/length(pop2)
    b_pop2<-(1-a_pop2)
    dif_pairs_fq<-c(dif_pairs_fq, ((a_pop1*b_pop2)+(a_pop2*b_pop1)))
  }
  dxy_fq<-sum(dif_pairs_fq)/(end_pos-start_pos)
  # pairwise theta within populations:
  # pop 1
  dif_pairs_theta_pop1<-c()
  for (i in seq(1,length(pop1))){
    for (j in seq(1,length(pop1))){
      if (j>i){
        dif_pairs_theta_pop1<-c(dif_pairs_theta_pop1, sum(!(gpop1[i,]==gpop1[j,])*1))
      }
    }
  }
  theta_pw_pop1<-(sum(dif_pairs_theta_pop1)/(end_pos-start_pos))/(((length(pop1)-1)*length(pop1))/2)
  # pop 2
  dif_pairs_theta_pop2<-c()
  for (i in seq(1,length(pop2))){
    for (j in seq(1,length(pop2))){
      if (j>i){
        dif_pairs_theta_pop2<-c(dif_pairs_theta_pop2, sum(!(gpop2[i,]==gpop2[j,])*1))
      }
    }
  }
  theta_pw_pop2<-(sum(dif_pairs_theta_pop2)/(end_pos-start_pos))/(((length(pop2)-1)*length(pop2))/2)
  # Watterson theta within populations:
  # pop1
  s_pop1<-c()
  # for (i in seq(1,length(gpop1[1,]))){
  #   s_pop1<-c(s_pop1, sum(gpop1[,i])!=length(pop1))
  # }
  for (i in seq(1,length(gpop1[1,]))){
    s_pop1<-c(s_pop1, !(sum(gpop1[,i]) %in% c(0,length(pop1))))
  }
  a1_pop1<-c()
  a2_pop1<-c()
  for (i in seq(1,length(pop1)-1)){
    a1_pop1<-c(a1_pop1, 1/i)
    a2_pop1<-c(a2_pop1, 1/(i^2))
  }
  a1_pop1<-sum(a1_pop1)
  a2_pop1<-sum(a2_pop1)
  b1_pop1<-(length(pop1)+1)/(3*(length(pop1)-1))
  b2_pop1<-(2*((length(pop1)^2)+length(pop1)+3))/(9*length(pop1)*(length(pop1)-1))
  c1_pop1<-b1_pop1-(1/a1_pop1)
  c2_pop1<-b2_pop1-((length(pop1)+2)/(a1_pop1*length(pop1)))+(a2_pop1/(a1_pop1^2))
  e1_pop1<-c1_pop1/a1_pop1
  e2_pop1<-c2_pop1/((a1_pop1^2)+a2_pop1)
  theta_w_pop1<-(sum(s_pop1)/(end_pos-start_pos))/a1_pop1
  #pop 2
  s_pop2<-c()
  # for (i in seq(1,length(gpop2[1,]))){
  #   s_pop2<-c(s_pop2, sum(gpop2[,i])!=length(pop2))
  # }
  for (i in seq(1,length(gpop2[1,]))){
    s_pop2<-c(s_pop2, !(sum(gpop2[,i]) %in% c(0,length(pop2))))
  }
  a1_pop2<-c()
  a2_pop2<-c()
  for (i in seq(1,length(pop2)-1)){
    a1_pop2<-c(a1_pop2, 1/i)
    a2_pop2<-c(a2_pop2, 1/(i^2))
  }
  a1_pop2<-sum(a1_pop2)
  a2_pop2<-sum(a2_pop2)
  b1_pop2<-(length(pop2)+1)/(3*(length(pop2)-1))
  b2_pop2<-(2*((length(pop2)^2)+length(pop2)+3))/(9*length(pop2)*(length(pop2)-1))
  c1_pop2<-b1_pop2-(1/a1_pop2)
  c2_pop2<-b2_pop2-((length(pop2)+2)/(a1_pop2*length(pop2)))+(a2_pop2/(a1_pop2^2))
  e1_pop2<-c1_pop2/a1_pop2
  e2_pop2<-c2_pop2/((a1_pop2^2)+a2_pop2)
  theta_w_pop2<-(sum(s_pop2)/(end_pos-start_pos))/a1_pop2
  # Tajima's D
  # pop 1
  theta_pw_pop1_unco<-(sum(dif_pairs_theta_pop1))/(((length(pop1)-1)*length(pop1))/2)
  theta_w_pop1_unco<-(sum(s_pop1))/a1_pop1
  Taj_D_pop1<-(theta_pw_pop1_unco- theta_w_pop1_unco)/(sqrt((e1_pop1*(sum(s_pop1)))+(e2_pop1*(sum(s_pop1))*(sum(s_pop1)-1))))
  # pop 2
  theta_pw_pop2_unco<-(sum(dif_pairs_theta_pop2))/(((length(pop2)-1)*length(pop2))/2)
  theta_w_pop2_unco<-(sum(s_pop2))/a1_pop2
  Taj_D_pop2<-(theta_pw_pop2_unco- theta_w_pop2_unco)/(sqrt((e1_pop2*(sum(s_pop2)))+(e2_pop2*(sum(s_pop2))*(sum(s_pop2)-1))))
  table_plot<-data.frame(PC1=ccm_pca$eigenvect[,1], PC3=ccm_pca$eigenvect[,3], sample=ccm_pca$sample)
  #
  #
#  # Dxy between parental JB22 non-mixed group and the same haplotype but mixed
#  if (length(pop_JB22_nomixed)>0 & length(pop_JB22_mixed)>0){
#  dif_pairs_nom_hap_JB22<-c()
#  for (i in seq(1,length(pop_JB22_nomixed))){
#    for (j in seq(1,length(pop_JB22_mixed))){
#      distanceij<-sum(!(gpop_JB22_nomixed[i,]==gpop_JB22_mixed[j,])*1)/(end_pos-start_pos)
#      dif_pairs_nom_hap_JB22<-c(dif_pairs_nom_hap_JB22, (distanceij*(1/(length(pop_JB22_nomixed)))*(1/(length(pop_JB22_mixed)))))
#    }
#  }
#  dxy_nom_hap_JB22<-((length(pop_JB22_nomixed)+length(pop_JB22_mixed))/(length(pop_JB22_nomixed)+length(pop_JB22_mixed)-1))*(sum(dif_pairs_nom_hap_JB22))
#  } else {
#    dxy_nom_hap_JB22<-"NA"
#  }
#  #pop_JB22_nomixed
#  if (length(pop_JB22_nomixed)>0){
#  dif_pairs_theta_pop_JB22_nomixed<-c()
#  for (i in seq(1,length(pop_JB22_nomixed))){
#    for (j in seq(1,length(pop_JB22_nomixed))){
#      if (j>i){
#        dif_pairs_theta_pop_JB22_nomixed<-c(dif_pairs_theta_pop_JB22_nomixed, sum(!(gpop_JB22_nomixed[i,]==gpop_JB22_nomixed[j,])*1))
#      }
#    }
#  }
#  theta_pw_pop_JB22_nomixed<-(sum(dif_pairs_theta_pop_JB22_nomixed)/(end_pos-start_pos))/(((length(pop_JB22_nomixed)-1)*length(pop_JB22_nomixed))/2)
#  } else {
#    theta_pw_pop_JB22_nomixed<-"NA"
#  }
#  #pop_JB22_mixed
#  if (length(pop_JB22_mixed)>0){
#  dif_pairs_theta_pop_JB22_mixed<-c()
#  for (i in seq(1,length(pop_JB22_mixed))){
#    for (j in seq(1,length(pop_JB22_mixed))){
#      if (j>i){
#        dif_pairs_theta_pop_JB22_mixed<-c(dif_pairs_theta_pop_JB22_mixed, sum(!(gpop_JB22_mixed[i,]==gpop_JB22_mixed[j,])*1))
#      }
#    }
#  }
#  theta_pw_pop_JB22_mixed<-(sum(dif_pairs_theta_pop_JB22_mixed)/(end_pos-start_pos))/(((length(pop_JB22_mixed)-1)*length(pop_JB22_mixed))/2)
#  } else {
#    theta_pw_pop_JB22_mixed<-"NA"
#  }
#  #
#  #
#  # Watterson theta within populations:
#  # pop_JB22_nomixed
#  if (length(pop_JB22_nomixed)>0){
#  s_pop_JB22_nomixed<-c()
#  for (i in seq(1,length(gpop_JB22_nomixed[1,]))){
#    s_pop_JB22_nomixed<-c(s_pop_JB22_nomixed, sum(gpop_JB22_nomixed[,i])!=length(pop_JB22_nomixed))
#  }
#  a1_pop_JB22_nomixed<-c()
#  a2_pop_JB22_nomixed<-c()
#  for (i in seq(1,length(pop_JB22_nomixed)-1)){
#    a1_pop_JB22_nomixed<-c(a1_pop_JB22_nomixed, 1/i)
#    a2_pop_JB22_nomixed<-c(a2_pop_JB22_nomixed, 1/(i^2))
#  }
#  a1_pop_JB22_nomixed<-sum(a1_pop_JB22_nomixed)
#  a2_pop_JB22_nomixed<-sum(a2_pop_JB22_nomixed)
#  b1_pop_JB22_nomixed<-(length(pop_JB22_nomixed)+1)/(3*(length(pop_JB22_nomixed)-1))
#  b2_pop_JB22_nomixed<-(2*((length(pop_JB22_nomixed)^2)+length(pop_JB22_nomixed)+3))/(9*length(pop_JB22_nomixed)*(length(pop_JB22_nomixed)-1))
#  c1_pop_JB22_nomixed<-b1_pop_JB22_nomixed-(1/a1_pop_JB22_nomixed)
#  c2_pop_JB22_nomixed<-b2_pop_JB22_nomixed-((length(pop_JB22_nomixed)+2)/(a1_pop_JB22_nomixed*length(pop_JB22_nomixed)))+(a2_pop_JB22_nomixed/(a1_pop_JB22_nomixed^2))
#  e1_pop_JB22_nomixed<-c1_pop_JB22_nomixed/a1_pop_JB22_nomixed
#  e2_pop_JB22_nomixed<-c2_pop_JB22_nomixed/((a1_pop_JB22_nomixed^2)+a2_pop_JB22_nomixed)
#  theta_w_pop_JB22_nomixed<-(sum(s_pop_JB22_nomixed)/(end_pos-start_pos))/a1_pop_JB22_nomixed
#  } else {
#    theta_w_pop_JB22_nomixed<-"NA"
#  }
#  #pop_JB22_mixed
#  if (length(pop_JB22_mixed)>0){
#  s_pop_JB22_mixed<-c()
#  for (i in seq(1,length(gpop_JB22_mixed[1,]))){
#    s_pop_JB22_mixed<-c(s_pop_JB22_mixed, sum(gpop_JB22_mixed[,i])!=length(pop_JB22_mixed))
#  }
#  a1_pop_JB22_mixed<-c()
#  a2_pop_JB22_mixed<-c()
#  for (i in seq(1,length(pop_JB22_mixed)-1)){
#    a1_pop_JB22_mixed<-c(a1_pop_JB22_mixed, 1/i)
#    a2_pop_JB22_mixed<-c(a2_pop_JB22_mixed, 1/(i^2))
#  }
#  a1_pop_JB22_mixed<-sum(a1_pop_JB22_mixed)
#  a2_pop_JB22_mixed<-sum(a2_pop_JB22_mixed)
#  b1_pop_JB22_mixed<-(length(pop_JB22_mixed)+1)/(3*(length(pop_JB22_mixed)-1))
#  b2_pop_JB22_mixed<-(2*((length(pop_JB22_mixed)^2)+length(pop_JB22_mixed)+3))/(9*length(pop_JB22_mixed)*(length(pop_JB22_mixed)-1))
#  c1_pop_JB22_mixed<-b1_pop_JB22_mixed-(1/a1_pop_JB22_mixed)
#  c2_pop_JB22_mixed<-b2_pop_JB22_mixed-((length(pop_JB22_mixed)+2)/(a1_pop_JB22_mixed*length(pop_JB22_mixed)))+(a2_pop_JB22_mixed/(a1_pop_JB22_mixed^2))
#  e1_pop_JB22_mixed<-c1_pop_JB22_mixed/a1_pop_JB22_mixed
#  e2_pop_JB22_mixed<-c2_pop_JB22_mixed/((a1_pop_JB22_mixed^2)+a2_pop_JB22_mixed)
#  theta_w_pop_JB22_mixed<-(sum(s_pop_JB22_mixed)/(end_pos-start_pos))/a1_pop_JB22_mixed
#  } else {
#    theta_w_pop_JB22_mixed<-"NA"
#  }
  #uncorrected_pw_theta
  #pop1
  variant_sites_pop1<-c()
  for (i in seq(1,length(gpop1[1,]))){
    if (!(sum(gpop1[,i]) %in% c(0,length(pop1)))){
      variant_sites_pop1<-c(variant_sites_pop1, i)
    }
  }
  dif_pairs_theta_pop1_uncorrected<-c()
  for (i in seq(1,length(pop1))){
    for (j in seq(1,length(pop1))){
      if (j>i){
        dif_pairs_theta_pop1_uncorrected<-c(dif_pairs_theta_pop1_uncorrected, 
          sum(!(gpop1[i,variant_sites_pop1]==gpop1[j,variant_sites_pop1])*1))
      }
    }
  }
  theta_pw_pop1_uncorrected<-(sum(dif_pairs_theta_pop1_uncorrected))/((((length(pop1)-1)*length(pop1))/2)*length(variant_sites_pop1))
  # pop 2
  variant_sites_pop2<-c()
  for (i in seq(1,length(gpop2[1,]))){
    if (!(sum(gpop2[,i]) %in% c(0,length(pop2)))){
      variant_sites_pop2<-c(variant_sites_pop2, i)
    }
  }
  dif_pairs_theta_pop2_uncorrected<-c()
  for (i in seq(1,length(pop2))){
    for (j in seq(1,length(pop2))){
      if (j>i){
        dif_pairs_theta_pop2_uncorrected<-c(dif_pairs_theta_pop2_uncorrected, 
          sum(!(gpop2[i,variant_sites_pop2]==gpop2[j,variant_sites_pop2])*1))
      }
    }
  }
  theta_pw_pop2_uncorrected<-(sum(dif_pairs_theta_pop2_uncorrected))/((((length(pop2)-1)*length(pop2))/2)*length(variant_sites_pop2))
  if (is.na(theta_pw_pop2) | is.na(theta_pw_pop1)) {
    theta_pw_pop1_uncorrected_pol<-NA
    theta_pw_pop2_uncorrected_pol<-NA
  } else if (theta_pw_pop2>theta_pw_pop1) {
    theta_pw_pop1_uncorrected_pol<-theta_pw_pop1_uncorrected
    theta_pw_pop2_uncorrected_pol<-theta_pw_pop2_uncorrected
  } else {
    theta_pw_pop1_uncorrected_pol<-theta_pw_pop2_uncorrected
    theta_pw_pop2_uncorrected_pol<-theta_pw_pop1_uncorrected
  }
#    ggplot(table_plot, aes(x=PC1, y=PC3, colour=sample, label=sample)) +
#    geom_point() + 
#    geom_text_repel(aes(label = sample), box.padding = unit(0.1, "lines")) +
#    theme(legend.position="none") +
#    labs(x = paste0("PC1 ",round(ccm_pca$varprop[1], digits=3)), 
#         y = paste0("PC3 ",round(ccm_pca$varprop[2], digits=3)), 
#         title = file_name ) +
#    ggsave(paste0("./PCA_byWindow_", number_snps, "_", overlapping_windows, "/", file_name, "_PC1_PC3", ".png"), width = 15, height = 10)
  varprop2<-data.frame(chromosome_name=chromosome, window_number=window, 
                       start_pos=start(new_vcf[1]), 
                       end_pos=end(new_vcf[length(new_vcf)]), 
                       pc1_varcom=ccm_pca$varprop[1], 
                       pc2_varcom=ccm_pca$varprop[2], 
                       pc3_varcom=ccm_pca$varprop[3], 
                       pc4_varcom=ccm_pca$varprop[4], 
                       pc5_varcom=ccm_pca$varprop[5], 
                       pop1N=length(pop1), 
                       pop2N=length(pop2), 
                       dxy_adj=dxy_adj, 
                       dxy_nom=dxy_nom, 
                       dxy_fq=dxy_fq, 
                       theta_pw_pop1=theta_pw_pop1, 
                       theta_pw_pop2=theta_pw_pop2, 
                       theta_pw_pop1_unco_pol=theta_pw_pop1_uncorrected_pol, 
                       theta_pw_pop2_unco_pol=theta_pw_pop1_uncorrected_pol, 
                       theta_w_pop1=theta_w_pop1, 
                       theta_w_pop2=theta_w_pop2, 
                       Taj_D_pop1=Taj_D_pop1, 
                       Taj_D_pop2=Taj_D_pop2, 
                       JB22_pop1=(paste0("JB22_",chromosome,"_ILL") %in% pop1))
#                       pop_JB22_nomixedN=length(pop_JB22_nomixed), 
#                       pop_JB22_mixedN=length(pop_JB22_mixed), 
#                       dxy_nom_hap_JB22=dxy_nom_hap_JB22,
#                       theta_pw_pop_JB22_nomixed=theta_pw_pop_JB22_nomixed,
#                       theta_pw_pop_JB22_mixed=theta_pw_pop_JB22_mixed,
#                       theta_w_pop_JB22_nomixed=theta_w_pop_JB22_nomixed,
#                       theta_w_pop_JB22_mixed=theta_w_pop_JB22_mixed, 
#                       JB869_pop_JB22_nomixed=("JB869" %in% pop_JB22_nomixed), 
#                       pop_JB22_nomixed_samples=paste(pop_JB22_nomixed, collapse = '_'))
  PC1_component<-ccm_pca$eigenvect[,1]
  PC2_component<-ccm_pca$eigenvect[,2]
  ###
  #polarization of haplotypes base on reference genome
  # if (((ccm_pca$eigenvect[ccm_pca$sample=="JB22",1]-min(PC1_component))/(max(PC1_component)-min(PC1_component)))>0.25){
  #   Nor_PC1_pol<-(1-((PC1_component-min(PC1_component))/(max(PC1_component)-min(PC1_component))))
  # }
  ###
  # polarization of haplotypes base on W theta differences
  if (is.na(theta_pw_pop2) | is.na(theta_pw_pop1)) {
    Nor_PC1_pol<-NA
  } else if (theta_pw_pop1<theta_pw_pop2){
    Nor_PC1_pol<-(PC1_component-min(PC1_component))/(max(PC1_component)-min(PC1_component))
  } else {
    Nor_PC1_pol<-(1-((PC1_component-min(PC1_component))/(max(PC1_component)-min(PC1_component))))
  }
  PC1_distribution2<-data.frame(chromosome_name=rep(chromosome, length(ccm_pca$eigenvect[,1])), 
    window_number=rep(window, length(ccm_pca$eigenvect[,1])), 
    start_pos=rep(start(new_vcf[1]), length(ccm_pca$eigenvect[,1])), 
    end_pos=rep(end(new_vcf[length(new_vcf)]), length(ccm_pca$eigenvect[,1])), 
    sample=ccm_pca$sample, 
    PC1=PC1_component, 
    PC2=ccm_pca$eigenvect[,2], 
    max_PC1=max(PC1_component), 
    min_PC1=min(PC1_component), 
    Nor_PC1=(PC1_component-min(PC1_component))/(max(PC1_component)-min(PC1_component)), 
    Nor_PC1_pol, 
    max_PC2=max(ccm_pca$eigenvect[,2]), 
    min_PC2=min(ccm_pca$eigenvect[,2]), 
    Nor_PC2=(PC2_component-min(PC2_component))/(max(PC2_component)-min(PC2_component)))
#  varprop<-rbind(varprop, varprop2)
#  PC1_distribution<-rbind(PC1_distribution, PC1_distribution2)
  write.table(varprop2, paste0("t_varprop_", chromosome, "_", number_snps, "_", overlapping_windows, "_", round_any(window, 1000), ".txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE, append=TRUE)
  write.table(PC1_distribution2, paste0("t_total_PC1_prop_", chromosome, "_", number_snps, "_", overlapping_windows, "_", round_any(window, 1000), ".txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE, append=TRUE)
  }
  snpgdsClose(genofile)
  start<-start+overlapping_windows
  end<-end+overlapping_windows
  #file.remove(paste0(file_name, ".vcf"))
}

#total_varprop<-varprop[-1,]
#total_PC1_distribution<-PC1_distribution[-1,]
#
#write.table(total_varprop, paste0("varprop_", chromosome, "_", number_snps, "_", overlapping_windows, ".txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)
#write.table(total_PC1_distribution, paste0("total_PC1_prop_", chromosome, "_", number_snps, "_", overlapping_windows, ".txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

#unlink(paste0("temporal_", temporal_number, ".gds")), force=TRUE)



