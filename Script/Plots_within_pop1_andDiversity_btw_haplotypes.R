# interactive -A snic2017-1-601 -n 1 -C fat -t 24:00:00
# this link is to use in Milou
module load R/3.4.3
export R_LIBS_USER=/home/sergio/R/libraries/3.4
# source("https://bioconductor.org/biocLite.R")
# biocLite("gdsfmt",  lib="/home/sergio/R/libraries/3.4")
# biocLite("SNPRelate",  lib="/home/sergio/R/libraries/3.4")
# install.packages('ggplot2', lib="/home/sergio/R/libraries/3.4")
# biocLite("VariantAnnotation",  lib="/home/sergio/R/libraries/3.4")
# install.packages('ggrepel', lib="/home/sergio/R/libraries/3.4")
# install.packages('tidyr', lib="/home/sergio/R/libraries/3.4")
# install.packages('plyr', lib="/home/sergio/R/libraries/3.4")
# install.packages('dplyr', lib="/home/sergio/R/libraries/3.4")
# install.packages('devtools', lib="/home/sergio/R/libraries/3.4")

R

library("SNPRelate")
library("gdsfmt")
library("ggplot2")
library("VariantAnnotation")
library("ggrepel")
library("tidyr")
library("plyr")
library("dplyr")

library("Hmisc")
library("reshape2")
library("ggstance")
library("ggtree") # to merge tree files witha notations. For instance to merge a tree with the heatmap
library("ggpubr") # for multiple plots
library("pvclust") # Ward Hierarchical Clustering with Bootstrapped p values
library("digest")
library("devtools")


# wtf_reg<-read.table("wtf_genes.txt", T) %>% mutate(xmin=start_pos, xmax=end_pos, ymin=-Inf, ymax=Inf) %>%  dplyr::select(chromosome_name, xmin, xmax, ymin, ymax) 
# wtf_actreg<-read.table("wtf_active_genes.txt", T) %>% mutate(xmin=start_pos, xmax=end_pos, ymin=-Inf, ymax=Inf) %>%  dplyr::select(chromosome_name, xmin, xmax, ymin, ymax) 
# tf_reg<-read.table("tf2_genes.txt", T) %>% mutate(xmin=start_pos, xmax=end_pos, ymin=-Inf, ymax=Inf) %>%  dplyr::select(chromosome_name, xmin, xmax, ymin, ymax) 
# ltr_reg<-read.table("ltr_genes.txt", T) %>% mutate(xmin=start_pos, xmax=end_pos, ymin=-Inf, ymax=Inf) %>%  dplyr::select(chromosome_name, xmin, xmax, ymin, ymax) 


total_varprop<-read.table("all_varprop_200_100.txt",F)

names(total_varprop)<-c("chromosome_name", 
                    "window_number", 
                    "start_pos", "end_pos", 
                    "pc1_varcom",
                    "pc2_varcom",
                    "pc3_varcom",
                    "pc4_varcom",
                    "pc5_varcom",
                    "pop1N",
                    "pop2N",
                    "dxy_adj",
                    "dxy_nom",
                    "dxy_fq",
                    "theta_pw_pop1",
                    "theta_pw_pop2",
                    "theta_pw_pop1_unco_pol",
                    "theta_pw_pop2_unco_pol",
                    "theta_w_pop1",
                    "theta_w_pop2",
                    "Taj_D_pop1",
                    "Taj_D_pop2",
                    "JB22_pop1")


head(total_varprop)


# final plots:
window_size<-200
overlap<-100

total_varprop %>% 
  filter(pc1_varcom>0.5) %>% 
  group_by(chromosome_name) %>% 
  summarise(num_windows=n())



total_varprop %>% 
  ggplot(aes(pc1_varcom)) +
  geom_histogram() +
  facet_grid(chromosome_name~.)



# PC var components:
total_varprop %>% 
  ggplot(aes(x=start_pos, y=pc1_varcom)) + 
  geom_line(alpha=0.8, colour="red") + 
  geom_line(aes(x=start_pos, y=pc2_varcom), alpha=0.8, colour="blue") + 
  #geom_point(aes(x=start_pos, y=theta_pw_pop2, colour="pop2"), alpha=0.8) + 
  #ylim(c(0,0.01)) + 
  xlab("Position") + ylab("PC Var comp.") + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + 
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") 
# ggsave(paste0("plot_theta_pw_noMix_red_haplotype_", window_size, "_", overlap, ".png"), width = 18, height = 5, dpi = 400)



# Pi per population:
total_varprop %>% 
  mutate(theta_pw_pop1_pol=ifelse(theta_pw_pop1<theta_pw_pop2, theta_pw_pop1, theta_pw_pop2), 
    theta_pw_pop2_pol=ifelse(theta_pw_pop1<theta_pw_pop2, theta_pw_pop2, theta_pw_pop1)) %>% 
  #filter(pc1_varcom>0.5 & !(log10(theta_pw_pop_JB22_mixed/dxy_nom_hap_JB22) %in% c("NaN", "-Inf")) & !(log10(theta_pw_pop_JB22_nomixed/dxy_nom_hap_JB22) %in% c("NaN", "-Inf"))) %>% 
  #ggplot(aes(x=start_pos, y=theta_pw_pop2_pol/theta_pw_pop1_pol)) + 
  ggplot(aes(x=start_pos, y=theta_pw_pop1_pol)) + 
  geom_line(alpha=0.8, colour="red") + 
  geom_line(aes(x=start_pos, y=theta_pw_pop2_pol), alpha=0.8, colour="blue") + 
  #geom_point(aes(x=start_pos, y=theta_pw_pop2, colour="pop2"), alpha=0.8) + 
  ylim(c(0,0.01)) + 
  xlab("Position") + ylab("Pi") + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + 
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") 
# ggsave(paste0("plot_theta_pw_noMix_red_haplotype_", window_size, "_", overlap, ".png"), width = 18, height = 5, dpi = 400)



# Dxy between populations:
total_varprop %>% 
  mutate(theta_pw_pop1_pol=ifelse(theta_pw_pop1<theta_pw_pop2, theta_pw_pop1, theta_pw_pop2), 
  theta_pw_pop2_pol=ifelse(theta_pw_pop1<theta_pw_pop2, theta_pw_pop2, theta_pw_pop1)) %>% 
  ggplot(aes(x=start_pos, y=dxy_nom)) + 
  geom_line(alpha=0.8, colour="black") + 
  geom_line(aes(x=start_pos, y=theta_pw_pop1_pol), alpha=0.5, colour="red") + 
  geom_line(aes(x=start_pos, y=theta_pw_pop2_pol), alpha=0.5, colour="blue") + 
  #geom_point(aes(x=start_pos, y=theta_pw_pop2, colour="pop2"), alpha=0.8) + 
  ylim(c(0,0.01)) + 
  xlab("Position") + ylab("Dxy Or Pi") + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + 
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") #+
  #ggsave(paste0("plot_theta_pw_Dxy_", window_size, "_", overlap, ".png"), width = 18, height = 5, dpi = 400)





# individual plot per sample and PC1:

data<-read.table(paste0("all_total_PC1_prop_200_100.txt"), F)

names(data)<-c("chromosome_name",
  "window_number", 
  "start_pos", 
  "end_pos", 
  "sample", 
  "PC1", 
  "PC2", 
  "max_PC1", 
  "min_PC1", 
  "Nor_PC1", 
  "Nor_PC1_pol", 
  "max_PC2", 
  "min_PC2", 
  "Nor_PC2")


data$start_pos<-as.numeric(data$start_pos)
data$end_pos<-as.numeric(data$end_pos)

pc1_varpro<- total_varprop
data<-merge(data, pc1_varpro, by=c("chromosome_name", "start_pos", "end_pos", "window_number"), all.x=T)

intervals<-data %>% 
  dplyr::select(chromosome_name, start_pos, end_pos, sample, Nor_PC1_pol) %>% 
  group_by(chromosome_name, start_pos, end_pos) %>% 
  spread(key=sample, value=Nor_PC1_pol) %>% 
  dplyr::select(chromosome_name, start_pos, end_pos) 

intervals$end_ed=intervals$end_pos

head(intervals)

intervals$end_ed=intervals$end_pos

for ( chromosome in levels(intervals$chromosome_name)){
  print(chromosome)
  table<- intervals %>% filter(chromosome_name==chromosome) 
  min_pos<-intervals %>% filter(chromosome_name==chromosome) %>% ungroup() %>% dplyr::select(start_pos) %>% unlist() %>% min()
  max_pos<-intervals %>% filter(chromosome_name==chromosome) %>% ungroup() %>% dplyr::select(end_pos) %>% unlist() %>% max()
  start_pos_tem<-min_pos
  end_pos_tem<-min_pos
  for (pos in sort(table$start_pos)){
    if (pos>start_pos_tem){
      intervals[intervals$chromosome_name==chromosome & intervals$start_pos==start_pos_tem,4]<-pos-1
      start_pos_tem<-pos
    }
  }
}

head(data)
data <- data %>% merge(intervals, by=c("chromosome_name", "start_pos", "end_pos"))

for ( chromosome in levels(intervals$chromosome_name)){
  print(chromosome)
  print(max(intervals$end_pos[intervals$chromosome_name==chromosome]))
}

data_small_windows<-matrix(rep(c(rep(-1, 5555844/1000), rep(-1, 4521265/1000), rep(-1, 2435469/1000)), length(levels(data$sample))), ncol=length(levels(data$sample)))

data_small_windows <- as.data.frame(data_small_windows)
names(data_small_windows) <- levels(data$sample)

data_small_windows <- data_small_windows %>% mutate(chromosome=c(rep("I", 5555844/1000), rep("II",  4521265/1000), rep("III", 2435469/1000)), position_bin=c(seq(0,(5555844/1000)-1), seq(0, ( 4521265/1000)-1), seq(0, (2435469/1000)-1)))

data_small_windows$position_bin<-as.numeric(data_small_windows$position_bin)

mark_values <- function(a,b,c,d,e){ data_small_windows[data_small_windows$chromosome==a & data_small_windows$position_bin %in% c(round((b/1000)):round((c/1000))), d] <<- e}

### this line will use the raw PC1 component to plot in the heat map:
#table_tobin<-data %>% dplyr::select('chromosome_name','start_pos','end_ed', 'sample', 'Nor_PC1_pol')
### and this one will simplify the PC1 component into 0, 0.5 or 1
data$Nor_PC1_pol_sim<-0
data$Nor_PC1_pol_sim[data$Nor_PC1_pol>0.4 & data$Nor_PC1_pol<0.6]<- 0.5
data$Nor_PC1_pol_sim[data$Nor_PC1_pol>0.6]<- 1
table_tobin<-data %>% dplyr::select('chromosome_name','start_pos','end_ed', 'sample', 'Nor_PC1_pol_sim')
#%>% filter(pc1_varcom>0.5)


apply(table_tobin, 1 , function(x) mark_values(x[1],as.numeric(x[2]),as.numeric(x[3]),x[4],as.numeric(x[5])))
data_small_windows_withoutpop1dif<-data_small_windows


data %>% 
  ungroup() %>% 
  select(chromosome_name, start_pos, end_ed, sample, Nor_PC1_pol_sim) %>% 
  write.table("ancestralhap_57ILL_LR.txt", quote = F, sep = "\t", 
            row.names = F, col.names = T)
  

## Difference in Theta W between ancestral population
pc1_varpro  %>% 
  filter(!(log2(theta_w_pop1/theta_w_pop2) %in% c("NaN", "-Inf", "Inf")) & pc1_varcom>0.5) %>% 
  ggplot(aes(start_pos, log2(theta_w_pop1/theta_w_pop2), colour=JB22_pop1)) + geom_point() + 
  xlab("Position") + 
  ylab("Log2 (W. Theta pop1 / W. Theta pop2 ) ") + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific) + 
  geom_hline(yintercept = 0, alpha=0.3) + 
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour="grey80"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"), axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black"), strip.background = element_rect(colour=NA, fill=NA)) + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") 
#plot_diff_theta_pop1_pop2 + ggsave("plot_diff_theta_pop1_pop2.png", width = 15, height = 4, dpi = 400)

log2_wt1_wt2<-pc1_varpro  %>% filter(!(log2(theta_w_pop1/theta_w_pop2) %in% c("NaN", "-Inf", "Inf")) & pc1_varcom>0.5) %>% mutate(log2_t1_t2=abs(log2(theta_w_pop1/theta_w_pop2))) %>% dplyr::select(chromosome_name, log2_t1_t2) 

ggplot(log2_wt1_wt2, aes(log2_t1_t2)) + geom_histogram() + 
  xlab("Absolute value Log2(W.theta pop 1 / pop 2)") + 
  ylab("Number of windows") + 
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour="grey80"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"), axis.text.x = element_text(hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") 
#hist_plot_log2_wt1_wt2_perChromosome + ggsave("hist_plot_log2_wt1_wt2_perChromosome.png", width = 7, height = 4, dpi = 400)

ggplot(log2_wt1_wt2, aes(log2_t1_t2)) + geom_histogram() + 
  xlab("Absolute value Log2(W.theta pop 1 / pop 2)") + 
  ylab("Number of windows") + 
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour="grey80"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"), axis.text.x = element_text(hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) 
#hist_plot_log2_wt1_wt2 + ggsave("hist_plot_log2_wt1_wt2.png", width = 4, height = 4, dpi = 400)



## Difference in Pi between ancestral population
pc1_varpro  %>% filter(!(log2(theta_pw_pop1/theta_pw_pop2) %in% c("NaN", "-Inf", "Inf")) & pc1_varcom>0.5) %>% ggplot(aes(start_pos, log2(theta_pw_pop1/theta_pw_pop2), colour=JB22_pop1)) + geom_point() + 
  xlab("Position") + 
  ylab("Log2 (Pi pop1 / Pi pop2 ) ") + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific) + 
  geom_hline(yintercept = 0, alpha=0.3) + 
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour="grey80"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"), axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black"), strip.background = element_rect(colour=NA, fill=NA)) + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") 
#plot_diff_theta_pw_pop1_pop2 + ggsave("plot_diff_theta_pw_pop1_pop2.png", width = 15, height = 4, dpi = 400)

log2_pwt1_pwt2<-pc1_varpro  %>% filter(!(log2(theta_pw_pop1/theta_pw_pop2) %in% c("NaN", "-Inf", "Inf")) & pc1_varcom>0.5) %>% mutate(log2_t1_t2=abs(log2(theta_pw_pop1/theta_pw_pop2))) %>% dplyr::select(chromosome_name, log2_t1_t2) 

ggplot(log2_pwt1_pwt2, aes(log2_t1_t2)) + geom_histogram() + 
  xlab("Absolute value Log2(Pi pop 1 / Pi pop 2)") + 
  ylab("Number of windows") + 
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour="grey80"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"), axis.text.x = element_text(hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") 
#hist_plot_log2_pwt1_pwt2_perChromosome + ggsave("hist_plot_log2_pwt1_pwt2_perChromosome.png", width = 7, height = 4, dpi = 400)

ggplot(log2_pwt1_pwt2, aes(log2_t1_t2)) + geom_histogram() + 
  xlab("Absolute value Log2(Pi pop 1 / Pi pop 2)") + 
  ylab("Number of windows") + 
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour="grey80"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"), axis.text.x = element_text(hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) 
#hist_plot_log2_pwt1_pwt2 + ggsave("hist_plot_log2_pwt1_pwt2.png", width = 4, height = 4, dpi = 400)



plot_list = list()
for (i in seq(1,length(levels(data$sample)))){
  print(i)
  sample_name<-levels(data$sample)[i]
  pc1_plot<-data %>% 
  filter(sample==sample_name) %>% 
  ggplot(aes(start_pos, Nor_PC1_pol, colour=chromosome_name)) + geom_line() +
    xlab("Position") + 
    ylab(paste0("PC1 - ",sample_name)) + 
    scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    geom_hline(yintercept = 0.5, alpha=0.3) + 
    theme_classic() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"), legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black"), strip.background = element_rect(colour=NA, fill=NA)) + 
    facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") 
  pc1_plot + ggsave(paste0("./plot_PC1_perSample/plot_PC1_perSample_", window_size, "_", overlap, "_",sample_name,".png"), width = 18, height = 2, dpi = 400)
  plot_list[[i]] = pc1_plot
}



plot_list = list()
for (i in seq(1,length(levels(data$sample)))){
  print(i)
  sample_name<-levels(data$sample)[i]
  pc1_plot<-data %>% 
  filter(sample==sample_name & pc1_varcom>0.5) %>% 
  ggplot(aes(start_pos, Nor_PC1_pol, colour=chromosome_name)) + geom_line() +
    xlab("Position") + 
    ylab(paste0(sample_name)) + 
    scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    geom_hline(yintercept = 0.5, alpha=0.3) + 
    theme_classic() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "grey"), legend.position="none", 
      axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), 
      axis.title=element_text(size=10,face="bold", colour="black"),
      axis.text.y = element_text(size=9, colour="black"), 
      strip.background = element_rect(colour=NA, fill=NA)) + 
    facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") 
  #pc1_plot + ggsave(paste0("./plot_PC1_perSample/plot_minPCval_PC1_perSample_", window_size, "_", overlap, "_",sample_name,".png"), width = 18, height = 2, dpi = 400)
  plot_list[[i]] = pc1_plot
}



for (seccion in seq(1,90,10)){
start_plot<-seccion
end_plot<-seccion+9
seq_plots<-seq(start_plot, end_plot)
for (i in seq(1,length(seq_plots))) {
  print(i)
  if (i==1){
    assign(paste("panel_", i , sep=""),plot_list[[seq_plots[i]]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x=element_blank()))
  } else if (i %in% seq(2,length(seq_plots)-1)) {
    assign(paste("panel_", i , sep=""),plot_list[[seq_plots[i]]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x=element_blank(), strip.text.x = element_blank()))
  } else if (i==length(seq_plots)) {
    assign(paste("panel_", i , sep=""),plot_list[[seq_plots[i]]] + theme(legend.position="none", strip.text.x = element_blank()))
  }
}
ggarrange(panel_1, panel_2, panel_3, panel_4, panel_5, panel_6, panel_7, panel_8, panel_9, panel_10, ncol = 1, nrow = length(seq_plots), 
  align="v", heights = c(rep(0.2,length(seq_plots)-1), 0.3)) + 
  ggsave(paste0("plot_PC1_perSample_", window_size, "_", overlap, "_panels", start_plot ,".png"), width = 12, height = 13, dpi = 400)
}


start_plot<-91
end_plot<-94
seq_plots<-seq(start_plot, end_plot)
for (i in seq(1,length(seq_plots))) {
  print(i)
  if (i==1){
    assign(paste("panel_", i , sep=""),plot_list[[seq_plots[i]]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x=element_blank()))
  } else if (i %in% seq(2,length(seq_plots)-1)) {
    assign(paste("panel_", i , sep=""),plot_list[[seq_plots[i]]] + theme(legend.position="none", axis.text.x = element_blank(), axis.title.x=element_blank(), strip.text.x = element_blank()))
  } else if (i==length(seq_plots)) {
    assign(paste("panel_", i , sep=""),plot_list[[seq_plots[i]]] + theme(legend.position="none", strip.text.x = element_blank()))
  }
}
ggarrange(panel_1, panel_2, panel_3, panel_4, ncol = 1, nrow = length(seq_plots), align="v", heights = c(rep(0.2,length(seq_plots)-1), 0.35)) + ggsave(paste0("plot_PC1_perSample_", window_size, "_", overlap, "_panels", start_plot ,".png"), width = 12, height = 5, dpi = 400)

#
#

head(data)
head(pc1_varpro)
# add bins column
head(data_small_windows)

# heat map plot_ final version:

head(data_small_windows)
# using only pop1 and pop2
dm <- melt(data_small_windows_withoutpop1dif,id.var=c("chromosome", "position_bin"))

# using raw PC1 values
# names(dm)<-c("chromosome_name", "start_pos", "sample", "Nor_PC1_pol")

# using simplified PC1 values (0,0.5 and 1)
names(dm)<-c("chromosome_name", "start_pos", "sample", "Nor_PC1_pol_sim")

# using raw PC1
# matrix_window <-data %>% dplyr::select(chromosome_name, window_number, start_pos, sample, Nor_PC1_pol) %>% mutate(windowID=paste0(chromosome_name, "_", window_number)) %>% dplyr::select(-chromosome_name, -window_number, -start_pos) %>% group_by(sample) %>% spread(windowID, Nor_PC1_pol) %>% ungroup() %>% dplyr::select(paste0(data$chromosome_name, "_", data$window_number)) %>% as.matrix()
# samples_ID<- data %>% dplyr::select(chromosome_name, window_number, start_pos, sample, Nor_PC1_pol) %>% mutate(windowID=paste0(chromosome_name, "_", window_number)) %>% dplyr::select(-chromosome_name, -window_number, -start_pos) %>% group_by(sample) %>% spread(windowID, Nor_PC1_pol) %>% ungroup() %>% dplyr::select(sample) %>% unlist() %>% as.vector()


# using simplified PC1
# matrix_window <-data %>% dplyr::select(chromosome_name, window_number, start_pos, sample, Nor_PC1_pol_sim) %>% mutate(windowID=paste0(chromosome_name, "_", window_number)) %>% dplyr::select(-chromosome_name, -window_number, -start_pos) %>% group_by(sample) %>% spread(windowID, Nor_PC1_pol_sim) %>% ungroup() %>% dplyr::select(paste0(data$chromosome_name, "_", data$window_number)) %>% as.matrix()
# samples_ID<- data %>% dplyr::select(chromosome_name, window_number, start_pos, sample, Nor_PC1_pol_sim) %>% mutate(windowID=paste0(chromosome_name, "_", window_number)) %>% dplyr::select(-chromosome_name, -window_number, -start_pos) %>% group_by(sample) %>% spread(windowID, Nor_PC1_pol_sim) %>% ungroup() %>% dplyr::select(sample) %>% unlist() %>% as.vector()

# another version just using "dm" intead of "data"
matrix_window <-dm %>% mutate(windowID=paste0(chromosome_name, "_", start_pos), Nor_PC1_pol_sim=as.factor(Nor_PC1_pol_sim)) %>% dplyr::select(-chromosome_name, -start_pos) %>% group_by(sample) %>% spread(windowID, Nor_PC1_pol_sim) %>% ungroup() %>% dplyr::select(paste0(dm$chromosome_name, "_", dm$start_pos)) %>% as.matrix()
samples_ID<- dm %>% mutate(windowID=paste0(chromosome_name, "_", start_pos)) %>% dplyr::select(-chromosome_name, -start_pos) %>% group_by(sample) %>% spread(windowID, Nor_PC1_pol_sim) %>% ungroup() %>% dplyr::select(sample) %>% unlist() %>% as.vector()

row.names(matrix_window)<- samples_ID



d <- dist(matrix_window, method = "euclidean")
H.fit <- hclust(d, method="ward")
H.fit$labels<-samples_ID
plot(H.fit)

dm$sample <- factor(dm$sample, levels = H.fit$labels[H.fit$order])

write(H.fit$labels[H.fit$order], file = "order_sample.txt")
order_samples<-scan("order_sample.txt", what="", sep="\n") 

dm$sample <- factor(dm$sample, levels = order_samples)


# plots _ variance components for PC1 and PC2
pc1Var_plot<-ggplot(pc1_varpro, aes(start_pos, pc1_varcom)) + geom_line() + ylim(c(0,1)) + labs(x=NULL) + ylab("PC1") + scale_x_continuous(breaks=seq(500000,5500000,500000), expand = c(0, 0)) + geom_hline(yintercept = c(0.25, 0.5, 0.75), colour="grey80", alpha=0.5) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")

pc2Var_plot<-ggplot(pc1_varpro, aes(start_pos, pc2_varcom)) + geom_line() + ylim(c(0,1)) + xlab("Position (bp)") + ylab("PC2") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand = c(0, 0)) + geom_hline(yintercept = c(0.25, 0.5, 0.75), colour="grey80", alpha=0.5) + theme_classic() + theme(axis.line = element_line(colour = "black"), strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")
ggarrange(pc1Var_plot, pc2Var_plot, ncol = 1, nrow = 2, align="v", heights = c(0.2, 0.25)) + ggsave(paste0("plot_merge_VarCom_", "_", window_size, "_", overlap, ".png"), width = 10, height = 6, dpi = 400)



# plot comparison Theta W withing pop1 and pop2
ggplot(pc1_varpro[is.finite(log10(pc1_varpro$theta_w_pop1/pc1_varpro$theta_w_pop2)),], aes(start_pos, log10(theta_pw_pop1/theta_pw_pop2), colour=JB22_pop1)) + geom_point(alpha=0.7) + geom_hline(yintercept = 0) +  labs(x=NULL) + ylab("Diff. log10(Theta_pw pop1 / pop2)") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand = c(0.05, 1)) +  theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), strip.background = element_blank(), strip.text.x = element_blank(), axis.text.y = element_text(colour="black"), legend.position="bottom") + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") #+ ggsave(paste0("plot_com_pop1_pop2_theta_w_", window_size, "_", overlap, ".png"), width = 8, height = 6, dpi = 600)


head(pc1_varpro)
# add columns of WTheta, PWTheta and T'sD polarizing pop1 and pop2 base on differences in Pi.
pc1_varpro<- pc1_varpro %>% mutate(theta_w_pop1_pol= (((theta_pw_pop1<theta_pw_pop2)*1)*theta_w_pop1)+((!((theta_pw_pop1<theta_pw_pop2))*1)*theta_w_pop2), 
                                   theta_w_pop2_pol= (((theta_pw_pop1<theta_pw_pop2)*1)*theta_w_pop2)+((!((theta_pw_pop1<theta_pw_pop2))*1)*theta_w_pop1), 
                                   theta_pw_pop1_pol= (((theta_pw_pop1<theta_pw_pop2)*1)*theta_pw_pop1)+((!((theta_pw_pop1<theta_pw_pop2))*1)*theta_pw_pop2), 
                                   theta_pw_pop2_pol= (((theta_pw_pop1<theta_pw_pop2)*1)*theta_pw_pop2)+((!((theta_pw_pop1<theta_pw_pop2))*1)*theta_pw_pop1), 
                                   Taj_D_pop1_pol= (((theta_pw_pop1<theta_pw_pop2)*1)*Taj_D_pop1)+((!((theta_pw_pop1<theta_pw_pop2))*1)*Taj_D_pop2), 
                                   Taj_D_pop2_pol= (((theta_pw_pop1<theta_pw_pop2)*1)*Taj_D_pop2)+((!((theta_pw_pop1<theta_pw_pop2))*1)*Taj_D_pop1), 
                                   N_pop1_pol= (((theta_pw_pop1<theta_pw_pop2)*1)*pop1N)+((!((theta_pw_pop1<theta_pw_pop2))*1)*pop2N), 
                                   N_pop2_pol= (((theta_pw_pop1<theta_pw_pop2)*1)*pop2N)+((!((theta_pw_pop1<theta_pw_pop2))*1)*pop1N)) 


head(dm)
heatplot<-ggplot((dm %>% filter(Nor_PC1_pol>(-1))), aes(start_pos*1000, sample)) + geom_tile(aes(fill = Nor_PC1_pol)) + scale_fill_gradient(low = "steelblue", high = "darkred" , expand=c(0,0)) + labs(x=NULL) + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + theme_classic() + ylab("\nSample") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none", axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black"), strip.text.x = element_blank()) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") 

ggplot((dm %>% filter(Nor_PC1_pol>(-1))), aes(start_pos*1000, sample)) + geom_tile(aes(fill = Nor_PC1_pol)) + 
  scale_fill_gradient(low = "steelblue", high = "darkred" , expand=c(0,0)) + labs(x=NULL) + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + theme_classic() + ylab("\nSample") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none", axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black"), strip.text.x = element_blank()) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + 
  ggsave(paste0("plot_heatmap_", window_size, "_", overlap, ".png"), width = 10, height = 20, dpi = 400)



plot<-ggplot((dm %>% filter(Nor_PC1_pol_sim>(-1))), aes(start_pos*1000, sample)) + geom_tile(aes(fill = Nor_PC1_pol_sim)) + 
  scale_fill_gradient(low = "darkred" , high = "steelblue" , expand=c(0,0)) + 
  labs(x=NULL) + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + 
  theme_classic() + ylab("\nSample") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none", 
    axis.text.x=element_text(angle = 45, hjust = 1, size=12, colour="black"), 
    #axis.text.y = element_blank(), 
    axis.text.y = element_text(colour="black"), 
    strip.text.x = element_blank()) + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")  

plot +
  ggsave(paste0("plot_heatmap_sim_", window_size, "_", overlap, ".png"), 
    width = 8, height = 11, dpi = 450)

# svg, but the file is massive: ~ 300 Mb
plot +
  ggsave(paste0("plot_heatmap_sim_", window_size, "_", overlap, ".svg"), 
    width = 10, height = 10, dpi = 450)


# heat map using breaking pop1 (red population) into JB22 group and JB4 group:
dm$an_population<-"JB22_pop1"
dm$an_population[dm$Nor_PC1_pol_sim=="0.1"]<-"JB4_pop1"
dm$an_population[dm$Nor_PC1_pol_sim=="0.5"]<-"Het_win"
dm$an_population[dm$Nor_PC1_pol_sim=="1"]<-"pop2"
dm$an_population<-as.factor(dm$an_population)

ggplot(dm, aes(start_pos*1000, sample)) + 
  geom_tile(aes(fill = an_population)) + 
  scale_fill_manual(values=c("darkorchid3", "darkred", "darkorange1", "steelblue")) + 
  labs(x=NULL) + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + theme_classic() + ylab("\nSample") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none", axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black"), strip.text.x = element_blank()) + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")  + 
  ggsave(paste0("plot_heatmap_", window_size, "_", overlap, "_withJB22andJB4.png"), width = 10, height = 20, dpi = 400)

# proportion of ancestral populations per sample
No_bins_perChromosome<-dm %>% 
  dplyr::select(-Nor_PC1_pol_sim) %>% 
  group_by(sample, chromosome_name) %>% 
  summarise(number_bins = n()) 


sample_order_AP<-dm %>% 
  dplyr::select(-Nor_PC1_pol_sim) %>% 
  group_by(sample, an_population) %>% 
  summarise(number_bins = n()) %>% 
  mutate(number_bins_fq_total=number_bins/12511) %>% 
  filter(an_population=="JB22_pop1") %>% 
  arrange(number_bins_fq_total) %>% 
  ungroup() %>% 
  select(sample) %>% 
  unlist() %>% 
  as.vector()


dm %>% 
  dplyr::select(-Nor_PC1_pol_sim) %>% 
  group_by(sample, an_population) %>% 
  summarise(number_bins = n()) %>% 
  mutate(number_bins_fq_total=number_bins/12511) %>% 
  ungroup() %>% 
  mutate(sample=factor(sample, levels=sample_order_AP)) %>% 
  ggplot(aes(x=sample, y=number_bins_fq_total, fill=an_population)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c("darkorange1", "darkred", "steelblue"), expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0,1,0.5)) + 
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(size=14, colour="black")) + 
  #facet_grid(chromosome_name ~ .) + 
  ggsave(paste0("proportion_an_populations_", window_size, "_", overlap, ".png"), width = 24, height = 5, dpi = 400)






#plot by chromosomes:
dm %>% 
  dplyr::select(-Nor_PC1_pol_sim) %>% 
  group_by(sample, chromosome_name, an_population) %>% 
  summarise(number_bins = n()) %>% 
  mutate(number_bins_fq_perChromosome = ((chromosome_name=="I")*(number_bins/5563))+((chromosome_name=="II")*(number_bins/4537))+((chromosome_name=="III")*(number_bins/2440))) %>% 
  ggplot(aes(x=sample, y=number_bins_fq_perChromosome, fill=an_population)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c("darkorchid3", "darkred", "darkorange1", "steelblue"), expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0,1,0.5)) + 
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(size=14, colour="black")) + 
  facet_grid(chromosome_name ~ .) + 
  ggsave(paste0("proportion_an_populations_bychromosome", window_size, "_", overlap, ".png"), width = 24, height = 15, dpi = 400)


##

dm %>% 
  dplyr::select(-Nor_PC1_pol_sim) %>% 
  group_by(sample, an_population) %>% 
  summarise(number_bins = n()) %>% 
  mutate(number_bins_fq_total=number_bins/12540) %>% 
  write.table(paste0("proportions_AncPop_161Samples_", window_size, "_", overlap, ".txt"), 
  	quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)


# proportion of ancestral populations per sample - using only two ancestral populations:
No_bins_perChromosome<-dm %>% 
  dplyr::select(-Nor_PC1_pol_sim) %>% 
  group_by(sample, chromosome_name) %>% 
  summarise(number_bins = n()) 

dm_forFig2<-dm %>% 
  dplyr::select(-Nor_PC1_pol_sim) %>% 
  group_by(sample, chromosome_name, an_population) %>% 
  summarise(number_bins = n()) %>% 
  mutate(number_bins_fq_total=number_bins/12540)

order_samples_proportions<-dm_forFig2 %>% filter(an_population=="JB22_pop1") %>% 
  group_by(sample) %>% 
  summarise(total_pop1 = sum(number_bins_fq_total)) 

order_samples_proportions<-order_samples_proportions[order(order_samples_proportions$total_pop1),"sample"] %>% unlist() %>% as.vector()

dm_forFig2$sample<-factor(dm_forFig2$sample, levels=order_samples_proportions)


dm_forFig2 %>% filter(sample!="JB1207") %>% 
  ggplot(aes(x=sample, y=number_bins_fq_total, fill=an_population)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c("darkorchid3", "darkred", "steelblue"), expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0,1,0.5)) + 
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), 
    axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=15, colour="black"), 
    axis.text.y = element_text(size=16, colour="black"), 
    legend.position="none") + 
  #facet_grid(chromosome_name ~ .) + 
  ggsave(paste0("proportion_an_populations_2AncPop", window_size, "_", overlap, ".png"), width = 14, height = 3, dpi = 400)

#plot by chromosomes:
dm %>% 
  dplyr::select(-Nor_PC1_pol_sim) %>% 
  group_by(sample, chromosome_name, an_population) %>% 
  summarise(number_bins = n()) %>% 
  mutate(number_bins_fq_perChromosome = ((chromosome_name=="I")*(number_bins/5563))+((chromosome_name=="II")*(number_bins/4537))+((chromosome_name=="III")*(number_bins/2440))) %>% 
  ggplot(aes(x=sample, y=number_bins_fq_perChromosome, fill=an_population)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c("darkorchid3", "darkred", "steelblue"), expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0,1,0.5)) + 
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(size=14, colour="black")) + 
  facet_grid(chromosome_name ~ .) + 
  ggsave(paste0("proportion_an_populations_bychromosome_2AncPop_", window_size, "_", overlap, ".png"), width = 24, height = 15, dpi = 400)

dm %>% 
  dplyr::select(-Nor_PC1_pol_sim) %>% 
  group_by(sample, chromosome_name, an_population) %>% 
  summarise(number_bins = n()) %>% 
  mutate(number_bins_fq_perChromosome = ((chromosome_name=="I")*(number_bins/5563))+((chromosome_name=="II")*(number_bins/4537))+((chromosome_name=="III")*(number_bins/2440))) %>% ungroup() %>% 
  dplyr::select(-number_bins) %>% 
  spread(key=chromosome_name, value=number_bins_fq_perChromosome) %>% 
  filter(an_population=="pop2") %>% 
  ggplot(aes(x=I, y=II, colour="Chr.II")) +
  geom_point(size=3, alpha=0.7) +
  geom_point(aes(x=I, y=III, colour="Chr.III"), size=3, alpha=0.5) +
  geom_smooth() + 
  #geom_smooth(method=lm) + 
  geom_smooth(aes(x=I, y=III, colour="Chr.III")) + 
  #geom_smooth(aes(x=I, y=III, colour="Chr.III"), method=lm) + 
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), 
    axis.text.x=element_text(size=15, colour="black"), 
    axis.text.y = element_text(size=15, colour="black"), 
    legend.position="bottom") + 
  ggsave(paste0("proportion_an_populations_pop2_bychromosome_smooth", window_size, "_", overlap, ".png"), width = 10, height = 10, dpi = 400)

dm %>% 
  dplyr::select(-Nor_PC1_pol_sim) %>% 
  group_by(sample, chromosome_name, an_population) %>% 
  summarise(number_bins = n()) %>% 
  mutate(number_bins_fq_perChromosome = ((chromosome_name=="I")*(number_bins/5563))+((chromosome_name=="II")*(number_bins/4537))+((chromosome_name=="III")*(number_bins/2440))) %>% ungroup() %>% 
  dplyr::select(-number_bins) %>% 
  spread(key=an_population, value=number_bins_fq_perChromosome) %>%   
  ggplot(aes(x=JB22_pop1, y=JB4_pop1, colour=chromosome_name)) +
  geom_point(size=3, alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, alpha=0.4) +
  geom_abline(intercept = 1, slope = -1, alpha=0.4) + 
  #scale_x_continuous(breaks=seq(0,1,0.2)) + 
  #scale_y_continuous(breaks=seq(0,1,0.2)) + 
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(size=9, colour="black"), axis.text.y = element_text(size=9, colour="black"), legend.position="bottom") + 
  ggsave(paste0("proportion_an_populations_pop1_sub_bychromosome", window_size, "_", overlap, ".png"), width = 10, height = 10, dpi = 400)


hybrid_samples<-dm %>% 
  dplyr::select(-Nor_PC1_pol_sim) %>% 
  group_by(sample, chromosome_name, an_population) %>% 
  summarise(number_bins = n()) %>% 
  mutate(number_bins_fq_total=number_bins/12540) %>% 
  filter(an_population=="pop2" & number_bins_fq_total>0.1) %>% 
  ungroup() %>% dplyr::select(sample) %>% unlist() %>% as.vector()


data %>% filter(pc1_varcom>0.5 & sample %in% hybrid_samples & !(sample %in% c("JB1169", "JB1207")) & Nor_PC1_pol_sim!=0.5) %>% 
  dplyr::select(sample, chromosome_name, start_pos, Nor_PC1_pol_sim) %>% 
  group_by(chromosome_name, start_pos) %>% 
  summarise(proportion_pop2 = mean(Nor_PC1_pol_sim)) %>% 
  ggplot(aes(start_pos, proportion_pop2)) + 
  #geom_rect(data=wtf_reg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="grey50", alpha=0.2, inherit.aes = FALSE) + 
  #geom_rect(data=wtf_actreg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="red", alpha=0.5, inherit.aes = FALSE) + 
  geom_hline(yintercept = 1, colour="gray80") + 
  geom_point(alpha=0.5) + 
  geom_line(alpha=0.5, colour="blue") + 
  ylab("Prop. pop2") + 
  xlab("Position") + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) +
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) +
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") +
  ggsave(paste0("proportion_an_populations_pop2_aloneGenome_withoutWTF", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)

proportion_an_populations_aloneGenome <- data %>% filter(sample %in% hybrid_samples & !(sample %in% c("JB1169", "JB1207")) & Nor_PC1_pol_sim!=0.5) %>% 
  dplyr::select(sample, chromosome_name, start_pos, Nor_PC1_pol_sim) %>% 
  group_by(chromosome_name, start_pos) %>% 
  summarise(proportion_pop2 = mean(Nor_PC1_pol_sim)) 



ggplot(proportion_an_populations_aloneGenome, aes(start_pos, proportion_pop2)) + 
  geom_rect(data=wtf_reg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="grey50", alpha=0.2, inherit.aes = FALSE) + 
  geom_rect(data=wtf_actreg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="red", alpha=0.5, inherit.aes = FALSE) + 
  geom_hline(yintercept = 1, colour="gray80") + 
  geom_point(alpha=0.5) + 
  geom_line(alpha=0.5, colour="blue") + 
  ylab("Prop. pop2") + 
  xlab("Position") + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) +
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) +
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") +
  ggsave(paste0("proportion_an_populations_pop2_aloneGenome_withWTF", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)


# now using simplified PC1
heatplot<-ggplot((dm %>% filter(Nor_PC1_pol_sim>(-1))), aes(start_pos*1000, sample)) + geom_tile(aes(fill = Nor_PC1_pol_sim)) + scale_fill_gradient(low = "darkred", high = "steelblue", expand=c(0,0)) + labs(x=NULL) + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + theme_classic() + ylab("\nSample") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none", axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black"), strip.text.x = element_blank()) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") 

ggplot((dm %>% filter(Nor_PC1_pol_sim>(-1))), aes(start_pos*1000, sample)) + geom_tile(aes(fill = Nor_PC1_pol_sim)) + scale_fill_gradient(low = "darkred", high = "steelblue", expand=c(0,0)) + labs(x=NULL) + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + theme_classic() + ylab("\nSample") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none", axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black"), strip.text.x = element_blank()) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_heatmap_sim_", window_size, "_", overlap, ".png"), width = 8, height = 10, dpi = 400)

dxy_plot<-ggplot(pc1_varpro[is.finite(pc1_varpro$dxy_nom),], aes(start_pos, dxy_nom)) + geom_line()  + labs(x=NULL) + ylab("\nDxy") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + ylim(c(0,0.015)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")

ggplot(pc1_varpro[is.finite(pc1_varpro$dxy_nom),], aes(start_pos, dxy_nom)) + geom_line()  + labs(x=NULL) + ylab("\nDxy") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_dxy_pop1Andpop2_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)


plot_theta_pw_pop1_pol<- ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$theta_pw_pop1_pol),], aes(start_pos, theta_pw_pop1_pol)) + geom_line(colour="darkred") + labs(x=NULL) + ylab("Pi\npop1") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) +   ylim(c(0,0.004)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")

ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$theta_pw_pop1_pol),], aes(start_pos, theta_pw_pop1_pol)) + geom_line(colour="darkred") + labs(x=NULL) + ylab("Pi\npop1") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_theta_pw_pop1_pol_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)

plot_theta_pw_pop2_pol<- ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$theta_pw_pop2_pol),], aes(start_pos, theta_pw_pop2_pol)) + geom_line(colour="steelblue") + labs(x=NULL) + ylab("Pi\npop2") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + ylim(c(0,0.010)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x = element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")


ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$theta_pw_pop2_pol),], aes(start_pos, theta_pw_pop2_pol)) + geom_line(colour="steelblue") + labs(x=NULL) + ylab("Pi\npop2") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_theta_pw_pop2_pol_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)


plot_theta_pw_overDxy <-pc1_varpro %>% filter(pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$dxy_nom) & is.finite(pc1_varpro$theta_pw_pop1_pol) & is.finite(pc1_varpro$theta_pw_pop2_pol)) %>% 
  ggplot(aes(start_pos, theta_pw_pop1_pol/dxy_nom, colour="pop1")) + geom_line()  + 
  geom_line(aes(start_pos, theta_pw_pop2_pol/dxy_nom, colour="pop2")) + 
                        labs(x=NULL) + 
                        ylab("\nPi/Dxy") + 
                        scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + 
  scale_color_manual(values = c("darkred", "steelblue")) + 
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), axis.text.x=element_blank(), axis.text.y = element_text(colour="black"), strip.text.x = element_blank(), legend.position="none") + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")


pc1_varpro %>% filter(pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$dxy_nom) & is.finite(pc1_varpro$theta_pw_pop1_pol) & is.finite(pc1_varpro$theta_pw_pop2_pol)) %>% 
  ggplot(aes(start_pos, theta_pw_pop1_pol/dxy_nom, colour="pop1")) + geom_line()  + 
  geom_line(aes(start_pos, theta_pw_pop2_pol/dxy_nom, colour="pop2")) + 
  labs(x=NULL) + 
  ylab("\nPi/Dxy") + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + 
  scale_color_manual(values = c("darkred", "steelblue")) + 
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black"), strip.text.x = element_blank(), legend.position="none") + 
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_theta_pw_overDxy_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)

ggarrange(dxy_plot, plot_theta_pw_pop1_pol, plot_theta_pw_pop2_pol, plot_theta_pw_overDxy, heatplot, ncol = 1, nrow = 5, align="v", heights = c(2.5, 1.5, 1.5, 1.5, 18)) + ggsave(paste0("plot_merge_heatmap_VarCom_", "_", window_size, "_", overlap, ".png"), width = 10, height = 12, dpi = 400)



## Pair-wise Theta 

head(pc1_varpro)

plot_theta_pw_pop1_pol<- ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$theta_pw_pop1_pol),], aes(start_pos, theta_pw_pop1_pol)) + geom_line(colour="darkred") + labs(x=NULL) + ylab("Theta_PW\npop1") + ylim(c(0,0.004)) + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) +  theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")

ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$theta_pw_pop1_pol),], aes(start_pos, theta_pw_pop1_pol)) + geom_line(colour="darkred") + labs(x=NULL) + ylab("Theta_PW\npop1") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + ylim(c(0,0.020)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_theta_pw_pop1_pol_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)

plot_theta_pw_pop2_pol<- ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$theta_pw_pop2_pol),], aes(start_pos, theta_pw_pop2_pol)) + geom_line(colour="steelblue") + labs(x=NULL) + ylab("Theta_PW\npop2") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + ylim(c(0,0.010)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x = element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")


ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$theta_pw_pop2_pol),], aes(start_pos, theta_pw_pop2_pol)) + geom_line(colour="steelblue") + labs(x=NULL) + ylab("Theta_PW\npop2") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + ylim(c(0,0.020)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_theta_pw_pop2_pol_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)

# pi, W. Theta and Taj's D merging all samples and only in hybrid samples (this exclude "pure" red samples, for instance JB22, and samples with high blue content)


total_varprop_allSamples_I<-read.table("./varprop_allSamples_I_200_100.txt",T)
total_varprop_allSamples_II<-read.table("./varprop_allSamples_II_200_100.txt",T)
total_varprop_allSamples_III<-read.table("./varprop_allSamples_III_200_100.txt",T)
total_varprop_allSamples <-rbind(total_varprop_allSamples_I, total_varprop_allSamples_II, total_varprop_allSamples_III)

total_varprop_allSamples_pc1_varpro <- pc1_varpro %>% dplyr::select(chromosome_name, window_number, start_pos, end_pos, pc1_varcom) %>% merge(proportion_an_populations_aloneGenome, by=c("chromosome_name", "start_pos")) %>%  merge(total_varprop_allSamples, by=c("chromosome_name", "window_number"))

total_varprop_allSamples_I<-read.table("./second_structure/varprop_allSamples_I_200_100.txt",T)
total_varprop_allSamples_II<-read.table("./second_structure/varprop_allSamples_II_200_100.txt",T)
total_varprop_allSamples_III<-read.table("./second_structure/varprop_allSamples_III_200_100.txt",T)
total_varprop_allSamples <-rbind(total_varprop_allSamples_I, total_varprop_allSamples_II, total_varprop_allSamples_III)

total_varprop_allSamples_pc1_varpro <- pc1_varpro %>% dplyr::select(chromosome_name, window_number, start_pos, end_pos, pc1_varcom) %>% merge(proportion_an_populations_aloneGenome, by=c("chromosome_name", "start_pos")) %>%  merge(total_varprop_allSamples, by=c("chromosome_name", "window_number"))

ggplot(total_varprop_allSamples_pc1_varpro %>% filter(pc1_varcom>0.6), aes(start_pos, proportion_pop2)) + 
  #geom_rect(data=wtf_reg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="grey50", alpha=0.2, inherit.aes = FALSE) + 
  #geom_rect(data=wtf_actreg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="red", alpha=0.5, inherit.aes = FALSE) + 
  #geom_rect(data=large_SV, aes(xmin=xmin, xmax=xmax, ymin=3, ymax=3.1), color="orange", alpha=0.5, inherit.aes = FALSE) + 
  geom_hline(yintercept = 1, colour="gray80") + 
  geom_point(size=2, alpha=0.5) + 
  geom_line(alpha=0.4, colour="black") + 
  ylab("Proportion Pop2") + 
  xlab("Position") + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) +
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour="gray80"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) +
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("proportion_an_populations_pop2_aloneGenome_", window_size, "_", overlap, ".png"), width = 12, height = 4, dpi = 400)


ggplot((total_varprop_allSamples_pc1_varpro %>% filter(pc1_varcom>0.6 & start_pos<5400000)), aes(proportion_pop2, Taj_D_hybrids)) + 
  geom_point(alpha=0.5) + 
  xlim(c(0,1)) +
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour="gray80"), panel.grid.minor = element_line(colour="gray80"), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) +
  facet_grid(. ~ chromosome_name) +
  ggsave(paste0("proportion_an_populations_pop2_Vs_TajD_", window_size, "_", overlap, ".png"), width = 12, height = 4, dpi = 400)




ggplot(total_varprop_allSamples_pc1_varpro %>% filter(pc1_varcom>0.6), aes(start_pos, Taj_D_hybrids, colour=proportion_pop2)) + 
  #geom_rect(data=wtf_reg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="orange", alpha=0.2, inherit.aes = FALSE) + 
  #geom_rect(data=wtf_actreg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="red", alpha=0.5, inherit.aes = FALSE) + 
  #geom_rect(data=large_SV, aes(xmin=xmin, xmax=xmax, ymin=3, ymax=3.1), color="orange", alpha=0.5, inherit.aes = FALSE) + 
  #geom_rect(data=tf_reg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="orange", alpha=0.5, inherit.aes = FALSE) + 
  geom_rect(data=ltr_reg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="orange", alpha=0.5, inherit.aes = FALSE) + 
  geom_hline(yintercept = 0, colour="gray80") + 
  geom_point(size=2, alpha=0.5) + 
  scale_colour_gradient2(midpoint=0.5, low="darkred", mid="white", high="steelblue", space ="Lab" , expand=c(0,0)) + 
  geom_line(alpha=0.7, colour="black") + 
  #geom_line(aes(start_pos, proportion_pop2), alpha=0.5, colour="red") + 
  ylab("Taj_D_hybrids") + 
  xlab("Position") + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific) +
  scale_y_continuous(breaks=seq(-5,5,1)) +
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour="gray85"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black"), legend.position="none") +
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") +
  ggsave(paste0("TajD_aloneGenome_withproportion_an_populations_pop2_withLTR", window_size, "_", overlap, ".png"), width = 12, height = 4, dpi = 400)




ggplot(total_varprop_allSamples_pc1_varpro %>% filter(pc1_varcom>0.6), aes(start_pos, theta_pw_hybrids, colour=Taj_D_hybrids)) + 
  #geom_rect(data=wtf_reg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="grey50", alpha=0.2, inherit.aes = FALSE) + 
  #geom_rect(data=wtf_actreg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="red", alpha=0.5, inherit.aes = FALSE) + 
  #geom_rect(data=large_SV, aes(xmin=xmin, xmax=xmax, ymin=3, ymax=3.1), color="orange", alpha=0.5, inherit.aes = FALSE) + 
  #geom_hline(yintercept = 1, colour="gray80") + 
  geom_point(size=2) + 
  geom_line(alpha=0.4, colour="black") + 
  scale_colour_gradient2(midpoint=0, low="darkred", mid="white", high="steelblue", space ="Lab" , expand=c(0,0)) + 
  ylim(c(0,0.006)) +
  ylab("Pi") + 
  xlab("Position") + 
  scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) +
  theme_classic() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour="gray80"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) +
  facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + 
  ggsave(paste0("Pi_aloneGenome_HybridSamples_", window_size, "_", overlap, ".png"), width = 12, height = 4, dpi = 400)




# 17 - 04 - 2018
# estimating t (number of genetations) from coalescesce between ancestral haplotypes:
# from Dxy and pair-wise theta

pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$theta_pw_pop1_pol),] %>% head()

theta_pw_pop1_pol
theta_pw_pop2_pol
dxy_nom
# mutation rate: 2x10-10 mutations site-1 generation-1 
# Ne: 12e+06
u<-2e-10
Ne<-12e+06

pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$theta_pw_pop1_pol),] %>% mutate(time_dxy=(dxy_nom-((theta_pw_pop1_pol+theta_pw_pop2_pol)/2))/(2*u)) %>% ggplot(aes(start_pos, time_dxy)) + geom_line(colour="black") + labs(x=NULL) + ylab("time (asexual generations)") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + ylim(c(0,30000000)) +theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_time_dxy_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)

pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$theta_pw_pop1_pol),] %>% mutate(time_dxy=(dxy_nom-((theta_pw_pop1_pol+theta_pw_pop2_pol)/2))/(2*u)) %>% select(time_dxy) %>% filter(time_dxy<30000000) %>% unlist() %>% as.vector() %>% mean()
#
#
#
## Tajima's D

head(pc1_varpro)

plot_Taj_D_pop1_pol<- ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$Taj_D_pop1_pol),], aes(start_pos, Taj_D_pop1_pol)) + geom_line(colour="darkred") + labs(x=NULL) + ylab("Taj_D\npop1") + ylim(c(-2.5,2.5)) + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) +  theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")

ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$Taj_D_pop1_pol),], aes(start_pos, Taj_D_pop1_pol)) + geom_line(colour="darkred") + labs(x=NULL) + ylab("Taj_D\npop1") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + ylim(c(-2.5,2.5)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_Taj_D_pop1_pol_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)

plot_Taj_D_pop2_pol<- ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$Taj_D_pop2_pol),], aes(start_pos, Taj_D_pop2_pol)) + geom_line(colour="steelblue") + labs(x=NULL) + ylab("Taj_D\npop2") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + ylim(c(-2.5,2.5)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x = element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")


ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$Taj_D_pop2_pol),], aes(start_pos, Taj_D_pop2_pol)) + geom_line(colour="steelblue") + labs(x=NULL) + ylab("Taj_D\npop2") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + ylim(c(-2.5,2.5)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_Taj_D_pop2_pol_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)

head(pc1_varpro)

ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$Taj_D_pop2_pol) & is.finite(pc1_varpro$Taj_D_pop1_pol),], aes(start_pos, Taj_D_pop1_pol, colour="Pop1")) + geom_line(colour="darkred") + 
  geom_line(aes(start_pos, Taj_D_pop2_pol),colour="steelblue") +
  labs(x=NULL) + ylab("Taj_D") + scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + theme_classic() + theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") #+ ggsave(paste0("plot_Taj_D_pop2_pol_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)

ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$Taj_D_pop2_pol) & is.finite(pc1_varpro$Taj_D_pop1_pol),], aes(Taj_D_pop1_pol, Taj_D_pop2_pol)) + geom_point(colour="black", alpha=0.6) + 
  geom_hline(yintercept = 0, colour="grey40") + 
  geom_vline(xintercept=0, colour="grey40") +
  theme_classic() + theme(axis.line = element_line(colour = "black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + facet_grid(. ~ chromosome_name) + ggsave(paste0("plot_Taj_D_pop1_Vs_pop2_pol_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)

ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$Taj_D_pop2_pol) & is.finite(pc1_varpro$Taj_D_pop1_pol),], aes(chromosome_name, Taj_D_pop1_pol)) + geom_boxplot(fill="darkred") + 
  ylim(c(-2.7,2.7)) +
  geom_hline(yintercept = 0, colour="grey40") + 
  #geom_vline(xintercept=0, colour="grey40") +
  theme_classic() + theme(axis.line = element_line(colour = "black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + ggsave(paste0("plot_Taj_D_pop1_perChro_", window_size, "_", overlap, ".png"), width = 5, height = 5, dpi = 400)


ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$Taj_D_pop2_pol) & is.finite(pc1_varpro$Taj_D_pop1_pol),], aes(chromosome_name, Taj_D_pop2_pol)) + geom_boxplot(fill="steelblue") +
  ylim(c(-2.7,2.7)) +
  geom_hline(yintercept = 0, colour="grey40") + 
  #geom_vline(xintercept=0, colour="grey40") +
  theme_classic() + theme(axis.line = element_line(colour = "black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + ggsave(paste0("plot_Taj_D_pop2_perChro_", window_size, "_", overlap, ".png"), width = 5, height = 5, dpi = 400)

##
##
# ggplot(pc1_varpro[pc1_varpro$pc1_varcom>0.5 & is.finite(pc1_varpro$Taj_D_pop2_pol) & is.finite(pc1_varpro$Taj_D_pop1_pol) & pc1_varpro$chromosome_name=="III",], 
#   aes(start_pos, Taj_D_pop1_pol, colour="Pop1")) + 
#   geom_rect(data=wtf_reg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="grey50", alpha=0.2, inherit.aes = FALSE) + 
#   geom_rect(data=wtf_actreg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="green", alpha=0.2, inherit.aes = FALSE) + 
#   geom_line(colour="darkred") + 
#   geom_line(aes(start_pos, Taj_D_pop2_pol), colour="steelblue") +
#   labs(x=NULL) + ylab("Taj_D") + 
#   scale_x_continuous(breaks=seq(500000,5500000,500000), labels = scales::scientific, expand=c(0,0)) + 
#   theme_classic() + theme(axis.line = element_line(colour="black"), axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), strip.background = element_blank(), axis.text.y = element_text(colour="black")) + 
#   facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") #+ ggsave(paste0("plot_Taj_D_pop2_pol_", window_size, "_", overlap, ".png"), width = 12, height = 5, dpi = 400)





#### inference of number of join points between haplotypes:

head(dm)

##
#

head(data)
head(data_small_windows)
head(dm)
names(dm)<-c("chromosome_name", "start_pos", "sample", "Nor_PC1_pol")


dm <- dm %>% mutate(join_points=rep(0,length(dm[,1])))
dm$chromosome_name<-as.factor(dm$chromosome_name)

for ( chromosome in levels(dm$chromosome_name)){
  print(chromosome)
  for (sampleID in levels(dm$sample)){
    print(sampleID)
    table<- dm %>% filter(chromosome_name==chromosome & sample==sampleID)
    haplotype<-table[table$start_pos==min(sort(table$start_pos)),4]
    if (haplotype<0.5){
      haplotype<-0
    } else {
      haplotype<-1
    }
    for (pos in sort(table$start_pos)){
      pc_value<-(table[table$start_pos==pos,4]>0.5)*1
      if (pc_value!=haplotype){
        dm[dm$chromosome_name==chromosome & dm$sample==sampleID & dm$start_pos==pos,5]<-1
        haplotype<-pc_value
      }
    }
  }
}


write.table(dm, paste0("haplotype_breaks_", window_size, "_", overlap, ".txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

dm<-read.table(paste0("haplotype_breaks_", window_size, "_", overlap, ".txt"),T)

head(dm)
# dm %>% filter(sample=="JB1174" ) %>% head()
# 
# dm %>% filter(chromosome_name=="I" & sample=="JB1174" & start_pos %in% seq(81,88))
# 
# dm %>% filter(sample=="JB1205") %>% 
#   ggplot(aes(start_pos*1000, join_points)) +
#   geom_point(alpha=0.5) +
#   facet_grid(. ~ chromosome_name, scale="free_x", space="free_x")
# 
# dm %>% filter(sample=="JB22" & join_points==1) 
# dm %>% filter(chromosome_name=="II" & sample=="JB22" & start_pos %in% seq(1683,1690))
# 
# dm %>% filter(join_points==1 & sample!="JB1207") %>% group_by(sample) %>% 
#   summarise(count = n()) %>% ggplot(aes(count)) + 
#   geom_histogram() + 
#   scale_x_continuous(breaks=seq(0,200,10))
# 
# data.frame(sample="JB869", count=26)


head(dm)
summary_table_joins_number<-dm %>% filter(join_points==1) %>% 
  group_by(sample) %>% 
  summarise(number_joins = n())

summary_table_joins_number[summary_table_joins_number$sample=="JB1207", 2]<-"NA"
# for 161 samples I also excluded JB1169
summary_table_joins_number[summary_table_joins_number$sample=="JB1169", 2]<-"NA"

ggplot(summary_table_joins_number, aes(as.factor(sample), as.numeric(number_joins))) + 
  geom_bar(stat="identity") + 
  #ylim(c(0,150)) + 
  xlab("Sample") + ylab("Num. Joins") + 
  #scale_y_continuous(breaks=seq(0,200,20)) + 
  coord_flip() + 
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), axis.text=element_text(hjust = 1, size=9, colour="black")) + 
  ggsave("number_joins_perSample.png", width = 2, height = 10, dpi = 400)


summary_table_joins_number %>% 
  filter(!(sample %in% hybrid_samples)) %>% 
  mutate(number_joins=as.numeric(number_joins)) %>% 
  ungroup() %>% 
  select(number_joins) %>% 
  unlist() %>% 
  as.vector() %>% 
  mean()

summary_table_joins_number %>% 
  filter(sample %in% hybrid_samples) %>% 
  filter(sample!="JB1207") %>% 
  mutate(number_joins=as.numeric(number_joins)-26) %>% 
  ggplot(aes(as.factor(sample), as.numeric(number_joins))) + 
  geom_bar(stat="identity") + 
  #ylim(c(0,150)) + 
  xlab("Sample") + ylab("Num. Joins") + 
  #scale_y_continuous(breaks=seq(0,200,20)) + 
  coord_flip() + 
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), axis.text=element_text(hjust = 1, size=9, colour="black")) + 
  ggsave("number_joins_perSample_hybrids.png", width = 2, height = 10, dpi = 400)
  
summary_table_joins_number %>% 
  filter(sample %in% hybrid_samples) %>% 
  filter(sample!="JB1207") %>% 
  mutate(number_joins=as.numeric(number_joins)-26) %>% 
  ggplot(aes(number_joins)) + 
  geom_histogram() +
  #ylim(c(0,150)) + 
  #scale_y_continuous(breaks=seq(0,200,20)) + 
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"), axis.text=element_text(hjust = 1, size=9, colour="black")) + 
  ggsave("number_joins_hist_hybrids.png", width = 2, height = 10, dpi = 400)
  

summary_table_joins_number %>% 
  filter(sample %in% hybrid_samples) %>% 
  filter(sample!="JB1207") %>% 
  mutate(number_joins=as.numeric(number_joins)-26) %>% 
  ungroup() %>% 
  select(number_joins) %>% 
  unlist() %>% 
  as.vector() %>% 
  mean()



summary_table_joins_number_bychromosome<-dm %>% filter(!(sample %in% c("JB1207","JB1169"))) %>% 
  group_by(sample, chromosome_name) %>% 
  summarise(mean_number_joins = mean(join_points))

summary_table_joins_number_bychromosome %>% filter(!(sample %in% c("JB1207","JB1169"))) %>% 
  ggplot(aes(sample, mean_number_joins, fill=chromosome_name)) +
  geom_bar(stat="identity", position="dodge") +
  xlab("Sample") + ylab("Mean Num. Joins\nper Window") + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y=element_text(hjust = 1, size=9, colour="black"), 
        #axis.title.y = element_text(size=9),
        legend.position="bottom") +
  guides(fill=guide_legend(title="Chr.")) +
  ggsave("Mean_number_joins_Chromosome_perSample.png", width = 12, height = 3, dpi = 400)

dm %>% filter(!(sample %in% c("JB1207","JB1169"))) %>% group_by(sample, chromosome_name) %>% 
  summarise(number_joins = sum(join_points)) %>% 
  filter(sample!="JB1207") %>% 
  ggplot(aes(sample, number_joins, fill=chromosome_name)) +
  geom_bar(stat="identity", position="dodge") +
  xlab("Sample") + ylab("Num. Joins") + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y=element_text(hjust = 1, size=9, colour="black"), 
        #axis.title.y = element_text(size=9),
        legend.position="bottom") +
  guides(fill=guide_legend(title="Chr.")) +
  facet_grid(chromosome_name ~ .) + 
  ggsave("number_joins_Chromosome_perSample.png", width = 12, height = 9, dpi = 400)


dm %>% group_by(sample, chromosome_name) %>% 
  summarise(number_joins = sum(join_points)) %>% 
  filter(sample!="JB1207") %>% 
  ggplot(aes(chromosome_name, number_joins)) +
  geom_boxplot() +
  xlab("Chromosome") + ylab("Num. Joins") + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y=element_text(hjust = 1, size=9, colour="black"), 
        #axis.title.y = element_text(size=9),
        legend.position="bottom") +
  #guides(fill=guide_legend(title="Chr.")) +
  #facet_grid(chromosome_name ~ .) + 
  ggsave("number_joins_Chromosome.png", width = 12, height = 9, dpi = 400)


dm %>% group_by(sample, chromosome_name) %>% 
  summarise(number_joins = sum(join_points)) %>% 
  filter(sample!="JB1207") %>% 
  filter(sample %in% hybrid_samples) %>% 
  ggplot(aes(chromosome_name, number_joins)) +
  geom_boxplot() +
  xlab("Chromosome") + ylab("Num. Joins") + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y=element_text(hjust = 1, size=9, colour="black"), 
        #axis.title.y = element_text(size=9),
        legend.position="bottom")


head(summary_table_joins_number_bychromosome)
summary_table_joins_number_bychromosome %>% filter(sample!="JB1207") %>% 
  ggplot(aes(chromosome_name, mean_number_joins)) + 
  geom_boxplot()

###

# The following Equations are taken from: Janzen T., Nolte A.W., and Traulsen A. 2018.
######## ######## ######## ######## ######## ######## ######## ######## 
######## Function to calculate K, following equation 12
######## ######## ######## ######## ######## ######## ######## ######## 
calc_K <- function(N = Inf, R = Inf, H_0 = 0.5, C = 1) {
  K <- H_0*C*2*N*R/(2*N*C+R);
  if(is.infinite(N)) {
    K <- H_0*R
  }
  if(is.infinite(R)) {
    K <- H_0*C*2*N
  } 
  if(is.infinite(N) && is.infinite(R)) {
    K <- Inf
  }
  return(K)
}
######## ######## ######## ######## ######## ######## ######## ######## 
######## Function to calculate J_t, following equation 13
######## ######## ######## ######## ######## ######## ######## ######## 
calculate_J <- function(N = Inf, R = Inf, H_0 = 0.5, C = 1, maxT = 100) {
  ######## N: population size
  ######## R: number of evenly spaced recombination sites
  ######## H_0: initial heterozygosity 2pq (at t = 0)
  ######## C: size of the chromosome in Morgan
  ######## maxT: maximum number of timesteps
  ########
  ######## returns: number of junctions in t = [0, maxT], including ZERO!
  t <- 0:maxT
  if(is.infinite(N) && is.infinite(R)) {
    # If both N and R are infinite, R gives 
    # numerical problems using equation 12 
    # to calculate K, so instead we use
    # equation 1
    jt <- H_0 * C * t;
    return(jt)
  }
  K <- calc_K(N, R, H_0, C)
  jt <- K - K *(1-H_0*C/K)^t
  return(jt) #jt in [0,maxT]
}
######## ######## ######## ######## ######## ######## ######## ######## 
######## Function to calculate the relative error in estimating t
######## following equation 33 in the Supplement.
######## ######## ######## ######## ######## ######## ######## ######## 
calculate_error <- function(J = NA, N = Inf, R = Inf, 
                            H_0 = 0.5, C = 1, t = 1, 
                            relative = TRUE) {
  # the flag relative determines whether we want the error
  # relative to K, or in absolute generations (relative = FALSE)
  K <- calc_K(N, R, H_0, C)
  u <- 1 - 1/(2*N) - C/R
  error <- log(u^t-1/K) / (log(u)*t) - 1
  if(relative)  return(error)
  if(!relative) return(t*error)
}
######## ######## ######## ######## ######## ######## ######## ######## 
######## Function to calculate the time since the onset 
######## of hybridization following equation 14
######## ######## ######## ######## ######## ######## ######## ######## 
calculate_time <- function(J = NA, N = Inf, R = Inf, H_0 = 0.5, C = 1) {
  if(is.na(J)) {
    cat("ERROR! did you forget to provide J?")
    return()
  }
  if(is.infinite(N) && is.infinite(R)) {
    cat("both N and R are infinite\n")
    cat("can not estimate t\n")
  }
  K <- calc_K(N, R, H_0, C)
  u <- 1 - 1/(2*N) - C/R
  t <- log(1-J/K) / (log(u))
  return(t)
}
######## ######## ######## ######## ######## ######## ######## ######## 
######## Function to calculate the maximum accurate time
######## following equation 15
######## ######## ######## ######## ######## ######## ######## ######## 
calculate_MAT <- function(N = Inf, R = Inf, H_0 = 0.5, C = 1) {
  if(is.infinite(N) && is.infinite(R)) {
    cat("both N and R are infinite\n")
    cat("can not estimate MAT\n")
  }
  K <- calc_K(N, R, H_0, C)
  u <- 1 - 1/(2*N) - C/R
  MAT = log(1/K)/log(u)
  return(MAT)
}
######## ######## ######## ######## ######## ######## ######## ######## 
######## Demonstration code of the previous functions
######## change to demonstrate = TRUE to evaluate
######## ######## ######## ######## ######## ######## ######## ######## 

head(dm)
#summary_table_joins_number_bychromosome<-
cal_time_f <- function(j_vector=c(1,1,1),n=1000000,r_vector=5500000,h0=0.5,c=1){
  cal_time<-c()
  for (line in seq(1,length(j_vector))){
    if (length(h0)==1){
      if (length(c)==1){
        cal_time<-c(cal_time, calculate_time(J = j_vector[line], N = n, R=r_vector[line], H_0=h0, C=c))
      } else {
        cal_time<-c(cal_time, calculate_time(J = j_vector[line], N = n, R=r_vector[line], H_0=h0, C=c[line]))
      }
    } else {
      if (length(c)==1){
        cal_time<-c(cal_time, calculate_time(J = j_vector[line], N = n, R=r_vector[line], H_0=h0[line], C=c))
      } else {
        cal_time<-c(cal_time, calculate_time(J = j_vector[line], N = n, R=r_vector[line], H_0=h0[line], C=c[line]))
      }
    }
  }
  return(cal_time)
}

c_per_chromosome<-data.frame(chromosome_name=levels(dm$chromosome_name), c_chromosome=c(18.8, 14.8, 10.8))
c_per_chromosome_low<-data.frame(chromosome_name=levels(dm$chromosome_name), c_chromosome_low=c(18.8, 14.8, 10.8)*17/45)


head(dm)
head(c_per_chromosome)
str(c_per_chromosome)

time_estimate <- merge(dm, c_per_chromosome, by="chromosome_name") 
time_estimate <- merge(time_estimate, c_per_chromosome_low, by="chromosome_name") %>% 
  filter(!(sample %in% c("JB1207","JB1169"))) %>% 
  mutate(haplotype_red=1*(Nor_PC1_pol<0.5), 
              haplotype_blue=1*(Nor_PC1_pol>0.5)) %>% 
  group_by(sample, chromosome_name) %>% 
  summarise(number_joins = sum(join_points), 
            chr_length=max(start_pos),#*1000, 
            num_red=sum(haplotype_red), 
            num_blue=sum(haplotype_blue), 
            c_per_chromosome=c_chromosome[1], 
            c_per_chromosome_low=c_chromosome_low[1]) %>% 
  mutate(red_ratio=num_red/(num_red+num_blue), 
         blue_ratio=num_blue/(num_red+num_blue), 
         heterogenicity=2*red_ratio*blue_ratio, 
         time_0.1=cal_time_f(j_vector=number_joins, r_vector=chr_length, h0=0.01, c=c_per_chromosome), 
         time_0.5=cal_time_f(j_vector=number_joins, r_vector=chr_length, h0=0.5, c=c_per_chromosome), 
         time_mean=cal_time_f(j_vector=number_joins, r_vector=chr_length, h0=0.298, c=c_per_chromosome), #h0=0.264 with 161 samples
         time_H=cal_time_f(j_vector=number_joins, r_vector=chr_length, h0=heterogenicity, c=c_per_chromosome), 
         time_mean_lower=cal_time_f(j_vector=number_joins, r_vector=chr_length, h0=0.298, c=c_per_chromosome_low)) #h0=0.264 with 161 samples

mean(time_estimate$heterogenicity)

# plot of hybridization estimate per chromosome
# in this plot, samples with heterogenicity lower than 0.1 were excluded (this correspond to minor haplotype frequency > ~0.05)
time_estimate %>% filter(heterogenicity>0.1) %>% 
  ggplot(aes(time_H)) +
  geom_histogram() + 
  #geom_density() +
  xlab("Time of Hybridization\n(Sexual generations)") + 
  ylab("Count") + 
  scale_y_continuous(breaks=seq(1,15,2)) + 
  #ylim(c(0,10)) +
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y=element_text(hjust = 1, size=9, colour="black"), 
        #axis.title.y = element_text(size=9),
        legend.position="bottom") +
  facet_grid(chromosome_name ~ .) +
  ggsave("time_estimate_Chromosome.png", width = 12, height = 9, dpi = 400)

time_estimate %>% 
  ggplot(aes(time_mean)) +
  geom_histogram() + 
  #geom_density() +
  xlab("Time of Hybridization\n(Sexual generations)") + 
  ylab("Count") + 
  scale_x_continuous(breaks=seq(0,1000,5)) + 
  #ylim(c(0,10)) +
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y=element_text(hjust = 1, size=9, colour="black"), 
        #axis.title.y = element_text(size=9),
        legend.position="bottom") +
  facet_grid(chromosome_name ~ .) +
  ggsave("time_estimate_mean_H0_Chromosome.png", width = 5, height = 5, dpi = 400)

  mean(time_estimate$time_mean)

time_estimate %>% filter(heterogenicity>0.1) %>% 
  ggplot(aes(time_mean_lower)) +
  geom_histogram() + 
  #geom_density() +
  xlab("Time of Hybridization\n(Sexual generations)") + 
  ylab("Count") + 
  #scale_x_continuous(breaks=seq(0,1000,5)) + 
  #ylim(c(0,10)) +
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y=element_text(hjust = 1, size=9, colour="black"), 
        #axis.title.y = element_text(size=9),
        legend.position="bottom") +
  facet_grid(chromosome_name ~ .) +
  ggsave("time_estimate_mean_H0_Chromosome_lowerC.png", width = 5, height = 5, dpi = 400)


time_estimate %>% filter(heterogenicity>0.1) %>% 
  group_by(sample) %>% 
  summarise(mean_time = mean(time_mean_lower)) %>% 
  ungroup() %>% 
  select(mean_time) %>%
  unlist() %>% as.vector() %>% mean() # %>% var() %>% sqrt()

#mean 21.5
#var 20.1
#sd 4.49



## time using the max number of joins:

time_estimate[which]
time_estimate[which(time_estimate$number_joins==max(time_estimate$number_joins)),"time_H"]

sample<-c()
number_joins<-c()
heterogenicity<-c()
time_H<-c()
for (sampleID in levels(time_estimate$sample)){
  if (sampleID!="JB1207"){
  data<-time_estimate[time_estimate$sample==sampleID,]
  print(data)
  sample<-c(sample, data[which(data$number_joins==max(data$number_joins)),][1,"sample"] %>% unlist() %>% as.vector())
  number_joins<-c(number_joins, data[which(data$number_joins==max(data$number_joins)),][1,"number_joins"] %>% unlist() %>% as.vector())
  heterogenicity<-c(heterogenicity, data[which(data$number_joins==max(data$number_joins)),][1,"heterogenicity"] %>% unlist() %>% as.vector())
  time_H<-c(time_H, data[which(data$number_joins==max(data$number_joins)),][1,"time_H"] %>% unlist() %>% as.vector())
  }
}

summary_time_estimate<-data.frame(sample, number_joins, heterogenicity, time_H) 
names(summary_time_estimate)<-c("sampleID", "num_joins", "heterogenicitySample", "time_HSample")
summary_time_estimate %>% 
  filter(heterogenicitySample>0.1) %>% 
  ggplot(aes(time_HSample)) +
  geom_histogram() +
  xlab("Time of Hybridization\n(Sexual generations)") + 
  ylab("Count") + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), 
        axis.text.y=element_text(hjust = 1, size=9, colour="black"), 
        #axis.title.y = element_text(size=9),
        legend.position="none") +
  ggsave("time_estimate_maxNJ.png", width = 4, height = 3, dpi = 400)

summary_time_estimate %>% filter(heterogenicitySample>0.1) %>% 
  dplyr::select(time_HSample) %>% unlist() %>% mean()

summary_time_estimate %>% filter(heterogenicitySample>0.1) %>% 
  dplyr::select(time_HSample) %>% unlist() %>% var() %>% sqrt()


### Correlation between number of joins and coleccion date:

collection_inf<-read.table("collection_info.txt",T)

time_estimate_withDate<-merge(time_estimate, collection_inf, by="sample")

time_estimate_withDate$date_collected <- as.character(time_estimate_withDate$date_collected)
time_estimate_withDate$date_collected[time_estimate_withDate$date_collected=="ND"]<-NA
time_estimate_withDate$date_collected <- as.numeric(time_estimate_withDate$date_collected)

# 
# time_estimate_withDate %>% filter(date_collected!="NA") %>% 
# ggplot(aes(date_collected, number_joins)) + geom_point() + facet_grid(chromosome_name ~ continent_collected) 
# 
# 
# no_hyb_samples<-time_estimate_withDate %>% filter(chromosome_name=="I" & number_joins<60) %>% select(sample) %>% unlist() %>% as.vector()
# 
# time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND" & !(sample %in% no_hyb_samples)) %>% 
# ggplot(aes(date_collected, number_joins)) + geom_smooth(method='lm',formula=y~x) + 
#   geom_point() + facet_grid(chromosome_name ~ ., scale="free") 
# 
# table_cor_sampling<-time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND") 
# sp <- ggscatter(table_cor_sampling, x = "date_collected", y = "number_joins",
#    color = "chromosome_name", palette = "jco",
#    add = "reg.line", conf.int = TRUE)
# sp + stat_cor(aes(color = chromosome_name), label.x = 1980)
# 
# table_cor_sampling<-time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND" & !(sample %in% no_hyb_samples)) 
# sp <- ggscatter(table_cor_sampling, x = "date_collected", y = "number_joins",
#    color = "chromosome_name", palette = "jco",
#    add = "reg.line", conf.int = TRUE)
# sp + stat_cor(aes(color = chromosome_name), label.x = 1980)
# 
# no_hyb_samplesI<-time_estimate_withDate %>% filter(chromosome_name=="I" & time_mean<11) %>% select(sample) %>% unlist() %>% as.vector()
# no_hyb_samplesII<-time_estimate_withDate %>% filter(chromosome_name=="II" & time_mean<9) %>% select(sample) %>% unlist() %>% as.vector()
# no_hyb_samplesIII<-time_estimate_withDate %>% filter(chromosome_name=="III" & time_mean<4) %>% select(sample) %>% unlist() %>% as.vector()
# 
# no_hyb_samples<-c(no_hyb_samplesI, no_hyb_samplesII, no_hyb_samplesIII)
# 
# table_cor_sampling<-time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND" & !(sample %in% no_hyb_samples)) %>% 
# ggscatter(x = "date_collected", y = "time_mean",
#    color = "chromosome_name", palette = "jco",
#    add = "reg.line", conf.int = TRUE)
# table_cor_sampling + stat_cor(aes(color = chromosome_name), label.x = 1980)

table_cor_sampling<-time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND") %>% 
ggscatter(x = "date_collected", y = "time_mean",
   color = "chromosome_name", palette = "jco",
   add = "reg.line", conf.int = TRUE)
table_cor_sampling + stat_cor(aes(color = chromosome_name), label.x = 1980) + ggsave("plot_cor_date_hibTime.png")


table_cor_sampling<-time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND") %>% 
ggscatter(x = "date_collected", y = "time_mean_lower",
   color = "chromosome_name", palette = "jco",
   add = "reg.line", conf.int = TRUE)
table_cor_sampling + stat_cor(aes(color = chromosome_name), label.x = 1980) + ggsave("plot_cor_date_hibTime_lowerC.png")



cor_table<-time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND" & chromosome_name=="I") 
summary(lm(time_mean ~ date_collected , data = cor_table))
-(coef(lm(time_mean ~ date_collected , data = cor_table))[1])/coef(lm(time_mean ~ date_collected , data = cor_table))[2]

cor_table<-time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND" & chromosome_name=="II") 
summary(lm(time_mean ~ date_collected , data = cor_table))
-(coef(lm(time_mean ~ date_collected , data = cor_table))[1])/coef(lm(time_mean ~ date_collected , data = cor_table))[2]

cor_table<-time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND" & chromosome_name=="III") 
summary(lm(time_mean ~ date_collected , data = cor_table))
-(coef(lm(time_mean ~ date_collected , data = cor_table))[1])/coef(lm(time_mean ~ date_collected , data = cor_table))[2]



cor_table<-time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND" & chromosome_name=="I") 
summary(lm(time_mean_lower ~ date_collected , data = cor_table))
-(coef(lm(time_mean_lower ~ date_collected , data = cor_table))[1])/coef(lm(time_mean_lower ~ date_collected , data = cor_table))[2]

cor_table<-time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND" & chromosome_name=="II") 
summary(lm(time_mean_lower ~ date_collected , data = cor_table))
-(coef(lm(time_mean_lower ~ date_collected , data = cor_table))[1])/coef(lm(time_mean_lower ~ date_collected , data = cor_table))[2]

cor_table<-time_estimate_withDate %>% filter(date_collected!="NA" & continent_collected!="ND" & chromosome_name=="III") 
summary(lm(time_mean_lower ~ date_collected , data = cor_table))
-(coef(lm(time_mean_lower ~ date_collected , data = cor_table))[1])/coef(lm(time_mean_lower ~ date_collected , data = cor_table))[2]









# producing PLINK file for heritability analysis:

head(matrix_window)

ancestralHap_matrix<-matrix_window
for (win in seq(1,dim(ancestralHap_matrix)[2])){
  values<-ancestralHap_matrix[,win] %>% as.vector()
  ancestralHap_matrix[,win]<-paste0("pop",(values>0.5)*2+(values<0.5)*1)
}

windowID<-data %>% dplyr::select(chromosome_name, window_number) %>% mutate(windowID=paste0(chromosome_name, "_", window_number)) %>% unique()
windowID<-windowID[with(windowID, order(chromosome_name, window_number)),] %>% ungroup()  %>% dplyr::select(windowID) %>% unlist() %>% as.vector() 
windowID_rep<-c()
for (win in windowID){
  windowID_rep<-c(windowID_rep, rep(win,2))
}

ancestralHap_matrix<-ancestralHap_matrix[,windowID_rep]

plink_ancestral_hap<-data.frame(names1=rownames(ancestralHap_matrix), names2=rownames(ancestralHap_matrix), rep1=rep(0,length(rownames(ancestralHap_matrix))), rep2=rep(0,length(rownames(ancestralHap_matrix))), rep3=rep(0,length(rownames(ancestralHap_matrix))), rep4=rep(0,length(rownames(ancestralHap_matrix))),  ancestralHap_matrix) 

write.table(plink_ancestral_hap, "ancestralhap.plink.ped", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

plink_ancestral_hap_map<-data %>% dplyr::select(chromosome_name, start_pos, end_pos) %>% mutate(chr_code=(chromosome_name=="I")*1+(chromosome_name=="II")*2+(chromosome_name=="III")*3, variant=paste0(chromosome_name,":",start_pos,"..", end_pos), other=0) %>%  unique() 

plink_ancestral_hap_map<-plink_ancestral_hap_map[with(plink_ancestral_hap_map, order(chromosome_name, start_pos)),] %>% ungroup()  %>% dplyr::select(chr_code, variant, other, start_pos) 

write.table(plink_ancestral_hap_map, "ancestralhap.plink.map", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# plink_ancestral_hap_bim<-data %>% dplyr::select(chromosome_name, start_pos, end_pos) %>% mutate(chr_code=(chromosome_name=="I")*1+(chromosome_name=="II")*2+(chromosome_name=="III")*3, variant=paste0(chromosome_name,":",start_pos,"..", end_pos), other=0, genotype1="pop1", genotype2="pop2") %>%  unique() 
# plink_ancestral_hap_bim<-plink_ancestral_hap_bim[with(plink_ancestral_hap_bim, order(chromosome_name, start_pos)),] %>% ungroup()  %>% dplyr::select(chr_code, variant, other, start_pos, genotype1, genotype2) 
# write.table(plink_ancestral_hap_bim, "ancestral_hap.bim", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# plink_ancestral_hap_fam<-data.frame(names1=rownames(ancestralHap_matrix), names2=rownames(ancestralHap_matrix), rep1=rep(0,length(rownames(ancestralHap_matrix))), rep2=rep(0,length(rownames(ancestralHap_matrix))), rep3=rep(0,length(rownames(ancestralHap_matrix))), rep4=rep((-9),length(rownames(ancestralHap_matrix)))) 
# write.table(plink_ancestral_hap_fam, "ancestral_hap.fam", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# plink_ancestral_hap_nosex<-data.frame(names1=rownames(ancestralHap_matrix), names2=rownames(ancestralHap_matrix))
# write.table(plink_ancestral_hap_nosex, "ancestral_hap.nosex", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE) 


### inpuut file for comparison with SV calls:





#
#
# produce bootstrap support for cluster
library(pvclust)
pvclust_fit <- pvclust(t(matrix_window), method.hclust="ward", method.dist="euclidean")
png(filename="pvclust_rearrangements_tree_95.png", width = 1400, height = 1000, pointsize = 20)
plot(pvclust_fit) # dendogram with p values
pvrect(pvclust_fit, alpha=.95) # add rectangles around groups highly supported by the data
dev.off()

####

# tree with groups by chromosome:

d <- dist(matrix_window, method = "euclidean")
H.fit <- hclust(d, method="ward")
H.fit$labels<-samples_ID
plot(H.fit)


tree_ed<-as.phylo(H.fit) 
tree_ed2<-ggtree(tree_ed) %>% rotate(63) %>% rotate(67) %>% rotate(68) %>% rotate(75) %>% rotate(91) %>% rotate(92) %>% rotate(71) %>% rotate(89)  %>% rotate(90)
tree_ed2 + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab() + ggsave(paste0("plot_rearrangements_rawTree_", window_size, "_", overlap, ".png"), width = 18, height = 12, dpi = 600)

samples_ID
H.fit$labels[H.fit$order]

groups_rearrangemets<-read.table("groups_rearrangemets.txt", T) 
row.names(groups_rearrangemets) <- NULL

groups_rearrangemets_heatmap<-as.data.frame(groups_rearrangemets[,2:5])
row.names(groups_rearrangemets_heatmap)<-groups_rearrangemets[,1]
groups_rearrangemets_heatmap$Group<-as.factor(groups_rearrangemets_heatmap$Group)
groups_rearrangemets_heatmap$Chr_I<-as.factor(groups_rearrangemets_heatmap$Chr_I)
groups_rearrangemets_heatmap$Chr_II<-as.factor(groups_rearrangemets_heatmap$Chr_II)
groups_rearrangemets_heatmap$Chr_III<-as.factor(groups_rearrangemets_heatmap$Chr_III)


# mycols <- colourPicker(1)
mycols <-c("lightgrey", "darkorange1", "darkgreen", "firebrick", "royalblue4", "darkgoldenrod", "cadetblue", "deeppink4", "firebrick1")

tree_ed3<- tree_ed2 + geom_tiplab() + geom_tippoint() + theme_tree2()  
gheatmap(tree_ed3, groups_rearrangemets_heatmap,offset = 13, width=0.2, font.size=4, colnames_angle=-45, hjust=0, color="black") + scale_fill_manual(breaks=c(as.character(0:8)), values=mycols) + theme(legend.position="none") + ggsave(paste0("plot_rearrangements_Tree_withGroups", window_size, "_", overlap, ".png"), width = 7, height = 13, dpi = 600)




#### checking for max values in axes by sample:

data_small_windows<-matrix(rep(c(rep(-1, 5563146/5000), rep(-1, 4537806/5000), rep(-1, 2440870/5000)), length(levels(data$sample))), ncol=length(levels(data$sample)))

data_small_windows <- as.data.frame(data_small_windows)
names(data_small_windows) <- levels(data$sample)

data_small_windows <- data_small_windows %>% mutate(chromosome=c(rep("I", 5563146/5000), rep("II", 4537806/5000), rep("III", 2440870/5000)), position_bin=c(seq(0,(5563146/5000)-1), seq(0, (4537806/5000)-1), seq(0, (2440870/5000)-1)))

data_small_windows$position_bin<-as.numeric(data_small_windows$position_bin)

mark_values <- function(a,b,c,d,e){ data_small_windows[data_small_windows$chromosome==a & data_small_windows$position_bin %in% c(round((b/5000)):round((c/5000))), d] <<- e}

apply(data[,c('chromosome_name','start_pos','end_ed', 'sample', 'Nor_PC1')], 1 , function(x) mark_values(x[1],as.numeric(x[2]),as.numeric(x[3]),x[4],as.numeric(x[5])))


head(data)
# add bins column
dm <- melt(data_small_windows,id.var=c("chromosome", "position_bin"))
names(dm)<-c("chromosome_name", "start_pos", "sample", "Nor_PC1")
dm <- dm %>% mutate(max_values=((Nor_PC1==0)*(-1))+((Nor_PC1==1)*1)) %>% filter(max_values!=0)
dm$max_values[dm$max_values==(-1)]<-0
# dm$max_values[dm$max_values==(-1)]<-"darkred"
# dm$max_values[dm$max_values==(1)]<-"steelblue"
dm$sample <- factor(dm$sample, levels = H.fit$labels[H.fit$order])
ggplot(dm, aes(start_pos*5000, sample)) + geom_tile(aes(fill = max_values)) + scale_fill_gradient(low = "darkred", high = "steelblue") + scale_x_continuous(breaks=seq(500000,5500000,500000), expand=c(0,0)) + theme_classic() + labs(x=NULL) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_blank(), legend.position="none") + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_maxPC1_heatmap_0_", "_", window_size, "_", overlap, ".png"), width = 18, height = 12, dpi = 800)


dm <- melt(data_small_windows,id.var=c("chromosome", "position_bin"))
names(dm)<-c("chromosome_name", "start_pos", "sample", "Nor_PC1")
dm <- dm %>% mutate(max_values=((Nor_PC1<0.02)*(-1))+((Nor_PC1>0.98)*1)) %>% filter(max_values!=0)
dm$max_values[dm$max_values==(-1)]<-0
dm$sample <- factor(dm$sample, levels = H.fit$labels[H.fit$order])
# dm$max_values[dm$max_values==(-1)]<-"darkred"
# dm$max_values[dm$max_values==(1)]<-"steelblue"
dm$sample <- factor(dm$sample, levels = H.fit$labels[H.fit$order])
ggplot(dm, aes(start_pos*5000, sample)) + geom_tile(aes(fill = max_values)) + scale_fill_gradient(low = "darkred", high = "steelblue") + scale_x_continuous(breaks=seq(500000,5500000,500000), expand=c(0,0)) + theme_classic() + labs(x=NULL) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_blank(), legend.position="none") + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plot_maxPC1_heatmap_2_", "_", window_size, "_", overlap, ".png"), width = 18, height = 12, dpi = 800)



# bar plot per sample:
data$sample <- factor(data$sample, levels = H.fit$labels[H.fit$order])
data %>% mutate(max_values=abs(((Nor_PC1<0.02)*(-1))+((Nor_PC1>0.98)*1))) %>% group_by(sample) %>% summarise(count_maxPC1=sum(max_values)) %>% ggplot(aes(x=sample, y=count_maxPC1)) + geom_bar(stat="identity") + theme_classic() +  coord_flip() + ggsave(paste0("plot_barplot_maxPC1_perSample_", "_", window_size, "_", overlap, ".png"), width = 5, height = 12, dpi = 600)







#### heatmaps and clustering by chromosome:

for (chromosome_used in c("I", "II", "III")){
  
  matrix_window <- data %>% filter(chromosome_name==chromosome_used) %>% dplyr::select(chromosome_name, window_number, start_pos, sample, Nor_PC1_pol) %>% mutate(windowID=paste0(chromosome_name, "_", window_number)) %>% dplyr::select(-chromosome_name, -window_number, -start_pos) %>% group_by(sample) %>% spread(windowID, Nor_PC1_pol) %>% ungroup() %>% dplyr::select(paste0(data$chromosome_name[data$chromosome_name==chromosome_used], "_", data$window_number[data$chromosome_name==chromosome_used])) %>% as.matrix()
  
  samples_ID<- data %>% filter(chromosome_name==chromosome_used) %>% dplyr::select(chromosome_name, window_number, start_pos, sample, Nor_PC1_pol) %>% mutate(windowID=paste0(chromosome_name, "_", window_number)) %>% dplyr::select(-chromosome_name, -window_number, -start_pos) %>% group_by(sample) %>% spread(windowID, Nor_PC1_pol) %>% ungroup() %>% dplyr::select(sample) %>% unlist() %>% as.vector()
  
  row.names(matrix_window)<- samples_ID
  
  d <- dist(matrix_window, method = "euclidean")
  H.fit <- hclust(d, method="ward")
  H.fit$labels<-samples_ID
  plot(H.fit)
  
  pvclust_fit <- pvclust(t(matrix_window), method.hclust="ward", method.dist="euclidean")
  png(filename=paste0("pvclust_rearrangements_chromosome_", chromosome_used ,"_tree_95.png"), width = 1400, height = 1000, pointsize = 20)
  plot(pvclust_fit) # dendogram with p values
  pvrect(pvclust_fit, alpha=.95) # add rectangles around groups highly supported by the data
  dev.off()
  
  dm$sample <- factor(dm$sample, levels = H.fit$labels[H.fit$order])
  
  ggplot(dm[dm$chromosome_name==chromosome_used,], aes(start_pos*5000, sample)) + geom_tile(aes(fill = Nor_PC1_pol)) + scale_fill_gradient(low = "darkred", high = "steelblue", expand=c(0,0)) + scale_x_continuous(breaks=seq(500000,5500000,500000), expand=c(0,0)) + theme_classic() + labs(x=NULL) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_blank(), legend.position="none") + facet_grid(. ~ chromosome_name, scale="free_x", space="free_x") + ggsave(paste0("plotheatmap_chromosome_", chromosome_used, "_", window_size, "_", overlap, ".png"), width = 8, height = 6, dpi = 800)
  
}





####
# 17 - 12 - 2018
# linkage disequilibrium analysis (analysis between windows):
####

# I only use hybrid samples for this section:
data_hybrid_sim<-data %>% 
	filter(sample %in% hybrid_samples) %>% 
  filter(pc1_varcom>0.4) %>%
	mutate(windowID=paste(chromosome_name, window_number, start_pos, end_pos, sep="__")) %>% 
	ungroup() %>% 
	select(windowID, sample, Nor_PC1_pol_sim)


write.table(data_hybrid_sim, "data_hybrid_sim.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)
data_hybrid_sim<-read.table("data_hybrid_sim.txt", T)

D_cal<-function(win_1="I__168__1186804", win_2="I__168__1186804"){
	estimated_val<-data_hybrid_sim %>% 
	filter(windowID %in% c(win_1, win_2)) %>% 
	spread(windowID, Nor_PC1_pol_sim) %>% 
	filter(get(win_1)!=0.5 & get(win_2)!=0.5) %>% 
	mutate(ind_p1q1=(get(win_1)==1 & get(win_2)==1)*1) %>% 
	ungroup() %>% 
	summarise(total=n(), 
		Np1=sum(get(win_1)==1), Nq1=sum(get(win_2)==1), Np1q1=sum(ind_p1q1), 
		p1=Np1/total, q1=Nq1/total, p1q1=Np1q1/total, 
		D=p1q1-(p1*q1), 
		Dp=ifelse(D<0, D/(max(c((-p1*q1),(-(1-p1)*(1-q1))))), 
			D/(min(c((p1*(1-q1)),((1-p1)*(q1)))))), 
		r2=(D^2)/(p1*(1-p1)*q1*(1-q1))) %>% 
	unlist() %>% 
	as.vector()
	return(paste(estimated_val, collapse ="___"))
}

list_windows<-unique(data_hybrid_sim$windowID)
LD_table<-data.frame(t(combn(list_windows, 2)))

colnames(LD_table)<-c("window1", "window2")
LD_table$window1<-as.vector(LD_table$window1)
LD_table$window2<-as.vector(LD_table$window2)


out <- data.frame(matrix(c(NA,NA,NA),ncol=3))
colnames(out) <- c("window1","window2","results")
write.table(out, "LD_results.txt", append=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
for (line in seq(1,dim(LD_table)[1])) {
	out[1,1]<-LD_table[line,1]
	out[1,2]<-LD_table[line,2]
	out[1,3]<-D_cal(LD_table[line,1], LD_table[line,2])
	write.table(out, "LD_results.txt", append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
}

LD_results<-read.table("LD_results.txt", T)


sbatch /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/10_HGAP_assembly_pacbio/scripts/35_LD.sh I
sbatch /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/10_HGAP_assembly_pacbio/scripts/35_LD.sh II
sbatch /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/10_HGAP_assembly_pacbio/scripts/35_LD.sh III

sbatch /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/10_HGAP_assembly_pacbio/scripts/35_LD_2chr.sh I II
sbatch /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/10_HGAP_assembly_pacbio/scripts/35_LD_2chr.sh I III
sbatch /home/sergio/private/Uppsala/Analyses/04_Genomic_analyses/10_HGAP_assembly_pacbio/scripts/35_LD_2chr.sh II III




# 11 - 06 - 2018
### 2d SFS using SNP data:

# This script was run in the directory:
# /proj/uppstore2017159/b2014286_nobackup/private/pac_bio/03_denovo_assembly/natural_strains_group1/alignment/07_genome_alignments_MI/fineSTRUCTURE/byWindow/different_WindowSize


module load R/3.4.3 MariaDB/10.2.11
R_LIBS_USER=/home/sergio/R/libraries_Rackham/3.4
R


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

#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges",  lib="/home/sergio/R/libraries_Rackham/3.4")


ancestralhap<-read.table("ancestralhap.txt", T) 

vcf_file_name<-"chromosome_III_all57_clean_withGAPS.recode.vcf"
temporal_number<-runif(1,0,100)
vcf <- readVcf(vcf_file_name)
snpgdsVCF2GDS("./chromosome_III_all57_clean_withGAPS.recode.vcf", paste0("temporal_", temporal_number, ".gds"), method="biallelic.only")
snpgdsSummary(paste0("temporal_", temporal_number, ".gds"))
genofile <- snpgdsOpen(paste0("temporal_", temporal_number, ".gds"))
genotypes <- t(snpgdsGetGeno(genofile))
# read.gdsn(index.gdsn(genofile, "sample.id"))
# read.gdsn(index.gdsn(genofile, "genotype"))
# read.gdsn(index.gdsn(genofile, "snp.position"))
# read.gdsn(index.gdsn(genofile, "snp.allele"))

colnames(genotypes)<-read.gdsn(index.gdsn(genofile, "sample.id"))
genotypes_III<-as.data.frame(genotypes) %>% mutate(ID=paste0("III_",read.gdsn(index.gdsn(genofile, "snp.position"))), chr="III", pos=read.gdsn(index.gdsn(genofile, "snp.position"))) %>% gather("sample", "GT", 1:57) 
ancestralhap_57samples<-ancestralhap %>% filter(sample %in% genotypes_III$sample) %>% mutate(Nor_PC1_pol=(Nor_PC1_pol==2)*1, windowID=as.factor(paste0(chromosome_name,"_", start_pos,"_", end_ed)))
SNP_list_genotypes_III<-genotypes_III %>% select(chr, pos) %>% mutate(ID=paste0(chr, "_", pos)) %>% unique() 
windowID_list<-ancestralhap_57samples %>% select(chromosome_name,start_pos,end_ed,windowID) %>% filter(chromosome_name=="III") %>% unique() 
isnps <- with(SNP_list_genotypes_III, IRanges(pos, width=1, names=ID))
igenes <- with(windowID_list, IRanges(start_pos, end_ed, names=windowID))
olaps <- findOverlaps(isnps, igenes)
genotypes_III_andGT<-cbind(SNP_list_genotypes_III[queryHits(olaps),], windowID_list[subjectHits(olaps),]) %>% select(ID, windowID) %>% merge(genotypes_III, by="ID") 
genotypes_III_andGT<-genotypes_III_andGT %>% merge(ancestralhap_57samples %>% select(windowID, sample, Nor_PC1_pol), by=c("windowID", "sample"))
snpgdsClose(genofile)



vcf_file_name<-"chromosome_II_all57_clean_withGAPS.recode.vcf"
temporal_number<-runif(1,0,100)
vcf <- readVcf(vcf_file_name)
snpgdsVCF2GDS("./chromosome_II_all57_clean_withGAPS.recode.vcf", paste0("temporal_", temporal_number, ".gds"), method="biallelic.only")
snpgdsSummary(paste0("temporal_", temporal_number, ".gds"))
genofile <- snpgdsOpen(paste0("temporal_", temporal_number, ".gds"))
genotypes <- t(snpgdsGetGeno(genofile))
# read.gdsn(index.gdsn(genofile, "sample.id"))
# read.gdsn(index.gdsn(genofile, "genotype"))
# read.gdsn(index.gdsn(genofile, "snp.position"))
# read.gdsn(index.gdsn(genofile, "snp.allele"))

colnames(genotypes)<-read.gdsn(index.gdsn(genofile, "sample.id"))
genotypes_II<-as.data.frame(genotypes) %>% mutate(ID=paste0("II_",read.gdsn(index.gdsn(genofile, "snp.position"))), chr="II", pos=read.gdsn(index.gdsn(genofile, "snp.position"))) %>% gather("sample", "GT", 1:57) 
ancestralhap_57samples<-ancestralhap %>% filter(sample %in% genotypes_II$sample) %>% mutate(Nor_PC1_pol=(Nor_PC1_pol==2)*1, windowID=as.factor(paste0(chromosome_name,"_", start_pos,"_", end_ed)))
SNP_list_genotypes_II<-genotypes_II %>% select(chr, pos) %>% mutate(ID=paste0(chr, "_", pos)) %>% unique() 
windowID_list<-ancestralhap_57samples %>% select(chromosome_name,start_pos,end_ed,windowID) %>% filter(chromosome_name=="II") %>% unique() 
isnps <- with(SNP_list_genotypes_II, IRanges(pos, width=1, names=ID))
igenes <- with(windowID_list, IRanges(start_pos, end_ed, names=windowID))
olaps <- findOverlaps(isnps, igenes)
genotypes_II_andGT<-cbind(SNP_list_genotypes_II[queryHits(olaps),], windowID_list[subjectHits(olaps),]) %>% select(ID, windowID) %>% merge(genotypes_II, by="ID") 
genotypes_II_andGT<-genotypes_II_andGT %>% merge(ancestralhap_57samples %>% select(windowID, sample, Nor_PC1_pol), by=c("windowID", "sample"))
snpgdsClose(genofile)







vcf_file_name<-"chromosome_I_all57_clean_withGAPS.recode.vcf"
temporal_number<-runif(1,0,100)
vcf <- readVcf(vcf_file_name)
snpgdsVCF2GDS("./chromosome_I_all57_clean_withGAPS.recode.vcf", paste0("temporal_", temporal_number, ".gds"), method="biallelic.only")
snpgdsSummary(paste0("temporal_", temporal_number, ".gds"))
genofile <- snpgdsOpen(paste0("temporal_", temporal_number, ".gds"))
genotypes <- t(snpgdsGetGeno(genofile))
# read.gdsn(index.gdsn(genofile, "sample.id"))
# read.gdsn(index.gdsn(genofile, "genotype"))
# read.gdsn(index.gdsn(genofile, "snp.position"))
# read.gdsn(index.gdsn(genofile, "snp.allele"))

colnames(genotypes)<-read.gdsn(index.gdsn(genofile, "sample.id"))
genotypes_I<-as.data.frame(genotypes) %>% mutate(ID=paste0("I_",read.gdsn(index.gdsn(genofile, "snp.position"))), chr="I", pos=read.gdsn(index.gdsn(genofile, "snp.position"))) %>% gather("sample", "GT", 1:57) 
ancestralhap_57samples<-ancestralhap %>% filter(sample %in% genotypes_I$sample) %>% mutate(Nor_PC1_pol=(Nor_PC1_pol==2)*1, windowID=as.factor(paste0(chromosome_name,"_", start_pos,"_", end_ed)))
SNP_list_genotypes_I<-genotypes_I %>% select(chr, pos) %>% mutate(ID=paste0(chr, "_", pos)) %>% unique() 
windowID_list<-ancestralhap_57samples %>% select(chromosome_name,start_pos,end_ed,windowID) %>% filter(chromosome_name=="I") %>% unique() 
isnps <- with(SNP_list_genotypes_I, IRanges(pos, width=1, names=ID))
igenes <- with(windowID_list, IRanges(start_pos, end_ed, names=windowID))
olaps <- findOverlaps(isnps, igenes)
genotypes_I_andGT<-cbind(SNP_list_genotypes_I[queryHits(olaps),], windowID_list[subjectHits(olaps),]) %>% select(ID, windowID) %>% merge(genotypes_I, by="ID") 
genotypes_I_andGT<-genotypes_I_andGT %>% merge(ancestralhap_57samples %>% select(windowID, sample, Nor_PC1_pol), by=c("windowID", "sample"))
snpgdsClose(genofile)


genotypes_all_andGT<-rbind(genotypes_I_andGT, genotypes_II_andGT, genotypes_III_andGT)

write.table(genotypes_all_andGT, "genotypes_all_andGT.txt", quote = F, row.names = F, col.names = T)
genotypes_all_andGT<-read.table("genotypes_all_andGT.txt",T)

genotypes_all_andGT$GT1_in0<-(genotypes_all_andGT$GT==1 & genotypes_all_andGT$Nor_PC1_pol==0)*1
genotypes_all_andGT$GT1_in1<-(genotypes_all_andGT$GT==1 & genotypes_all_andGT$Nor_PC1_pol==1)*1
genotypes_all_andGT$GT0_in0<-(genotypes_all_andGT$GT==0 & genotypes_all_andGT$Nor_PC1_pol==0)*1
genotypes_all_andGT$GT0_in1<-(genotypes_all_andGT$GT==0 & genotypes_all_andGT$Nor_PC1_pol==1)*1



genotypes_all_andGT_pop<-genotypes_all_andGT %>% 
  select(-windowID) %>% 
  group_by(chr, pos, ID) %>% 
    summarise(total = n(), NAnc0=sum(Nor_PC1_pol==0), NAnc1=sum(Nor_PC1_pol==1), 
    NAnc0_in0=sum(GT0_in0), NAnc0_in1=sum(GT0_in1), 
    NAnc1_in0=sum(GT1_in0), NAnc1_in1=sum(GT1_in1),) %>% 
  mutate(fq_NAnc0_in0=NAnc0_in0/NAnc0, fq_NAnc0_in1=NAnc0_in1/NAnc1, 
       fq_NAnc1_in0=NAnc1_in0/NAnc0, fq_NAnc1_in1=NAnc1_in1/NAnc1) 

genotypes_all_andGT_pop[is.na(genotypes_all_andGT_pop)] <- 0

genotypes_all_andGT_pop<-genotypes_all_andGT_pop %>% 
  mutate(folded_fq_pop0=(((fq_NAnc1_in0+fq_NAnc1_in1)>=1)*1*fq_NAnc0_in0)+(((fq_NAnc1_in0+fq_NAnc1_in1)<1)*1*(fq_NAnc1_in0)), 
    folded_fq_pop1=(((fq_NAnc1_in0+fq_NAnc1_in1)>=1)*1*fq_NAnc0_in1)+(((fq_NAnc1_in0+fq_NAnc1_in1)<1)*1*(fq_NAnc1_in1))) 

# genotypes_all_andGT_pop<-genotypes_all_andGT_pop %>% 
#   mutate(folded_fq_pop0=((fq_NAnc1_in0<=0.5)*1*fq_NAnc1_in0)+((fq_NAnc1_in0>0.5)*1*(1-fq_NAnc1_in0)), 
#     folded_fq_pop1=((fq_NAnc1_in1<=0.5)*1*fq_NAnc1_in1)+((fq_NAnc1_in1>0.5)*1*(1-fq_NAnc1_in1))) 


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
folded_2dSFS_data$labels[folded_2dSFS_data$labels<0.001]<-NA

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
  xlab("Allele Fq. Pop 0") + 
  ylab("Allele Fq. Pop 1") + 
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
      legend.text = element_text(colour="black", size = 14)) + 
  ggsave(paste0("2dSFS_folded_SNPs_", window_size, "_", overlap, ".png"), 
      width = 6, height = 6, dpi = 400)



genotypes_all_andGT_pop_min5_chr_Total<-genotypes_all_andGT_pop_min5 %>% 
  group_by(chr) %>% 
  summarise(total_Var = n())

genotypes_all_andGT_pop_min5 %>% 
  merge(genotypes_all_andGT_pop_min5_chr_Total, by="chr") %>% 
  select(chr, folded_fq_pop0_group, folded_fq_pop1_group, total_Var) %>% 
  group_by(chr, folded_fq_pop0_group, folded_fq_pop1_group) %>%
  summarise(sfs_2d = n()/total_Var[1]) %>% 
  ggplot(aes(folded_fq_pop0_group, folded_fq_pop1_group)) + 
  geom_tile(aes(fill = sfs_2d)) + 
  geom_text(aes(label= (formatC(sfs_2d*100, format = "e", digits = 1)))) + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme_classic() + 
  xlab("Allele Fq. Pop 0") + 
  ylab("Allele Fq. Pop 1") + 
  scale_x_continuous(breaks=seq(0,1,by=0.1)) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none", axis.text.x=element_text(angle = 45, hjust = 1, size=9, colour="black"), axis.text.y = element_text(colour="black")) +
  facet_grid(. ~ chr) #+ ggsave(paste0("2dSFS_folded_SNPs_byChr_", window_size, "_", overlap, ".png"), width = 25, height = 7, dpi = 400)




