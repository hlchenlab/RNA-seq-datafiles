##Packages

library(ggplot2)
library(cowplot)
library(circlize)
library(viridis)
library(dplyr)

##Plotting SnpEff stats:

### Variant Total across samples

input = list.files("GATK", full.names = T, pattern = "SNP.ann.stats.csv")
VT_SNP_TYPE <- lapply(1:length(input), function(x) {res =read.csv(input[[x]], header =T, skip = 26, nrows = 1, blank.lines.skip =F)
res$Sample = gsub(".SNP.ann.stats.csv","",basename(input[[x]]))
return(res)}) %>% do.call(rbind.data.frame,.)

input = list.files("GATK", full.names = T, pattern = "INDEL.ann.stats.csv")
VT_INDEL_TYPE <- lapply(1:length(input), function(x) {res =read.csv(input[[x]], header =T, skip = 26, nrows = 2, blank.lines.skip =F)
res$Sample = gsub(".INDEL.ann.stats.csv","",basename(input[[x]]))
return(res)}) %>% do.call(rbind.data.frame,.)

VT_Type =rbind.data.frame(VT_SNP_TYPE, VT_INDEL_TYPE)
p1 =VT_Type %>% ggplot(aes(x= Sample, y = Count, fill =Type)) +
  geom_col() +
  scale_fill_viridis_d(option = "A") +
  theme_bw() +
  coord_flip()

## Count by genomic region
input = list.files("GATK", full.names = T, pattern = "SNP.ann.stats.csv")
VT_SNP_TYPE <- lapply(1:length(input), function(x) {res =read.csv(input[[x]], header =T, skip = 48, nrows = 14, blank.lines.skip =F)
res$Sample = gsub(".SNP.ann.stats.csv","",basename(input[[x]]))
return(res)}) %>% do.call(rbind.data.frame,.)

input = list.files("GATK", full.names = T, pattern = "INDEL.ann.stats.csv")
VT_INDEL_TYPE <- lapply(1:length(input), function(x) {res =read.csv(input[[x]], header =T, skip = 45, nrows = 18, blank.lines.skip =F)
res$Sample = gsub(".INDEL.ann.stats.csv","",basename(input[[x]]))
return(res)}) %>% do.call(rbind.data.frame,.)

load("Effecttypes.rdata")
VT_SNP_TYPE$Percent = gsub("%","",VT_SNP_TYPE$Percent) %>% as.numeric()
VT_SNP_TYPE$Type =factor(VT_SNP_TYPE$Type, levels =unique(VT_SNP_TYPE[order(VT_SNP_TYPE$Percent, decreasing = T),]$Type))
VT_SNP_TYPE$Group ="SNP"
VT_INDEL_TYPE$Percent = gsub("%","",VT_INDEL_TYPE$Percent) %>% as.numeric()
VT_INDEL_TYPE$Type =factor(VT_INDEL_TYPE$Type, levels =unique(VT_INDEL_TYPE[order(VT_INDEL_TYPE$Percent, decreasing = T),]$Type))
VT_INDEL_TYPE$Group ="INDEL"
VT_Type = rbind.data.frame(VT_SNP_TYPE, VT_INDEL_TYPE)


p2= ggplot(data =VT_Type)+
  geom_col(aes(x= Sample, y = Percent, fill =Type)) +
  scale_fill_viridis_d(option = "C",guide = guide_legend(ncol =2)) +
  theme_bw() +
  coord_flip() +
theme(legend.position = "bottom") +
  facet_grid(~Group)

## Changes by Impact-type 

## Count by genomic region
input = list.files("GATK", full.names = T, pattern = "SNP.ann.stats.csv")
VT_SNP_TYPE <- lapply(1:length(input), function(x) {res =read.csv(input[[x]], header =T, skip = 31, nrows = 4, blank.lines.skip =F)
res$Sample = gsub(".SNP.ann.stats.csv","",basename(input[[x]]))
return(res)}) %>% do.call(rbind.data.frame,.)

input = list.files("GATK", full.names = T, pattern = "INDEL.ann.stats.csv")
VT_INDEL_TYPE <- lapply(1:length(input), function(x) {res =read.csv(input[[x]], header =T, skip = 32, nrows = 3,  blank.lines.skip =F)
res$Sample = gsub(".INDEL.ann.stats.csv","",basename(input[[x]]))
return(res)}) %>% do.call(rbind.data.frame,.)

VT_SNP_TYPE$Group ="SNP"
VT_INDEL_TYPE$Group ="INDEL"
VT_Type = rbind.data.frame(VT_SNP_TYPE, VT_INDEL_TYPE)
VT_Type$Percent = gsub("%","",VT_Type$Percent) %>% as.numeric()
VT_Type$Type =factor(VT_Type$Type, levels = unique(VT_Type$Type)[c(4,2,3,1)])

p3 = ggplot(data =VT_Type,aes(x= Sample, y = Count, fill =Type))+
  geom_bar(position ="dodge", stat = "identity", width = 0.7) +
  scale_fill_viridis_d(option = "D", guide =guide_legend(reverse = T)) +
  theme_bw() +
  facet_grid(~Group)+
  coord_flip()

pa =plot_grid(p1, p3, nrow = 2, labels = LETTERS)

pdf("SNPeff_stats.pdf", height = 7, width = 10)
plot_grid(pa, p2, ncol =2, labels = c("",LETTERS[3]))
dev.off()
