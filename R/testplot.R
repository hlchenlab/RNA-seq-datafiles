##Packages

library(ggplot2)
library(cowplot)
library(circlize)
library(viridis)
library(dplyr)

plotlist =list()

for (i in 1:length(modules)) {
  try( { plotlist[[i]] =plot(qc_plot_collection(qc_out, modules[i]) +theme_bw() +
                                     theme(strip.text = element_text(size = 5),
                                           axis.title = element_text(size = 7, family = "sans"), 
                                           axis.text = element_text(size = 6, family = 'sans'),
                                           legend.text = element_text(size = 6, family = 'sans')))},
   function(e){plotlist[[i]] = NULL})
}
  
save(plotlist, file ="plotlist.rdata")
plot_grid(plotlist = plotlist[c(3:4)])

##Plotting Duplication stats

input = list.files("RMDups", full.names = T, pattern = "_dedup_metrics.txt")
dup_metrics <- lapply(1:length(input), function(x) {res =read.delim(input[[x]], header =T, skip = 6, nrows = 1, blank.lines.skip =F)
res$Sample = gsub("_dedup_metrics.txt","",basename(input[[x]]))
return(res)}) %>% do.call(rbind.data.frame,.)
ggplot(data =dup_metrics)+
  geom_col(aes(x = Sample, y =PERCENT_DUPLICATION), fill = ifelse(grepl('MOCK',dup_metrics$Sample), "steelblue2", "tan2"))


input = list.files("RMDups", full.names = T, pattern = "_dedup_metrics.txt")
dup_metrics <- lapply(1:length(input), function(x) {res =read.delim(input[[x]], header =T, skip = 9)
res$Sample = gsub("_dedup_metrics.txt","",basename(input[[x]]))
return(res)}) %>% do.call(rbind.data.frame,.)


