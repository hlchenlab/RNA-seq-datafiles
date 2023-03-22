
#Packages

if( !require("EnhancedVolcano")){
  BiocManager::install("EnhancedVolcano")
}
if( !require("DESeq2")){
  BiocManager::install("DESeq2")
}
if( !require("cowplot")){
  install.packages("cowplot")
}
if( !require("dplyr")){
  BiocManager::install("dplyr")
}
if( !require("IsoformSwitchAnalyzeR")){
  BiocManager::install("IsoformSwitchAnalyzeR", dependencies =T)
}
if( !require("data.table")){
  BiocManager::install("data.table", dependencies = T)
}
if( !require("ggplot2")){
  BiocManager::install("ggplot2")
}
if( !require("ggsignif")){
  BiocManager::install("ggsignif")
}
if( !require("GenomicFeatures")){
  BiocManager::install("GenomicFeatures")
}
if( !require("ranger")){
  BiocManager::install("ranger")
}
if( !require("ComplexHeatmap")){
  BiocManager::install("ComplexHeatmap")
}
if( !require("gridExtra")){
  install.packages("gridExtra")
}
if( !require("lemon")){
  install.packages("lemon")
}
if( !require("viridis")){
  BiocManager::install("viridis")
}
if( !require("circlize")){
  BiocManager::install("circlize")
}


#### File formatting (Only when needed): ####
library(doMC)     # Only available on UNIX Version of R
registerDoMC(8)
setwd(paste('/home/hlchen/Conor/Studies/',gseid,'/',sep=''))
#runs <- as.vector(read.delim(paste('SRR_Acc_List.txt',sep=''),header=F)[,1])
samples <- list.files('fastq', pattern = "\\.gz$") %>% sapply(function(x)strsplit(x,'[._]')[[1]][1]) %>% as.vector
pr <- duplicated(samples) %>% samples[.]
runs <- pr

#if(!dir.exists('fastq')) {system(paste('mkdir fastq'))}
setwd('fastq/')
for(i in 1:length(runs)) {
  input1 = paste0(runs[i],"_1.fq.gz" )
  input2 = paste0(runs[i],"_2.fq.gz" )
  output1 = paste0(runs[i],'_1.fastq.gz')
  output2 = paste0(runs[i],'_2.fastq.gz')
  system(paste('mv', input1, output1))
  system(paste('mv', input2, output2))
}

library(dplyr)
setwd('/home/hlchen/Conor/Studies/IAV')
studies <- list.dirs(recursive=F) %>% sapply(function(x) strsplit(x,'/',T)[[1]][2]) %>% as.vector
for (i in studies) {
  
  files <- list.files(paste0(i,"/Alignments"), full.names = T,recursive =T, pattern = "SJ.out.tab$")
  files <- lapply(files, function(x){ g =gsub("SJ.out.tab","",basename(x))
  x =read.delim(x, header = F)
  x[,c(ncol(x)+1)] = g
  x[,c(ncol(x)+1)] =i
  return(x)})
  files =do.call(rbind.data.frame, files) 
  load(paste0(i,"/metadata.rdata"))
  m =match(files$V10 , metadata$Run)
  f.a =!is.na(m)
  f.t = m[f.a]
  files$Condition =0
  files[f.a,]$Condition = metadata[f.t,]$Condition
}

ggplot(files) +
  geom_col(aes(x =Condition))

#### Make index for STAR BAMfiles ####

library(dplyr)
setwd('/home/hlchen/Conor/Studies/IAV')
studies <- list.dirs(recursive=F) %>% sapply(function(x) strsplit(x,'/',T)[[1]][2]) %>% as.vector
for(i in studies){
  samtools = "/software/samtools/1.14/bin/samtools"
  bams <- list.files(paste0(i,"/TwoPass/"), pattern = ".bam$", full.names = T)
  for (j in bams) {
    if(!file.exists(paste0(j,".bai"))){
      cmd <- paste(samtools, "index -@ 10", j)
      system(cmd)
    }
  }
}

#### Step 1: Meta data file construction: ####

setwd('/home/hlchen/Conor/Studies/IAV/')
homedir = getwd()
studies = list.files()

logname = "Log.final.out"
filedir = paste(i,"TwoPass", sep="/") # Specify file path
f1 = list.files(filedir, full.names = T, pattern = logname)

# for Loop formula: Getting file paths to get SRR and SRP numbers

Ns = list()
i = 1

for( n in f1 ){
  N = list()
  logfile = n
  
  if( file.exists(logfile) ) {
    logstats =read.delim(logfile, header = F)
    N=logstats
  } 
  
  Ns[[i]] = N
  
  i = i + 1
  
}
data = f1
mat = matrix( unlist(strsplit( as.character(data) , "/" ) ),ncol = 3, by=T )
mat[,3] = gsub("Log.final.out", "", mat[,3]  )
sraids = mat[,c(1,3)]

filenames = paste( sraids[,1], "/TwoPass/", sraids[,2], "Log.final.out", sep="")
projectid = sraids[,1]
sampleid = sraids[,2]

#  Extracting the appropriate STATs from alignment log files.

filt = c(8, 9, 24, 29, 31)
f2 <- sapply(1:length(f1), function(x) Ns[[x]][filt,2] )
f2 = t(f2)
f2 <- gsub("%", "", as.matrix(f2))
mat = as.data.frame(cbind(sampleid, f2))
colnames(mat) = c("Run", "Uniquely mapped reads number", "Unique_Mapping", "Multi-Mapping", "Unmapped-Missmatch", "Unmapped-Too_Short")
mat[,2:6] <- lapply(mat[,2:6], function(x) as.numeric(as.character(x)))

load(paste0(i,'/DE.rdata'))
meta = colData(dds)[,c(1,7)]
meta = read.csv(paste0(i,'/SraRunTable.txt')) 
metadata = merge.data.frame(meta,mat, by.x = 'Run', by.y='Run' )
meta2$Run = meta2$Experiment
meta2 = meta2 %>% unique
meta = meta[,-c(3,5:13)]
metadata[is.na(metadata$infection),]$hours_infected = 0
metadata$Condition = 'Infected'
metadata[grep('SM', metadata$Run),]$Condition = 'Mock' 
meta$Transduced = 'No'
meta[grep("trasnduced", meta$source_name),]$Transduced = 'Yes'

metadata$treatment = gsub(' IAV', '',metadata$treatment)
metadata$Replicate = metadata$source_name %>% sapply(function(x)strsplit(x,'_rep')[[1]][2]) %>% as.vector()
metadata$source_name = metadata$source_name %>% sapply(function(x)strsplit(x,'_rep')[[1]][1]) %>% as.vector()
metadata$Condition = 'Infected'
metadata[grep("mock", metadata$infection),]$Condition = 'Mock'
metadata[grep('COVID', metadata$disease_state),]$Condition = 'Infected'
metadata[grep('IFNb', metadata$treatment),]$mutant = '1'
meta$source_name = gsub("treated A549 cells trasnduced with a vector expressing human", "", meta$source_name)
metadata$Condition ="Infected"
metadata[grep('MOCK', metadata$Run),]$Condition = 'Mock'
meta$Celltype ='A549'
meta[grep('postmortem', meta$source_name),]$Celltype = 'Biopsy'
meta[grep('NHBE', meta$source_name),]$Celltype = 'NHBE'
colnames(metadata)[2] = c('Condition')

metadata$bamReads = paste0(homedir,"/",i,"/TwoPass/",metadata$Run,"Aligned.sortedByCoord.out.bam")
# Step 1: Count data file construction

files = matrix(unlist(strsplit(as.character(list.files(paste0(homedir,"/",i,"/TwoPass/"),pattern = "ReadsPerGene.out.tab")), "Reads")),ncol=2, byro=T)[,1]
genecounts = "ReadsPerGene.out.tab"
logname = "Log.final.out"
load("/home/conor93/References/gene_annotations_v39.Rdata") # get version 28.Only if not using v29 Gencode
dir = paste(getwd(),i, sep ="/")
filedir = paste(dir, "TwoPass","", sep="/") # Specify file path

#### Step 2: Counting Loop formula: #####

Ns = list()
i = 1
for( n in files ){
  N = list()
  countfile = paste0(filedir,n, genecounts)
  logfile = paste0(filedir,n, logname)
  
  if( file.exists(countfile) ) {
    print(countfile)
    counts =  read.table(countfile)
    
    log1 =read.table(logfile, sep="\t", nrows=6)
    log2 =read.table(logfile, sep="\t", skip=8, nrows=14)
    log3 =read.table(logfile, sep="\t", skip=23, nrows=4)
    log4 =read.table(logfile, sep="\t", skip=28, nrows=3)
    
    N$mapinfo = rbind(log1,log2,log3,log4)
    N$unmapped =  counts[1,]
    N$multimapping = counts[2,]
    N$noFeature =   counts[3,]
    N$ambiguous = counts[4,]
    N$length = dim(counts)[1]-4
    N$genes = counts[ (1:N$length)+4,1]
    N$counts1 = counts[ (1:N$length)+4,2]
    N$counts2 = counts[ (1:N$length)+4,3]
    N$counts3 = counts[ (1:N$length)+4,4]
    print(paste("counts 1:", sum(N$counts1),"counts 2:",sum(N$counts2),"counts 3:", sum(N$counts3)))
    
    if(sum(N$counts1) > sum(N$counts2) & sum(N$counts1) > sum(N$counts3)) {
      print(paste('count1 columns has most reads for',n ))
      ff = N$counts1
    } else if (sum(N$counts2) > sum(N$counts1) & sum(N$counts2) > sum(N$counts3)) {
      print(paste('count2 columns has most reads for',n ))
      ff = N$counts2
    } else {
      print(paste('count3 columns has most reads for',n ))
      ff = N$counts3
    }
    ff = N$counts3
    
  } else {
    N$counts1 = rep(0, length(attr$ensemblID) ) ## Change if appropriate counts are in col 3/4
  }
  if( i > 1  ){
    counts_exp = cbind(counts_exp, ff) ## Change if appropriate counts are in col 3/4
  } else {
    counts_exp = N$counts1 ## Change if appropriate counts are in col 3/4
  }
  Ns[[i]] = N
  print(i)
  i = i + 1
}

# Step 2a: Format Row Names to remove the decimal from ENSG names and Save:
genes = list()
colnames(counts_exp) = as.character(files) 
if( !require("stringr")){
  BiocManager::install("stringr")
}
for (i in 1:length(N$genes)) {
  if(str_detect(N$genes[i],".")){
    genes[i] = matrix(unlist(strsplit( as.character(N$genes[i]), "\\.") ), ncol=2, byro=T)[,1]
  } else {
    genes[i] = as.character(N$genes[i])
  }
}
genes = do.call(rbind.data.frame,genes)
rownames(counts_exp) = as.character(genes[,1])

save(Ns, counts_exp, file = paste0(studies[14],'/countdata.Rdata'))
save(metadata, file = paste0(studies[14],"/metadata.rdata"))

#### Step 3: Differential Expression: ####

countdata = DESeq2::counts(dds)
countdata <- countdata[,as.character(coldata$Run)]

coldata = metadata
countdata = counts_exp[,coldata$Run]

keep <- rowSums(countdata) > 5*ncol(countdata)
countdata <- countdata[keep,]
all(coldata$Run == colnames(countdata)) # Must be "TRUE"

dds <- DESeqDataSetFromMatrix(countData = countdata, 
                              colData = coldata, design = ~Condition)
dds <- DESeq(dds)
res = results(dds, contrast = c('Condition','Infected','Mock'))

save(dds,file = paste0(studies[4],'/DE.rdata'))

#### Plot 1: Volplot ####

setwd(paste0("E:/Studies/COVID"))
studies = list.files()[-14]
reslist = list()
n =1

for (i in studies) {
  load(paste0(i,"/DE.rdata"))
  vsd <- vst(dds, blind = FALSE)
  res = results(dds, contrast = c('Condition','Infected','Mock'))
  reslist[[n]] <- as.data.frame(res)
  reslist[[n]]$Studies <- i
  reslist[[n]]$genes = rownames(reslist[[n]])
  rownames(reslist[[n]]) =NULL
  n =n+1
}
res = do.call(rbind.data.frame, reslist)
res[is.na(res$pvalue),]$pvalue =1
res$Studies = factor(res$Studies, levels = unique(res$Studies)[c(12:14,1:11,15)])
legend =c("Not significant","-1.5 > log2FC > 1.5","p <= 0.05","Significant")

p1 =EnhancedVolcano(res,
                    lab = res$genes,
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    selectLab = unique(res$genes[!grepl("ENSG", res$genes)]),
                    xlab = expression(paste('log'[2],'(FC)')),
                    ylab = expression(paste('-log'[10],'(P)')),
                    #pCutoff = 0.05,
                    FCcutoff = 1.5,
                    pointSize = 0.1,
                    labSize  = 1.5,
                    max.overlaps = 50,
                    drawConnectors = TRUE,
                    endsConnectors = "last",
                    widthConnectors = 0.3,
                    min.segment.length = 0.01,
                    title = NULL,
                    caption = NULL,
                    cutoffLineType = "dashed",
                    subtitle = NULL,
                    col=c('black', 'green3', 'royalblue1', 'coral1'),
                    legendLabels = legend,
                    legendIconSize = 2,
                    colAlpha = 0.5,
                    border = 'full',
                    legendPosition = 'none') +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.text = element_text(size =6, family = "sans"),
        legend.margin = margin(0, 0, 0, 0,unit = "cm"),
        legend.spacing.x = unit("1", "mm"),
        legend.box.spacing =unit("1", "mm"),
        plot.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text = element_text(size =6 ,family = "sans"),
        axis.title = element_text(size =7, family = "sans"),
        strip.text.x = element_text(size = 7, family = "sans")) +
  facet_wrap(vars(Studies), nrow = 3)

#### Plot 2: Study Description table ####

df = data.frame(Studies = factor(studies, levels = rev(studies)))
df$Celltype = c("Primary","Primary","Primary","Primary","Immortalized",
                "Immortalized","Immortalized",'Primary','Immortalized',
                'Immortalized','Primary','Immortalized','Immortalized',
                'Primary','Immortalized')

df$Celltype = c('Immortalized','Primary','Immortalized','Primary','Immortalized',
                'Primary','Primary','Primary','Primary','Primary','Primary','Immortalized',
                'Primary')

p2 =ggplot(data =df, aes(x =Studies, fill = Celltype)) +
  geom_bar(width = 0.8) +
  coord_flip() +
  theme_bw() +
  ylab("Celltype")+
  xlab("Studies")+
  theme(axis.text.y = element_text(size =6, family = "sans"),
        axis.title.y = element_text(size =7, family = "sans"),
        axis.text.x = element_text(colour = 'white', size =6),
        axis.title.x = element_text(size =7, family = "sans"),
        axis.ticks.x = element_line(colour = "white"),
        legend.title = element_text(size = 7, family ='sans', face = 'bold'),
        legend.key.size = unit("2","mm"),
        legend.text = element_text(size =6, family = "sans"),
        plot.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_manual('Celltype Used\nin study', 
                    values = c("Immortalized" = "mediumvioletred",
                               "Primary" = "cornflowerblue"))

#### Plot 3: Barplot of DE totals per study ####

ids = colnames(genesets.up)
dfup = cbind.data.frame(colSums(genesets.up, na.rm = TRUE),rep("Upregulated", length(colSums(genesets.up, na.rm = TRUE))), names(colSums(genesets.up, na.rm = TRUE)))
dfdown = cbind.data.frame(colSums(genesets.down, na.rm = TRUE),rep("Downregulated", length(colSums(genesets.down, na.rm = TRUE))), names(colSums(genesets.down, na.rm = TRUE)))
colnames(dfup) = c('Totals', 'Class','Study')
colnames(dfdown) = c('Totals', 'Class','Study')
df = rbind.data.frame(dfup,dfdown)
df$Totals = as.numeric(df$Totals)
for (i in 1:nrow(df)) {
  if(df[i,]$Class == "Downregulated"){
    df[i,1] = df[i,1]*-1
  }
}
df$Class = factor(df$Class, levels = c("Upregulated", "Downregulated"))

p3 = ggplot(df,aes(y=factor(Study,levels = rev(colnames(genesets.up)[c(12:14,1:11,15)])),x=Totals,fill=Class)) +
  geom_col(width = 0.5)+
  labs(x = "Differentially Expressed genes", y = "Studies") +
  theme_classic() + 
  scale_fill_manual("Class", values = c("Downregulated" = "brown2","Upregulated" = "darkolivegreen3")) +
  geom_text(
    aes(x = Totals+(Totals/abs(Totals)*500),label = abs(Totals),group = Study), size = 2) +
  theme_bw() +
  theme(axis.text = element_text(size =6, family = "sans"),
        axis.title = element_text(size =7, family = "sans"),
        legend.title = element_text(size = 7, family ='sans', face = 'bold'),
        legend.key.size = unit("2","mm"),
        legend.text = element_text(size =6, family = "sans"),
        plot.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

#### Plot 4: Sample totals per study ####

reslist = list()
n =1

for (i in studies) {
  load(paste0(i,"/DE.rdata"))
  df = as.data.frame(colData(dds)) %>% .[,c('Run','Condition')] 
  df = as.data.frame(table(df$Condition))
  colnames(df) = c("Condition","Total")
  df$Studies = i
  reslist[[n]] =df
  n =n+1
}

res2 = do.call(rbind.data.frame, reslist)
res2$Studies = factor(res2$Studies, levels = rev(unique(res2$Studies)[c(12:14,1:11,15)]))
res2$Condition = factor(res2$Condition, levels = rev(unique(res2$Condition)))

p4 = ggplot(data = res2, aes(x = Studies, y = Total, fill =Condition)) +
  geom_bar(position ="dodge", stat = "identity", width = 0.7) +
  labs(x = "Studies", y ="Samples (N)") +
  theme_bw() +
  coord_flip() +
  theme(legend.title = element_text(face = 'bold', size=7, family = "sans"),
        legend.key.size = unit("2", "mm"),
        legend.text = element_text(size =6, family = "sans"),
        plot.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title = element_text(size = 7, family = "sans"), 
        axis.text = element_text(size = 6, family = 'sans')) +
  scale_fill_manual('Condition', values = c("Infected" = "tan2",
                                            "Mock" = "steelblue2"),guide = guide_legend(reverse = T)) 

#### Plot 5: Barplot of Unique reads per study ####

df =list()
n =1
for (i in studies) {
  load(paste0(i,"/DE.rdata"))
  meta = colData(dds)
  df[[n]] = as.data.frame(meta[,c('Condition','Unique_Mapping')])
  df[[n]]$Studies = i
  n = n+1
}
df1 =do.call(rbind.data.frame, df)
df1$Studies =factor(df1$Studies, levels = rev(unique(df1$Studies)[c(12:14,1:11,15)]))

p5 = ggplot(data =df1,aes(x =Studies, y = Unique_Mapping/100, colour =Condition)) +
  stat_summary(fun=mean, 
               fun.max = function(z){quantile(z,0.25)},
               fun.min = function(z){quantile(z,0.75)},
               size=0.1)+
  coord_flip() +
  ylab("Unique-Mapped (%)") +
  theme_bw() +
  theme(axis.title = element_text(size = 7, family = "sans"), 
        axis.text = element_text(size = 6, family = 'sans'),
        legend.position = "none",
        plot.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_manual('Condition', values = c("Infected" = "tan2",
                                              "Mock" = "steelblue2"),guide = guide_legend(reverse = T)) 


legend2 = g_legend(p2 + theme(legend.background = element_blank()))
legend3 = g_legend(p3 + theme(legend.background = element_blank()))
legend4 = g_legend(p4 + theme(legend.background = element_blank()))
blk = plot.new()

legends <- plot_grid(blk,legend2,legend3,blk, ncol = 1, align = 'v', rel_heights = c(1.6,1,1,1.6)) +theme(plot.background = element_blank())
f1b = plot_grid(p2 + theme(legend.position = "hidden"), 
                p3+ theme(legend.position = "hidden",
                          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                legends, ncol = 3, rel_widths = c(0.25,0.5, 0.2))
legends <- plot_grid(blk,legend4,blk, ncol = 1, align = 'v', rel_heights = c(1.6,1,1.6)) +theme(plot.background = element_blank())
f1c = plot_grid(p4 + theme(legend.position = "hidden"), 
                p5 + theme(legend.position = "hidden",
                           axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                legends, ncol = 3, rel_widths = c(0.8,0.4,0.3), align = 'h' )


f1bc = plot_grid(f1b, f1c, rel_widths = c(1,0.85), labels = c(letters[2:3]), label_fontfamily = "Helvetica", label_size = 8)

pdf("E:/ProjectWork/Figures/Figure 1b.pdf", height = 10,width =7, paper = "A4",useDingbats = F)
plot_grid(p1, f1bc,g, nrow = 3, 
          labels = c(letters[1],"", letters[4]), 
          rel_heights = c(1,0.6,1.2),
          label_fontfamily = "Helvetica", 
          label_size = 8 )
dev.off()


#### Step 4: Recurrence: ####

load("E:/gene_annotations_v39.Rdata") # get version 28.Only if not using v29 Gencode
setwd("E:/Studies/IAV")
ids = list.files()
n = length(ids)
attr2 = attr
attr2$Chromosome.scaffold.name = paste0("chr", attr2$Chromosome.scaffold.name)
attr2[is.na(match(attr2$Chromosome.scaffold.name,paste0("chr", c(1:22,"X", "Y", "MT")))),]$Chromosome.scaffold.name = "Scaffold"
attr2 = attr2[attr2$Chromosome.scaffold.name %in% paste0("chr", c(1:22)),]
ngenes = matrix(unlist(strsplit( as.character(unique(attr2$Gene.stable.ID.version )), "\\.") ) , ncol=2, byro=T)[,1]
N = length(ngenes)
genesets.up = matrix(0, ncol=n, nrow= N ) 
genesets.down = matrix(0, ncol=n, nrow= N ) 

# Loop function for recurrence analysis across all data sets.

for( i in 1:n ){
  
  id = ids[i]
  filename = paste0(id,"/DE.rdata") 
  load(filename)
  res = results(dds, contrast = c('Condition','Infected','Mock'))
  
  m = match(ngenes, rownames(res) )  
  f.a = !is.na(m)
  f.r = m[f.a] 
  
  genes.up = (res$log2FoldChange >= 1 & res$padj < 0.05 ) *1
  genes.down = (res$log2FoldChange <= -1 & res$padj < 0.05)*1  
  
  genesets.up[f.a,i] = genes.up[f.r]
  genesets.down[f.a,i] = genes.down[f.r]
}

colnames(genesets.up) = ids
rownames(genesets.up) = ngenes
colnames(genesets.down) = ids
rownames(genesets.down) = ngenes
genesets.up = genesets.up %>% as.data.frame()
genesets.down = genesets.down %>% as.data.frame()
genesets.up[is.na(genesets.up)] = 0
genesets.down[is.na(genesets.down)] = 0

save(genesets.up, genesets.down, file ="E:/ProjectWork/Rfiles/recurrence_genesets.rdata")
#### Plot 7: Recurrence plot ####

df.fdr = as.data.frame(rep(0, 4))
colnames(df.fdr) = 'Virus'
df.fdr$`Total` = rep(0,4)
df.fdr$`Differential Expression` = rep(0,4)
df.fdr$counts = rep(0,4)

res = list()
load("E:/ProjectWork/Rfiles/recurrence_genesets.rdata")
n=1
genesets = genesets.up
recur = rowSums(genesets, na.rm = TRUE)
fdrs = calc_fdrs_recur(genesets)
df.fdr[n,1] = 'Influenza'
df.fdr$`Differential Expression`[n] = 'Upregulated'
df.fdr$`Total`[n] = fdrs$sig
df.fdr$counts[n] = fdrs$Pt
df = data.frame(recur[recur>0])
res[[n]] = data.frame(table(df), row.names = NULL)
colnames(res[[n]]) = c("Recurrence","counts")
res[[n]]$Virus = 'Influenza'
res[[n]]$`Differential Expression` = 'Upregulated'

n=2
genesets = genesets.down
recur = rowSums(genesets, na.rm = TRUE)
fdrs = calc_fdrs_recur(genesets)
df.fdr[n,1] = 'Influenza'
df.fdr$`Differential Expression`[n] = 'Downregulated'
df.fdr$`Total`[n] = fdrs$sig
df.fdr$counts[n] = fdrs$Pt
df = data.frame(recur[recur>0])
res[[n]] = data.frame(table(df), row.names = NULL)
colnames(res[[n]]) = c("Recurrence","counts")
res[[n]]$Virus = 'Influenza'
res[[n]]$`Differential Expression` = 'Downregulated'
res[[n]]$counts = res[[n]]$counts*-1

load("E:/ProjectWork/Rfiles/COVID_recurrence_genesets.rdata")
n=3
genesets = genesets.up
recur = rowSums(genesets, na.rm = TRUE)
fdrs = calc_fdrs_recur(genesets)
df.fdr[n,1] = 'SARS-CoV-2'
df.fdr$`Differential Expression`[n] = 'Upregulated'
df.fdr$`Total`[n] = fdrs$sig
df.fdr$counts[n] = fdrs$Pt
df = data.frame(recur[recur>0])
res[[n]] = data.frame(table(df), row.names = NULL)
colnames(res[[n]]) = c("Recurrence","counts")
res[[n]]$Virus = 'SARS-CoV-2'
res[[n]]$`Differential Expression` = 'Upregulated'

n=4
genesets = genesets.down
recur = rowSums(genesets, na.rm = TRUE)
fdrs = calc_fdrs_recur(genesets)
df.fdr[n,1] = 'SARS-CoV-2'
df.fdr$`Differential Expression`[n] = 'Downregulated'
df.fdr$`Total`[n] = fdrs$sig
df.fdr$counts[n] = fdrs$Pt
df = data.frame(recur[recur>0])
res[[n]] = data.frame(table(df), row.names = NULL)
colnames(res[[n]]) = c("Recurrence","counts")
res[[n]]$Virus = 'SARS-CoV-2'
res[[n]]$`Differential Expression` = 'Downregulated'
res[[n]]$counts = res[[n]]$counts*-1

res1 = do.call(rbind.data.frame, res)
res1$`Differential Expression` = factor(res1$`Differential Expression`, levels = c('Upregulated','Downregulated'))
colnames(df.fdr)[4] = 'Recurrence'
df.fdr$counts = rep(0, 4)
res1$Recurrence = as.numeric(res1$Recurrence)
df.fdr$ymin =c(0,-Inf,0, -Inf)
df.fdr$ymax =c(Inf,0, Inf, 0)
df.fdr$xmax =rep(Inf, 4)
df.fdr$yend = c(Inf, -Inf, Inf, -Inf)

p7 =ggplot(data = res1, aes(x = Recurrence, y = counts,fill = `Differential Expression`)) +  
  geom_rect(data = df.fdr, aes(xmin = Recurrence-.5,
                               xmax = xmax,
                               ymin = ymin, 
                               ymax = ymax),fill="lightblue", show.legend = T) +
  geom_segment(data = df.fdr, aes(x=Recurrence-.5,
                                  y = 0,
                                  xend = Recurrence-.5, 
                                  yend = yend, color = 'FDR <= 0.05'), 
               linetype="dashed",
               size = 0.5) +
  geom_col() +
  geom_text(aes(y = counts+(counts/abs(counts)*400),label = abs(counts),group = Recurrence, angle =90), 
            size = 2) +
  labs(x = "Recurrence", y = "Number of Genes" ) +
  theme_bw() +
  scale_fill_manual('Class', values = c("Downregulated" = "brown2",
                                        "Upregulated" = "darkolivegreen3"),
                    guide = guide_legend(reverse = T)) +
  scale_x_continuous(breaks = 1:max(res1$Recurrence)) +
  facet_grid(~ Virus, scales = 'free', space ="free") +
  scale_color_manual('Recurrence\nThreshold', values = c("FDR <= 0.05" = "red")) +
  theme(legend.direction = 'vertical',
        legend.position = 'right',
        legend.title = element_text(family = 'sans', face = 'bold', size = 7),
        legend.key.size = unit("2","mm"),
        legend.spacing.x = unit("1", "mm"),
        legend.box.spacing =unit("1", "mm"),
        legend.text = element_text(size =6, family = "sans"),
        axis.title = element_text(size = 7, family = "sans"), 
        strip.text.x = element_text(size = 7, family = "sans"), 
        axis.text = element_text(size = 6, family = 'sans'),
        plot.background = element_blank())

save(res1, df.fdr, file ="E:/ProjectWork/Rfiles/RecurrencePlotfiles.Rdata")

#### Plot 8: DE Prior v Recurrence ####

load("E:/DE_prior.rdata")
load("E:/gene_annotations_v39.Rdata")
load("E:/ProjectWork/Rfiles/recurrence_genesets.rdata")
df =list()
genesets = genesets.up
recur = rowSums(genesets, na.rm = TRUE)
gmat = as.data.frame(matrix(0, ncol = 5, nrow = nrow(DE_list)))
colnames(gmat) = c('Recurrence', 'DE prior', 'Name', 'Class', "Virus")
gmat$Name = rownames(DE_list)
m = match(gmat$Name, names(recur))
f.a = !is.na(m)
f.t = m[f.a]
gmat[f.a,1] = recur[f.t]
gmat[,2] = (DE_list[,4])
gmat = as.data.frame(gmat[gmat[,1] >= 1,])
gmat = as.data.frame(gmat)
gmat$Class = 'Upregulated'
gmat$Virus ="Influenza"
df[[1]] =gmat
genesets = genesets.down
recur = rowSums(genesets, na.rm = TRUE)
gmat = as.data.frame(matrix(0, ncol = 5, nrow = nrow(DE_list)))
colnames(gmat) = c('Recurrence', 'DE prior', 'Name', 'Class', "Virus")
gmat$Name = rownames(DE_list)
m = match(gmat$Name, names(recur))
f.a = !is.na(m)
f.t = m[f.a]
gmat[f.a,1] = recur[f.t]
gmat[,2] = (DE_list[,4])
gmat = as.data.frame(gmat[gmat[,1] >= 1,])
gmat = as.data.frame(gmat)
gmat$Class = 'Downregulated'
gmat$Virus ="Influenza"
df[[2]] = gmat

load("E:/ProjectWork/Rfiles/COVID_recurrence_genesets.rdata")
genesets = genesets.up
recur = rowSums(genesets, na.rm = TRUE)
gmat = as.data.frame(matrix(0, ncol = 5, nrow = nrow(DE_list)))
colnames(gmat) = c('Recurrence', 'DE prior', 'Name', 'Class', "Virus")
gmat$Name = rownames(DE_list)
m = match(gmat$Name, names(recur))
f.a = !is.na(m)
f.t = m[f.a]
gmat[f.a,1] = recur[f.t]
gmat[,2] = (DE_list[,4])
gmat = as.data.frame(gmat[gmat[,1] >= 1,])
gmat = as.data.frame(gmat)
gmat$Class = 'Upregulated'
gmat$Virus ="SARS-CoV-2"
df[[3]] =gmat
genesets = genesets.down
recur = rowSums(genesets, na.rm = TRUE)
gmat = as.data.frame(matrix(0, ncol = 5, nrow = nrow(DE_list)))
colnames(gmat) = c('Recurrence', 'DE prior', 'Name', 'Class', "Virus")
gmat$Name = rownames(DE_list)
m = match(gmat$Name, names(recur))
f.a = !is.na(m)
f.t = m[f.a]
gmat[f.a,1] = recur[f.t]
gmat[,2] = (DE_list[,4])
gmat = as.data.frame(gmat[gmat[,1] >= 1,])
gmat = as.data.frame(gmat)
gmat$Class = 'Downregulated'
gmat$Virus ="SARS-CoV-2"
df[[4]] = gmat

df = do.call(rbind.data.frame, df)
df$Class = factor(df$Class, levels = c("Upregulated","Downregulated"))
p8 <- ggplot() + 
  geom_boxplot(data = df, aes(x= factor(Recurrence, levels = 0:max(Recurrence)), 
                              y = `DE prior`, fill =Class), 
               width = 0.7, outlier.shape = NA,position = position_dodge(preserve = "single")) +
  xlab("Recurrence") +
  scale_fill_manual(values = c("Downregulated" = "brown2","Upregulated" = "darkolivegreen3")) +
  facet_grid(~ Virus, scales = 'free_x', space = "free") +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.text = element_text(size = 6, family = "sans"),
        legend.title = element_text(size = 7, family = "sans", face = 'bold'),
        legend.key.size = unit("3","mm"),
        legend.spacing.x = unit("1", "mm"),
        legend.box.spacing =unit("1", "mm"),
        axis.text = element_text(size = 6, family = "sans"),
        axis.title = element_text(size = 7, family = "sans"),
        strip.text.x = element_text(size = 7, family = "sans"))

#### Plot 9: Normalized expression recurrent genes ####

setwd("E:/Studies/IAV")
studies = list.files()
load("E:/ProjectWork/Rfiles/recurrence_genesets.rdata")
#or 
setwd("E:/Studies/COVID")
studies = list.files()[-14]
load("E:/ProjectWork/Rfiles/COVID_recurrence_genesets.rdata")

reslist = list()
genes = rownames(genesets.up[rowSums(genesets.up) > 2,])
n =1
test =list()
for (i in studies) {
  load(paste0(i,"/DE.rdata"))
  m = match(genes, rownames(dds))
  f.a =!is.na(m)
  genes = genes[f.a]
  tcounts <- t(log2((counts(dds[genes, ], normalized=TRUE, replaced=FALSE)+1))) %>%
    merge(colData(dds), ., by="row.names") %>%
    tidyr::gather(gene, expression, (ncol(.)-length(genes)+1):ncol(.))
  tcounts = tcounts %>% 
    dplyr::select(Row.names, Condition, gene, expression)
  colnames(tcounts) = c("SampleID",'Condition',"Gene","expression")
  tcounts$expression = as.numeric(tcounts$expression)
  tcounts$Studies = i
  test[[n]] =tcounts
  n = n+1
}
test = do.call(rbind.data.frame, test)
test$Condition = factor(test$Condition, levels = c("Mock","Infected"))
test$Studies = factor(test$Studies, levels = unique(test$Studies)[c(12:14,1:11,15)])
test$Class = 'Upregulated'
reslist[[1]] = test
genes = rownames(genesets.down[rowSums(genesets.down) > 3,])
n =1
test =list()
for (i in studies) {
  load(paste0(i,"/DE.rdata"))
  m = match(genes, rownames(dds))
  f.a =!is.na(m)
  genes = genes[f.a]
  tcounts <- t(log2((counts(dds[genes, ], normalized=TRUE, replaced=FALSE)+1))) %>%
    merge(colData(dds), ., by="row.names") %>%
    tidyr::gather(gene, expression, (ncol(.)-length(genes)+1):ncol(.))
  tcounts = tcounts %>% 
    dplyr::select(Row.names, Condition, gene, expression)
  colnames(tcounts) = c("SampleID",'Condition',"Gene","expression")
  tcounts$expression = as.numeric(tcounts$expression)
  tcounts$Studies = i
  test[[n]] =tcounts
  n = n+1
}
test = do.call(rbind.data.frame, test)
test$Condition = factor(test$Condition, levels = c("Mock","Infected"))
test$Studies = factor(test$Studies, levels = unique(test$Studies)[c(12:14,1:11,15)])
test$Class = 'Downregulated'
reslist[[2]] = test

reslist = do.call(rbind.data.frame, reslist)
reslist$Class = factor(reslist$Class, levels = c("Upregulated","Downregulated"))

p11 = ggplot(reslist, aes(Condition, expression,fill = Condition)) + 
  geom_violin(color =NA) + 
  geom_boxplot(width = 0.2,fill ="white", outlier.shape = NA) +
  labs(x="Conditions", 
       y="Normalized Gene Expression") +
  facet_grid(cols = vars(Class), rows = vars(Studies),switch = "y") +
  theme_bw()+
  coord_flip() +
  geom_signif(map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "ns" = 1) , textsize=4, 
              comparisons = list(c("Mock", "Infected")),
              test = "wilcox.test",
              margin_top = -0.1,
              vjust = 0.2,hjust =0.4 ) +
  scale_fill_manual(values = c("Infected" = "tan2",
                               "Mock" = "steelblue2"), 
                    name="Conditions", guide = guide_legend(reverse = T)) +
  theme(legend.text =element_text(family = 'sans', size =6),
        legend.key.size = unit("2","mm"),
        legend.title = element_text(family = 'sans', size = 7, face ="bold"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(family = 'sans', size = 7),
        axis.text.x = element_text(family = 'sans', size = 6),
        strip.text.x = element_text(family = 'sans', size = 7),
        strip.text.y.left = element_text(family = 'sans', size = 7, angle = 0),
        strip.background = element_blank())

blk = plot.new()
legend = plot_grid(blk, legend1, blk, ncol = 1, rel_heights = c(0.5,1,0.5))

f2a = plot_grid(p7+theme(legend.position = 'bottom', legend.direction = 'vertical'),
                p8 +theme(legend.position = 'bottom', legend.direction = 'vertical'),
                p9 +theme(legend.position = 'bottom'), rel_heights = c(1,1,1.3), 
                labels = c(letters[1:3]), 
                label_fontfamily = "Helvetica",
                ncol =1,nrow =3, label_size = 8 )
f2b = plot_grid(p10 , 
                p11 ,
                labels = c(letters[4:5]), rel_width = 4,               
                label_fontfamily = "Helvetica", 
                nrow = 2,ncol =1,label_size = 8)

pdf("E:/ProjectWork/Figures/Figure 2.pdf", height = 10,width =7, paper = "A4",useDingbats = F)
plot_grid(f2a, f2b, ncol = 2,nrow = 1,rel_width = c(1,0.4),
          label_fontfamily = "Helvetica", 
          label_size = 8, vjust = 0.5 )
dev.off()

#### Alternative Splicing ####

#### StringTie (For transcript assembly and quantification) ####

#Generate GTF files per sample
sortedBAM <- list.files("TwoPass",pattern = ".bam$", recursive = T, full.names = T)
#sortedBAM <- paste0('Alignments/',myDesign$sampleID,"Aligned.sortedByCoord.out.bam")

threads <- 10
system(paste0('rm -r StringTie'))
if(!dir.exists("StringTie")){system(paste0("mkdir StringTie"))}
ref="/home/hlchen/Conor/references/human_wsn.gtf" #Or /home/hlchen/Conor/references/human_covid-19.gtf
stringtie <- "/home/conor93/stringtie-2.2.1.Linux_x86_64/stringtie"

for (i in 1:length(sortedBAM)) {
  label <- gsub("Aligned.sortedByCoord.out.bam","", basename(sortedBAM[i]))
  if(!dir.exists(paste0("StringTie/",label))) system(paste("mkdir", paste0("StringTie/",label)))
  params <- paste("-G",ref,
                  "-l", label,
                  "-p", threads,
                  "-o", paste0("StringTie/",label,"/", label,".gtf"))
  cmd <- paste(stringtie, params, sortedBAM[i])
  system(cmd)                
}

temp <- as.data.frame(list.files("StringTie",pattern = '.gtf$', recursive = T, full.names = T))[-c(1,2),]

#### Or for split studies ####
#temp <- as.data.frame(paste0(myDesign$batch,"/StringTie/",myDesign$sampleID, "/",myDesign$sampleID,".gtf"))

write.table(temp, row.names = F, col.names = F, quote = F, file ="StringTie/gtflist.txt")
cmd <- paste(stringtie, "--merge",
             "-G",ref,
             "-p",threads,
             "-o","StringTie/merged.gtf",
             "-i StringTie/gtflist.txt")
system(cmd)
system(paste("rm StringTie/gtflist.txt"))

#Evaluate merged GTF

gffcompare ="/software/gffcompare/0.10.6/gffcompare"

#cmd <- paste(gffcompare, "-r /home/hlchen/Conor/references/human_covid-19.gtf -T",
#            "-o StringTie/merged StringTie/merged.gtf")
#Or
cmd <- paste(gffcompare, "-r /home/hlchen/Conor/references/human_wsn.gtf -T",
             "-o StringTie/merged StringTie/merged.gtf")
system(cmd)

#Re-quantification step
for (i in 1:length(sortedBAM)) {
  stringtie <- "/home/conor93/stringtie-2.2.1.Linux_x86_64/stringtie"
  label <- gsub("Aligned.sortedByCoord.out.bam","", basename(sortedBAM[i]))
  if(!dir.exists(paste0("StringTie/",label))) system(paste("mkdir", paste0("StringTie/",label)))
  params <- paste("-eB -G StringTie/merged.annotated.gtf",
                  "-l", label,
                  "-p", threads,
                  "-o", paste0("StringTie/",label,"/", label,".gtf"))
  cmd <- paste(stringtie, params, sortedBAM[i])
  system(cmd)   
}

#### Fix Issues with NA gene (Influenza studies only)####

files = list.files("StringTie", pattern ="t_data.ctab$", full.names = T, recursive = T)

for (i in 1:length(files)) {
  ff = read.table(files[i], header = T)
  f.a = is.na(ff$chr)
  ff[f.a,]$chr = "Nea"
  ff[f.a,]$t_name = "Nea"
  ff[f.a,]$gene_name = "Nea"
  write.table(ff,col.names = T,sep ="\t", row.names = F,quote = F, file = files[i])
  print(files[i])
}

file = read.delim("StringTie/merged.annotated.gtf", header = F)
f.a =is.na(file[,1])
file[f.a,1] = 'Nea'
file[f.a,9] = gsub("NA","Nea",file[f.a,9])
write.table(file,col.names =F,sep ="\t",quote = F, row.names = F, file = "StringTie/merged.annotated.gtf")

#### IsoformSwitchAnalyzeR ####

library(DESeq2)
library(dplyr)
library(IsoformSwitchAnalyzeR)
library(BSgenome.Hsapiens.UCSC.hg38)

load("DE.rdata")
metadata <- colData(dds)
myDesign <- data.frame(sampleID = metadata$Run, condition = metadata$Condition)

total = as.numeric(unlist(lapply(list.files("TwoPass", full.names = T,pattern = "Log.final.out$"),function(x) read.delim(x, header =F)[6,2])))

samples = paste0("StringTie/",myDesign$sampleID, "/t_data.ctab")
names(samples) = myDesign$sampleID
quant <- importIsoformExpression(sampleVector  = samples, 
                                 readLength = 150,
                                 addIsofomIdAsColumn =TRUE)


aSwitchList <- importRdata(isoformCountMatrix   = quant$counts,
                           isoformRepExpression = quant$abundance,
                           designMatrix         = myDesign,
                           isoformExonAnnoation = "StringTie/merged.annotated.gtf",
                           removeNonConvensionalChr = TRUE,
                           showProgress = TRUE)

#ORF Analysis:

aSwitchList <- addORFfromGTF(aSwitchList, 
                             pathToGTF = "/home/hlchen/Conor/references/human_covid-19.gtf",
                             removeNonConvensionalChr = TRUE,
                             overwriteExistingORF=TRUE)
aSwitchList <- analyzeORF(aSwitchList, 
                          genomeObject = BSgenome.Hsapiens.UCSC.hg38)

filtSwitchList <- preFilter(aSwitchList,
                            reduceToSwitchingGenes = FALSE)

#Isoform Switching:

SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = filtSwitchList,
  reduceToSwitchingGenes=FALSE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  showProgress = TRUE)

#Get Sequences:
library(IsoformSwitchAnalyzeR)
system('rm -r Sequences')
dir.create("Sequences")
SwitchListAnalyzed <- extractSequence(
  SwitchListAnalyzed, 
  writeToFile=T, # to avoid output when running this example data
  alsoSplitFastaFile = TRUE,
  pathToOutput ="Sequences/"
)
save(aSwitchList, filtSwitchList, SwitchListAnalyzed, file ="SwitchLists.rdata")

#Attach CPAT results:
#format CPAT output
library(IsoformSwitchAnalyzeR)
load("SwitchLists.rdata")

cpout <- read.delim("CPAT/iso.ORF_prob.tsv")
cpout <- cpout[,c(1,2,7,8,9,10)]
colnames(cpout) <- c("Sequence.Name", "RNA.size", "ORF.size", "Ficket.Score", "Hexamer.Score", "Coding.Probability")
cpout <- cbind(Data.ID =cpout$Sequence.Name,cpout)
cpout$Sequence.Name <- as.character(matrix(unlist(strsplit(cpout$Sequence.Name,'_ORF')), ncol = 2, byrow = T)[,1])
cpout$Coding.Label <- ifelse(cpout$Coding.Probability >= 0.364,"yes", "no")
write.table(cpout, row.names = F, col.names =T, quote = F, sep ="\t", file ="CPAT/iso_results.txt")

SwitchListAnalyzed <- analyzeCPAT(SwitchListAnalyzed,
                                  pathToCPATresultFile = "CPAT/iso_results.txt",
                                  codingCutoff = 0.364,
                                  removeNoncodinORFs = FALSE)

#Attach PFAM results:

SwitchListAnalyzed <- analyzePFAM(SwitchListAnalyzed,
                                  pathToPFAMresultFile = as.character(list.files("Pfam", full.names = T)))

#Attach IUPred2A results:

SwitchListAnalyzed <- analyzeIUPred2A(SwitchListAnalyzed,
                                      pathToIUPred2AresultFile = as.character(list.files("IUPred2A",full.names = T)))

#Attach SignalP results:

SwitchListAnalyzed <- analyzeSignalP(SwitchListAnalyzed,
                                     pathToSignalPresultFile = as.character(list.files("SignalP",full.names = T)))

#Predicting Alternative Splicing:

SwitchListAnalyzed <- analyzeAlternativeSplicing(SwitchListAnalyzed, 
                                                 onlySwitchingGenes = F)

#Define consequential switches

SwitchListAnalyzed <- analyzeSwitchConsequences(
  SwitchListAnalyzed,
  consequencesToAnalyze=c('tss','tts','last_exon', 'isoform_seq_similarity','isoform_length',
                          'exon_number','intron_structure','intron_retention','coding_potential',
                          'ORF_seq_similarity','ORF_genomic','ORF_length','5_utr_seq_similarity',
                          '5_utr_length','3_utr_seq_similarity','3_utr_length','NMD_status',
                          'domains_identified','domain_length','IDR_identified','IDR_length','IDR_type',
                          'signal_peptide_identified'),
  showProgress=TRUE)

save(aSwitchList, filtSwitchList, SwitchListAnalyzed, file ="SwitchLists.rdata")

#### Post Analysis: ####

attr <- read.delim("E:/ensembl_v105.txt")
load("E:/ProjectWork/Rfiles/Influenza_clusters.rdata")
load("E:/ProjectWork/Rfiles/COVID_clusters.rdata")

##Edit gen_ids of unknown matched transcripts to ensure gene ids are correct

#Extract coordinates of novels/know transcripts as BED file
studies = list.files()
library(dplyr)
homedir =getwd()
for (i in studies) {
  setwd(paste0(homedir,"/",i))
  gtf = read.delim("StringTie/merged.annotated.gtf", header =F)
  gtf = gtf[gtf$V3 == "transcript",]
  gtf$class_code = gtf$V9 %>% sapply(function(x)strsplit(x,'class_code ')[[1]][2]) %>% as.vector()
  gtf$class_code = gtf$class_code %>% sapply(function(x)strsplit(x,';')[[1]][1]) %>% as.vector()
  gtf$transcript_id <- gtf$V9 %>% sapply(function(x)strsplit(x,'; gene_id')[[1]][1]) %>% as.vector()
  gtf$transcript_id = gsub("transcript_id ","",gtf$transcript_id)
  gtf$ref_gene = '.'
  
  filt =grepl("ref_gene_id",gtf$V9) 
  gtf[filt,]$ref_gene = gtf[filt,]$V9 %>% sapply(function(x) strsplit(x,'ref_gene_id ')[[1]][2]) %>% as.character()
  gtf[filt,]$ref_gene = gtf[filt,]$ref_gene %>% sapply(function(x) strsplit(x,'[.;]')[[1]][1]) %>% as.character()
  filt2 = (!grepl("ref_gene_id",gtf$V9) & grepl("gene_name",gtf$V9))  
  gtf[filt2,]$ref_gene = gtf[filt2,]$V9 %>% sapply(function(x) strsplit(x,'gene_name ')[[1]][2]) %>% as.character()
  gtf[filt2,]$ref_gene = gtf[filt2,]$ref_gene %>% sapply(function(x) strsplit(x,';')[[1]][1]) %>% as.character()
  
  bed = gtf[,c(1,4,5,11,6,7,10,12)]
  write.table(bed, col.names = F, row.names = F, quote = F, sep ="\t", file ="StringTie/merged.annotated.bed")
}

#Ensure gene ids are correct, convert gene names to ensembl ids
setwd("E:/Studies/COVID")
studies <- list.files()[-14]

for (i in studies) {
  bed = read.delim(paste0(i,"/merged.annotated.bed"), header =F, stringsAsFactors = F)
  m = match(bed$V8, attr$Gene.name)
  f.a =!is.na(m)
  f.t =m[f.a]
  bed[f.a,]$V8 = attr[f.t,]$Gene.stable.ID
  load(paste0(i,"/SwitchLists.rdata"))
  
  m = match(SwitchListAnalyzed$isoformFeatures$isoform_id, bed$V4)
  f.a =!is.na(m)
  f.t =m[f.a]
  SwitchListAnalyzed$isoformFeatures$gene_id[f.a] = bed$V8[f.t]
  
  #Edit gen_ids of known matched transcripts to ensure gene ids are correct
  
  #Change ids in isofeatures table:
  
  m = match(SwitchListAnalyzed$isoformFeatures$isoform_id, attr$Transcript.stable.ID.version)
  f.a =!is.na(m)
  f.t =m[f.a]
  SwitchListAnalyzed$isoformFeatures$gene_id[f.a] = attr$Gene.stable.ID[f.t]
  m = match(SwitchListAnalyzed$isoformFeatures$gene_id, attr$Gene.stable.ID)
  f.a =!is.na(m)
  f.t = m[f.a]
  SwitchListAnalyzed$isoformFeatures$gene_name[f.a] = attr$Gene.name[f.t]
  SwitchListAnalyzed$isoformFeatures$gene_biotype[f.a] = attr$Gene.type[f.t]
  m = match(SwitchListAnalyzed$isoformFeatures$isoform_id, attr$Transcript.stable.ID.version)
  f.a =!is.na(m)
  f.t =m[f.a]
  SwitchListAnalyzed$isoformFeatures$iso_biotype[f.a] = attr$Transcript.type[f.t]
  save(aSwitchList, filtSwitchList, SwitchListAnalyzed, file =paste0(i,"/SwitchLists_annotated.rdata"))
}

#### Extract Isoforms across studies as dataframe ####

setwd("E:/Studies/IAV")
n =1
isos = list()
studies = list.files()[-14]
for (k in studies) {
  load(paste0(k,"/SwitchLists_annotated.rdata"))
  totals =list()
  isos[[n]] = as.data.frame(SwitchListAnalyzed$isoformFeatures)
  isos[[n]]$Studies = k
  load(paste0(k,"/DE.rdata"))
  res = results(dds, contrast = c('Condition','Infected','Mock'))
  m = match(isos[[n]]$gene_id, rownames(res))
  f.a = !is.na(m)
  f.t = m[f.a]
  isos[[n]][f.a,]$gene_log2_fold_change = res[f.t,]$log2FoldChange
  isos[[n]][f.a,]$gene_q_value = res[f.t,]$padj
  
  ##Assign class codes to isos
  bed = read.delim(paste0(k,"/merged.annotated.bed"), header =F, stringsAsFactors = F)
  m = match(bed$V8, attr$Gene.name)
  f.a =!is.na(m)
  f.t =m[f.a]
  bed$V8[f.a] = attr$Gene.stable.ID[f.t]
  
  m = match(isos[[n]]$isoform_id, bed$V4)
  f.a =!is.na(m)
  f.t =m[f.a]
  isos[[n]]$`Class Code` = 0
  isos[[n]][f.a,]$`Class Code` =bed[f.t,]$V7
  n =n+1
}

isos = do.call(rbind.data.frame, isos)
m = match(isos$gene_id, clusters.up$V1)
f.a = !is.na(m)
f.t =m[f.a]
isos$Clusters = 'Not recurrent'
isos[f.a,]$Clusters = clusters.up[f.t,]$unmergedLabels.mod
#isos[is.na(isos$isoform_switch_q_value),]$isoform_switch_q_value = 1
isos$Studies = factor(isos$Studies, levels = unique(isos$Studies)[c(12:14,1:11,15)])
isos$Annotated ="Not Annotated"
isos$Annotated[isos$`Class Code` %in% c("c","=")] = "Annotated"
isos$`Significant\nIsoform\nSwitching` = "Not significant"

isos[is.na(isos$gene_q_value),]$gene_q_value = 1
save(isos, file ="E:/ProjectWork/Rfiles/SwitchPlot_isoformFeatures.rdata")
#### Filtering for significant isoforms ####

DEfilt = abs(isos$gene_log2_fold_change) >=1 & isos$gene_q_value <= 0.05
isos_sig = isos[DEfilt,]
isos_sig = isos_sig[!is.na(isos_sig$Clusters),]
isos_sig$Ranking =0
filt1 = isos_sig$Clusters != "Not recurrent"
isos_sig$Recurrence = 'Not recurrent'
isos_sig[filt1,]$Recurrence ="Recurrent"
m = match(isos_sig$gene_id, clusters.down$V1)
f.a =!is.na(m)
isos_sig$Recurrence[f.a] = "Recurrent"

isos_sig$Class = "Upregulated"
isos_sig[isos_sig$gene_log2_fold_change < 0,]$Class = 'Downregulated'
isos_sig$iso_log2_fold_change = -isos_sig$iso_log2_fold_change
isos_sig$dIF =-isos_sig$dIF

#### Plot gene totals with at least 1 significant isoform across studies ####
isos_sig$Clusters = 'Not recurrent'
m = match(isos_sig$gene_id, clusters.up$V1)
f.a = !is.na(m)
f.t =m[f.a]
isos_sig[f.a,]$Clusters = paste0(clusters.up[f.t,]$unmergedLabels.mod,"-Up")
m = match(isos_sig$gene_id, clusters.down$V1)
f.a = !is.na(m)
f.t =m[f.a]
isos_sig[f.a,]$Clusters = paste0(clusters.down[f.t,]$unmergedLabels.mod,"-Down")
df =isos_sig[isos_sig$isoform_switch_q_value <= 0.05,c("Studies","Clusters","gene_id", "Annotated")] %>% group_by(Studies, Clusters, Annotated) %>% dplyr::summarise(Total =n())
df$Studies =factor(df$Studies, levels = studies[c(12:14,1:11,15)])
df = df[df$Clusters != "Not recurrent",]
df$Clusters = factor(df$Clusters, levels = rev(c(paste0(1:4,"-Up"),paste0(1:2,"-Down"))))

ggplot(data = df,aes(x =Clusters,y =Total,fill =Annotated)) +
  geom_col(position  = "stack")+
  labs(x = "Totals of significant isoform switches (FDR <= 0.05)", y = "Studies") +
  scale_fill_manual(values = c("royalblue","red2")) +
  theme_bw() +
  coord_flip() +
  facet_wrap(vars(Studies)) +
  theme(axis.text = element_text(size =6, family = "sans"),
        axis.title = element_text(size =7, family = "sans"),
        legend.spacing.x = unit("1", "mm"),
        legend.box.spacing =unit("1", "mm"),
        legend.title = element_text(size = 7, family ='sans', face = 'bold'),
        legend.key.size = unit(2,"mm"),
        legend.text = element_text(size =6, family = "sans"),
        plot.background = element_blank())

#### Plot gene totals with at least 1 significant isoform across clusters ####

m = match(isos_sig$gene_id, clusters.up$V1)
f.a = !is.na(m)
f.t =m[f.a]
isos_sig$Clusters = 'Not recurrent'
isos_sig[f.a,]$Clusters = paste0(clusters.up[f.t,]$unmergedLabels.mod,"-Up")
df1 =unique(isos_sig[isos_sig$isoform_switch_q_value <= 0.05 & isos_sig$Clusters != "Not recurrent",c("Studies","gene_id", "Clusters")]) %>% group_by(Studies, Clusters) %>% dplyr::summarise(Total =n())
df1$Class = "Upregulated"
m = match(isos_sig$gene_id, clusters.down$V1)
f.a = !is.na(m)
f.t =m[f.a]
isos_sig$Clusters = 'Not recurrent'
isos_sig[f.a,]$Clusters = paste0(clusters.down[f.t,]$unmergedLabels.mod,"-Down")
df2 =unique(isos_sig[isos_sig$isoform_switch_q_value <= 0.05 & isos_sig$Clusters != "Not recurrent",c("Studies","gene_id", "Clusters")]) %>% group_by(Studies, Clusters) %>% dplyr::summarise(Total =n())
df2$Class = "Downregulated"
df1 = rbind.data.frame(df1, df2)

df1$Class = factor(df1$Class, levels =c("Upregulated", "Downregulated"))
df1$Studies =factor(df1$Studies, levels =rev( studies[c(12:14,1:11,15)]))
df1$Clusters = factor(df1$Clusters, levels = c(paste0(1:4,"-Up"),paste0(1:2,"-Down")))
p13 = ggplot(data = df1,aes(y =Studies,x =Total,fill =Clusters, group =Class)) +
  geom_col(position  = "stack")+
  labs(x = "Differentially Expressed genes with at least one\nsignificant isoform switch (FDR <= 0.05)", y = "Studies") +
  scale_fill_manual(values = brewer.pal(6, "Set3")) +
  theme_bw() +
  theme(axis.text = element_text(size =6, family = "sans"),
        axis.title = element_text(size =7, family = "sans"),
        legend.spacing.x = unit("1", "mm"),
        legend.box.spacing =unit("1", "mm"),
        legend.title = element_text(size = 7, family ='sans', face = 'bold'),
        legend.key.size = unit(2,"mm"),
        legend.text = element_text(size =6, family = "sans"),
        plot.background = element_blank()) +
  facet_grid(~Class)


#### Plot dIF vs Gene logFC to determine if differential expression and significant isoform switching are mutually exclusive - No ####
isos[!DEfilt,]$Annotated = "Not DE"
p14 =ggplot(data=isos_sig, aes(x=gene_log2_fold_change, y=-dIF,color=Annotated)) +
  geom_point(size =0.8) +
  facet_wrap(vars(Studies),nrow =3) +
  geom_hline(yintercept = c(-.1,.1), linetype='dashed') +
  geom_vline(xintercept = c(-1,1), linetype='dashed') +
  scale_color_manual('Annotated', values = c("royalblue","red2", 'grey')) +
  labs(x='Gene log2 fold change', y=expression(paste(Delta,'IF'))) +
  theme_bw() +
  theme(axis.text = element_text(size =6, family = "sans"),
        axis.title = element_text(size =7, family = "sans"),
        legend.title = element_text(size = 7, family ='sans', face = 'bold'),
        legend.key.size = unit("2","mm"),
        legend.text = element_text(size =6, family = "sans"),
        plot.background = element_blank(),
        strip.text = element_text(family = 'sans', size = 7))

ggplot(data=isos_sig[isos_sig$isoform_switch_q_value <= 0.05,], aes(x=Annotated, y=dIF,fill=Annotated)) +
  geom_violin(color =NA) + 
  geom_boxplot(width = 0.2,fill ="white", outlier.shape = NA) +
  labs(x="Annotated", 
       y=expression(paste(Delta,'IF'))) +
  theme_bw()+
  coord_flip()+
  geom_signif(map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "ns" = 1) , textsize=4, 
              comparisons = list(c("Annotated", "Not Annotated")),
              test = "t.test",
              vjust = 0, margin_top = -0.02, tip_length = 0.01) +
  scale_fill_manual(values = c("Annotated" = "red2",
                               "Not Annotated" = "royalblue", "Not DE" = "grey"), 
                    name="Annotated", guide = guide_legend(reverse = T)) +
  facet_grid(rows = vars(Studies),switch = "y") +
  theme(legend.text =element_text(family = 'sans', size =6),
        legend.key.size = unit("2","mm"),
        legend.title = element_text(family = 'sans', size = 7, face ="bold"),
        axis.text.x = element_text(family = 'sans', size = 6),
        axis.title.x = element_text(family = 'sans', size = 7),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(family = 'sans', size = 7, angle = 0),
        strip.background = element_blank())

#### Plot distributuion of dIF across studies ####

p15 <- ggplot(data =isos_sig,aes(y =-log10(isoform_switch_q_value), x =dIF, color =`Recurrence`)) +
  geom_point(size=0.8) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  geom_hline(yintercept = 1.3, linetype='dashed') + # default cutoff
  facet_wrap(vars(Studies),nrow =3, scales = 'free_y') +
  labs(x=expression(paste(Delta,'IF')), y =expression(paste('-log'[10],'(Isoform Switch Q Value)'))) +
  theme_bw() +
  scale_color_manual(values = c('Recurrent' = "brown2", 'Not recurrent' ="royalblue"))+
  theme(axis.text = element_text(size =6, family = "sans"),
        axis.title = element_text(size =7, family = "sans"),
        legend.title = element_text(size = 7, family ='sans', face = 'bold'),
        legend.key.size = unit("2","mm"),
        legend.text = element_text(size =6, family = "sans"),
        plot.background = element_blank(),
        strip.text = element_text(family = 'sans', size = 7))


#### Make dataframe of splice events/consequences across studies ####

n=1
splicelist = list()
for (k in studies) {
  load(paste0(k,"/SwitchLists_annotated.rdata"))
  tryCatch({ 
    totals =
      extractSplicingEnrichment(
        SwitchListAnalyzed,
        alpha = 0.05,
        dIFcutoff = 0,
        plot =F,
        returnResult =T, 
        returnSummary = F)
  }, error = function(e) {})
  totals$Studies = k
  splicelist[[n]] = totals
  n = n+1
}  

splicelist =do.call(rbind.data.frame, splicelist)

n=1
conseqlist = list()
for (k in studies) {
  load(paste0(k,"/SwitchLists_annotated.rdata"))
  tryCatch({ 
    totals =
      extractConsequenceEnrichment(
        SwitchListAnalyzed,
        alpha = 0.05,
        dIFcutoff = 0,
        plot =F,
        returnResult =T, 
        returnSummary = F)
  }, error = function(e) {})
  totals$Studies = k
  bed = read.delim(paste0(k,"/merged.annotated.bed"), header =F, stringsAsFactors = F)
  m = match(bed$V8, attr$Gene.name)
  f.a =!is.na(m)
  f.t =m[f.a]
  bed[f.a,]$V8 = attr[f.t,]$Gene.stable.ID.version
  m = match(totals$gene_id, bed$V4)
  f.a =!is.na(m)
  f.t =m[f.a]
  totals[f.a,]$gene_id = bed[f.t,]$V8
  conseqlist[[n]] = totals
  n = n+1
}  

conseqlist =do.call(rbind.data.frame, conseqlist)
m = match(conseqlist$isoformUpregulated, attr$Transcript.stable.ID.version)
f.a =!is.na(m)
f.t =m[f.a]
conseqlist[f.a,]$gene_id = attr[f.t,]$Gene.stable.ID
save(splicelist,conseqlist, file ="E:/ProjectWork/Rfiles/Splicelist_COVID.rdata")

#### Plot total distribution of AS events across clusters ####

clusters = clusters.up

splicelist$Clusters = "Not recurrent"
m = match(splicelist$gene_id, clusters$V1)
f.a =!is.na(m)
f.t =m[f.a]
splicelist[f.a,]$Clusters = clusters[f.t,]$unmergedLabels.mod

# filter out isoforms without AS changes
df1 = splicelist[splicelist$Clusters != "Not recurrent" & splicelist$ASchange != 'No change',] %>% group_by(Studies, Clusters, AStype, ASchange) %>% summarise(Total =n()) %>% 
  mutate(Props = Total/sum(Total))
df1$id = paste(df1$Clusters,df1$AStype, df1$ASchange, sep ="_")

mat = matrix(0, nrow = length(studies), ncol = length(unique(df1$id)))  
rownames(mat) = studies
colnames(mat) = unique(df1$id)
for (i in studies) {
  test = df1[df1$Studies == i,]
  mat[i,test$id] = test$Total
}

cc =as.data.frame(df1[,c("Clusters","AStype","ASchange","id"),drop =F]) %>% unique()
rownames(cc) =cc$id
cc =cc[,-4, drop =F]
cc_col = unique(clusters[,c(3)])
names(cc_col) = unique(clusters$unmergedLabels.mod)
cc_col2 = brewer.pal(n = length(unique(cc$AStype)), name = "Dark2")
names(cc_col2) = unique(cc$AStype)
cc_col3 = brewer.pal(n = length(unique(cc$ASchange)), name = "Dark2")[-3]
names(cc_col3) = unique(cc$ASchange)
ha =columnAnnotation(df = cc, 
                     Totals = anno_barplot(colSums(mat), axis_param = list(gp =gpar(fontsize = 6))),
                     col = list(Clusters = cc_col,                                         
                                ASchange =cc_col3,
                                AStype =cc_col2 
                     ),
                     annotation_height =unit(c(10,3,3,3), "mm"),
                     annotation_name_gp = gpar(fontsize = 6),
                     annotation_legend_param = list(
                       labels_gp = gpar(fontsize = 6),
                       border =T,
                       grid_height = unit(2, "mm"),
                       grid_width = unit(2, "mm"),
                       title_gp =gpar(fontsize =7,fontface ="bold")
                     )
)

hb = rowAnnotation(Totals =anno_barplot(rowSums(mat),axis_param = list(gp =gpar(fontsize = 6))),
                   annotation_name_gp = gpar(fontsize = 6)
)
colfun =colorRamp2(c(0,100, 500, 1000, 1500, 2000),brewer.pal(6, "Reds"))
ht_down =Heatmap(mat,
                 col = colfun,
                 top_annotation = ha,
                 right_annotation = hb,
                 column_split = cc$Clusters,
                 column_title = NULL,
                 cluster_column_slices = F,
                 cluster_columns = T, 
                 row_names_gp = gpar(fontsize =6),
                 show_row_dend = F,
                 show_column_dend = F, 
                 row_names_side = "left", 
                 show_column_names = F,
                 heatmap_legend_param = list(
                   at =c(0, 500, 1000, 1500, 2000),
                   border =T,
                   title = "Event totals",
                   labels_gp = gpar(fontsize = 6),
                   grid_width = unit(3,"mm"),
                   legend_height = unit(0.6,"in"),
                   title_gp =gpar(fontsize =7,fontface ="bold")
                 ))

ht_up = grid.grabExpr(draw(ht_up, padding = unit(c(1, 1, 1, 1), "mm"), merge_legends =T))
ht_down = grid.grabExpr(draw(ht_down, padding = unit(c(1, 1, 1, 1), "mm"), merge_legends =T))
plot_grid(ht_up, ht_down,plot.new(), nrow = 1, rel_widths = c(1.5,1,0.5))

ht_COVID = grid.grabExpr(draw(ht_COVID, padding = unit(c(1, 1, 1, 1), "mm"), merge_legends =T))
ht_COVID_down = grid.grabExpr(draw(ht_COVID_down, padding = unit(c(1, 1, 1, 1), "mm"), merge_legends =T))

p16 = plot_grid(ht_up, ht_down, labels = letters, 
                label_size = 8, label_fontfamily = "sans" )


#### Consequence investigation ####

clusters = clusters.up

splicelist$Clusters = "Not recurrent"
m = match(splicelist$gene_id, clusters$V1)
f.a =!is.na(m)
f.t =m[f.a]
splicelist[f.a,]$Clusters = clusters[f.t,]$unmergedLabels.mod
conseqlist$Clusters = "Not recurrent"
m = match(conseqlist$gene_id, clusters$V1)
f.a =!is.na(m)
f.t =m[f.a]
conseqlist[f.a,]$Clusters = clusters[f.t,]$unmergedLabels.mod

con_list = list()
n =1
for (i in studies) {
  df1 = splicelist[splicelist$Clusters != "Not recurrent" & splicelist$ASchange != 'No change' & splicelist$Studies == i,]
  df1$id = paste(df1$isoformUpregulated,df1$isoformDownregulated, sep ="_")
  df2 = conseqlist[conseqlist$Clusters != "Not recurrent" & conseqlist$Studies == i,]
  df2$id = paste(df2$isoformUpregulated,df2$isoformDownregulated, sep ="_")
  m = match(df2$id, unique(df1$id))
  f.a =!is.na(m)
  con_list[[n]] = df2[f.a,]
  n = n+1
}
df_conq =do.call(rbind.data.frame, con_list)
df_conq = df_conq %>% group_by(Studies, Clusters, switchConsequence) %>% summarise(Total = n()) %>% mutate(Props = Total/sum(Total))


#### Plot dIF distribution of significant isoforms across

ggplot(data =isos_sig,aes(x =IF1, y =IF2, color =Clusters)) +
  geom_point(size=0.8) +
  geom_abline(slope = 1) + # default cutoff
  facet_wrap(vars(Studies),nrow =3)

ggplot(data =isos_sig,aes(x =iso_log2_fold_change, y =dIF, color =Recurrence)) +
  geom_point(size=0.8) +
  facet_wrap(vars(Studies),nrow =3) +
  labs(y='dIF', x='Isoform Log2 Fold Change') +
  theme_bw() 

#### Plot 1: Cluster Expression ####
load("E:/ProjectWork/Rfiles/COVID_clusters.rdata")
attr <- read.delim("E:/ensembl_v105.txt")

setwd("E:/Studies/IAV")
studies <- list.files()
test1 =list()
clusters = clusters.up
genes = as.character(clusters$V1)
n =1
for (i in studies) {
  load(paste0(i,"/DE.rdata"))
  m = match(genes, rownames(dds))
  f.a =!is.na(m)
  genes = genes[f.a]
  tcounts <- t(log2((counts(dds[genes, ], normalized=TRUE, replaced=FALSE)+1))) %>%
    merge(colData(dds), ., by="row.names") %>%
    tidyr::gather(gene, expression, (ncol(.)-length(genes)+1):ncol(.))
  tcounts = tcounts %>% 
    dplyr::select(Row.names, Condition, gene, expression)
  test = list()
  ##Averaging genes in clusters per sample:
  for (k in 1:length(unique(clusters$unmergedLabels.mod))) {
    samps = tcounts[tcounts$gene %in% clusters[clusters$unmergedLabels.mod == k,]$V1,]
    trs = list()
    for (j in 1:length(unique(samps$Row.names))) {
      bd = samps[samps$Row.names == unique(samps$Row.names)[j],]
      md = mean(bd$expression)
      trs[[j]] = cbind(as.character(unique(bd$Row.names)), as.character(unique(bd$Condition)),k,md)
    }
    test[[k]] = do.call(rbind.data.frame,trs)
  }
  test = do.call(rbind.data.frame, test)
  colnames(test) = c("SampleID",'Condition', 'Cluster',"expression")
  test$expression = as.numeric(test$expression)
  test$Studies = i
  test1[[n]] =test
  n = n+1
}
test1 = do.call(rbind.data.frame, test1)
test1$Condition = factor(test1$Condition, levels = c("Mock","Infected"))
test1$Cluster = factor(test1$Cluster, levels = c(1:max(test1$Cluster)))
test1$Studies = factor(test1$Studies, levels = unique(test1$Studies)[c(12:14,1:11,15)])

#Or with tcounts
m = match(tcounts$gene, clusters$V1)
f.a =!is.na(m)
f.t =m[f.a]
tcounts$Cluster = 0
tcounts[f.a,]$Cluster = clusters[f.t,]$unmergedLabels.mod

p1 = ggplot(test1, aes(Condition, expression, fill = Condition)) + 
  geom_violin(color =NA) + 
  geom_boxplot(width = 0.2,fill ="white", outlier.shape = NA) +
  labs(x="Conditions", 
       y="Normalized Gene Expression") +
  facet_grid(rows= vars(Cluster)) +
  theme_bw()+
  coord_flip() +
  geom_signif(map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "ns" = 1) , textsize=4, 
              comparisons = list(c("Mock", "Infected")),
              test = "wilcox.test",
              vjust = 0.2,hjust =0.4 ) +
  scale_fill_manual(values = c("Infected" = "tan2",
                               "Mock" = "steelblue2"), 
                    name="Conditions", guide = guide_legend(reverse = T)) +
  theme(legend.text =element_text(family = 'sans', size =6),
        legend.key.size = unit("2","mm"),
        legend.title = element_text(family = 'sans', size = 7, face ="bold"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(family = 'sans', size = 7),
        axis.text.x = element_text(family = 'sans', size = 6),
        strip.text.x = element_text(family = 'sans', size = 7),
        strip.text.y.left = element_text(family = 'sans', size = 7, angle = 0),
        strip.background = element_blank())
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- unique(clusters$V3)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
g1 =g

pdf("E:/ProjectWork/Figures/Splicing/IAV_Cluster+expression.pdf", width = 5, height = 5)
grid::grid.draw(g1)
dev.off()

#### Plot 2: Isoforms types identified across all samples per cluster across studies ####

dfg = list()
for (i in 1:max(clusters$unmergedLabels.mod)) {
  m = match(attr$Gene.stable.ID, clusters[clusters$unmergedLabels.mod ==i,]$V1)
  f.a =!is.na(m)
  cc = table(attr[f.a,]$Gene.stable.ID)
  f = data.frame(Genetype =c('SingleIso','MultiIso'), 
                 Totals = c(length(cc[cc==1]),length(cc[cc>1])),
                 Cluster = i,
                 Virus = "IAV")
  dfg[[i]] =f
}
dfg1 = do.call(rbind.data.frame, dfg)
ggplot(dfg1, aes(Genetype, y =Totals))+
  geom_col()+
  facet_wrap(~Cluster)

studies <- list.files()[c(12:14,1:11,15)]
load("E:/ProjectWork/Rfiles/recurrence_genesets.Rdata")
recur = rowSums(genesets.up, na.rm = T)
dfs = list()
k =1
df1 = as.data.frame(matrix(0, nrow = length(unique(attr$Gene.stable.ID)), ncol = length(studies)))
rownames(df1) = unique(attr$Gene.stable.ID)
colnames(df1) = studies

## Global Splicing
n =1
ll =list()
for (j in studies) {
  load(paste0(j,"/SwitchLists_annotated.rdata"))
  bed = read.delim(paste0(j,"/merged.annotated.bed"), header =F, stringsAsFactors = F)
  
  #Subset to only expressed transcripts of interest
  filt = abs(SwitchListAnalyzed$isoformFeatures$dIF) >= 0.1 & SwitchListAnalyzed$isoformFeatures$isoform_switch_q_value <= 0.05
  m = match(SwitchListAnalyzed$isoformFeatures[filt,]$isoform_id, bed$V4)
  f.a =!is.na(m)
  f.t =m[f.a]
  bed =bed[f.t,]
  bed = bed[bed$V1 %in% paste0('chr',1:22),]
  bed$V7 = gsub("=", 'Constitutive', bed$V7)
  bed$V7 = gsub("c", 'Constitutive', bed$V7)
  tryCatch({bed[bed$V7 != "Constitutive",]$V7 = "Alternative"}, error = function(e){})    
  bed = bed[bed$V8 != ".",]
  m = match(bed$V8, attr$Gene.name)
  f.a =!is.na(m)
  f.t =m[f.a]
  bed$V8[f.a] = attr$Gene.stable.ID[f.t]
  
  if(isFALSE(bed$V7 == "Alternative")){
    known = unique(bed[bed1$V7 == "Constitutive",]$V8)
    Splicing = c("Alternative", "Constitutive")
    df = as.data.frame(Splicing)
    df$Proportion <- as.numeric(c(0/length(unique(bed$V8)), length(known)/length(unique(bed$V8))))
    df$Total <- as.numeric(c(0, length(known)))
    df$Studies = j
  } else {
    found = unique(bed[bed$V7 == "Alternative",]$V8)
    known = unique(bed[bed$V7 == "Constitutive",]$V8)
    m = match(known,found)
    f.a =is.na(m)
    known = known[f.a] %>% unique()
    Splicing = c("Alternative", "Constitutive")
    df = as.data.frame(Splicing)
    df$Proportion <- as.numeric(c(length(found)/length(unique(bed$V8)), length(known)/length(unique(bed$V8))))
    df$Total <- as.numeric(c(length(found), length(known)))
    df$Studies = j
  } 
  dfs[[k]] = df
  
  #Customized Matrix for alternatively spliced genes per study
  df1[found, k] =1
  k = k+1
  ll[[n]] = SwitchListAnalyzed$designMatrix
  ll[[n]]$batch = j
  ll[[n]] =ll[[n]][,c('sampleID','condition','batch')]
  n =n+1
}
df = do.call(rbind.data.frame, dfs)
colnames(df)[1] ="Splicing"
df2 = df
df2$Studies = factor(df2$Studies, levels = unique(df2$Studies)[c(12:14,1:11,15)])
df2[grep("Constitutive", df2$Splicing),]$Total = df2[grep("Constitutive", df2$Splicing),]$Total*-1

avgs <- df %>% dplyr::group_by(Splicing) %>% dplyr::summarise(Averages = mean(Proportion, na.rm =T))
avgs[grep("Alternative", avgs$Splicing),]$Averages = 1 -avgs[grep("Alternative", avgs$Splicing),]$Averages
avgs =as.data.frame(avgs)

p1 <- ggplot(df2, aes(x = Studies, fill = Splicing)) +
  geom_bar(aes(y = Proportion),position="stack", stat="identity") +
  scale_fill_hue() +
  theme_classic() +
  geom_hline(data =avgs, mapping =aes(yintercept =Averages), size = 1, linetype =2, color ='black') +
  geom_text(data =avgs, aes(x = 7.5, y = Averages+0.05, label =paste("Average =",round(Averages[1], digits = 2)))) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1,vjust = 1, size = 6))

p2 <- ggplot(df2, aes(x = Studies, y = Totals, fill = Splicing)) +
  geom_col(aes(y=Total)) +
  scale_fill_hue() +
  theme_classic() +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1,vjust = 1, size = 6))

pdf("E:/ProjectWork/Figures/Splicing/IAV_Splicing_proportions-totals.pdf", width = 9, height = 5)
plot_grid(p1, p2, nrow = 1)
dev.off()

## Frequency of Alt splicing across studies

altrecur = rowSums(df1, na.rm = T)
res = data.frame(altrecur[altrecur>0])
res = data.frame(table(res), row.names = NULL)
colnames(res) = c("Recurrence","counts")
res$counts = as.numeric(res$counts)

clustdf = list()
for(i in 1:max(clusters$unmergedLabels.mod)){
  m = match(names(altrecur), clusters[clusters$unmergedLabels.mod == i,]$V1)
  f.a =!is.na(m)
  clustrec = altrecur[f.a][altrecur[f.a] >0]
  clustrec = data.frame(table(clustrec), row.names = NULL)
  colnames(clustrec) = c("Recurrence","counts")
  clustrec$Cluster =i
  clustdf[[i]] = clustrec
}
clustdf = do.call(rbind.data.frame, clustdf)
clustdf$Proportion = sapply(1:nrow(clustdf), function(x) clustdf[x,]$counts/sum(clustdf[clustdf$Recurrence == clustdf[x,]$Recurrence,]$counts))
clustdf$Cluster = factor(clustdf$Cluster, levels = 1:max(clustdf$Cluster))
clustdf$Recurrence = factor(as.numeric(clustdf$Recurrence), levels = 1:max(as.numeric(clustdf$Recurrence)))
p3a <- ggplot(data = res) +
  geom_col(aes(x = Recurrence, y = counts),color = "black", fill = "coral3") +
  labs(x = "Recurrence", y = "Genes with at least one isoform\nderived from Alternative Splicing" ) +
  theme_classic() +
  geom_text(aes(x =Recurrence,y = counts+80,label = abs(counts),group = Recurrence), 
            size = 3) 
cols = c(unique(clusters$V3))
names(cols) = c(unique(clusters$unmergedLabels.mod))
p3b <- ggplot(data = clustdf) +
  geom_col(aes(x = Recurrence, y = counts, fill = Cluster),color = "black") +
  labs(x = "Recurrence", y = "Genes totals in clusters\nwith Alternative Splicing" ) +
  theme_classic() +
  scale_fill_manual(values =cols)
p3c <- ggplot(data = clustdf) +
  geom_col(aes(x = Recurrence, y =Proportion, fill = Cluster),color = "black") +
  labs(x = "Recurrence", y = "Gene proportions in clusters\nwith Alternative Splicing" ) +
  theme_classic() +
  scale_fill_manual(values =cols)+
  theme(legend.position = "none",
        axis.title.x = element_blank())
p3 = p3b + annotation_custom(ggplotGrob(p3c), xmin = 5, xmax = 10, ymin = 60, ymax =150) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

attr2 = attr
attr2$Chromosome.scaffold.name = paste0("chr", attr2$Chromosome.scaffold.name)
attr2[is.na(match(attr2$Chromosome.scaffold.name,paste0("chr", c(1:22,"X", "Y", "MT")))),]$Chromosome.scaffold.name = "Scaffold"
res.list <- list()
n=1
for(i in studies){
  res.list[[n]] = data.frame(unique(attr2$Chromosome.scaffold.name))
  colnames(res.list[[n]]) = "Chromosome"
  chrmdf = data.frame(table(attr2[attr$Gene.stable.ID %in% rownames(df1[df1[,i]==1,]),]$Chromosome.scaffold.name))
  res.list[[n]]$Total =0
  m = match(res.list[[n]]$Chromosome, chrmdf$Var1)
  f.a =!is.na(m)
  f.t =m[f.a]
  res.list[[n]][f.a,]$Total = chrmdf[f.t,]$Freq
  res.list[[n]]$Studies =i
  res.list[[n]]$Proportion = res.list[[n]]$Total/sum(res.list[[n]]$Total)
  n = n+1
}
res.list = do.call(rbind.data.frame, res.list)
res.list$Chromosome =factor(res.list$Chromosome, levels = c(paste0("chr", c(1:22,"X", "Y", "MT")), "Scaffold"))
p4 <- ggplot(res.list, aes(x =Chromosome, fill = Chromosome))+
  geom_bar(aes(y = Total), stat="identity") +
  facet_wrap(~Studies, ncol =3) +
  scale_fill_discrete() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x.bottom = element_text(angle = 45, hjust = 1,vjust = 1, size = 7))


pdf("E:/ProjectWork/Figures/Splicing/IAV_Splicing_totals across Chrom per Recurrence.pdf", height = 12, width = 8, paper ="a4", useDingbats = T)
p4
dev.off()

##Overlap -VennDiagram
library(VennDiagram)
res = data.frame(altrecur[altrecur>0])
alt.rec = rownames(res)
recur = rowSums(genesets.up, na.rm =T)
m = match(names(recur), names(altrecur))
f.a =!is.na(m)
f.t =m[f.a]
venn.plot <- venn.diagram(
  x = list(
    Alt = alt.rec %>% unlist() , 
    Recurrence = inf_clusters$V1 %>% unlist()
  ),
  filename ="E:/ProjectWork/Figures/Splicing/IAV_Venn2.tiff",
  disable.logging =T,
  category.names = c("AS genes" , "Recurrence"),
  col=c("#440154ff", '#21908dff'),
  fill = c("red3",'green3'),
  scaled = TRUE,
  ext.text = TRUE,
  ext.line.lwd = 2,
  ext.dist = c(-.15,-0.09),
  ext.length = c(0.9,0.75),
  ext.pos = c(-45,235),
  fontfamily = "sans",
  fontface = "bold",
  cex = 1,
  cat.cex = 1,
  cat.default.pos = "outer",
  cat.fontface ="bold",
  cat.fontfamily = "sans",
  cat.just = list(c(1,2),c(0,2)),
  cat.col = c("red3", 'green3'),
  margin =0.1)

##Splicing across clusters
dfs =list()
k=1
for (j in studies) {
  load(paste0(j,"/SwitchLists.rdata"))
  bed = read.delim(paste0(j,"/merged.annotated.bed"), header =F, stringsAsFactors = F)
  
  #Subset to only expressed transcripts of interest
  
  filt = abs(SwitchListAnalyzed$isoformFeatures$iso_log2_fold_change)>= 0.1 & abs(SwitchListAnalyzed$isoformFeatures$dIF) >= 0.1 & SwitchListAnalyzed$isoformFeatures$isoform_switch_q_value <= 0.05
  m = match(SwitchListAnalyzed$isoformFeatures[filt,]$isoform_id, bed$V4)
  f.a =!is.na(m)
  f.t = m[f.a]
  bed =bed[f.t,]
  m = match(bed$V8, attr$Gene.name)
  f.a =!is.na(m)
  f.t =m[f.a]
  bed$V8[f.a] = attr$Gene.stable.ID[f.t]
  m = match(bed$V8, clusters$V1)
  f.a =!is.na(m)
  f.t =m[f.a]
  bed$cluster = 0
  bed$cluster[f.a] = clusters$unmergedLabels.mod[f.t]
  df =list()
  for (i in 1:max(clusters$unmergedLabels.mod)) {
    bed1 = bed[bed$cluster == i,]
    bed1$V7 = gsub("=", 'Constitutive', bed1$V7)
    bed1$V7 = gsub("c", 'Constitutive', bed1$V7)
    tryCatch({bed1[bed1$V7 != "Constitutive",]$V7 = "Alternative"}, error = function(e){})    
    if(isFALSE(bed1$V7 == "Alternative")){
      known = unique(bed1[bed1$V7 == "Constitutive",]$V8)
      Splicing = c("Alternative", "Constitutive")
      df[[i]] = as.data.frame(Splicing)
      df[[i]]$Proportion <- as.numeric(c(0/length(unique(bed1$V8)), length(known)/length(unique(bed1$V8))))
      df[[i]]$Total <- as.numeric(c(0, length(known)))
      df[[i]]$Cluster <- i
    } else {
      found = unique(bed1[bed1$V7 == "Alternative",]$V8)
      known = unique(bed1[bed1$V7 == "Constitutive",]$V8)
      m = match(known,found)
      f.a =is.na(m)
      known = known[f.a] %>% unique()
      Splicing = c("Alternative", "Constitutive")
      df[[i]] = as.data.frame(Splicing)
      df[[i]]$Proportion <- as.numeric(c(length(found)/length(unique(bed1$V8)), length(known)/length(unique(bed1$V8))))
      df[[i]]$Total <- as.numeric(c(length(found), length(known)))
      df[[i]]$Cluster <- i
    }
  }
  df = do.call(rbind.data.frame, df)
  df$Studies = j
  dfs[[k]] = df
  k = k+1
}
df = do.call(rbind.data.frame, dfs)
df2 = df
df2$Studies = factor(df2$Studies, levels = unique(df2$Studies)[c(12:14,1:11,15)])
df2[grep("Constitutive", df2$Splicing),]$Total = df2[grep("Constitutive", df2$Splicing),]$Total*-1

avgs <- df %>% group_by(Splicing,Cluster) %>% summarise(Averages = mean(Proportion, na.rm =T))
avgs =as.data.frame(avgs)

p5 <- ggplot(df2, aes(x = Studies, fill = Splicing)) +
  geom_bar(aes(y = Proportion),position="stack", stat="identity") +
  scale_fill_hue() +
  theme_bw() +
  facet_grid(cols  =vars(Cluster)) +
  geom_hline(data =avgs[avgs$Splicing == 'Constitutive',], mapping =aes(yintercept =Averages), size = 1, linetype =2, color ='black') +
  geom_text(data =avgs[avgs$Splicing == 'Constitutive',], 
            aes(x = 7.5, y = Averages+0.05, label =paste("Mean =",round(Averages, digits = 2)))) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1,vjust = 1.1, size = 6))

g <- ggplot_gtable(ggplot_build(p5))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- unique(clusters$V3)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
g2 =g

p6 =ggplot(df2, aes(x = Studies, fill = Splicing)) +
  geom_col(aes(y=Total)) +
  scale_fill_hue() +
  theme_bw() +
  facet_grid(cols  =vars(Cluster)) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1,vjust = 1.1, size = 6))

g <- ggplot_gtable(ggplot_build(p6))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- unique(clusters$V3)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
g3 =g

grid::grid.draw(g3)

plot_grid(g2,g3,nrow = 1)

#### Plot 3: Boxplot of IsoformSwitch totals per cluster per study ####

setwd("E:/Studies/IAV")
df = list()
n=1
isos = list()
clusters = inf_clusters
for (k in studies) {
  load(paste0(k,"/SwitchLists_annotated.rdata"))
  totals =list()
  #m = match(SwitchListAnalyzed$isoformFeatures$gene_id, clusters$V1)
  #f.a = !is.na(m)
  #exampleSwitchListAnalyzedSubset <- subsetSwitchAnalyzeRlist(SwitchListAnalyzed, f.a)
  isos[[n]] = as.data.frame(SwitchListAnalyzed$isoformFeatures)
  for (i in 1:max(clusters$unmergedLabels.mod)) {
    m = match(SwitchListAnalyzed$isoformFeatures$gene_id, clusters[clusters$unmergedLabels.mod == i,]$V1)
    f.a = !is.na(m)
    exampleSwitchListAnalyzedSubset <- subsetSwitchAnalyzeRlist(
      SwitchListAnalyzed, f.a)
    totalConQ <- length(na.omit(unique(exampleSwitchListAnalyzedSubset$isoformFeatures[exampleSwitchListAnalyzedSubset$isoformFeatures$isoform_switch_q_value <= 0.05 & abs(exampleSwitchListAnalyzedSubset$isoformFeatures$dIF) >= 0.1,]$gene_id)))
    totals[[i]] <- data.frame(matrix(0, nrow = 1, ncol = 2))
    totals[[i]][,1] = totalConQ
    totals[[i]][,2] = totalConQ/length(clusters[clusters$unmergedLabels.mod ==i,]$V1)
  }
  totals = do.call(rbind.data.frame, totals)
  colnames(totals) =c("Total","Proportion")
  totals$Clusters = 1:nrow(totals)
  df[[n]] = totals
  #df[[n]]$Proportion <- df[[n]]$Total/length(unique(SwitchListAnalyzed$isoformFeatures[SwitchListAnalyzed$isoformFeatures$switchConsequencesGene == TRUE,]$gene_id))
  df[[n]]$Studies = k
  isos[[n]]$Studies = k
  load(paste0(k,"/DE.rdata"))
  m = match(isos[[n]]$gene_id, rownames(res))
  f.a = !is.na(m)
  f.t = m[f.a]
  isos[[n]][f.a,]$gene_log2_fold_change = res[f.t,]$log2FoldChange
  isos[[n]][f.a,]$gene_q_value = res[f.t,]$padj
  
  ##Assign class codes to isos
  bed = read.delim(paste0(k,"/merged.annotated.bed"), header =F, stringsAsFactors = F)
  m = match(bed$V8, attr$Gene.name)
  f.a =!is.na(m)
  f.t =m[f.a]
  bed$V8[f.a] = attr$Gene.stable.ID[f.t]
  
  m = match(isos[[n]]$isoform_id, bed$V4)
  f.a =!is.na(m)
  f.t =m[f.a]
  isos[[n]]$`Class Code` = 0
  isos[[n]][f.a,]$`Class Code` =bed[f.t,]$V7
  n =n+1
}

isos = do.call(rbind.data.frame, isos)
m = match(isos$gene_id, clusters$V1)
f.a = !is.na(m)
f.t =m[f.a]
isos$Clusters = 0
isos[f.a,]$Clusters = clusters[f.t,]$unmergedLabels.mod
isos[is.na(isos$isoform_switch_q_value),]$isoform_switch_q_value = 1
isos$Studies = factor(isos$Studies, levels = unique(isos$Studies)[c(12:14,1:11,15)])

filt1 = abs(isos$dIF) >= 0.1 & isos$isoform_switch_q_value <= 0.05 & abs(isos$iso_log2_fold_change) >= 1
filt2 = abs(isos$dIF) >= 0.1 & isos$isoform_switch_q_value <= 0.05 & abs(isos$iso_log2_fold_change) < 1
isos$`Significant Isoform Switching` = 'Not Significant'
isos[filt1,]$`Significant Isoform Switching` ="FDR <= 0.05 +\nLog2FC + dIF"
isos[filt2,]$`Significant Isoform Switching` ="FDR <= 0.05\n+ dIF"

q =ggplot(data =isos_sig,aes(x =iso_log2_fold_change, y =dIF, color =`Significant Isoform Switching`)) +
  geom_point(size=0.8) +
  geom_vline(xintercept = c(-1,1), linetype='dashed') + # default cutoff
  geom_hline(yintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_grid(~Studies) +
  labs(y='dIF', x='Gene Log2 Fold Change') +
  theme_bw() +
  scale_color_manual('Signficant\nIsoform\nSwitching', values = rev(c('gray','red', 'blue'))) 

filt = abs(isos$dIF) >= 0.1 & isos$isoform_switch_q_value <= 0.05
dfs = isos[filt,] %>% group_by(Clusters, Studies,`Class Code`) %>% summarise(Total = length(isoform_id),
                                                                             Proportion = length(isoform_id)) %>% as.data.frame()
for (i in 1:nrow(dfs)) {
  st = dfs[i,]$`Class Code`
  cl = dfs[i,]$Clusters
  total = sum(dfs[dfs$`Class Code` == st & dfs$Clusters == cl,]$Total)
  dfs[i,]$Proportion = dfs[i,]$Proportion/total
}                                                             

ggplot(data =psi,aes(x = `ΔPSI....`, y =-log10(T.test.p.value), fill =Event.Type)) +
  geom_point() + 
  facet_wrap(~Clusters) +
  labs(y='dIF', x='Gene Log2 Fold Change') +
  theme_bw() 


g <- ggplot_gtable(ggplot_build(q))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- unique(clusters$V3)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("E:/ProjectWork/Figures/Splicing/IAV_clusters_DEI_volplots.pdf", height = 13, width = 8, paper ="A4")
grid::grid.draw(g)
dev.off()


df = do.call(rbind.data.frame, df)
df$Studies = factor(df$Studies, levels = unique(df$Studies)[c(12:14,1:11,15)])
df$Clusters = factor(df$Clusters)
p8 <- ggplot(df, aes(x = Clusters,y = Total,group =Clusters, fill = Clusters)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, size =1) +
  theme_bw() +
  scale_fill_manual(values = cols, 
                    name="Clusters") +
  geom_signif(map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05) , textsize=5, 
              comparisons = list(c("1", "3"),c("2","3"),c("3","4"),c("3","5")),
              show.legend = T, 
              test = "t.test",
              vjust = 0.4,
              step_increase = 0.1 )

#### Plot 4: Splicing Recurrence ####

n=1
splicelist = list()
for (k in studies) {
  load(paste0(k,"/SwitchLists_annotated.rdata"))
  totals =list()
  for (i in 1:max(clusters$unmergedLabels.mod)) {
    m = match(SwitchListAnalyzed$isoformFeatures$gene_id, clusters[clusters$unmergedLabels.mod == i,]$V1)
    f.a = !is.na(m)
    exampleSwitchListAnalyzedSubset <- subsetSwitchAnalyzeRlist(
      SwitchListAnalyzed, f.a)
    ##Barplots of gene totals per splicing event (e.g ATSS, ES)
    tryCatch({ 
      totals[[i]] =
        extractSplicingEnrichment(
          exampleSwitchListAnalyzedSubset,
          asFractionTotal = T,
          alpha = 0.05,
          dIFcutoff = 0,
          plot =F,
          plotGenes=F,
          returnResult =T)
      totals[[i]]$Clusters = i
    }, error = function(e) {})
    
  }
  totals = do.call(rbind.data.frame, totals)
  if(nrow(totals) > 0) totals$Studies = k
  splicelist[[n]] = totals
  n = n+1
}  

splicelist = do.call(rbind.data.frame, splicelist)
splicelist$Studies = factor(splicelist$Studies, levels = unique(splicelist$Studies)[c(12:14,1:11,15)])
clust.up = genesets.up[clusters$V1,]
m = match(splicelist$Studies, colnames(genesets.up))
f.a =!is.na(m)
f.t =m[f.a]  

splicelist$`DE totals per study` = 0
splicelist[f.a,]$`DE totals per study` = colSums(genesets.up, na.rm = T)[f.t]
splicelist$`DE totals per Cluster` = 0
splicelist[f.a,]$`DE totals per Cluster` = colSums(clust.up, na.rm = T)[f.t]
m = match(splicelist$Clusters, names(table(clusters$unmergedLabels.mod)))
f.a =!is.na(m)
f.t =m[f.a]  
splicelist$`Totals per Cluster` =0
splicelist[f.a,]$`Totals per Cluster` = table(clusters$unmergedLabels.mod)[f.t]

p6 <-ggplot(splicelist, aes(x= AStype,y = nrGenesWithConsequences, group =AStype, fill = AStype)) +
  geom_boxplot(width =0.5,outlier.shape = NA)+
  geom_jitter(width = 0.2, size =0.8) +
  labs(y = "Number of significant genes\n(with at least one event)", x = 'Alternative transcription event') +
  scale_fill_viridis_d() +
  facet_grid(cols =vars(Clusters)) +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "none")

g <- ggplot_gtable(ggplot_build(p6))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- unique(clusters$V3)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

g4cov =g
#### Plot 5: Consequence Summary ####

n=1
ConSeqlist = list()
for (k in studies) {
  load(paste0(k,"/SwitchLists.rdata"))
  totals =list()
  for (i in 1:max(clusters$unmergedLabels.mod)) {
    m = match(SwitchListAnalyzed$isoformFeatures$gene_id, clusters[clusters$unmergedLabels.mod == i,]$V1)
    f.a = !is.na(m)
    exampleSwitchListAnalyzedSubset <- subsetSwitchAnalyzeRlist(
      SwitchListAnalyzed, f.a)
    ##Barplots of gene totals per splicing event (e.g ATSS, ES)
    tryCatch({ 
      totals[[i]] =
        extractConsequenceSummary(
          exampleSwitchListAnalyzedSubset,
          plotGenes=F,
          
          plot = F,
          returnResult =T)
      totals[[i]]$Clusters = i
      totals[[i]]$Studies = k
    }, error = function(e) {})
    
  }
  totals = do.call(rbind.data.frame, totals)
  ConSeqlist[[n]] = totals
  n = n+1
}  

ConSeqlist = do.call(rbind.data.frame, ConSeqlist)

cc =list()
for (j in 1:length(unique(ConSeqlist$Studies))) {
  gg = list()
  dd = list()
  for (k in 1:max(clusters$unmergedLabels.mod)) {
    filt = ConSeqlist$Clusters == k & ConSeqlist$Studies == unique(ConSeqlist$Studies)[j]
    gg[[k]] = ConSeqlist[filt,]
    m = match(unique(ConSeqlist$switchConsequence), gg[[k]]$switchConsequence)
    f.a =is.na(m)
    tryCatch({switchfeat = unique(ConSeqlist$switchConsequence)[f.a]
    dd[[k]] = data.frame(Comparison ="Infected vs Mock", featureCompared ="Unk", switchConsequence = switchfeat, 
                         nrGenesWithConsequences =0,nrIsoWithConsequences =0, Clusters = k, Studies = unique(ConSeqlist$Studies)[j])}, error =function(e){})
  }
  gg = do.call(rbind.data.frame, gg)
  dd = do.call(rbind.data.frame, dd)
  cc[[j]] = rbind(gg, dd)
}

cc = do.call(rbind.data.frame, cc)
ConSeqlist2 = ConSeqlist
ConSeqlist = cc

ConSeqlist$Studies = factor(ConSeqlist$Studies, levels = unique(ConSeqlist$Studies)[c(12:14,1:11,15)])
ConSeqlist = ConSeqlist[ConSeqlist$switchConsequence != "Any consequence",]

avgs = ConSeqlist[ConSeqlist$Clusters == 1,] %>% group_by(switchConsequence) %>% dplyr::summarise(Averages = median(nrIsoWithConsequences, na.rm =T)) 
avgs = as.data.frame(avgs[order(avgs$Averages),])[,1] %>% as.character()
ConSeqlist$switchConsequence = factor(ConSeqlist$switchConsequence, levels =avgs)

p7 <- ggplot(data =ConSeqlist, aes(x = switchConsequence,y = nrGenesWithConsequences , fill = switchConsequence)) +
  geom_boxplot(width =0.5,outlier.shape = NA)+
  geom_jitter(width = 0.2, size =0.8) +
  scale_fill_viridis_d() +
  labs(y = "Number of significant genes\n(with at least one event)", x = 'Consequence of Isoform Switch') +
  facet_grid(rows =vars(Clusters)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1,vjust = 1, size = 8),
        legend.position = "none")
g <- ggplot_gtable(ggplot_build(p7))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- unique(clusters$V3)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

g5 =g

#### Plot 6: Isoform Switch frequency across studies ####

load("E:/ProjectWork/Rfiles/COVID_recurrence_genesets.Rdata")
recur = rowSums(genesets.up, na.rm = T)
studies = list.files(recursive = T, pattern = "SwitchLists.rdata")
n=1
df1 = as.data.frame(matrix(0, nrow = length(unique(attr$Gene.stable.ID)), ncol = length(studies)))
rownames(df1) = unique(attr$Gene.stable.ID)
colnames(df1) = dirname(studies)
for (i in studies) {
  load(paste0(i,"/SwitchLists.rdata"))
  set = SwitchListAnalyzed$isoformFeatures[SwitchListAnalyzed$isoformFeatures$switchConsequencesGene == T,]$gene_id %>% unique %>% na.omit()
  m = match(set, attr$Gene.name)
  f.a = !is.na(m)
  f.t = m[f.a]
  set[f.a] =attr[f.t,]$Gene.stable.ID
  set = unique(set)
  m = match(set, attr$Gene.stable.ID)
  f.a =!is.na(m)
  set =set[f.a]
  df1[set,n] =1
  n =n+1
}

IsoSwitchRecur =rowSums(df1, na.rm = T)
merge(recur, IsoSwitchRecur)
plot(recur, rowSums(df1, na.rm = T))


switchPlot(SwitchListAnalyzed, gene = "ENSG00000204677")


df1 = isos[,c("IF1", "IF2","Clusters")] %>% reshape2::melt(value.name = "Fraction")
colnames(df1)[2] ="Condition"

ggplot(data=df1, aes(x=Condition, y = -log(Fraction))) +
  geom_boxplot() +
  facet_wrap(~Clusters)
theme_bw()
