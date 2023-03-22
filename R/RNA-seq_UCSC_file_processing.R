##Date: 16/2/23
##Written by: Conor Cremin
##Title:UCSC Analysis for RNA-Seq

###################Description:#######################

#Function built to utilize the makeUCSCfile function from HOMER. A metadata file with columns specifying
#Run, BamReads, SampleID, Condition, Antibody (must specifiy the 'input' samples in ChIP-seq data), UCSColor

######################################################

#Packages:
if( !require("dplyr")){
  BiocManager::install("dplyr")
  library(dplyr)
} 
if( !require("Rsamtools")){
  BiocManager::install("Rsamtools")
  library(Rsamtools)
}
if( !require("DESeq2")){
  BiocManager::install("DESeq2")
  library(DESeq2)
}
if( !require("doMC")){
  BiocManager::install("doMC")
  library(doMC)     # Only available on UNIX Version of R
} 

old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path,"/software/samtools/1.11/bin/samtools","/home/conor93/HOMER/bin/", sep =":"))

setwd("/home/hlchen/Conor/Studies/IAV")

studies = list.files()

#Functions:

Study_screen <- function(study, bamdir ="TwoPass", n = 6){
  load("DE.rdata")
  stdout = SummarizedExperiment::colData(dds)[,c('Run', 'Condition')]
  samples = paste0(bamdir,"/",stdout$Run,"Aligned.sortedByCoord.out.bam")
  registerDoMC(n)
  t =  foreach(i = 1:length(samples)) %dopar% {
    testPairedEndBam(samples[i])
  }
  t = unlist(t)
  stdout$Paired = t
  stdout$BamReads = samples
  return(stdout)
}

make_tag_dir <- function(stdout, n =6){
  tagfun <- '/home/conor93/HOMER/bin/makeTagDirectory'
  samtools <- "/software/samtools/1.11/bin/samtools"
  registerDoMC(n)
  for(i in 1:nrow(stdout)){
    cmd1 = paste(samtools, "view -h -@",n,stdout[i,]$BamReads,"-o",paste0(gsub(".bam","",stdout[i,]$BamReads),".sam"))
    system(cmd1)
    if(isTRUE(stdout$Paired[i])){
      params = paste("-sspe -flip -genome hg38 -tbp 1")
    } else {
      params = paste("-genome hg38 -tbp 1")
    }
    if(!dir.exists("UCSC")){dir.create('UCSC')}
    cmd2 = paste(tagfun, paste0("UCSC/",stdout[i,]$Run),params,paste0(gsub(".bam","",stdout[i,]$BamReads),".sam"))
    system(cmd2)
    unlink(paste0(gsub(".bam","",stdout[i,]$BamReads),".sam"))
  }
}

scan_tag_dir <- function(dir, chr = paste0('chr',1:22, '.tags.tsv')){
  unlink(list.files(dir,full.names = T, pattern = '.tags.tsv')[!list.files(dir,pattern = '.tags.tsv') %in% chr])
}

make_Bedgraph <-function(stdout, homer_style = "rnaseq"){
  bgfunc <- "/home/conor93/HOMER/bin/makeUCSCfile"
  if (!dir.exists("UCSC/Bedgraph")) dir.create('UCSC/Bedgraph')
  for(i in 1:nrow(stdout)){
    if(isTRUE(stdout$Paired[i])){
      params = paste("-o",paste0("UCSC/Bedgraph/",stdout[i,]$Run), "-pseudo 1 -log -fsize 20e7","-style",homer_style, 
                     "-color", as.character(cols[stdout[i,]$Condition]) , 
                     "-strand separate -name", paste0(stdout[i,]$Condition,"-",i))
    } else {
      params = paste("-o",paste0("UCSC/Bedgraph/",stdout[i,]$Run), "-pseudo 1 -log -fsize 20e7","-style",homer_style, 
                     "-color", as.character(cols[stdout[i,]$Condition]) , 
                     "-strand both -name", paste0(stdout[i,]$Condition,"-",i))
    }
    cmd1 = paste(bgfunc, paste0("UCSC/",stdout[i,]$Run),params)
    system(cmd1)
  }
}

cols =c('Infected' = "210,180,140", "Mock" = "70,130,180")

for (study in studies) {
  setwd(study)
  stdout = Study_screen(study, bamdir ="TwoPass", n = 6)
  # make_tag_dir(stdout, n =6)
  for (dir in paste0("UCSC/",stdout$Run)) {
    scan_tag_dir(dir, chr = paste0('chr',1:22, '.tags.tsv'))
  }
  make_Bedgraph(stdout, homer_style = "rnaseq")
  setwd("..")
}

scan_motifs <- function(file){
  func <- "/home/conor93/HOMER/bin/findMotifsGenome.pl"
  if(!dir.create("Motifs")) dir.create("Motifs")
  params =paste('hg38 Motifs -h -p 6')
  cmd =paste(func, file, params)
  system(cmd)
}

