##Title: STARscript
##Author: Conor Cremin
##Company: University of Hong Kong

##Description:
###Requires sample directories to have the same format

##Functions

environ_set_up <- function(){
  if(!require(dplyr)){
    print(paste("dplyr pkg not available, installing now and re-run"))
    BiocManager::install("dplyr")
  }
  if(!require(data.table)){
    print(paste("data.table pkg not available installing and re-run"))
    install.packages("data.table")
  }
  if(!require(doMC)){
    print(paste("dplyr pkg not available, installing now and re-run"))
    install.packages("doMC")
    library(doMC)     # Only available on UNIX Version of R
  }
}

mapSingleRead <- function(inputdir ="fastq", run, star_outdir, IDX){
  print(paste('Mapping',run))
  fastqFile <- list.files(inputdir, full.names = T, pattern =run)
  fileNamePrefix <- paste0(star_outdir,"/",run)
  STARcmd <- paste(STAR,STARParameters, 
                   '--outFileNamePrefix', fileNamePrefix,
                   "--readFilesIn", fastqFile)
  system(STARcmd)
}

mapPairedRead <- function(inputdir ='fastq',run,star_outdir,IDX){
  print(paste('Mapping',run))
  fileNamePrefix <- paste0(star_outdir,"/",run)
  fastqFile1 <- list.files(inputdir, full.names = T, pattern =run)[1]
  fastqFile2 <- list.files(inputdir, full.names = T, pattern =run)[2]
  STARcmd <- paste(STAR, STARParameters,
                   '--outFileNamePrefix', fileNamePrefix,
                   '--readFilesIn', fastqFile1,fastqFile2)
  system(STARcmd)
}

sample_selection <- function(inputdir ="fastq"){
  samples <- tryCatch(
    {list.files(inputdir, full.names =F,pattern = ".gz") %>% 
    sapply(function(x)strsplit(x,'[._]')[[1]][1]) %>% as.vector
    }, error = function(e) {
      print("sample_selection: Error in extracted sample names, consider renaming if (.-) present in filenames")
      stop()
    })
  pr <- duplicated(samples) %>% samples[.]
  sr <- setdiff(samples%>%unique%>%as.vector,pr)
  print(paste('Single read:',length(sr),'-------','Paired reads:',length(pr)))
  if(length(sr) > 0 & length(pr)==0) {
    names(sr) =rep("single",length(sr))
    return(sr)
  } else if(length(pr) > 0 & length(sr) == 0) {
    names(pr) =rep("paired", length(pr))
    return(pr)
  } else{
    names(sr) =rep("single",length(sr))
    names(pr) =rep("paired", length(pr))
    runs = c(pr,sr)
    return(runs)
  }
  if(length(sr) == 0 & length(pr)==0) print(paste("sample_selection error: Cannot identify if fastq files are single/paired."))
}

mapRNAStudies <- function(inputdir ='fastq',study, IDX,star_outdir,parallel_threads = 1){
  paste('Mapping study:',study)
  if(!dir.exists(star_outdir)) dir.create(star_outdir)
  registerDoMC(parallel_threads)
  runs = sample_selection(inputdir)    
  foreach (j = 1:length(runs)) %dopar% {
    if(names(runs)[j] =="paired"){
      mapPairedRead(inputdir,runs[j],star_outdir,IDX)
      
    } else if(names(runs)[j] =="single"){
      mapSingleRead(inputdir,runs[j],star_outdir,IDX)
    }
    samidxcmd = paste(samtools, "index", paste0(star_outdir,"/",runs[j],"Aligned.sortedByCoord.out.bam"), 
                      paste0(star_outdir,"/",runs[j],"Aligned.sortedByCoord.out.bam.bai"))
    system(samidxcmd)
  }
  if (length(runs) == 0){ print(paste("mapRNAStudies:","for",study,", No Alignments IDs identified -Check fastq filenames!"))
  stop()
  }
}

rmdupPCR <- function(star_outdir,rmdupdir){
  input =list.files(star_outdir, full.names = T, pattern = ".RG.bam$")[8]
  if(!dir.exists(rmdupdir)) dir.create(rmdupdir)
  for(j in 1:length(input)){
    cmd1 <- paste(gatk, "MarkDuplicatesSpark",
                  "-I",input[j],
                  "-O", paste0(rmdupdir,"/",gsub("Aligned.RG", ".rmdup",basename(input[j]))),
                  "-M", paste0(rmdupdir,"/",gsub("Aligned.RG.bam", "_dedup_metrics.txt",basename(input[j]))),
                  rmdup_params)
    system(cmd1)
  }
}

Add_ReadGroups <- function(star_outdir,STARthreads){
  registerDoMC(STARthreads)
  foreach (k = list.files(star_outdir,pattern ="Aligned.sortedByCoord.out.bam$",full.names = T))%dopar% {
    cmd <- paste(java, "-jar", picard, "AddOrReplaceReadGroups",paste0("O=",gsub("Aligned.sortedByCoord.out.bam","Aligned.RG.bam", k)), 
                 paste0("I=",k),paste("RGLB=lib1 RGID=4 RGPL=ILLUMINA RGPU=unit1 RGSM=20"))
                 
    system(cmd)
  }
}

SplitBAM <- function(split_indir,split_outdir, split_parallel_threads){
  input =list.files(split_indir, full.names = T, pattern = split_suffix)
  if(!dir.exists(split_outdir)) dir.create(split_outdir)
  registerDoMC(split_parallel_threads)
  if(Do_chromosome_subset == "unmapped"){
    foreach(j = 1:length(input)) %dopar%{
      cmd1 <- paste(samtools, "view",subsetparams, input[j], ">",
                    paste0(split_outdir,"/",gsub("rmdup","unmapped",basename(input[j]))))
      system(cmd1)
    }
  } else if(Do_chromosome_subset == "ByChromosome"){
    foreach(j = 1:length(input)) %dopar%{
      cmd1 <- paste(bamtools, "split -in",input[j], subsetparams, "-refPrefix", 
                    paste0(split_outdir,"/",gsub(split_suffix,"",basename(input[j]))))
      
      system(cmd1)
    }
  }
}
  

