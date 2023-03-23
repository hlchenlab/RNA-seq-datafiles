##Title: Varient calling
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


GATK_CigarSplitter <- function(GATK_input, GATK_outdir, GATK_suffix,GATK_parallel_threads){
  input =list.files(GATK_input, full.names = T, pattern = GATK_suffix)
  if(!dir.exists(GATK_outdir)) {
    dir.create(GATK_outdir)
  } else{
      print("GATK Directories already exists, will override entries!!!!!")
    unlink(GATK_outdir,recursive=TRUE)
    dir.create(GATK_outdir)
  }
  if(!file.exists(paste(ref_fasta,".fai"))){
    cmd1 = paste(samtools,"faidx", ref_fasta)
    system(cmd1)
  } else{
    print(paste("GATK_CigarSplitter: Check path for",ref_fasta,"as SAMtools can't make index"))
  }
  if(!file.exists(paste0(gsub(".fa","",ref_fasta),".dict"))){
    cmd2 = paste(java,"-jar", picard,"CreateSequenceDictionary", gatkparams1, "-O",paste0(gsub(".fa","",ref_fasta),".dict"))
    system(cmd2)
  } else{
    print(paste("GATK_CigarSplitter: Check path for",ref_fasta,"as GATK can't make .dict file"))
  }
  registerDoMC(GATK_parallel_threads)
  foreach(j = 1:length(input)) %dopar% {
    cmd3 <- paste(gatk, "SplitNCigarReads",
                  "-I",input[j],
                  "-O", paste0(GATK_outdir,"/",gsub(".rmdup",".GATKfilt",basename(input[j]))),
                  gatkparams1)
    system(cmd3)
  }
}

GATK_BaseRecalibratorion <- function(GATK_outdir, GATK_suffix){
  input =list.files(GATK_outdir, full.names = T, pattern = GATK_suffix)
  registerDoMC(GATK_parallel_threads)
  foreach(j = 1:length(input)) %dopar% {
    cmd1 <- paste(gatk, "BaseRecalibrator",
                 "-I",input[j],
                 "-O", paste0(gsub(GATK_suffix,"",input[j]),"tab"),
                 gatkparams2)
    system(cmd1)
  }
}

GATK_ApplyBQSR <- function(GATK_outdir, GATK_suffix1, GATK_suffix2, plot ="yes"){
  input =list.files(GATK_outdir, full.names = T, pattern = GATK_suffix1)
  if(!require(gsalib)){
    install.packages("gsalib")
    library(gsalib)
  }
  registerDoMC(GATK_parallel_threads)
  foreach(j = 1:length(input)) %dopar% {
    if(!file.exists(gsub(GATK_suffix1,GATK_suffix2,input[j]))){
      cmd1 <- paste(gatk, "ApplyBQSR",
                    "-I",input[j],
                    "--bqsr-recal-file", paste0(gsub(GATK_suffix1,"",input[j]),"tab"),
                    "-O",gsub(GATK_suffix1,GATK_suffix2,input[j]),
                    gatkparams1)
      system(cmd1)
    }
    if(plot=="yes"){
      cmd2 <- paste(gatk, "AnalyzeCovariates",
                    "-bqsr",paste0(gsub(GATK_suffix1,"",input[j]),"tab"),
                    "-plots",paste0(gsub(GATK_suffix1,"",input[j]),"pdf"))
      system(cmd2)
    }
  }
}

GATK_HaploTypeCaller <- function(GATK_outdir, GATK_suffix2){
  input =list.files(GATK_outdir, full.names = T, pattern = GATK_suffix2)
  registerDoMC(GATK_parallel_threads)
  foreach(j = 1:length(input)) %dopar% {
    if(!file.exists(GATK_suffix2,"vcf.gz",input[j])){
      cmd1 <- paste(gatk, "HaplotypeCaller",
                    "-I",input[j],
                    "-O",gsub(GATK_suffix2,"vcf.gz",input[j]),
                    gatkparams3)
      system(cmd1)
    }
  }
}

GATK_VariantFilteration <- function(GATK_outdir){
  input =list.files(GATK_outdir, full.names = T, pattern = GATK_suffix3)
  registerDoMC(GATK_parallel_threads)
  foreach(j = 1:length(input)) %dopar% {
    if(!file.exists(GATK_suffix2,"vcf.gz",input[j])){
      cmd1 <- paste(gatk, "VariantFiltration",
                    "-V",input[j],
                    "-O",gsub(GATK_suffix3,GATK_suffix4,input[j]),
                    gatkparams4)
      system(cmd1)
    }
  }
}

GATK_SPLIT <- function(GATK_outdir, GATK_suffix4 ){
  input =list.files(GATK_outdir, full.names = T, pattern = paste0(GATK_suffix4,"$"))
  registerDoMC(GATK_parallel_threads)
  foreach(j = 1:length(input)) %dopar% {
    cmd1 <- paste(gatk, "SelectVariants",
                  "-V", input[j],
                  "--select-type-to-include SNP",
                  "-O", gsub(GATK_suffix4, "SNP.vcf.gz", input[j]))
    system(cmd1)
    cmd2 <- paste(gatk, "SelectVariants",
                  "-V", input[j],
                  "--select-type-to-include INDEL",
                  "-O", gsub(GATK_suffix4, "INDEL.vcf.gz", input[j]))
    system(cmd2)
  }
}

GATK_snpEff <- function(GATK_outdir, GATK_suffix4 ){
  input =list.files(GATK_outdir, full.names = T, pattern = paste0(GATK_suffix4,"$"))
  registerDoMC(GATK_parallel_threads)
  foreach(j = 1:length(input)) %dopar% {
    cmd1 =paste(java,"-jar", paste0(snpEff,"/snpEff.jar ann -v -csvStats ",gsub(GATK_suffix4, "SNP.ann.stats.csv", input[j]),
                                    " -o gatk GRCh37.75"), gsub(GATK_suffix4, "SNP.vcf.gz", input[j]) ,
                "| /software/pindel/htslib/bgzip >", gsub(GATK_suffix4, "SNP.ann.vcf.gz", input[j]))
    system(cmd1)
    system(paste("/software/pindel/htslib/tabix -p vcf", gsub(GATK_suffix4, "SNP.ann.vcf.gz", input[j])))
    cmd2 <- paste(java, "-jar",picard,"CollectVariantCallingMetrics",
                  "-DBSNP /home/hlchen/Conor/references/GATK_references/Homo_sapiens_assembly38.chr1_22.dbsnp138.vcf",
                  "-I",gsub(GATK_suffix4, "SNP.ann.vcf.gz", input[j]), 
                  "-O", gsub(GATK_suffix4, "SNP.ann.eval", input[j]))
    system(cmd2)
    cmd3 <- paste(gatk, "VariantAnnotator", gatkparams1,
                  "--variant", gsub(GATK_suffix4, "SNP.vcf.gz", input[j]),
                  "--resource", gsub(GATK_suffix4, "SNP.ann.vcf.gz", input[j]),
                  "-L",gsub(GATK_suffix4, "SNP.vcf.gz", input[j]),
                  '-O', gsub(GATK_suffix4, "SNP.GATK.ann.vcf.gz", input[j]))
    system(cmd3)
    ##For INDELs
    cmd4 =paste(java,"-jar", paste0(snpEff,"/snpEff.jar ann -v -csvStats ",gsub(GATK_suffix4, "INDEL.ann.stats.csv", input[j]),
                                    " -o gatk GRCh37.75"), gsub(GATK_suffix4, "INDEL.vcf.gz", input[j]) ,
                "| /software/pindel/htslib/bgzip >", gsub(GATK_suffix4, "INDEL.ann.vcf.gz", input[j]))
    system(cmd4)
    system(paste("/software/pindel/htslib/tabix -p vcf", gsub(GATK_suffix4, "INDEL.ann.vcf.gz", input[j])))
    cmd5 <- paste(java, "-jar",picard,"CollectVariantCallingMetrics",
                  "-DBSNP /home/hlchen/Conor/references/GATK_references/Mills_and_1000G_gold_standard.chr1_22.indels.hg38.vcf",
                  "-I",gsub(GATK_suffix4, "INDEL.ann.vcf.gz", input[j]), 
                  "-O", gsub(GATK_suffix4, "INDEL.ann.eval.txt", input[j]))
    system(cmd5)
    cmd6<- paste(gatk, "VariantAnnotator", gatkparams1,
                  "--variant", gsub(GATK_suffix4, "INDEL.vcf.gz", input[j]),
                  "--resource", gsub(GATK_suffix4, "INDEL.ann.vcf.gz", input[j]),
                  "-L",gsub(GATK_suffix4, "INDEL.vcf.gz", input[j]),
                  '-O', gsub(GATK_suffix4, "INDEL.GATK.ann.vcf.gz", input[j]))
    system(cmd6)
    
  }
}
