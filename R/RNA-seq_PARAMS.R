#######################################################

#################RNA-Seq Analysis modules (yes/no)#############

Do_FasQC_1 ="yes"
#Or
Do_FastQC_Trim = "no"

Do_STAR ="yes"
Do_STAR_TWO_PASS ="yes" #Do_STAR must be "yes" for two pass method

RM_duplicates_from_BAMs = "yes"
Do_chromosome_subset ="no"

Do_variant_call ="yes"
Do_DESeq2 = "yes"

###################Essential Paths####################
java <- "/software/java/13.0.2/bin/java"
samtools <- '/software/samtools/1.14/bin/samtools'
picard <- "/software/Picard/2.25.2/picard.jar"
gatk <- "/software/gatk/4.2.5.0/gatk"
STAR <- '/software/STAR/2.7.9a/bin/STAR'
fastqc.path = "/software/FastQC/0.11.9/fastqc"
bamtools ="/software/bamtools/2.5.1/bin/bamtools"
snpEff <- "/home/conor93/snpEff"
zip7 <- "/software/7-Zip/16.02/bin/7z"
work_dir <- "RNA-Seq-datafiles/R"

home_dir ="/home/hlchen/Conor/Studies/IAV" 
studies = "GSE156060"
ref_fasta = "/home/hlchen/Conor/references/human_wsn.fa"
IDX = "/home/hlchen/Conor/references/Gencode_v39/STAR_Human"
inputdir = "fastq"
################RNA-seq Parameters#####################

if(Do_FasQC_1 =="yes"){ 
  source(paste0(work_dir,"/FastQC.R"))
###################FastQC#############################
  fastqc_threads = 8
  fastqc_out ="QC"
##############Running operation#######################
  setwd(home_dir)
  
  environ_set_up()
  
  for (i in studies) {
    setwd(i)
    if(Do_FasQC_1 == "yes"){
      FastQC_analysis(inputdir ='fastq',fastqc_out, fastqc_threads,report_out ="yes") 
    }
    if(Do_FastQC_Trim == "yes"){
      Trim_Fastq(inputdir ='fastq', trimoutdir = 'FastqTrim')
    }
    setwd("..")
  }
}
######################################################

###############Post-Alignment-Processing##################

if(RM_duplicates_from_BAMs == "yes"){
  source(paste0(work_dir,"/STARscript.R"))
  ###############Deduplication of BAMs######################
  
  rmdup_threads = 5
  temp = "tempdir"
  rmdupdir ="RMDups"
  rmdup_params <- paste("--remove-sequencing-duplicates true",
                        paste0("--conf \'spark.executor.cores=",rmdup_threads,"\'"),
                        paste0("--conf \'spark.local.dir=",temp,"\'"),
                        paste0("--conf \'spark.executor.instances=10\'"),
                        paste0("--conf \'spark.executor.memory=8G\'"))
  
} 
##########################################################

#####################Subsetting BAM files ################

split_indir = rmdupdir
split_outdir = "Splitfiles"
split_suffix =".bam$"
split_parallel_threads =10

if(Do_chromosome_subset == "unmapped"){ 
  subsetparams =paste("-b -f 4")
} else if(Do_chromosome_subset == "ByChromosome"){
  subsetparams =paste("-reference")
}
##########################################################

if(Do_STAR =="yes"){
  source(paste0(work_dir,"/STARscript.R"))

#####Star Alignment Parameters (Change if necessary)####
  STARthreads = 10
  star_outdir ="BAM"
  BAM_output ="BAM SortedByCoordinate"
  
  STARParameters = paste('--readFilesCommand gunzip -c',
                         "--runThreadN", STARthreads,
                         '--genomeDir', IDX,
                         '--outSAMtype',BAM_output,
                         '--quantMode GeneCounts')

  if(Do_STAR_TWO_PASS == "yes"){
    star_outdir ="TwoPass"
    star_add_params <- paste('--outSAMstrandField intronMotif --outSAMunmapped Within --twopassMode Basic  --chimSegmentMin 12', 
                             "--chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10",
                             "--alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5", 
                             "--chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20",
                             "--chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right", 
                             "--alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30")
    STARParameters =paste(STARParameters, star_add_params)
  }
  
  ##############Running operation#######################
  setwd(home_dir)
  
  environ_set_up()
  
  for (i in studies) {
    setwd(i)
    if(Do_STAR == "yes"){
      mapRNAStudies(inputdir,study=i, IDX, star_outdir, parallel_threads =1) 
    }
    if(RM_duplicates_from_BAMs == "yes"){
      dir.create(temp)
      Add_ReadGroups(star_outdir,STARthreads)
      rmdupPCR(star_outdir,rmdupdir)
      unlink(temp, recursive=TRUE)
      print(paste("for", i , "Reads deduplicated and Read groups assigned successfully!"))
    }
    if(Do_chromosome_subset =="unmapped" | Do_chromosome_subset == "ByChromosome"){
      SplitBAM(split_indir,split_outdir)
    }
    setwd("..")
  }
}
##########################################################

##################Variant Calling#########################

GATK_parallel_threads =6
GATK_outdir ="GATK"
GATK_suffix1 ="GATKfilt.bam$"
GATK_suffix2 ="BQSR.bam"
GATK_suffix3 = "vcf.gz$"
GATK_suffix4 = "filt.vcf.gz"
GATK_suffix5 = "ann.vcf.gz"
GATK_Intervals = c("chr10")
gatkparams1 = paste("-R", ref_fasta)

##May need to format knowsite vcf files to only relevant chromosomes


gatkparams2 = paste("-R", ref_fasta, "-known-sites /home/hlchen/Conor/references/GATK_references/Homo_sapiens_assembly38.dbsnp138.vcf",
                    "-known-sites /home/hlchen/Conor/references/GATK_references/Mills_and_1000G_gold_standard.indels.hg38.vcf")
gatkparams3 = paste("-R", ref_fasta,paste("-L",GATK_Intervals))
gatkparams4 <- paste("-R", ref_fasta, "-window 35 -cluster 3 -filter-name \"Filter1\" -filter-expression \"QD < 2.0\" -filter-name \"Filter2\" -filter-expression \"FS > 30.0\"")


if(Do_variant_call =="yes"){
  if(RM_duplicates_from_BAMs == "yes"){
    GATK_input = rmdupdir
  } else {
    GATK_input = star_outdir
  }
  
  ##############Running operation#######################
  setwd(home_dir)
  
  environ_set_up()
  
  for (i in studies) {
    setwd(i)
    if(Do_variant_call =="yes"){
      GATK_CigarSplitter(GATK_input, GATK_outdir,GATK_suffix, GATK_parallel_threads)
      GATK_BaseRecalibratorion(GATK_outdir, GATK_suffix)
      GATK_ApplyBQSR(GATK_outdir, GATK_suffix1,GATK_suffix2)
      GATK_HaploTypeCaller(GATK_outdir, GATK_suffix2)
      GATK_SPLIT(GATK_outdir, GATK_suffix4 )
      GATK_snpEff(GATK_outdir, GATK_suffix4 )
    }
    setwd(home_dir)
  }
}



