##Title: TrimGalore and FastQC
##Author: Conor Cremin
##Company: University of Hong Kong

##Description:
###Requires sample directories to have the same format

environ_set_up <- function(){
  if(!require(dplyr)){
    print(paste("dplyr pkg not available, installing now and re-run"))
    BiocManager::install("dplyr")
    require(dplyr)
  }
  if(!require(fastqcr)){
    print(paste("fastqcr pkg not available installing and re-run"))
    install.packages("fastqcr")
    require(fastqcr)
  }
}

FastQC_analysis <- function(inputdir,fastqc_out, fastqc_threads,report_out ="yes"){
  if(!dir.exists('QC')) system('mkdir QC')
  fastqc(fq.dir =inputdir, 
         qc.dir = fastqc_out,
         threads = fastqc_threads, 
         fastqc.path = fastqc.path)
  if(report_out == "yes"){
    qc <- qc_aggregate(fastqc_out)
    save(qc, file = paste0(fastqc_out,"/Multi_sample_report.rdata"))
    tryCatch(modules = as.data.frame(qc_fails(qc, "module")[,"module"])[,1] %>% as.character(),
             function(e){
               print("FastQC_analysis: Unable to find failed modules for plotting")
             })
    qc_out <-  qc_read_collection(list.files(fastqc_out, full.names = T, pattern = ".zip$")[c(1:6,13:18)], 
                                  modules = modules,
                                  sample_names = gsub("_fastqc.zip","",list.files(fastqc_out, pattern = ".zip$"))[c(1:6,13:18)])
    save(qc_out, modules, file = paste0(fastqc_out,"/FastqcPlot.rdata"))
    for (i in modules) {
      pdf(paste0(getwd(),"/",i,".pdf"), height = 7, width = 7)
      plot(qc_plot_collection(qc_out, i) +theme_bw() +
             theme(strip.text = element_text(size = 5)))
      dev.off()
    }
  }
}

###
Trim_Fastq <- function(inputdir,trimoutdir){
  library(dplyr)
  trim = '/software/TrimGalore/0.4.5/trim_galore'
  cutapt ='/software/cutadapt/3.4/bin/cutadapt'
  files = list.files(inputdir, full.names = T, pattern = '.gz') %>% unlist() %>% cat() %>% capture.output()
  if(!dir.exists('QC')) system('mkdir QC')
  if(!dir.exists(trimoutdir)) dir.create(trimoutdir)
  if(names(i) =="paired"){
    params = paste('-j',threads,
                   '--path_to_cutadapt', cutapt,
                   '--paired', 
                   '--fastqc_args', dQuote("--outdir /QC/"), 
                   '--gzip -o',trimoutdir)
  } else {
    params = paste('-j',threads,
                   '--path_to_cutadapt', cutapt,
                   '-o',trimoutdir,
                   '--fastqc_args', dQuote("--outdir /QC/"), 
                   '--gzip')
  }
  cmd = paste(trim, params, files)
  print(paste('Running:', cmd))
  system(paste(cmd))
}

##Running operation
source("Alignment-params.txt")
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

