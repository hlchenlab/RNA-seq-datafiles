

##Running Code:

##Base Call

if(DO_Base_Call == "yes"){
  for (i in studies){
    cmd = paste("bcl2fastq -r 10 -p 10 -w 10 -i", i)
    system(cmd)
  }
  print("BCL2Fastq: Fastq files generated")
}

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

##FastQC
source("Alignment-params.txt")

setwd(home_dir)

environ_set_up()

for (i in studies) {
  setwd(i)
  if(Do_FasQC_1 == "yes"){
    FastQC_analysis(inputdir ='fastq',fastqc_out, fastqc_threads,report_out ="yes") 
  }
  setwd("..")
}

