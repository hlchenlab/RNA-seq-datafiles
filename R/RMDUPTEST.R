rmdupPCR <- function(star_outdir,rmdupdir){
  input =list.files(star_outdir, full.names = T, pattern = ".out.bam$")
  if(!dir.exists(rmdupdir)) dir.create(rmdupdir)
  for(j in 1:length(input)){
    cmd1 <- paste(gatk, "MarkDuplicatesSpark",
                  "-I",input[j],
                  "-O", paste0(rmdupdir,"/",gsub("Aligned.sortedByCoord.out", ".rmdup",basename(input[j]))),
                  "-M", paste0(rmdupdir,"/",gsub("Aligned.sortedByCoord.out.bam", "_dedup_metrics.txt",basename(input[j]))),
                  rmdup_params)
    system(cmd1)
    cmd2 <- paste(java,"-jar",picard, "CollectAlignmentSummaryMetrics",
                  '-I',paste0(rmdupdir,"/",gsub("Aligned.sortedByCoord.out", ".rmdup",basename(input[j]))),
                  "-R",ref_fasta,
                  "-O", paste0(rmdupdir,"/",gsub("Aligned.sortedByCoord.out.bam", "_alignment_metrics.txt",basename(input[j]))))
    system(cmd2)
    cmd3 <- paste(java,"-jar",picard, "CollectInsertSizeMetrics",
                  '-I',paste0(rmdupdir,"/",gsub("Aligned.sortedByCoord.out", ".rmdup",basename(input[j]))),
                  "-H",paste0(rmdupdir,"/",gsub("Aligned.sortedByCoord.out", "insert_size_histogram.pdf",basename(input[j]))),
                  "-O", paste0(rmdupdir,"/",gsub("Aligned.sortedByCoord.out.bam", "_insert_metrics.txt",basename(input[j]))))
    system(cmd3)
    cmd4 <- paste(samtools, "depth -a", paste0(rmdupdir,"/",gsub("Aligned.sortedByCoord.out", ".rmdup",basename(input[j]))), ">",
                  paste0(rmdupdir,"/",gsub("Aligned.sortedByCoord.out.bam", "_depth_out.txt",basename(input[j]))))
    system(cmd4)
  }
}

