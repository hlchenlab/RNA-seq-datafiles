

if( !require("ggplot2")){
  BiocManager::install("ggplot2")
}
if (!require("ChIPseeker")){
  BiocManager::install("ChIPseeker")
}
if (!require("dplyr")){
  BiocManager::install("dplyr")
}
if (!require("GenomicRanges")){
  BiocManager::install("GenomicRanges")
}
if (!require("clusterProfiler")){
  BiocManager::install("clusterProfiler", ask = F)
}
if (!require("ReactomePA")){
  BiocManager::install("ReactomePA", ask = F)
}
if (!require("org.Hs.eg.db")){
  BiocManager::install("org.Hs.eg.db", ask = F)
}

genes = lapply(1:ncol(genesets.up), function(i) {x =rownames(genesets.down[genesets.down[,i]==1,])
x =stack(mget(x, org.Hs.egENSEMBL2EG, ifnotfound = NA))[,1]
x =na.omit(x)
x = data.frame(Entrez = x, Studies =colnames(genesets.down)[i], `Differential Expression` = "Downregulated")
return(x)})
genes = do.call(rbind.data.frame, genes)

genes.down = genes
genes = lapply(1:ncol(genesets.up), function(i) {x =rownames(genesets.up[genesets.up[,i]==1,])
x =stack(mget(x, org.Hs.egENSEMBL2EG, ifnotfound = NA))[,1]
x =na.omit(x)
x = data.frame(Entrez = x, Studies =colnames(genesets.up)[i], `Differential Expression` = "Upregulated")
return(x)})
genes = do.call(rbind.data.frame, genes)

genes.up = genes
genes = rbind.data.frame(genes.up, genes.down)
genes$Studies =factor(genes$Studies, levels = studies[c(12:14,1:11,15)])

## PSI enriched genes
genes = stack(mget(psi1$ensemblID, org.Hs.egENSEMBL2EG, ifnotfound = NA))
m = match(psi1$ensemblID, genes$ind)
f.a =!is.na(m)
f.t =m[f.a]
psi1$Entrez =NA
psi1[f.a,]$Entrez = genes[f.t,]$values
genes = psi1[psi1$T.test.p.value <= 0.05,c("Entrez","Studies","pos")]
genes$Studies =factor(genes$Studies, levels = studies)


compGO <- compareCluster(Entrez~`Differential Expression`,
                         data = genes,
                         fun           = 'enrichGO',
                         pvalueCutoff  = 0.1,
                         OrgDb='org.Hs.eg.db',
                         keyType ="ENTREZID",
                         ont ="BP",
                         minGSSize = 10,
                         maxGSSize = 500)

compGO_simp =simplify(compGO,
                      cutoff=-0.1, 
                      by="p.adjust", 
                      select_fun=min, 
                      measure = "Wang",
                      semData = NULL)

compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = 'enrichKEGG',
                           pvalueCutoff  = 0.1,
                           minGSSize = 10,
                           maxGSSize = 500)

compRA <- compareCluster(geneCluster   = genes,
                         fun = 'enrichPathway',
                         pvalueCutoff  = 0.1,
                         minGSSize = 10,
                         maxGSSize = 500)

compDO <- compareCluster(geneCluster   = genes,
                         fun = 'enrichDO',
                         pvalueCutoff  = 0.1,
                         minGSSize = 10,
                         maxGSSize = 500)

keyterms =list("inflammatory response","leukocyte differentiation", "positive regulation of cell death","I-kappaB kinase/NF-kappaB signaling",
               "cellular response to interferon-gamma","cellular response to type I interferon","tumor necrosis factor-mediated signaling pathway",
               "chemokine-mediated signaling pathway", "mitotic nuclear division", "cytoplasmic translation","meiotic cell cycle","regulation of ion transport",
               "cilium movement", "leukocyte migration","negative regulation of viral process","hematopoietic or lymphoid organ development",
               "sensory perception of light stimulus", "cytoplasmic translation")


p6 =dotplot(compGO, showCategory =1) +
  theme_bw() +
  xlab("Studies") +
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
  coord_flip() +
  scale_size_continuous(name ="GeneRatio",range = c(0,5)) +
  guides(colour =guide_colorbar(title = "P.adjust", barwidth = unit(3, "mm"), barheight = unit(15, "mm"))) +
  theme(legend.title = element_text(size =7, family = "sans", face = "bold"),
        legend.text = element_text(size =6, family = "sans"),
        legend.spacing.x = unit("1", "mm"),
        legend.box.spacing =unit("1", "mm"),
        legend.box.just = "left",
        plot.title = element_text(size =7, family = "sans", face = "bold", hjust =0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.y = element_text(size =6 ,family = "sans"),
        axis.text.x = element_text(size =6 ,family = "sans", vjust = 1, angle = 45, hjust = 1),
        axis.title = element_text(size =7, family = "sans"),
        strip.text.x = element_text(size = 7, family = "sans")) +
  facet_grid(rows = vars(pos), space ="free", scales = 'free_y')

p6$data$pos = factor(p6$data$pos, levels = c("Upregulated","Downregulated"))
lvls = levels(p6$data$Description)
p6$data$Description = as.character(p6$data$Description)
df = table(p6$data[p6$data$Description %in% lvls,]$Description)
df =df[order(df, decreasing = T)]
p6$data$Description = factor(as.character(p6$data$Description), levels = names(df))
p6$data$Cluster = factor(matrix(unlist(strsplit(as.character(p6$data$Cluster), "[.]")), byrow = T, ncol = 2)[,1], exclude = "", levels = rev(c(paste0(studies))))
p6$guides$colour$frame.colour = "black"
  p6$guides$colour$ticks.colour = "black"
    
  g <- ggplot_gtable(ggplot_build(p6))
  
  stripr <- which(grepl('strip', g$layout$name))
  
  fills <- c("darkolivegreen3","brown2")
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  
  
  #### GSEA across clusters ####
  
  
  genes = lapply(1:max(clusters.up$unmergedLabels.mod), function(i) {x =clusters.up[clusters.up$unmergedLabels.mod == i,]$V1
  x =stack(mget(x, org.Hs.egENSEMBL2EG, ifnotfound = NA))[,1]
  x =na.omit(x)
  x = data.frame(Entrez = x, Clusters =i)
  return(x)})
  genes = do.call(rbind.data.frame, genes)
  genes.up = genes
  
  genes = lapply(1:max(clusters.down$unmergedLabels.mod), function(i) {x =clusters.down[clusters.down$unmergedLabels.mod == i,]$V1
  x =stack(mget(x, org.Hs.egENSEMBL2EG, ifnotfound = NA))[,1]
  x =na.omit(x)
  x = data.frame(Entrez = x, Clusters =i)
  return(x)})
  genes = do.call(rbind.data.frame, genes)
  genes.down = genes
  
  genes = rbind.data.frame(genes.up, genes.down)
  
  compGO <- compareCluster(Entrez~Clusters,
                           data = genes,
                           fun           = 'enrichGO',
                           pvalueCutoff  = 0.1,
                           OrgDb='org.Hs.eg.db',
                           keyType ="ENTREZID",
                           ont ="BP",
                           minGSSize = 20,
                           maxGSSize = 1000)
  
  p6 <- dotplot(compGO, showCategory =10) +
    theme_bw() +
    xlab("Clusters") +
    scale_y_discrete(guide = guide_axis(check.overlap = TRUE), labels = function(x) stringr::str_wrap(x, width = 30))+
    coord_flip() +
    scale_size_continuous(name ="GeneRatio",range = c(0,5)) +
    guides(colour =guide_colorbar(title = "P.adjust", barwidth = unit(4, "mm"), barheight = unit(15, "mm"))) +
    theme(legend.title = element_text(size =7, family = "sans", face = "bold"),
          legend.text = element_text(size =6, family = "sans"),
          legend.spacing.x = unit("1", "mm"),
          legend.box.spacing =unit("1", "mm"),
          legend.box.just = "left",
          plot.title = element_text(size =7, family = "sans", face = "bold", hjust =0.5),
          axis.text.y = element_text(size =6 ,family = "sans"),
          axis.text.x = element_text(size =6 ,family = "sans", vjust = 1.05, angle = 45, hjust = 1),
          axis.title = element_text(size =7, family = "sans"),
          strip.text.x = element_text(size = 7, family = "sans"))
  
  p6$data$Cluster = factor(matrix(unlist(strsplit(as.character(p6$data$Cluster), "\n")), byrow = T, ncol = 2)[,1], exclude = "", levels = rev(1:max(clusters$unmergedLabels.mod)))
  p6$guides$colour$frame.colour = "black"
    p6$guides$colour$ticks.colour = "black"
      
    