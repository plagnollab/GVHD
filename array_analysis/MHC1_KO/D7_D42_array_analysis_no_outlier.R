library(oligo)
library(limma)
##library(affy)

all.files <- list.celfiles('/cluster/project8/vyp/Winship_GVHD/claire/data_files/Teresa_microarray', full.names = TRUE, recursive = TRUE)
annotations <- read.csv('/cluster/project8/vyp/Winship_GVHD/oldFiles/data_Hannah_Shorrock/MoGene-2_0-st-v1.design-time.20120706.transcript.csv', skip = 9)

#for (choice in c("TM006wt_vs_TM006ko")){
for (choice in c("TM008wt_vs_TM008ko")){

  if (choice== 'TM008wt_vs_TM008ko'){
    selected.files <- subset(all.files,grepl(pattern = "TM008_wt|TM008_ko",all.files))
    design <- model.matrix(~factor(c(2,2,1,1,1)))
  }
  
  if (choice== 'TM006wt_vs_TM006ko'){
    selected.files <- subset(all.files,grepl(pattern = "TM006_wt|TM006_ko",all.files))
    design <- model.matrix(~factor(c(2,2,2,1,1,1)))
  }
  
## need to specify wild type and knock out in each case then run analysis

  set.probes <- 'core'  ##extended, full are other options
  
  affyexpression <- read.celfiles(selected.files)
  summaries <- rma(affyexpression, target = set.probes)

  ##Library - package mogene20stcdf not installed
  new.expr <- justRMA(filenames = basename(selected.files), celfile.path = "/cluster/project8/vyp/Winship_GVHD/claire/data_files/Teresa_microarray/")
  print(exprs(summaries)[head(subset(tab.annotated, P.Value == 1)$probeset_id),])
  
  fit <- lmFit(summaries,design)
  ebayes <- eBayes(fit)
  lod <- -log10(ebayes[["p.value"]][,2])
  mtstat<- ebayes[["t"]][,2]
  tab <- topTable(ebayes, coef=2, adjust="fdr", genelist = ebayes$genes, n = nrow(ebayes))
  IDS <- rownames(tab)
  IDS <- as.vector(IDS)
  tab$new.col <- IDS
  ##print(IDS)
  colnames(tab)[7] <- c('probeset_id')
                                        #  probes <- rownames(tab)
                                        #  print(probes)
                                        #  print(p.adj)
                                        #  test <- IDS[p.adj<0.05]
                                        #  print(test)
  tab.annotated <- merge(tab, annotations, by = 'probeset_id', all.x = TRUE)
  tab.annotated <- tab.annotated[ order(tab.annotated$P.Value, decreasing = FALSE), ]
  tab.annotated <- subset(tab.annotated, seqname != '---')
  ##print(tab.annotated)
  significant <- subset(tab.annotated,P.Value < 5.0e-03)
                                        # write.csv(x= tab, file = "Teresa_test1.csv",sep = "\n")
  write.csv(x= tab.annotated, file = paste('array_analysis/MHC1_KO/results/', choice,'.csv',sep =""))
  esetSel <- as.matrix(significant[1:100,]) 

  #pdf("array_analysis/MHC1_KO/figs/heatmap_DE.pdf",width = 9, height = 5)
  #par(las = 2)
  #heatmap(exprs(esetSel))
  #dev.off()
}
