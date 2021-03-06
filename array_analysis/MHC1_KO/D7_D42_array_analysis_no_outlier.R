library(oligo)
library(limma)
library(lattice)

all.files <- list.celfiles('/cluster/project8/vyp/Winship_GVHD/claire/data_files/Teresa_microarray', full.names = TRUE, recursive = TRUE)
annotations <- read.csv('/cluster/project8/vyp/Winship_GVHD/oldFiles/data_Hannah_Shorrock/MoGene-2_0-st-v1.design-time.20120706.transcript.csv', skip = 9)

set.probes <- 'core'  ##extended, full are other options

affyexpression <- read.celfiles(all.files)
summaries <- rma(affyexpression, target = set.probes)

pheno.data.main <- pData(summaries)


#for (choice in c("TM008wt_vs_TM008ko")){
for (choice in c("TM008wt_vs_TM008ko", "TM006wt_vs_TM006ko")){

  if (choice== 'TM008wt_vs_TM008ko'){
    loc.summaries <- summaries[, !grepl(pattern = "TM008_ko4", row.names(pheno.data.main)) & grepl(pattern = "TM008", row.names(pheno.data.main) )]  ## this is the subsetted file, I remove the outlier here
  }
  
  if (choice== 'TM006wt_vs_TM006ko'){
    loc.summaries <- summaries[, grepl(pattern = "TM006", row.names(pheno.data.main) )]  ## this is the subsetted file
  }

  design1 <- model.matrix(~factor(ifelse( grepl(pattern = "wt", row.names(pData(loc.summaries))), 1, 2)))
  fit1 <- lmFit(loc.summaries,design1)
  
  ##design0 <- model.matrix(~rep(1, times = nrow(pData(loc.summaries))))[,1]
  ##fit0 <- lmFit(loc.summaries,design0)

  ebayes <- eBayes(fit1, proportion=0.55)
  lod <- -log10(ebayes[["p.value"]][,2])
  mtstat<- ebayes[["t"]][,2]
  tab <- topTable(ebayes, coef=2, adjust="fdr", genelist = ebayes$genes, n = nrow(ebayes))
  IDS <- rownames(tab)
  IDS <- as.vector(IDS)
  tab$new.col <- IDS
  ##print(IDS)
  colnames(tab)[7] <- c('probeset_id')
  tab.annotated <- merge(tab, annotations, by = 'probeset_id', all.x = TRUE)
  tab.annotated <- tab.annotated[ order(tab.annotated$P.Value, decreasing = FALSE), ]
  tab.annotated <- subset(tab.annotated, seqname != '---')
  ##print(tab.annotated)
  significant <- subset(tab.annotated,P.Value < 5.0e-03)
                                        # write.csv(x= tab, file = "Teresa_test1.csv",sep = "\n")
  write.csv(x= tab.annotated, file = paste('/cluster/project8/vyp/Winship_GVHD/claire/GVHD/array_analysis/MHC1_KO/results/', choice,'.csv',sep =""))
}
