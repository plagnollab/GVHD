library(oligo)
library(limma)


all.files <- list.celfiles("/cluster/project8/vyp/Winship_GVHD/rawData/Langerhans_cells", full.names = TRUE, recursive = TRUE)
annotations <- read.csv('/cluster/project8/vyp/Winship_GVHD/oldFiles/data_Hannah_Shorrock/MoGene-2_0-st-v1.design-time.20120706.transcript.csv', skip = 9)


for (choice in c("BMT_male_vs_BMT_female", "BMT_female_vs_naive_male", "BMT_male_vs_naive_male")){
  if (choice== 'BMT_male_vs_BMT_female'){
    selected.files <- subset(all.files,grepl(pattern = "E9.|F9.",all.files))
  }
  if (choice== 'BMT_female_vs_naive_male'){
    selected.files <- subset(all.files,grepl(pattern = "F9.|G9.",all.files))
  }
  
  
  if (choice== 'BMT_male_vs_naive_male'){
    selected.files <- subset(all.files,grepl(pattern = "E9.|G9.",all.files))
  }
  
  set.probes <- 'core' 
  
  affyexpression <- read.celfiles(selected.files)
  summaries <- rma(affyexpression, target = set.probes)
  design <- model.matrix(~factor(c(1,1,1,2,2,2)))
  fit <- lmFit(summaries,design)
  ebayes <- eBayes(fit)
  
  lod <- -log10(ebayes[["p.value"]][,2])
  mtstat<- ebayes[["t"]][,2]
  tab <- topTable(ebayes, coef=2, adjust="fdr", genelist = ebayes$genes, n = nrow(ebayes))
  IDS <- rownames(tab)
  IDS <- as.vector(IDS)
  tab$new.col <- IDS
  colnames(tab)[7] <- c('probeset_id')
  
  tab.annotated <- merge(tab, annotations, by = 'probeset_id', all.x = TRUE)
  tab.annotated <- tab.annotated[ order(tab.annotated$P.Value, decreasing = FALSE), ]
  tab.annotated <- subset(tab.annotated, seqname != '---')
  significant <- subset(tab.annotated,P.Value < 5.0e-02)
  
  write.csv(x= tab.annotated, file = paste('array_analysis/syn_allo_bmt/results/', choice,'.csv',sep =""))
  write.csv(x= significant, file = paste('array_analysis/syn_allo_bmt/results/', choice,'_significant.csv',sep =""))
  
}


