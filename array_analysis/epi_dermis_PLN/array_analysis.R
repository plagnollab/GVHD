#As a first step, a good comparison would be  (1) epidermis versus epidermis DT (2)  dermis versus dermis DT and (3) peripheral LN versus peripheral LN +DT.

##1- look at different groupings
##2- different probe sets
##3- add batch as covariate
##4- redo PCA analysis within R



library(oligo)
library(limma)

annotations <- read.csv('/cluster/project8/vyp/Winship_GVHD/oldFiles/data_Hannah_Shorrock/Affymetrix/Mouse_exon_array_2.0.ST/MoGene-2_0-st-v1.design-time.20120706.transcript.csv', skip = 9)

all.files <- list.celfiles('/cluster/project8/vyp/Winship_GVHD/rawData/MataHari_model/MH_CEL_files/', full.names = TRUE, recursive = TRUE)


#choice <- 'epidermis_vs_epidermisDT'
choice <- 'dermis_vs_dermisDT'
#choice <- 'PLN_vs_PLNDT'

set.probes <- 'core'  ##extended, full are other options



for (choice in c('dermis_vs_dermisDT', 'PLN_vs_PLNDT')) {
  
  if (choice == 'epidermis_vs_epidermisDT') {
    selected.files <- sort(subset(all.files, grepl(pattern = 'D9|A9', all.files)))
    design <- model.matrix(~factor(c(1,1,1,2,2,2)))
  }
  
  if (choice == 'dermis_vs_dermisDT') {
    selected.files <- sort(subset(all.files, grepl(pattern = 'D10|A10', all.files)))
    design <- model.matrix(~factor(c(1,1,1,2,2,2)))
  }

  if (choice == 'PLN_vs_PLNDT') {
    selected.files <- sort(subset(all.files, grepl(pattern = 'D11|A11', all.files)))
    design <- model.matrix(~factor(c(1,1,1,2,2,2)))
  }

  
  affyExpressionFS <- read.celfiles(selected.files)
  Summaries <- rma(affyExpressionFS, target= set.probes)
  fit <- lmFit(Summaries, design)
  ebayes <- eBayes(fit)
  lod <- -log10(ebayes[["p.value"]][,2])
  mtstat<- ebayes[["t"]][,2]
  tab <- topTable(ebayes, coef=2, adjust="fdr", n = nrow(ebayes$genes))
  names(tab)[1] <- 'probeset_id'
  
  tab.annotated <- merge(tab, annotations, by = 'probeset_id', all.x = TRUE)
  tab.annotated <- tab.annotated[ order(tab.annotated$P.Value, decreasing = FALSE), ]
  tab.annotated <- subset(tab.annotated, seqname != '---')
  
  write.csv(x = tab.annotated, file = paste("/cluster/project8/vyp/Winship_GVHD/claire/GVHD/array_analysis/epi_dermis_PLN/results/", choice, '.csv', sep = ''))
}



