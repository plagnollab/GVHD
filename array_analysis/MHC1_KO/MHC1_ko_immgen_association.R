
for (choice in c("D7WT_VS_D7KO", "D42WT_VS_D42KO")){

if (choice == "D7WT_VS_D7KO"){
  ext <- "_ALLANALYSIS_TM008"
}

if (choice == "D42WT_VS_D42KO"){
    ext <- "_ALLANALYSIS_TM006"
  }

if (!file.exists("array_analysis/MHC1_KO/results/immgen_association")) dir.create("array_analysis/MHC1_KO/results/immgen_association")

data <- read.table("/cluster/project8/vyp/Winship_GVHD/claire/data_files/biomart_annotations_mmusculus_gene_ensembl.tab", sep = "\t", comment.char = "", header = TRUE)

my.table <- data.frame (gene.name = unique(tolower(data$external_gene_name)),
	 DE = NA)

########## list of DE genes
DE <- read.csv(paste("array_analysis/MHC1_KO/results/DE_", choice, ext,  ".csv", sep = ""), sep = ",", header = TRUE)
#DE$gene_assignment <- gsub(DE$gene_assignment, pattern = ";.*", replacement = "")   ## optoinal line to remove multiple mapping probes
DE.gene.list <- unique(strsplit(paste(DE$gene_assignment, collapse = " // "), split = " // ")[[1]])
DE.gene.list <- gsub(DE.gene.list, pattern = "^ ", DE.gene.list, replacement = "")
DE.gene.list <- tolower(DE.gene.list [ DE.gene.list != "" ])

########## list of non-DE genes
nonDE <- read.csv(paste("array_analysis/MHC1_KO/results/non_DE_", choice, ext, ".csv", sep = ""), sep = ",", header = TRUE)
#nonDE$gene_assignment <- gsub(nonDE$gene_assignment, pattern = ";.*", replacement = "")   ## optoinal line to remove multiple mapping probes
nonDE.gene.list <- unique(strsplit(paste(nonDE$gene_assignment, collapse = " // "), split = " // ")[[1]])
nonDE.gene.list <- gsub(nonDE.gene.list, pattern = "^ ", nonDE.gene.list, replacement = "")
nonDE.gene.list <- tolower(nonDE.gene.list [ DE.gene.list != "" ])



my.table <- subset(my.table, my.table$gene.name %in% c(DE.gene.list, nonDE.gene.list))
my.table$DE <- my.table$gene.name %in% DE.gene.list

modules <- read.table("/cluster/project8/vyp/Winship_GVHD/claire/data_files/gene_assignment.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
all.modules <- unique(modules$Fine.module)
table.results <- data.frame(module = all.modules, pvalue = NA, neg_log_pvalue = NA)

for (i in 1:nrow(table.results)) {
   genes.in.module <- tolower(subset(modules, Fine.module == table.results$module[ i ], Gene, drop = TRUE))
  my.table$in.module <- my.table$gene.name %in% genes.in.module
  table1 <- table(my.table$DE, my.table$in.module)
  my.test <- fisher.test(table(my.table$DE, my.table$in.module))
  #print(my.test)
    message("Module ",table.results$module[ i ] , " has a P-value of ", my.test$p.value) 
  if (table.results$module[i] == 3){
    print(table1)
  }
    
  table.results$pvalue[ i ] <- my.test$p.value
  table.results$neg_log_pvalue[i] <- -log10(my.test$p.value)
}
test1 = cbind(table.results$module, table.results$neg_log_pvalue)

write.table(table.results, file = paste("array_analysis/MHC1_KO/results/immgen_association/", choice, "_ImmGen_fine__module_association_Fishers_exact_test_Pvalues.txt", sep = ""), sep = "\t")
colnames(modules)[1] <- "gene.name"
modules$gene.name <- tolower(modules$gene.name)
final.table = merge(my.table, modules, by = "gene.name", all.x = TRUE)
final.table$module = NULL
write.table(final.table, file = paste("array_analysis/MHC1_KO/results/immgen_association/", choice, "_final_table_fine.txt", sep = ""), sep = "\t")

########## plot #############
#########       #############


xmin = 0
xmax = 334
ymin = 0
ymax = 12

pdf(paste("array_analysis/MHC1_KO/results/immgen_association/", choice, "_fine_module_association_graph.pdf", sep = ""))

plot(test1, type = "n", main = paste("Fisher's Exact Test p-values for association \n of the ", choice, " differential gene \n expression dataset with the ImmGen fine modules", sep = ""), xlab = "Fine module number", ylab = "-log10(p-value)", xlim = c(xmin, xmax), ylim = c(ymin, ymax), xaxs = "i", yaxs = "i")
axis(side=1, at=seq(0, 334, by = 40))
axis(side=2, at=seq(0, 12, by = 2))
lines(test1, type = "h", col = "blue")

dev.off()

}
