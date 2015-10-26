for (choice in c("BMT_male_vs_BMT_female", "BMT_female_vs_naive_male", "BMT_male_vs_naive_male")){

data <- read.table("/cluster/project8/vyp/Winship_GVHD/claire/data_files/biomart_annotations_mmusculus_gene_ensembl.tab", sep = "\t", comment.char = "", header = TRUE)


my.table <- data.frame (gene.name = unique(tolower(data$external_gene_name)),
                        DE = NA)

########## list of DE genes
DE <- read.csv(paste("/cluster/project8/vyp/Winship_GVHD/claire/results/syn_allo_bmt/results/DE_", choice, ".csv", sep = ""), sep = ",", header = TRUE)
#DE$gene_assignment <- gsub(DE$gene_assignment, pattern = ";.*", replacement = "") 
DE.gene.list <- unique(strsplit(paste(DE$gene_assignment, collapse = " // "), split = " // ")[[1]])
DE.gene.list <- gsub(DE.gene.list, pattern = " ", DE.gene.list, replacement = "")
DE.gene.list <- tolower(DE.gene.list [ DE.gene.list != "" ])

########## list of non-DE genes
nonDE <- read.csv(paste("/cluster/project8/vyp/Winship_GVHD/claire/results/syn_allo_bmt/results/non_DE_", choice,".csv", sep = ""), sep = ",", header = TRUE)
#nonDE$gene_assignment <- gsub(nonDE$gene_assignment, pattern = ";.*", replacement = "") 
nonDE.gene.list <- unique(strsplit(paste(nonDE$gene_assignment, collapse = " // "), split = " // ")[[1]])
nonDE.gene.list <- gsub(nonDE.gene.list, pattern = " ", nonDE.gene.list, replacement = "")
nonDE.gene.list <- tolower(nonDE.gene.list [ DE.gene.list != "" ])


my.table <- subset(my.table, my.table$gene.name %in% c(DE.gene.list, nonDE.gene.list))
my.table$DE <- my.table$gene.name %in% DE.gene.list

modules <- read.table("/cluster/project8/vyp/Winship_GVHD/claire/data_files/gene_assignment.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
all.modules <- unique(modules$Coarse.module)
table.results <- data.frame(module = all.modules, pvalue = NA)

for (i in 1:nrow(table.results)) {
  genes.in.module <- tolower(subset(modules, Coarse.module == table.results$module[ i ], Gene, drop = TRUE))
  #print(genes.in.module)
  my.table$in.module <- my.table$gene.name %in% genes.in.module
  #print(my.table$DE)
  table1 <- table(my.table$DE, my.table$in.module)
  #print(table1)
  my.test <- fisher.test(table(my.table$DE, my.table$in.module))
  #print(my.test)
  message("Module ",table.results$module[ i ] , " has a P-value of ", my.test$p.value) 
 # if (table.results$module[i] == 52){
  #  print(table1)
  # }
  table.results$pvalue[ i ] <- my.test$p.value
  table.results$neg_log_pvalue[i] <- -log10(my.test$p.value)
}

test1 = cbind(table.results$module, table.results$neg_log_pvalue)

write.table(table.results, file = paste("/cluster/project8/vyp/Winship_GVHD/claire/results/syn_allo_bmt/results/", choice,"_ImmGen_coarse__module_association_Fishers_exact_test_Pvalues.txt", sep = ""), sep = "\t")
colnames(modules)[1] <- "gene.name"
modules$gene.name <- tolower(modules$gene.name)
final.table = merge(my.table, modules, by = "gene.name", all.x = TRUE)
final.table$module = NULL
write.table(final.table, file = paste("/cluster/project8/vyp/Winship_GVHD/claire/results/syn_allo_bmt/results/", choice, "_final_table_coarse.txt", sep = ""), sep = "\t")


########## plot #############
#########       #############


xmin = 0
xmax = 82
ymin = 0
ymax = 18

pdf(paste("/cluster/project8/vyp/Winship_GVHD/claire/results/syn_allo_bmt/figs/", choice, "_coarse_module_association_graph.pdf", sep = ""))
plot(test1, type = "n", main = paste("Fisher's Exact Test p-values for association \n of the ", choice, " differential gene expression \n dataset with the ImmGen Coarse modules", sep = ""), xlab = "Coarse module number", ylab = "-log10(p-value)", xlim = c(xmin, xmax), ylim = c(ymin, ymax), xaxs = "i", yaxs = "i")
axis(side=1, at=seq(0, 82, by = 2))
axis(side=2, at=seq(0, 18, by = 2))
lines(test1, type = "h", col = "blue")
dev.off()
}

