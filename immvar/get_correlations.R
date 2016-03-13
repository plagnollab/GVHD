all.data <- read.csv("/cluster/project8/vyp/Winship_GVHD/claire/data_files/ImmGen_sample_data/all_celltypes.csv", stringsAsFactors = FALSE)
all.data[,2] <- gsub(pattern = "^ ", all.data[,2], replacement = "")





######## some magic below with that dplyr package
library(dplyr)
collapsed.data <- all.data[, -c(1,3)] %>% group_by(GeneSymbol) %>% summarise_each(funs(mean))  ## note the grouping by gene symbol then summarize each column by mean
collapsed.data.matrix <- as.matrix(as.data.frame(collapsed.data)[, -1])
dimnames(collapsed.data.matrix)[[1]] <- collapsed.data$GeneSymbol



####### now let's start looking at that turquoise module
turquoise.genes <- scan("WGCNA_turquoise_mod_genelist.csv", sep = ",", what = character())
missing.genes <- turquoise.genes[ ! turquoise.genes %in% dimnames(collapsed.data.matrix)[[1]] ] ## a few quirks, gene names that are not there. Not many...
turquoise.genes <- turquoise.genes[ ! turquoise.genes %in% missing.genes ]

  
turquoise.mod <- collapsed.data.matrix[ turquoise.genes, ]


### I want to pick a representative gene of the module, so I propose to pick the one showing the maximum absolute value correlation
full.cor.matrix <- cor(t(turquoise.mod))
my.mean.cor <- apply(full.cor.matrix, MAR = 1, FUN = function(x) {mean(abs(x), na.rm = TRUE)})  ## note abs is important here, I just want high correlation, plus or minus
my.representative.gene <- names(which.max(my.mean.cor))
my.cor.values <- cor(t(turquoise.mod), turquoise.mod[ my.representative.gene,])


## positively and negatively correlated genes wrt to the "representative" gene
positive.genes <- dimnames(my.cor.values)[[1]][ my.cor.values[,1] > 0 ]
negative.genes <- dimnames(my.cor.values)[[1]][ my.cor.values[,1] < 0 ]

cor.matrix.pos.genes <- cor(t(turquoise.mod[ positive.genes, ]))
cor.matrix.neg.genes <- cor(t(turquoise.mod[ negative.genes, ]))


##### now some final check, once the subsetting is done between positive and negative genes, we should really mostly get positive values
print(table(sign(cor.matrix.neg.genes)))
print(table(sign(cor.matrix.pos.genes)))


message("Number of positive genes: ", length(positive.genes))
message("Number of negative genes: ", length(negative.genes))
