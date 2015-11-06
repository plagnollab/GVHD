library(oligo)
library(limma)
library(lattice)

all.files <- list.celfiles('/cluster/project8/vyp/Winship_GVHD/rawData/Langerhans_cells/', full.names = TRUE, recursive = TRUE)

set.probes <- 'core'  ##extended, full are other options
affyexpression <- read.celfiles(all.files)
summaries <- rma(affyexpression, target = set.probes)

save(list = "summaries", file = "/cluster/project8/vyp/Winship_GVHD/claire/results/syn_allo_bmt/results/summary_core.RData")

my.pca <- prcomp(t(exprs(summaries)), retx = TRUE)$x

pdf("/cluster/project8/vyp/Winship_GVHD/claire/results/syn_allo_bmt/figs/PCA.pdf")
plot(x = my.pca[, 1], ##first PC
     y = my.pca[, 2])  ## second PC

text(x = my.pca[, 1], ##first PC
     y = my.pca[, 2], ##second PC
     label = gsub(row.names(my.pca), pattern = ".CEL", replacement = ""))  ## second PC

plot(x = my.pca[, 1], ##third PC
     y = my.pca[, 2])  ## fourth PC

text(x = my.pca[, 1], ##third PC
     y = my.pca[, 2], ##fourth PC
     label = gsub(row.names(my.pca), pattern = ".CEL", replacement = ""))  ## second PC

dev.off()


#### prettier ggplot version
library(ggplot2)
my.frame <- as.data.frame(my.pca[, 1:4])  ##ggplot wants data frame
my.frame$array <- gsub(row.names(my.pca), pattern = ".CEL", replacement = "")
my.frame$group = c("BMT_male", "BMT_male", "BMT_male", "BMT_female", "BMT_female", "BMT_female", "naive_male", "naive_male", "naive_male")

g <- ggplot(data = my.frame, aes(x = PC1, y = PC2)) + geom_point()
g <- g + ggplot2::geom_text(ggplot2::aes(x = PC1, y = PC2, label = array, vjust = -1, size = 2)) + ggplot2::scale_size(guide = "none")
#g <- g + geom_text(size = 0.3)
g <- g + geom_point(aes(colour = factor(my.frame$group)), size = 5)
g <- g + labs(colour='group')
ggsave(g, width = 10, file = "/cluster/project8/vyp/Winship_GVHD/claire/results/syn_allo_bmt/figs/PCA_prettier.pdf")




