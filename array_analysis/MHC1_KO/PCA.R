library(oligo)
library(limma)
library(lattice)

all.files <- list.celfiles('/cluster/project8/vyp/Winship_GVHD/claire/data_files/Teresa_microarray', full.names = TRUE, recursive = TRUE)

set.probes <- 'core'  ##extended, full are other options
affyexpression <- read.celfiles(all.files)
summaries <- rma(affyexpression, target = set.probes)

save(list = "summaries", file = "/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/results/summary_core.RData")

my.pca <- prcomp(t(exprs(summaries)), retx = TRUE)$x

pdf("/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/figs/PCA.pdf")
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

g <- ggplot(data = my.frame, aes(x = PC1, y = PC2)) + geom_point()
g <- g + geom_text(aes(x = PC1, y = PC2, label = array))
ggsave(g, file = "/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/figs/PCA_prettier.pdf")


### now we now that TM008_ko4 is an outlier, it makes the graph less interesting. Let us remove and redo

summaries.no.outlier <- summaries[, row.names(pData(summaries)) != "TM008_ko4.CEL"] ##here we remove TM008_ko4.CEL
my.pca <- prcomp(t(exprs(summaries.no.outlier)), retx = TRUE)$x  ## and we redo a PCA
my.frame <- as.data.frame(my.pca[, 1:4])  ##ggplot wants data frame
my.frame$array <- gsub(row.names(my.pca), pattern = ".CEL", replacement = "")

g <- ggplot(data = my.frame, ggplot2::aes(x = PC1, y = PC2)) + ggplot2::geom_point()
g <- g + ggplot2::geom_text(ggplot2::aes(x = PC1, y = PC2, label = array, vjust = -1, size = 2)) + ggplot2::scale_size(guide = "none")
g <- g + ggplot2::ggtitle("TM008_ko4.CEL removed")
ggsave(g, width = 10, height = 5, file = "/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/figs/PCA_prettier_outlier_removed.pdf")
