library(oligo)
library(limma)
library(lattice)

all.files <- list.celfiles('/cluster/project8/vyp/Winship_GVHD/rawData/MataHari_model/MH_CEL_files/', full.names = TRUE, recursive = TRUE)
selected.files <- sort(subset(all.files, grepl(pattern = 'D9|A9|D10|A10|D11|A11', all.files)))

set.probes <- 'core'  ##extended, full are other options
affyexpression <- read.celfiles(selected.files)
summaries <- rma(affyexpression, target = set.probes)

save(list = "summaries", file = "/cluster/project8/vyp/Winship_GVHD/claire/results/epi_dermis_PLN/results/summary_core.RData")

my.pca <- prcomp(t(exprs(summaries)), retx = TRUE)$x

pdf("/cluster/project8/vyp/Winship_GVHD/claire/results/epi_dermis_PLN/figs/PCA.pdf")
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
my.frame$group = c("epidermisDT", "epidermisDT", "epidermisDT","epidermis", "epidermis", "epidermis", "dermisDT", "dermisDT", "dermisDT","dermis", "dermis", "dermis", "PLN", "PLN", "PLN", "PLNDT", "PLNDT", "PLNDT")

g <- ggplot(data = my.frame, aes(x = PC1, y = PC2)) + geom_point()
g <- g + ggplot2::geom_text(ggplot2::aes(x = PC1, y = PC2, label = array, vjust = -1, size = 2)) + ggplot2::scale_size(guide = "none")
g <- g + geom_point(aes(colour = factor(my.frame$group)), size = 5)
g <- g + labs(colour='group')
ggsave(g, width = 10, file = "/cluster/project8/vyp/Winship_GVHD/claire/results/epi_dermis_PLN/figs/PCA_prettier.pdf")




