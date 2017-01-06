modules <- read.table("paper/support/module_assignment.tab", header = TRUE, stringsAsFactors = FALSE)


mod28.genes <- subset(modules, module <- id == "MH_28")[,1]

data <- read.csv("/cluster/project8/vyp/Winship_GVHD/claire/results/epi_dermis_PLN/results/epidermis_vs_epidermisDT.csv", stringsAsFactors = FALSE)


data$gene <- gsub(pattern = " //.*", replacement = "", gsub(data$gene_assignment, pattern = "^[.|A-Z|_|0-9]* // ", replacement = ""))

## restrict to the genes that have a module
data <- data[ data$gene %in% modules$gene, ]


data$module <- modules$module_id[ match(data$gene, table = modules$gene) ]

print(table(sign(data$logFC), data$module))


### a basic Fisher test to assess significance at P < 0.05


list.modules <-  sort(unique(modules$module_id))
list.modules <- list.modules[ ! list.modules == "grey" ]

my.res.table <- data.frame(module = list.modules, p.value.direction = NA, p.value.basic = NA, n.genes = NA, stringsAsFactors = FALSE)

for (i in 1:nrow(my.res.table)) {
  choice.module <- my.res.table$module[ i ]
  loc.module <- subset(data, module == choice.module)
  
  my.res.table$n.genes[ i ] <- nrow(loc.module)
  
  my.res.table$p.value.basic[ i ] <- fisher.test(table(data$module == choice.module, data$P.Value < 0.05))$p.value
  
#### now for each gene we find the chi square statistic associated with each gene
  loc.module$chisq <- qchisq(p = loc.module[, "P.Value"], df = 1, lower.tail = FALSE)

####### now we need to test two possible directions, and we have two different stats
  stat.all.positive <- sum(ifelse (loc.module[,"logFC"] > 0, loc.module$chisq, 0))
  stat.all.negative <- sum(ifelse (loc.module[,"logFC"] < 0, loc.module$chisq, 0))
  
  final.stat <- max(stat.all.positive, stat.all.negative)
  
###### now we need to assess the distribution of the statistic under the null
### it is a mixture model. 50% of the time there is no improvement, 50% chisq
### we start by assessning the number of time we get a non zero chisq stat
  proba.nb.nonzeros.under.null <- dbinom(size = nrow(loc.module), prob = 0.5, x = 1:nrow(loc.module) )
  proba.nb.nonzeros.under.null <- ifelse (proba.nb.nonzeros.under.null < 10^(-8), 0, proba.nb.nonzeros.under.null) ## avoid comp issues with small nbs

#### now we have a probability to get each nb of non-zero chisquares and for each of these values, we use chi2 in R to compare to the observed stat
  proba.under.null <- sum(pchisq(df = 1:nrow(loc.module), q = final.stat, lower.tail = FALSE)*proba.nb.nonzeros.under.null)
  my.res.table$p.value.direction[ i ] <- proba.under.null
}
write.csv(x = my.res.table, row.names = FALSE, quote = FALSE, file = "modules_epidermisPBSvsDT.csv")

#data <- data[ ! grepl(pattern = "7SK|[A-Z][0-9][0-9]+|^Olfr|^linc|^Gm[0-9]|U1|U4|U6|U7|^Mir|5S_rRNA|^[0-9]+4|^[0-9][0-9][0-9]|^LOC[0-9]|^AC[0-9]|^SNOR",
#data$gene), ]
#print(table(mod28.genes %in% data$gene))
#print(mod28.genes[ ! mod28.genes %in% data$gene ])
