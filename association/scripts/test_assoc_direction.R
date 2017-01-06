data <- read.csv("association/overall_DofE_table.txt", header = TRUE, stringsAsFactors = FALSE)

data$module.DofE.num <- as.numeric(factor(data$module_DofE))
assoc.pvalues <- data.frame( module = unique(data$module_id),
                            basic.assoc.pvalue = NA,
                            refined.assoc.pvalue = NA,
                            n.genes = NA,
                            stringsAsFactors = FALSE)

assoc.pvalues <- assoc.pvalues[ ! is.na(assoc.pvalues$module), ]

data <- data[ ! is.na(data$module_DofE) & !is.na(data$epidermis_vs_epidermisDT_p.value), ]

data$epidermis_vs_epidermisDT_logFC.signed <- ifelse (data$module_DofE == "-", - data$epidermis_vs_epidermisDT_logFC, data$epidermis_vs_epidermisDT_logFC)

### first a look whether DoE in one module is a good predictor of module DoE
#print( summary(glm(data = data, formula = "module.DofE.num ~ epidermis_vs_epidermisDT_logFC:module_id")))


pdf("modules_direction.pdf")
for (i in 1:nrow(assoc.pvalues)) {
#for (i in 1:5) {

  colour <- assoc.pvalues$module[ i ]
  message("Module ", i, " out of ", nrow(assoc.pvalues), ", colour ", colour)
  
  assoc.pvalues$n.genes[ i ] <- sum(data$module_id == assoc.pvalues$module[ i ], na.rm = TRUE)
  assoc.pvalues$basic.assoc.pvalue[ i ] <- fisher.test(table(data$module_id == assoc.pvalues$module[ i ], data$dermis_vs_dermisDT_p.value < 0.05))$p.value

  loc.module <- subset(data, module_id == colour)

  
  

  
  #### now for each gene we find the chi square statistic associated with each gene
  loc.module$chisq <- qchisq(p = loc.module$epidermis_vs_epidermisDT_p.value, df = 1, lower.tail = FALSE)

  ####### now we need to test two possible directions, and we have two different stats
  stat.all.positive <- sum(ifelse (loc.module$epidermis_vs_epidermisDT_logFC.signed > 0, loc.module$chisq, 0))
  stat.all.negative <- sum(ifelse (loc.module$epidermis_vs_epidermisDT_logFC.signed < 0, loc.module$chisq, 0))

  final.stat <- max(stat.all.positive, stat.all.negative)
  
  ###### now we need to assess the distribution of the statistic under the null
  ### it is a mixture model. 50% of the time there is no improvement, 50% chisq
  ### we start by assessning the number of time we get a non zero chisq stat
  proba.nb.nonzeros.under.null <- dbinom(size = nrow(loc.module), prob = 0.5, x = 1:nrow(loc.module) )
  proba.nb.nonzeros.under.null <- ifelse (proba.nb.nonzeros.under.null < 10^(-8), 0, proba.nb.nonzeros.under.null) ## just to avoid computational issues with small nbs

  #### now we have a probability to get each nb of non-zero chisquares and for each of these values, we can use the chi-squared function in R to compare to the observed stat
  proba.under.null <- sum(pchisq(df = 1:nrow(loc.module), q = final.stat, lower.tail = FALSE)*proba.nb.nonzeros.under.null)
  assoc.pvalues$refined.assoc.pvalue[ i ] <- proba.under.null

########
  plot(x = loc.module$epidermis_vs_epidermisDT_logFC, y = 1:nrow(loc.module), main =
       paste(colour, ", P = ", signif(assoc.pvalues$basic.assoc.pvalue[ i ], 3), " and refined P-value is ", signif(proba.under.null, 3)),
       pch = "+")
  abline(v = 0, col = "red")
  
}
dev.off()

print(subset(assoc.pvalues, refined.assoc.pvalue < 0.001 | basic.assoc.pvalue < 0.0001))
