data <- read.csv("association/overall_DofE_table.txt", header = TRUE, stringsAsFactors = FALSE)

data$module.DofE.num <- as.numeric(factor(data$module_DofE))
assoc.pvalues <- data.frame( module = unique(data$module_id),
                            basic.assoc.pvalue = NA,
                            refined.assoc.pvalue = NA,
                            n.genes = NA,
                            stringsAsFactors = FALSE)

assoc.pvalues <- assoc.pvalues[ ! is.na(assoc.pvalues$module), ]

data <- data[ ! is.na(data$module_DofE), ]

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
  plot(x = loc.module$epidermis_vs_epidermisDT_logFC, y = 1:nrow(loc.module), main = paste(colour, ", P = ", signif(assoc.pvalues$basic.assoc.pvalue[ i ], 3)), pch = "+")
  abline(v = 0, col = "red")
  
  
  ###### now trying something directional
  data$good.pvalue <- data$dermis_vs_dermisDT_p.value < 0.05 & sign(data$dermis_vs_dermisDT_logFC)  == 1 & data$module.DofE.num == 2
  first.pval.with.direction <- fisher.test(table(data$module_id == assoc.pvalues$module[ i ], data$good.pvalue))$p.value


  data$good.pvalue <- data$dermis_vs_dermisDT_p.value < 0.05 & sign(data$dermis_vs_dermisDT_logFC)  == -1 & data$module.DofE.num == 2
  second.pval.with.direction <- fisher.test(table(data$module_id == assoc.pvalues$module[ i ], data$good.pvalue < 0.05))$p.value

  assoc.pvalues$refined.assoc.pvalue[ i ] <- min(first.pval.with.direction, second.pval.with.direction)

  
}
dev.off()

 print(subset(assoc.pvalues, refined.assoc.pvalue < 0.001))
