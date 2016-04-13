library(ggplot2)
library(ggrepel)

data <- read.csv("/cluster/project8/vyp/Winship_GVHD/vincent/GVHD/association/overall_DofE_table.txt", header = TRUE, stringsAsFactors = FALSE)

for(choice in c( "epidermis_vs_epidermisDT", "dermis_vs_dermisDT","PLN_vs_PLNDT","D7","D42","BMT_male_vs_BMT_female","BMT_male_vs_naive_male","BMT_female_vs_naive_male"   )){
  data$module.DofE.num <- as.numeric(factor(data$module_DofE))
  assoc.pvalues <- data.frame( module = unique(data$module_id),
                               basic.assoc.pvalue = NA,
                               refined.assoc.pvalue = NA,
                               neg.log.pvalue = NA,
                               n.genes = NA,
                               stringsAsFactors = FALSE)
  
  assoc.pvalues <- assoc.pvalues[ ! is.na(assoc.pvalues$module), ]
  
  choice_p.value =paste(choice,"_p.value", sep="")
  choice_logFC =paste(choice,"_logFC", sep="")
  choice_p.value.signed =paste(choice,"_p.value.signed", sep="")
  choice_logFC.signed =paste(choice,"_logFC.signed", sep="")
  
  data <- data[ ! is.na(data$module_DofE) & !is.na(data[,choice_p.value]), ]
  
  data[,choice_logFC.signed] <- ifelse (data$module_DofE == "-", - data[,choice_logFC], data[,choice_logFC])
  
  ### first a look whether DoE in one module is a good predictor of module DoE
  #print( summary(glm(data = data, formula = "module.DofE.num ~ epidermis_vs_epidermisDT_logFC:module_id")))
  
  
  pdf(paste("/cluster/project8/vyp/Winship_GVHD/claire/results/module_association/", choice, "_modules_direction.pdf", sep=""))
  for (i in 1:nrow(assoc.pvalues)) {
    #for (i in 1:5) {
    
    colour <- assoc.pvalues$module[ i ]
    message("Module ", i, " out of ", nrow(assoc.pvalues), ", colour ", colour)
    
    assoc.pvalues$n.genes[ i ] <- sum(data$module_id == assoc.pvalues$module[ i ], na.rm = TRUE)
    assoc.pvalues$basic.assoc.pvalue[ i ] <- fisher.test(table(data$module_id == assoc.pvalues$module[ i ], data[,choice_p.value] < 0.05))$p.value
    
    loc.module <- subset(data, module_id == colour)
    
    #### now for each gene we find the chi square statistic associated with each gene
    loc.module$chisq <- qchisq(p = loc.module[,choice_p.value], df = 1, lower.tail = FALSE)
    
    ####### now we need to test two possible directions, and we have two different stats
    stat.all.positive <- sum(ifelse (loc.module[,choice_logFC.signed] > 0, loc.module$chisq, 0))
    stat.all.negative <- sum(ifelse (loc.module[,choice_logFC.signed] < 0, loc.module$chisq, 0))
    
    final.stat <- max(stat.all.positive, stat.all.negative)
    
    ###### now we need to assess the distribution of the statistic under the null
    ### it is a mixture model. 50% of the time there is no improvement, 50% chisq
    ### we start by assessning the number of time we get a non zero chisq stat
    proba.nb.nonzeros.under.null <- dbinom(size = nrow(loc.module), prob = 0.5, x = 1:nrow(loc.module) )
    proba.nb.nonzeros.under.null <- ifelse (proba.nb.nonzeros.under.null < 10^(-8), 0, proba.nb.nonzeros.under.null) ## just to avoid computational issues with small nbs
    
    #### now we have a probability to get each nb of non-zero chisquares and for each of these values, we can use the chi-squared function in R to compare to the observed stat
    proba.under.null <- sum(pchisq(df = 1:nrow(loc.module), q = final.stat, lower.tail = FALSE)*proba.nb.nonzeros.under.null)
    assoc.pvalues$refined.assoc.pvalue[ i ] <- proba.under.null
    assoc.pvalues$neg.log.pvalue[i] <- -log10(assoc.pvalues$refined.assoc.pvalue[ i ])
    
    ########
    plot(loc.module[,choice_logFC], 1:length(loc.module[,choice_logFC]), main =
           paste(colour, ", P = ", signif(assoc.pvalues$basic.assoc.pvalue[ i ], 3), " and refined P-value is ", signif(proba.under.null, 3)),
         pch = "+")
    abline(v = 0, col = "red")
    
  }
  dev.off()
  assoc.pvalues$log.n.genes <- log(assoc.pvalues$n.genes)
  write.table(assoc.pvalues, paste("/cluster/project8/vyp/Winship_GVHD/claire/results/module_association/", choice, "_association_table.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
  print(subset(assoc.pvalues, refined.assoc.pvalue < 0.001 | basic.assoc.pvalue < 0.0001))
  significant = subset(assoc.pvalues, refined.assoc.pvalue < 0.001 | basic.assoc.pvalue < 0.0001)
  plot_data = assoc.pvalues[c("module", "neg.log.pvalue")]
  
  xmin = 0
  xmax = nrow(plot_data)
  ymin = 0
  ymax = max(plot_data$neg.log.pvalue) +1
  
  module.list = as.vector(assoc.pvalues$module)
  pdf(file = paste("/cluster/project8/vyp/Winship_GVHD/claire/results/module_association/plots/", choice, "_WGCNA_modules_refined_association_graph.pdf", sep = ""), height = 12, width = 12)
  plot(plot_data$neg.log.pvalue, type = "n", main = paste("Refined Fishers Exact Test p-values for association of the ", choice, " differential gene expression \n dataset with the new WGCNA defined modules", sep= ""),xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "Module", ylab = "-log10(p-value)", xaxs = "i", yaxs = "i", xaxt = "n")
  # axis(side=1, at=seq(1,xmax,by = 2))
  axis(side=2, at=seq(0, ymax, by = 2))
  text(x = 1:nrow(plot_data), y = 0,
       srt = 45, adj= 1, xpd = TRUE,
       labels = plot_data$module,
       cex=0.7)
  lines(plot_data$neg.log.pvalue, type = "h", col = "blue")
  dev.off()
  x = ggplot(assoc.pvalues, aes(y = log.n.genes, x = neg.log.pvalue)) +
    ggtitle(paste(choice, "module association p values \n with respect to module size", sep = " ")) +
    geom_point() +
    geom_point() +geom_text_repel(data = subset(assoc.pvalues, neg.log.pvalue > 4), aes(neg.log.pvalue, label = module)) +
    theme(plot.title = element_text(family = "Helvetica", color="#666666", face="bold", size=20)) +
    theme(axis.title = element_text(family = "Helvetica", color="#666666", face=NULL, size=15)) 
  #geom_point() +geom_point(data = subset(assoc.pvalues, neg.log.pvalue > 4), aes(color = color_palette))
  ggsave(paste("/cluster/project8/vyp/Winship_GVHD/claire/results/module_association/plots/", choice, "_WGCNA_modules_association_vs_module_size_graph.pdf", sep = ""), plot = x, width = 11, height = 11)
  pdf(paste("/cluster/project8/vyp/Winship_GVHD/claire/results/module_association/", choice, "_significant_modules_direction.pdf", sep=""))
  for (i in 1:nrow(assoc.pvalues)) {
    colour <- assoc.pvalues$module[ i ]
    if (colour %in% c(significant$module)){
      print("matched")
      assoc.pvalues$n.genes[ i ] <- sum(data$module_id == assoc.pvalues$module[ i ], na.rm = TRUE)
      assoc.pvalues$basic.assoc.pvalue[ i ] <- fisher.test(table(data$module_id == assoc.pvalues$module[ i ], data[,choice_p.value] < 0.05))$p.value
      
      loc.module <- subset(data, module_id == colour)
      
      #### now for each gene we find the chi square statistic associated with each gene
      loc.module$chisq <- qchisq(p = loc.module[,choice_p.value], df = 1, lower.tail = FALSE)
      
      ####### now we need to test two possible directions, and we have two different stats
      stat.all.positive <- sum(ifelse (loc.module[,choice_logFC.signed] > 0, loc.module$chisq, 0))
      stat.all.negative <- sum(ifelse (loc.module[,choice_logFC.signed] < 0, loc.module$chisq, 0))
      
      final.stat <- max(stat.all.positive, stat.all.negative)
      
      proba.nb.nonzeros.under.null <- dbinom(size = nrow(loc.module), prob = 0.5, x = 1:nrow(loc.module) )
      proba.nb.nonzeros.under.null <- ifelse (proba.nb.nonzeros.under.null < 10^(-8), 0, proba.nb.nonzeros.under.null) ## just to avoid computational issues with small nbs
      
      #### now we have a probability to get each nb of non-zero chisquares and for each of these values, we can use the chi-squared function in R to compare to the observed stat
      proba.under.null <- sum(pchisq(df = 1:nrow(loc.module), q = final.stat, lower.tail = FALSE)*proba.nb.nonzeros.under.null)
      assoc.pvalues$refined.assoc.pvalue[ i ] <- proba.under.null
      
      ########
      plot(loc.module[,choice_logFC], 1:length(loc.module[,choice_logFC]), main =
             paste(colour, ", P = ", signif(assoc.pvalues$basic.assoc.pvalue[ i ], 3), " and refined P-value is ", signif(proba.under.null, 3)),
           pch = "+")
      abline(v = 0, col = "red")
      gene_data = loc.module[c("gene", choice_p.value, choice_logFC.signed)]
      colnames(gene_data)[2] = "p.value"
      colnames(gene_data)[3] = "logFC.signed"
      write.table(gene_data, paste("/cluster/project8/vyp/Winship_GVHD/claire/results/module_association/", choice, "_", colour, "_gene_data.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
    }
    
    #  dev.off()
  }
  dev.off()
}

