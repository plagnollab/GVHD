
groupbased.pvalue <- function(pvalues, sign.association, genomewide.prop.positive = 0.5) {
  if (length(pvalues) != length(sign.association)) stop("Argument lengths do not match")
  if (class(pvalues) != "numeric") stop("Pvalues must be numeric")
  
  ngenes <- length(pvalues)
  gene.chisq <- qchisq(p = pvalues, df = 1, lower.tail = FALSE)
  stat.all.positive <- sum(ifelse(sign.association > 0, gene.chisq, 0))
  stat.all.negative <- sum(ifelse (sign.association < 0, gene.chisq, 0))
  final.stat <- max(stat.all.positive, stat.all.negative)
  
  proba.nb.nonzeros.under.null <- dbinom(size = ngenes, prob = genomewide.prop.positive, x = 1:ngenes )
  proba.nb.nonzeros.under.null <- ifelse (proba.nb.nonzeros.under.null < 10^(-8), 0, proba.nb.nonzeros.under.null) ## just to avoid computational issues with small nbs
  
#### now we have a probability to get each nb of non-zero chisquares and for each of these values, we can use the chi-squared function in R to compare to the observed stat
  proba.under.null <- sum(pchisq(df = 1:ngenes, q = final.stat, lower.tail = FALSE)*proba.nb.nonzeros.under.null)
  return(proba.under.null)
}
