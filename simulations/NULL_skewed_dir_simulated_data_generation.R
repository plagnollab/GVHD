library(snpStats)

source("association/scripts/group_pval.R")
###############    testing the null hypothesis - uniform distribution ###################
increment_jump = c("-0.1") #, "-0.01")

ngenes <- 10000    
ntissues <- 16

output.pdf <- "results/null_simulations.pdf"
pdf(output.pdf, height = 4, width = 4)

for (proba.minus in 0.5) {
  
  ### generating "p-values" ###
  
  fraction.positive.tests <- 0.
  n.positive.tests <- 0

  ## we don't need the stuff below
  #sig_pvalues_sample = pmax(rnorm(n.positive.tests, mean = 2, sd = 1), 0) ## VP comment: the pmax caps the P-values between 0 and 1
  #sig_pvalues = 10^(-sig_pvalues_sample)
  #sig_pvalues <- runif(n.positive.tests)
  uniform_pvalues = runif(min = 0, max = 1, n = ngenes*ntissues) 
  comb_pvalues = uniform_pvalues
  
  randm_dir = sign(runif(min = -1, max = 1, n = ngenes*ntissues))  ## VP comment: I propose to set the background 1-1

  ## No need for the stuff below
  ##skewed_dir = sign(-(runif(min = 0, max = 1, n = n.positive.tests) - proba.minus)) ## VP comment: but modify the strength of the skew for positive results
  comb.directions <- randm_dir

  ## VP COMMENT: Claire, important: you need to randomize your results here, you set it up to have all the sig P-values at the end
  my.random <- sample(1:length(comb_pvalues), replace = FALSE)
  comb_pvalues <- comb_pvalues[ my.random ]
  comb.directions <- comb.directions[ my.random  ]

  pvalues <- matrix(nrow = ngenes, ncol = ntissues, data = as.numeric(comb_pvalues))
  DoE <- matrix(nrow = ngenes, ncol = ntissues, data = comb.directions)
  
  full.list <- c()
  for (tissue in 1:ntissues) {
    
    for (module in 1:50) {
      #print(module)
      selected.genes <- sample(1:ngenes, size = 200, replace = FALSE)
      #print(selected.genes)
      my.pvalue <- groupbased.pvalue(pvalues = pvalues[selected.genes, tissue], 
                                     sign.association = DoE[selected.genes, tissue], 
                                     genomewide.prop.positive = 0.5)  ## VP comment: I propose to test assuming a 1-1 distributin
      #print(my.pvalue)
      full.list <- c(full.list, my.pvalue)
    }  
  }



  qq.chisq(-2*log(full.list), df = 2, main = paste("Null hypothesis"))
  abline(0, 1, col = 'red')
}

dev.off()
message(output.pdf)
