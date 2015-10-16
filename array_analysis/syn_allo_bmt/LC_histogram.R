for (choice in c("BMT_male_vs_BMT_female", "BMT_female_vs_naive_male", "BMT_male_vs_naive_male")){

infile = read.csv(paste("/cluster/project8/vyp/Winship_GVHD/claire/GVHD/array_analysis/syn_allo_bmt/results/", choice, ".csv", sep = ""), header = TRUE)

significant = subset(infile, as.vector(infile$P.Value) < 0.05)

non_significant = subset(infile, as.vector(infile$P.Value) > 0.05)

write.csv(x= significant, file = paste('/cluster/project8/vyp/Winship_GVHD/claire/GVHD/array_analysis/syn_allo_bmt/results/DE_', choice, '.csv',sep =""))
write.csv(x= non_significant, file = paste('/cluster/project8/vyp/Winship_GVHD/claire/GVHD/array_analysis/syn_allo_bmt/results/non_DE_', choice, '.csv',sep =""))

########### plot of p values ########

pdf(paste("/cluster/project8/vyp/Winship_GVHD/claire/GVHD/array_analysis/syn_allo_bmt/figs/", choice, "_histogram.pdf", sep = ""))

hist(infile$P.Value, breaks=seq(0, 1, by=0.02), main = paste("Histogram of ", choice, " p-values", sep = ""), xlab = "p-value", col = "red")

dev.off()

}
