TM008_no_outlier = read.csv("/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/no_outlier/TM008wt_vs_TM008ko.csv", header = TRUE)
TM006_no_outlier = read.csv("/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/no_outlier/TM006wt_vs_TM006ko.csv", header = TRUE)

significant_TM006 = subset(TM006_no_outlier, as.vector(TM006_no_outlier$P.Value) < 0.05)
significant_TM008 = subset(TM008_no_outlier, as.vector(TM008_no_outlier$P.Value) < 0.05)

non_significant_TM006 = subset(TM006_no_outlier, as.vector(TM006_no_outlier$P.Value) > 0.05)
non_significant_TM008 = subset(TM008_no_outlier, as.vector(TM008_no_outlier$P.Value) > 0.05)

write.csv(x= significant_TM006, file = paste('/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/DE_D42WT_VS_D42KO_ALLANALYSIS_TM006.csv',sep =""))
write.csv(x= non_significant_TM006, file = paste('/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/non_DE_D42WT_VS_D42KO_ALLANALYSIS_TM006.csv',sep =""))
write.csv(x= significant_TM008, file = paste('/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/DE_D7WT_VS_D7KO_ALLANALYSIS_TM008.csv',sep =""))
write.csv(x= non_significant_TM008, file = paste('/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/non_DE_D7WT_VS_D7KO_ALLANALYSIS_TM008.csv',sep =""))

########### plot of p values ########

pdf("/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/figs/D7_histogram.pdf")

hist(TM008_no_outlier$P.Value, breaks=seq(0, 1, by=0.02), main = "Histogram of D7 (TM008) MHC1-KO vs WT p-values", xlab = "p-value", col = "red")
dev.off()

pdf("/cluster/project8/vyp/Winship_GVHD/claire/results/mhc1_ko/figs/D42_histogram.pdf")
TM006HIST = hist(TM006_no_outlier$P.Value, breaks=seq(0, 1, by=0.02), main = "Histogram of D42 (TM006) MHC1-KO vs WT p-values", xlab = "p-value", col = "red" )

dev.off()
