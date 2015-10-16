library(WGCNA);
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#RCD4datad in the female liver data set
EA = read.table("/cluster/project8/vyp/Winship_GVHD/claire/data_files/immvar_expression/GSE56033_GSM.ImmVarCD4.EA.PC12.txt", header = T)
AA = read.table("/cluster/project8/vyp/Winship_GVHD/claire/data_files/immvar_expression/GSE56033_GSM.ImmVarCD4.AA.PC12.txt", header = T)
EU = read.table("/cluster/project8/vyp/Winship_GVHD/claire/data_files/immvar_expression/GSE56033_GSM.ImmVarCD4.EU.PC20.txt", header = T)
EU$ID_REF = NULL
AA$ID_REF = NULL

expression_temp = cbind(EA, AA)
expression_all = cbind(expression_temp, EU)
CD4data = as.data.frame(expression_all)
# Take a quick look at what is in the data set:
dim(CD4data);

datExpr0 = as.data.frame(t(CD4data[,-c(1)]));
names(datExpr0) = CD4data$ID_REF;
rownames(datExpr0) = names(CD4data)[-c(1)];
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(25,9)
#png(file = "/cluster/project8/vyp/Winship_GVHD/claire/results/immvar/plots/sampleClustering.png", width = 15, height = 9);
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#     cex.axis = 1.5, cex.main = 2)
#abline(h = 42, col = "red")
#dev.off()

clust = cutreeStatic(sampleTree, cutHeight = 42)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, dataIsExpr = TRUE, powerVector = powers, verbose = 5)

pdf(file = "/cluster/project8/vyp/Winship_GVHD/claire/GVHD/immvar/results/plots/scale_free_topology_fit.png", width = 12, height = 9);
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
#plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
 #    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
  #   main = paste("Scale independence"));
#text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
 #    labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
#abline(h=0.97,col="red")
# Mean connectivity as a function of the soft-thresholding power
#plot(sft$fitIndices[,1], sft$fitIndices[,5],
 #    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
  #   main = paste("Mean connectivity"))
#text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()

##################### check scale free topology ################################
k=softConnectivity(datE=datExpr,power=3)
pdf(file = "/cluster/project8/vyp/Winship_GVHD/claire/GVHD/immvar/results/plots/check_scale_free_topology.pdf", width = 12, height = 9);
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
dev.off()

net = blockwiseModules(datExpr, power = 3,
  TOMType = "unsigned", minModuleSize = 20,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "TOM",
  verbose = 3)

table(net$colors)

pdf(file = "/cluster/project8/vyp/Winship_GVHD/claire/GVHD/immvar/results/plots/cluster_dendogram.pdf", width = 12, height = 9);
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "immvar-networkConstruction-auto.RData")

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#################### annotation ########################

annot = read.csv(file = "/cluster/project8/vyp/Winship_GVHD/claire/scripts/immvar/HuGene-1_0-st.csv", header = T);
dim(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$transcript_cluster_id)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.
gene.list <- unique(strsplit(paste(annot$gene_assignment, collapse = " /// "), split = " /// ")[[1]])
gene.list <- unique(strsplit(paste(gene.list, collapse = " // "), split = " // ")[[1]])
gene.list <- gsub(gene.list, pattern = "ENST.+", gene.list, replacement = "")
x = select(org.Hs.eg.db, gene.list, c("ENTREZID"), "ALIAS")
allLLIDs = x$ENTREZID[probes2annot];
#print(allLLIDs)

GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "human", nBestP = 10);
tab = GOenr$bestPTerms[[4]]$enrichment

geneInfo0 = data.frame(transcript_cluster_id = probes,
  geneSymbol = annot$gene_symbol[probes2annot],
  LocusLinkID = annot$LocusLinkID[probes2annot],
  moduleColor = moduleColors)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));

write.csv(geneInfo, file = "/cluster/project8/vyp/Winship_GVHD/claire/GVHD/immvar/results/geneInfo.csv")
