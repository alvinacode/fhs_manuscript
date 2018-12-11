###################################################
# Vectors/ Constants
###################################################

cons_names = c("EXP22BCLEX1", "EXP22BCL03A", "EXP22BCL04A", "EXP22BCLSB", "EXP22DEEX1", "EXP22DE03A", "EXP22DE04A", "EXP22DESB",
               "EXP22PEEX1", "EXP22PE03A", "EXP22PE04A", "EXP22PESB", "EXP26BCLEX1", "EXP26BCL03A", "EXP26BCL04A", "EXP26BCLSB",
               "EXP26DE04A", "EXP26DEEX1", "EXP26DESB", "EXP26PE03A", "EXP26PE04A", "EXP26PEEX1", "EXP26PESB", "EXP28BCLEX1",
               "EXP28BCL03A", "EXP28BCL04A", "EXP28BCLSB", "EXP28DE03A", "EXP28DE04A", "EXP28DEEX1", "EXP28DESB", "EXP28PE04A",
               "EXP28PEEX1", "EXP28PESB")


ord_samples = c("EXP22DESB", "EXP26DESB","EXP28DESB",
                "EXP22DEEX1", "EXP26DEEX1","EXP28DEEX1",
                "EXP22DE03A", "EXP28DE03A",
                "EXP22DE04A", "EXP26DE04A", "EXP28DE04A",
                "EXP22PESB", "EXP26PESB", "EXP28PESB",
                "EXP22PEEX1", "EXP26PEEX1", "EXP28PEEX1",
                "EXP22PE03A", "EXP26PE03A",
                "EXP22PE04A", "EXP26PE04A", "EXP28PE04A",
                "EXP22BCLSB", "EXP26BCLSB", "EXP28BCLSB",
                "EXP22BCLEX1", "EXP26BCLEX1", "EXP28BCLEX1",
                "EXP22BCL03A", "EXP26BCL03A", "EXP28BCL03A",
                "EXP22BCL04A", "EXP26BCL04A", "EXP28BCL04A")


me_barplots_col = c(rep(c("dodgerblue4"), times = 3),
                    rep(c("dodgerblue3"), times = 3),
                    rep(c("firebrick3"), times = 2),
                    rep(c("firebrick2"), times = 3),
                    rep(c("dodgerblue4"), times = 3),
                    rep(c("dodgerblue3"), times = 3),
                    rep(c("firebrick3"), times = 2),
                    rep(c("firebrick2"), times = 3),
                    rep(c("dodgerblue4"), times = 3),
                    rep(c("dodgerblue3"), times = 3),
                    rep(c("firebrick3"), times = 3),
                    rep(c("firebrick2"), times = 3))

###################################################



###################################################
# Tidying up inputs
###################################################

#count table
counts = read.delim("~/OneDrive/oxford/summer_internship/deseq2/counts.txt", row.names = 1)
colnames(counts) = cons_names #make names consistent
library(dplyr)
counts = select(counts, ord_samples) #put in order and take out exception

#sample info
colData = read.delim("~/OneDrive/oxford/summer_internship/deseq2/colData.txt", row.names = 1)

#filter out genes with no expression
library(edgeR)
isexpr = rowSums(cpm(counts) > 10) >= 0.5 * ncol(counts)
counts = counts[isexpr, ] #filter for expressed genes

#filter genes with low central tendency i.e. consistently low median expression, likely to represent noise
##see spread of expression
#expr_median = apply(counts, 1, median)
#hist(expr_median)
#plot(expr_median) #still a lot of genes with a median of zero
#library(DGCA)
#library(matrixStats)
#counts = filterGenes(counts,
#                     filterTypes = "central",
#                     filterCentralType = "median",
#                     filterCentralPercentile = 0.25) #filter to 11418 genes

#normalisation using variance stabilising transformation as recommended by WGCNA authors
library(DESeq2)
vsd = vst(as.matrix(counts), blind = FALSE)

#batch correction using ComBat as recommended by WGCNA authors
library(sva)
batch = colData$Experiment
modcombat = model.matrix( ~1, data = colData)
batch_cor = ComBat(dat = vsd, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
batch_cor = t(batch_cor) #transpose

#check for good genes
library(WGCNA)
gsg = goodSamplesGenes(batch_cor, verbose = 3)
gsg$allOK #returned TRUE

#find any sample outliers (shouldn't have any, as seen in earlier PCA plots)
sampleTree = hclust(dist(batch_cor), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "", xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) #no outliers
###################################################


###################################################
#WGCNA Analysis
###################################################

#choose soft-thresholding power (lowest power of the elbow of graph)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(batch_cor, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red") #use soft power of 7

#topological overlap matrix calculation
##restart R and open WGCNA package -- conflicts with other packages
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
softPower = 7
##I have used signed network as recommended by authors
adj = adjacency(batch_cor,type = "unsigned", power = softPower, corFnc = "bicor") #calculate the adjacency matrix, use unsigned because negative correlations are of interest as well, use bicor (biweight mid-correlation) because robust and recommended by authors
TOM = TOMsimilarityFromExpr(batch_cor,
                            corType = "bicor",
                            networkType = "unsigned",
                            TOMType = "unsigned",
                            nThreads = 4,
                            power = softPower) #turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
gene_names = substr(rownames(counts), 1, 15) #get rid of last few digits signifying ensembl version for conversion later on
colnames(TOM) = rownames(TOM) = gene_names #name by gene name
dissTOM = 1-TOM #calculate dissimilarity

#module detection
library(flashClust)
geneTree = flashClust(as.dist(dissTOM), method = "average") #hierarchical clustering of the genes based on the TOM dissimilarity

##plot gene clustering
pdf("wgcna.gene.tree.pdf")
plot(geneTree, xlab = "", sub = "",cex = 0.3) #plot
dev.off()

minModuleSize = 30 #set module size to be bigger (smaller modules are quite uninformative)
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize) #module identification using dynamic tree cut
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods) #put color as labels
table(dynamicColors)

##plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
plotDendroAndColors(geneTree, dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

##plot module eigengene per sample
MEList = moduleEigengenes(batch_cor, colors = dynamicColors) #calculate eigengenes
MEs = MEList$eigengenes
pdf("ME_barplots_unmerged.pdf",width = 10, height = 4)
par(las = 2, mar = c(8,5,4,1))
for (i in 1:length(unique(dynamicColors))){
  which.module = unique(dynamicColors)[i]
  ME = MEs[, paste("ME", which.module, sep = "")]
  barplot(ME, ylab = "eigengene expression", xlab="",
          main = which.module,
          col = me_barplots_col,
          names.arg = rownames(batch_cor), cex.names=0.75)}
dev.off()


#merge modules whose expression are very similar (too many modules!)

##merge highly coexpressed modules by calculating eigengenes and clustering them on their correlation
MEDiss = 1 - cor(MEs) #calculate dissimilarity of module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average") #cluster module eigengenes

##plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

##merge
MEDissThres = 0.1 #merge correlation of 0.9
abline(h = MEDissThres, col = "red")
merge = mergeCloseModules(batch_cor, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

##plot merged results
pdf("merged_cluster_dendogram.pdf", width = 10)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
table(mergedColors)

##plot module eigengene per sample
pdf("ME_barplots_merged.pdf",width = 10, height = 4)
par(las = 2, mar = c(8,5,4,1))
for (i in 1:length(unique(mergedColors))){
  which.module = unique(mergedColors)[i]
  ME = mergedMEs[, paste("ME", which.module, sep = "")]
  barplot(ME, ylab = "eigengene expression", xlab="",
          main = which.module,
          col = me_barplots_col,
          names.arg = rownames(batch_cor), cex.names=0.75)}
dev.off()

#calculate gene information
##prepare data
datExpr = as.data.frame(batch_cor) #put as data frame
nGenes = ncol(datExpr) #no. of genes
nSamples = nrow(datExpr) #no. of samples
modNames = substring(names(MEs), 3) #module names

##calculate gene module membership
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p")) #gene module membership
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #gene module membership p value
names(geneModuleMembership) = paste("MM", modNames, sep = "") #label
names(MMPvalue) = paste("p.MM", modNames, sep = "") #label


#import gene information
##create the starting data frame
annot = read.delim("~/OneDrive/oxford/summer_internship/wgcna/gene_annotation.txt")
genes2annot = match(rownames(counts), annot$GeneID)
sum(is.na(genes2annot))
geneInfo = data.frame(ensembl = annot$GeneID[genes2annot],
                      geneSymbol = annot$GeneName[genes2annot],
                      moduleColor = mergedColors)
write.csv(geneInfo, file = "geneInfo.csv")
#for(mod in 1: ncol(geneModuleMembership)) {
#  oldNames = names(geneInfo)
#  geneInfo = data.frame(geneInfo,
#                        geneModuleMembership,
#                        MMPvalue)
#  names(geneInfo) = c(oldNames,
#                      paste("MM.", modNames, sep = ""),
#                      paste("p.MM", modNames, sep = ""))
#} #to import module membership for each gene


##inspect
geneInfo[which(geneInfo$geneSymbol == "HNF4A"), ] #HNF4A is in blue module
geneInfo[which(geneInfo$geneSymbol == "HNF1A"), ] #HNF1A is in blue module


save(MEs, mergedMEs, dynamicColors, mergedColors, geneTree,
     file = "wgcna_construction.RData")
###################################################

###################################################
#WGCNA Gene Ontology Enrichment Analysis
###################################################
options(stringsAsFactors = FALSE)
library("anRichment")

#data prep
data = read.csv("~/OneDrive/Oxford/summer_internship/wgcna/geneInfo.csv", header = TRUE, rownames = 1)
moduleColor = data$moduleColor
entrez = convert2entrez(organism = "human", symbol = data$GeneName) #convert to entrez
table(is.finite(entrez)) #how many conversions are successful?

#make reference collection (GO and NCBI biosystems collection)
GOcollection = buildGOcollection(organism = "human")
biosysCollection = BioSystemsCollection("human") #KEGG, REACTOME, BIOCYC and Lipid Maps
combinedCollection = mergeCollections(GOcollection, biosysCollection)

#one-by-one analysis of interesting modules
active = entrez[moduleColor == "blue"]
all = entrez

GOenrichment_blue = enrichmentAnalysis(active = active,
                                       inactive = all,
                                       refCollection = combinedCollection,
                                       useBackground = "intersection",
                                       threshold = 0.01,
                                       thresholdType = "Bonferroni")


##explore results
summary = GOenrichment_blue$enrichmentTable
summary$overlapGenes = shortenStrings(summary$overlapGenes, maxLength = 70, split = "|");
head(summary)

#export
write.csv(summary$enrichmentTable, file = "wgcna_GOenrichmentTable_blue.csv",
          row.names = FALSE)
###################################################
