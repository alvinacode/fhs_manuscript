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
isexpr = rowSums(cpm(counts) > 1) >= 0.5 * ncol(counts)
counts = counts[isexpr, ] #filter for expressed genes

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


#merge modules whose expression are very similar

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
geneTreeMerge = merge$dendro

##plot merged results
pdf("merged_module_cluster.pdf")
plot(geneTreeMerge, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()

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

#calculate gene module membership
geneModuleMembership = as.data.frame(cor(datExpr, mergedMEs, use = "p")) #gene module membership
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #gene module membership p value
names(geneModuleMembership) = paste("MM", modNames, sep = "") #label
names(MMPvalue) = paste("p.MM", modNames, sep = "") #label

#relating modules to external clinical traits
##define numbers of genes and samples
datExpr = as.data.frame(batch_cor) #put as data frame
nGenes = ncol(datExpr) #no. of genes
nSamples = nrow(datExpr) #no. of samples
modNames = substring(names(mergedMEs), 3) #module names
mergedMEs2 = orderMEs(mergedMEs) #ordered

##create data frame of disease status
datTraits = as.data.frame(rownames(datExpr))
datTraitsDEPE = data.frame(bla = substring(datTraits[1:22,1], 8))
datTraitsBLC = data.frame(bla = substring(datTraits[23:34,1], 9))
datTraits = rbind(datTraitsDEPE, datTraitsBLC)
datTraits[,1] = gsub("EX1", "0", datTraits[,1]) #binary classification: WT = 0, MODY = 1
datTraits[,1] = gsub("SB", "0", datTraits[,1])
datTraits[,1] = gsub("03A", "1", datTraits[,1])
datTraits[,1] = gsub("04A", "1", datTraits[,1])
rownames(datTraits) = c() #no rownames for plot later on

##calculate correlation with MODY
moduleTraitCor = cor(mergedMEs2, (as.numeric(datTraits$bla)), use = "p") #pearson correlation of ME with MODY
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

##plot
pdf("moduleTraitRelationships.pdf", width = 10, height = 10)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "",
               yLabels = names(mergedMEs2),
               ySymbols = names(mergedMEs2),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


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

#analyse hub genes
##put in right format
nSets = 1
setLabels = c("HNF4A")
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(batch_cor))
names(multiExpr[[1]]$data) = colnames(batch_cor)
rownames(multiExpr[[1]]$data) = rownames(batch_cor)
checkSets(multiExpr, checkStructure = TRUE)

##calculate gene significance
GS1 = as.numeric(cor(batch_cor, as.numeric(datTraits$bla), use = "p"))
GeneSignificance = abs(GS1)

##plot
alldegrees1 = intramodularConnectivity(adj, mergedColors) #calculate the intramodular connectivity for each gene
colorlevels = unique(mergedColors)
sizeGrWindow(9, 6)
par(mfrow = c(2, as.integer(0.5 + length(colorlevels) / 2)))
par(mar = c(4,5,3,1))
pdf("temp.pdf")
for (i in c(1:length(colorlevels))) {
  which.module = colorlevels[[i]];
  restrict1 = (mergedColors == which.module);
  verboseScatterplot(alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col = mergedColors[restrict1],
                     main = which.module,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}
dev.off()

##filter for hub genes
datKME = signedKME(batch_cor, mergedMEs, outputColumnName = "MM.") #calcular module membership
filterGenes = abs(GS1) > 0.20 & abs(datKME$MM.violet) > 0.95
table(filterGenes)
string = dimnames(data.frame(datExpr))[[2]][filterGenes]
string = annot[which(annot$GeneID %in% string), 2]
string
write.table(string, "hubGenes_violet_GS20_KME95.txt", quote = F, sep = "\t", row.names = FALSE, col.names = FALSE)

#export to Cytoscape to visualise
##select module genes
modules = c("violet")
probes = names(datExpr)
inModule = is.finite(match(dynamicColors, modules))
modProbes = probes[inModule]
modGenes = annot$GeneName[match(modProbes, annot$GeneID)]

##select corresponding TOM
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-20", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-20", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2, #very high threshold because too many genes (chosen from barplot above)
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = dynamicColors[inModule])


save(MEs, mergedMEs, dynamicColors, mergedColors, geneTree,
     file = "wgcna_construction.RData")
###################################################

###################################################
#WGCNA Gene Ontology Enrichment Analysis
###################################################
options(stringsAsFactors = FALSE)
library("anRichment")

#data prep
data = read.csv("~/OneDrive/Oxford/summer_internship/wgcna/geneInfo.csv", header = TRUE, row.names = 1)
moduleColor = data$moduleColor
entrez = convert2entrez(organism = "human", symbol = data$geneSymbol) #convert to entrez
table(is.finite(entrez)) #how many conversions are successful?

#make reference collection (GO and NCBI biosystems collection)
load("~/OneDrive/oxford/summer_internship/wgcna/combinedCollection_GOBPonly.RData")

GOenrichment = enrichmentAnalysis(classLabels = moduleColor,
                                  identifiers = entrez,
                                  refCollection = combinedCollection,
                                  useBackground = "given",
                                  threshold = 0.05,
                                  nBestDataSets = 50,
                                  thresholdType = "Bonferroni",
                                  getOverlapEntrez = TRUE,
                                  getOverlapSymbols = TRUE,
                                  ignoreLabels = "grey")
collectGarbage()

write.csv(GOenrichment$enrichmentTable, file = "GOenrichment-enrichmentTable.csv", row.names = FALSE)
###################################################
