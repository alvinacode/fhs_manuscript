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
# Functions
###################################################
#function to compute the dispersion value that quantifies the correlation between 2 conditions
#for a pair of genes drawn from module c1 and module c2
#in case c1 = c2, the function quantifies the differential coexpression in c1
dispersionModule2Module = function(c1, c2, datC1, datC2, colorh1C1C2) {
  if(c1 == c2) {
    difCor = (cor(datC1[,which(colorh1C1C2 == c1)], method = "spearman") -
                cor(datC2[,which(colorh1C1C2 == c1)], method = "spearman"))^2
    n = length(which(colorh1C1C2 == c1))
    (1/((n^2 - n)/2) * (sum(difCor)/2))^(0.5)
  }
  else if(c1 != c2) {
    difCor = (cor(datC1[,which(colorh1C1C2 == c1)], datC1[,which(colorh1C1C2 == c2)], method = "spearman") -
                cor(datC2[, which(colorh1C1C2 == c1)], datC2[,which(colorh1C1C2 == c2)], method = "spearman"))^2
    n1 = length(which(colorh1C1C2 == c1))
    n2 = length(which(colorh1C1C2 == c2))
    (1/((n1 * n2)) * (sum(difCor)))^(0.5)
  }
}

#function to calculate the dispersion value of a module to module coexpression change on permuted data
permutationProcedureModule2Module = function(permutation, d, c1, c2, colorh1C1C2) {
  d1 = d[permutation,]
  d2 = d[-permutation,]
  dispersionModule2Module(c1, c2, d1, d2, colorh1C1C2)
}

#function to plot permutation matrix
plotMatrix = function(mat) {
  mat[which(row(mat) > col(mat))] = 1001
  image(mat,
        col = c(gray.colors(4), "white"),
        breaks = c(0,0.1,50,100,1000,1001),
        xaxt = 'n',
        yaxt = 'n',
        xlim = c(-0.2,1.2),
        ylim = c(-0.2,1.2),
        bty = 'n',
        asp = 1)
  text(0:(nrow(mat)-1) / (nrow(mat)-1),
       1.1,
       rownames(mat),
       cex = 1,
       col = rownames(mat))
  text(-0.15,
       0:(ncol(mat)-1) / (ncol(mat)-1),
       colnames(mat),
       cex = 1,
       col = colnames(mat))
  text(apply(matrix(0:(nrow(mat)-1) / (nrow(mat)-1)), 1, rep, ncol(mat)),
       rep(0:(ncol(mat)-1) / (ncol(mat)-1), nrow(mat)),
       as.numeric(t(mat)),
       col = "white",
       cex = 1.5)
}

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
###################################################

###################################################
# DiffCoEx 2010
###################################################
#make matrices
datC1 = batch_cor[c(7:11, 18:22, 29:34),] #mody samples
datC2 = batch_cor[c(1:6, 12:17, 23:28),] #wt samples

adjMatC1 = sign(cor(datC1, method = "spearman")) * (cor(datC1, method = "spearman"))^2 #make adjacency matrix for mutants
adjMatC2 = sign(cor(datC2, method = "spearman")) * (cor(datC2, method = "spearman"))^2 #make adjacency matrix for wt
diag(adjMatC1) = 0
diag(adjMatC2) = 0
collectGarbage()

#pick softpower threshold
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(datC1, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

sft = pickSoftThreshold(datC2, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
softPower = 10
dissTOMC1C2 = TOMdist((abs(adjMatC1 - adjMatC2)/ 2)^(softPower/ 2)) #compute TOM based on dissimilarity matrix
collectGarbage()

library(flashClust)
geneTreeC1C2 = flashClust(as.dist(dissTOMC1C2), method = "average") #hierarchical clustering using the Topological Overlap of the adjacency difference as input distance matrix

#plot the resulting clustering tree (dendrogram)
pdf("diffCoExpr_hierarchical_tree.pdf", height = 1000, width = 1000)
plot(geneTreeC1C2, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
dev.off()

#module detection
dynamicModsHybridC1C2 = cutreeDynamic(dendro = geneTreeC1C2,
                                      distM = dissTOMC1C2,
                                      method = "hybrid",
                                      cutHeight=.996,
                                      deepSplit = 2,
                                      pamRespectsDendro = F,
                                      minClusterSize = 30)

dynamicColorsHybridC1C2 = labels2colors(dynamicModsHybridC1C2)
colorh1C1C2 = dynamicColorsHybridC1C2
table(dynamicModsHybridC1C2)
table(dynamicColorsHybridC1C2)

#merge similar modules
#mergeC1C2 = mergeCloseModules(rbind(datC1, datC2),
#                              dynamicColorsHybridC1C2,
#                              cutHeight = 0.1) #0.9 correlation
#colorh1C1C2 = mergeC1C2$color
#mergedC1C2MEs = mergeC1C2$newMEs
#rownames(mergedC1C2MEs) = rownames(rbind(datC1, datC2))
#table(colorh1C1C2)

#mergedC1C2MEs = as.data.frame(t(mergedC1C2MEs))
#mergedC1C2MEs = as.data.frame(t(select(mergedC1C2MEs, ord_samples))) #put in order

#pdf("diffCoExpr_barplots_trial.pdf",width = 10, height = 4)
#par(las = 2, mar = c(8,5,4,1))
#for (i in 1:length(unique(colorh1C1C2))){
#  which.module = unique(colorh1C1C2)[i]
#  ME = mergedC1C2MEs[, paste("ME", which.module, sep = "")]
#  barplot(ME, ylab = "eigengene expression", xlab="",
#          main = which.module,
#          col = me_barplots_col,
#          names.arg = rownames(mergedC1C2MEs), cex.names = 0.75)}
#dev.off()

#significance testing using permutations
permutations = NULL
for(i in 1:1000) {
  permutations = rbind(permutations, sample(1:(nrow(datC1) + nrow(datC2)), nrow(datC1)))
} #generate a set of 1000 permuted indexes
d = rbind(scale(datC1), scale(datC2)) #scale the data in both conditions to mean 0 and variance 1

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
dispersion_matrix = matrix(nrow = length(unique(colorh1C1C2)) - 1, ncol = length(unique(colorh1C1C2)) - 1)
null_distrib = list()
i <- j <- 0
for (c1 in setdiff(unique(colorh1C1C2), "grey")) {
  i = i + 1
  j = 0
  null_distrib[[c1]] = list()
  for (c2 in setdiff(unique(colorh1C1C2), "grey")) {
    j = j + 1
    dispersion_matrix[i,j] = dispersionModule2Module(c1, c2, datC1, datC2, colorh1C1C2)
    null_distrib[[c1]][[c2]] = apply(permutations, 1, permutationProcedureModule2Module, d, c2, c1, colorh1C1C2)
  }
} #compute all pairwise module to module dispersion values, and generate a null distribution from permuted scaled data (very slow)

#create a summary matrix indicating for each module to module
#differential coexpression the number of permuted data yielding
#an equal or higher dispersion
permutationSummary = matrix(nrow = 8, ncol = 8) #n is number of modules
colnames(permutationSummary) = setdiff(unique(colorh1C1C2), "grey")
rownames(permutationSummary) = setdiff(unique(colorh1C1C2), "grey")
for(i in 1:8) {
  for(j in 1:8) {
    permutationSummary[i,j] = length(which(null_distrib[[i]][[j]] >= dispersion_matrix[i,j]))
  }
}

pdf("diffCoExpr_significance.pdf")
plotMatrix(permutationSummary) #plot the result
dev.off()

#visualising module eigengene per sample
MEListC1 = moduleEigengenes(datC1, colors = colorh1C1C2) #calculate eigengenes
MEsC1 = MEListC1$eigengenes
rownames(MEsC1) = rownames(datC1)
MEListC2 = moduleEigengenes(datC2, colors = colorh1C1C2)
MEsC2 = MEListC2$eigengenes
rownames(MEsC2) = rownames(datC2)
MEsC1C2 = as.data.frame(t(rbind(MEsC1, MEsC2)))
MEsC1C2 = as.data.frame(t(select(MEsC1C2, ord_samples))) #put in order

pdf("diffCoExpr_barplots_trial.pdf",width = 10, height = 4)
par(las = 2, mar = c(8,5,4,1))
for (i in 1:length(unique(colorh1C1C2))){
  which.module = unique(colorh1C1C2)[i]
  ME = MEsC1C2[, paste("ME", which.module, sep = "")]
  barplot(ME, ylab = "eigengene expression", xlab="",
          main = which.module,
          col = me_barplots_col,
          names.arg = rownames(MEsC1C2), cex.names=0.75)}
dev.off()

#export gene information
annot = read.delim("~/OneDrive/oxford/summer_internship/wgcna/gene_annotation.txt")
genes2annot = match(colnames(datC1), annot$GeneID)
sum(is.na(genes2annot))
geneInfo_diffCoExpr = data.frame(ensembl = annot$GeneID[genes2annot],
                                 geneSymbol = annot$GeneName[genes2annot],
                                 moduleColor = colorh1C1C2)
write.csv(geneInfo_diffCoExpr, file = "geneInfo_diffCoExpr.csv")

save(datC1, datC2, adjMatC1, adjMatC2, colorh1C1C2, geneTreeC1C2,
     null_distrib, dispersion_matrix,
     file = "diffCoExpr_construction.RData")
###################################################
