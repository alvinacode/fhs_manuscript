##########################################################
# Functions
##########################################################
order_by_donors.hnf4a = function(counts, donors. = hnf4a.donors) {
  nc = counts[, grepl("SB" , names(counts))]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in donors.[-1])  {
    i = counts[, grepl(s , names(counts))]
    nc = cbind(nc, i)
  }
  
  return(nc)
}

same_donor = function(x) {
  x = gsub("SB", "D1", x)
  x = gsub("EX1", "D1", x)
  x = gsub("03A", "D2", x)
  x = gsub("04A", "D2", x)
}

same_disease = function(x) {
  x = gsub("SB", "0", x)
  x = gsub("EX1", "0", x)
  x = gsub("03A", "1", x)
  x = gsub("04A", "1", x)
}
##########################################################

##########################################################
# Getting differentially expressed genes 
##########################################################
#count matrix
library(edgeR)
counts = read.delim("~/OneDrive/oxford/summer_internship/deseq2/counts_pe.txt", row.names = 1)
counts = counts[, -9] #remove outlier 
rownames(counts) = substr(rownames(counts), 1, 15)
colnames(counts)[2] = c("EXP22PE03A")
colnames(counts)[3] = c("EXP22PE04A")
##tidy up
hnf4a.donors = c("SB", "EX1", "03A", "04A")  #data comes from four 'donors' (2 indiduals, 4 cell clones)
counts = order_by_donors.hnf4a(counts) #order


#design matrix 
experiment = substr(colnames(counts), 1, 5)
cell_line = substr(colnames(counts), 8, nchar(colnames(counts)))
hnf4a_donor = substr(colnames(counts), 8, nchar(colnames(counts)))
donor = same_donor(hnf4a_donor)
disease = same_disease(hnf4a_donor)
##make metadata data frame 
metadata = data.frame(Run = colnames(counts),
                      Experiment = experiment, 
                      CellLine = cell_line,
                      Donor = donor,
                      Disease = disease) #form colData
rownames(metadata) = metadata$Run
metadata = metadata[, 2:5]

#filter
isexpr = rowSums(cpm(counts) > 1)>= 0.5 * ncol(counts)
counts = counts[isexpr, ] #filter for expressed genes

#normalisation
dge = DGEList(counts = counts, samples = metadata)
dge = calcNormFactors(dge) #TMM normalisation
design = model.matrix( ~ Disease + Experiment, metadata)
v = voom(dge, design, plot = F)

#correct for batch and same donor effect using linear mixed model
form = ~ Disease + (1|Donor) + (1|Experiment)
L = getContrast(v, form, metadata, "Disease1")
library(doParallel)
library(variancePartition)
cl <- makeCluster(4)
registerDoParallel(cl)
fit = dream(v, form, metadata, L)
fit_bayes = eBayes(fit) #fit empirical bayes for moderated t-statistics
deg = topTable(fit_bayes, p.value = 0.05, number = Inf)
pe_all = topTable(fit_bayes, p.value = 1, number = Inf)
write.csv(deg, "dream.deg.pe.07122018.csv")
##########################################################


##########################################################
# Filter for genes with HNF4A motif on promoter
##########################################################
#get gene symbol and ensembl 
deg$X = rownames(deg)
geneInfo = read.delim("~/OneDrive/oxford/summer_internship/differential_gene_expression/gene_annotation.txt")
geneInfo$GeneID = substr(geneInfo$GeneID, 1, 15)
deg$GeneSymbol = geneInfo[match(deg$X, geneInfo$GeneID), 2] #get gene symbol
library(anRichmentMethods)
deg$entrez = convert2entrez(organism = "human", symbol = deg$GeneSymbol) #get entrez id
table(is.finite(deg$entrez)) #how many conversions are successful? (15 unsuccessful)
deg = deg[,c(7:9, 1:6)] #reorder

#get promoter sequences of deg 
library(biomaRt)
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
pro_seq = getSequence(id = deg$entrez,
                      type = "entrezgene",
                      seqType = "coding_gene_flank",
                      upstream = 2000, 
                      mart = mart)
exportFASTA(pro_seq, "pe.pro.fasta.txt")

#find genes with HNF4A motif on promoter with FIMO (on WTCHG server)
# export PATH=${PATH}:/apps/well/meme/4.11.2_2/bin
# fimo <your_motifs.meme> <your_fasta_file.fasta>

#identify deg with HNF4A motif
fimo = read.delim("~/OneDrive/oxford/summer_internship/differential_gene_expression/pe.fimo.txt") #import fimo results 
colnames(fimo) = c("pattern.name", "seq.name", "start", "stop", "strand", "score", "p.val", "q.val", "matched.seq") #rename columns
fimo_unique = unique(fimo, by = "seq.name") #drop repeated genes
deg_hnf4a_t = deg[deg$entrez %in% fimo_unique$seq.name,] #reduced to 38 deg

#explore and export 
write.csv(deg_hnf4a_t, "pe.hnf4a.target.all.csv") #all (regardless of direction)

library(dplyr)
deg_hnf4a_t_down = filter(deg_hnf4a_t, logFC < 0) #downregulated genes
write.csv(deg_hnf4a_t_down, "pe.hnf4a.target.down.csv")

library(dplyr)
deg_hnf4a_t_up = filter(deg_hnf4a_t, logFC > 0) #upregulated genes
write.csv(deg_hnf4a_t_up, "pe.hnf4a.target.up.csv")

##########################################################
