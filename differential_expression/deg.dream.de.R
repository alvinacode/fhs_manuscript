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
counts = read.delim("~/OneDrive/oxford/summer_internship/deseq2/counts_de.txt", header = TRUE, row.names = 1)
colnames(counts)[2] = c("EXP22DE03A")
colnames(counts)[3] = c("EXP22DE04A")
rownames(counts) = substr(rownames(counts), 1, 15)
##make names stuff for easier reading later
#gene_name = read.delim("~/OneDrive/oxford/summer_internship/deseq2/gene_name.txt", header = TRUE, row.names = NULL)
#gene_name = cbind(rownames(counts), gene_name)
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
de_all = topTable(fit_bayes, p.value = 1, number = Inf)
write.csv(deg, "dream.deg.de.07122018.csv")
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
exportFASTA(pro_seq, "de.pro.fasta.txt")

#find genes with HNF4A motif on promoter with FIMO (on WTCHG server)
# export PATH=${PATH}:/apps/well/meme/4.11.2_2/bin
# fimo <your_motifs.meme> <your_fasta_file.fasta>

#identify deg with HNF4A motif
fimo = read.delim("~/OneDrive/oxford/summer_internship/differential_gene_expression/de.fimo.txt") #import fimo results 
colnames(fimo) = c("pattern.name", "seq.name", "start", "stop", "strand", "score", "p.val", "q.val", "matched.seq") #rename columns
fimo_unique = unique(fimo, by = "seq.name") #drop repeated genes
deg_hnf4a_t = deg[deg$entrez %in% fimo_unique$seq.name,] #reduced to 38 deg

#explore and export 
write.csv(deg_hnf4a_t, "de.hnf4a.target.all.csv") #all (regardless of direction)

library(dplyr)
deg_hnf4a_t_down = filter(deg_hnf4a_t, logFC < 0) #downregulated genes
write.csv(deg_hnf4a_t_down, "de.hnf4a.target.down.csv")

library(dplyr)
deg_hnf4a_t_up = filter(deg_hnf4a_t, logFC > 0) #upregulated genes
write.csv(deg_hnf4a_t_up, "de.hnf4a.target.up.csv")

##########################################################

##########################################################
# Gene Ontology Enrichment Analysis Preparation
##########################################################
#create background genes i.e. genes that are tested and then filtered for HNF4A motif
genes_tested = data.frame(X = rownames(batch_corrected)) #i.e. the 8080 genes passed through deg analysis 
genes_tested$GeneSymbol= geneInfo[match(genes_tested$X, geneInfo$GeneID), 2] #get gene symbol
genes_tested$entrez = convert2entrez(organism = "human", symbol = genes_tested$GeneSymbol) #get entrez id
table(is.finite(genes_tested$entrez)) #how many conversions are successful? (170 unsuccessful)

#get promoter sequences of tested genes 
universe_pro_seq = getSequence(id = genes_tested$entrez,
                               type = "entrezgene",
                               seqType = "coding_gene_flank",
                               upstream = 2000, 
                               mart = mart) #slow, probably a bad idea
exportFASTA(universe_pro_seq, "universe.pro.fasta.txt")

#FIMO analysis 

#identify tested genes with HN4A motif
universe_fimo = read.delim("~/OneDrive/oxford/summer_internship/differential_gene_expression/universe.fimo.txt") #import fimo results 
colnames(universe_fimo) = c("pattern.name", "seq.name", "start", "stop", "strand", "score", "p.val", "q.val", "matched.seq") #rename columns
universe_fimo_unique = unique(universe_fimo, by = "seq.name") #drop repeated genes
universe_hnf4a_t = genes_tested[genes_tested$entrez %in% universe_fimo_unique$seq.name,] #save for gene ontology 
save(universe_hnf4a_t, file = "background.go.enrichment.Rdata")
##########################################################
