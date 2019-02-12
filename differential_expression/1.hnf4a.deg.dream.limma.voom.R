setwd("~/OneDrive/oxford/summer_internship/differential_gene_expression")
library(edgeR)
library(limma)
library(variancePartition)
library(xlsx)
library(dplyr)
library(anRichmentMethods) 
library(biomaRt)
library(VennDiagram)
library(ggsci)

##########################################################
# Functions
##########################################################
order_by_donors.hnf4a = function(counts, donors. = hnf4a.donors) {
  nc = counts[, grepl("SB" , names(counts))]  # takes columns whose name matches x
  
  # all stages 
  
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

volcanoPlot = function(res) {
  with(res, plot(logFC, -log10(adj.P.Val), 
                 pch = 20, main = "Volcano plot", 
                 xlim = c(-2.5, 2))) # black points 
  
  with(subset(res, adj.P.Val < 0.05 ), 
       points(logFC, -log10(adj.P.Val), 
              pch = 20, col = "red")) # red points if adj. p-val < 0.05 
  
  with(subset(res, abs(logFC) > 1), 
       points(logFC, -log10(adj.P.Val), 
              pch = 20, col = "orange")) # orange points if abs(log2FC) > 1
  
  with(subset(res, adj.P.Val < 0.05 & abs(logFC) > 1), 
       points(logFC, -log10(adj.P.Val), 
              pch = 20, col = "green")) # green for significant DEG 
  
  # pink for DEG with HNF4A motif
}

##########################################################

##########################################################
# Prepare data (count matrix and sample info)
##########################################################

# put together count matrix 

counts = read.delim("~/OneDrive/oxford/summer_internship/counts/hnf4a_tidyCounts.txt", header = TRUE, row.names = 1, check.names = FALSE)
counts = counts[,-which(colnames(counts) %in% c("EXP28PE03A"))] # remove outlier

de = counts[,1:11]
pe = counts[,12:22]
blc = counts[,23:34]

count_list = list(de = as.matrix(de),
                  pe = as.matrix(pe), 
                  blc = as.matrix(blc))

# put together sample info
metadata = list()
for(i in 1:3) {
  Experiment = substr(colnames(count_list[[i]]), 1, 5)
  CellLine = substr(colnames(count_list[[i]]), 8, nchar(colnames(count_list[[i]])))
  Donor = substr(colnames(count_list[[i]]), 8, nchar(colnames(count_list[[i]])))
  Disease = substr(colnames(count_list[[i]]), 8, nchar(colnames(count_list[[i]])))
  metadata[[i]] = data.frame(Experiment,
                             CellLine,
                             Donor,
                             Disease)
  rownames(metadata[[i]]) = colnames(count_list[[i]])
}
names(metadata) = c("de", "pe", "blc")

# blc only needs editing because diff. no. of characters i.e. DE & PE (2) vs BLC (3)
metadata$blc$CellLine = as.character(metadata$blc$CellLine)
metadata$blc$CellLine = substr(metadata$blc$CellLine, 2, nchar(metadata$blc$CellLine))

metadata$blc$Donor = as.character(metadata$blc$Donor)
metadata$blc$Donor = substr(metadata$blc$Donor, 2, nchar(metadata$blc$Donor))

metadata$blc$Disease = as.character(metadata$blc$Disease)
metadata$blc$Disease = substr(metadata$blc$Disease, 2, nchar(metadata$blc$Disease))


for(i in 1:3) {
  metadata[[i]]$Donor = same_donor(metadata[[i]]$Donor) # samples derived from 2 individuals
  metadata[[i]]$Disease = same_disease(metadata[[i]]$Disease) # samples derived from WT and MODY
}

# filter counts 
for(i in 1:3) {
  isexpr = rowSums(cpm(count_list[[i]]) > 1) >= 3 # expressed at least in one sample group e.g. EX1
  count_list[[i]] = count_list[[i]][isexpr, ]
}

# normalisation 
dge = list()
v = list()
for(i in 1:3) {
  dge[[i]] = DGEList(counts = count_list[[i]], samples = metadata[[i]])
  dge[[i]] = calcNormFactors(dge[[i]]) # TMM normalisation 
  design = model.matrix(~ Disease + Experiment, metadata[[i]])
  v[[i]] = voom(dge[[i]], design, plot = F) # normalise reads taking acct of TMM factors
}

##########################################################

##########################################################
# dream + limma + voom 
##########################################################

# correct for batch and same donor effects 
library(doParallel) 
cl <- makeCluster(4)
registerDoParallel(cl)

form = ~ Disease + (1|Donor) + (1|Experiment)

fit_bayes = list()
for(i in 1:3) {
  L = getContrast(v[[i]], form, metadata[[i]], "Disease1")
  fit = dream(v[[i]], form, metadata[[i]], L)
  fit_bayes[[i]] = eBayes(fit) # fit empirical bayes for moderated t-statistics 
}

L = getContrast(v[[1]], form, metadata[[1]], "Disease1")
fit = dream(v[[1]], form, metadata[[1]], L)
fit_bayes[[1]] = eBayes(fit) # fit empirical bayes for moderated t-statistics 


##########################################################

##########################################################
# DEG results 
##########################################################

## NB that LFC in limma is actually log2FC i.e. LFC < 0 is actually log2FC < 0, which means it's actually actually logFC < 1 

results = list()
for(i in 1:3) {
  results[[i]] = data.frame(all_genes = topTable(fit_bayes[[i]], p.value = 1, number = Inf), # for plots
                            deg_genes = topTable(fit_bayes[[i]], p.value = 1, number = Inf)) # for further analysis  
}


write.xlsx(results[[1]], file = "dream.deg.xlsx", sheetName = "de_unfiltered")
write.xlsx(results[[2]], file = "dream.deg.xlsx", sheetName = "pe_unfiltered", append = TRUE)
write.xlsx(results[[3]], file = "dream.deg.xlsx", sheetName = "blc_unfiltered", append = TRUE)

deg = list()
for(i in 1:3) {
  deg[[i]] = topTable(fit_bayes[[i]], p.value = 1, number = Inf) # easier to analyse in this format
}
names(deg) = c("de", "pe", "blc")

##########################################################

##########################################################
# Filter for genes with HNF4A motif on promoter
##########################################################
geneInfo = read.delim("~/OneDrive/oxford/summer_internship/differential_gene_expression/gene_annotation.txt")
geneInfo$GeneID = substr(geneInfo$GeneID, 1, 15)

for(i in 1:3) {
  deg[[i]]$X = rownames(deg[[i]])
  deg[[i]]$GeneSymbol = geneInfo[match(deg[[i]]$X, geneInfo$GeneID), 2] # get gene symbol 
  deg[[i]]$entrez = convert2entrez(organism = "human", symbol = deg[[i]]$GeneSymbol) # get entrex id
}

# inspect and reorder 
deg$de = deg$de[,c(7:9, 1:6)] 
deg$pe = deg$pe[,c(7:9, 1:6)] 
deg$blc = deg$blc[,c(7:9, 1:6)] 

# get promoter sequences of deg 
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

pro_seq = list()
for(i in 1:3) {
  pro_seq[[i]] = getSequence(id = deg[[i]]$entrez,
                             type = "entrezgene",
                             seqType = "coding_gene_flank",
                             upstream = 2000,
                             mart = mart)
}

exportFASTA(pro_seq[[1]], "de.pro.fasta.txt")
exportFASTA(pro_seq[[2]], "de.pro.fasta.txt")
exportFASTA(pro_seq[[3]], "de.pro.fasta.txt")

# find genes with HNF4A motif on promoter with FIMO (on WTCHG server)

# import results

de_fimo = read.delim("de.fimo.txt")
pe_fimo = read.delim("pe.fimo.txt")
blc_fimo = read.delim("blc.fimo.txt")

fimo = list(de_fimo, pe_fimo, blc_fimo)
deg_hnf4a_t = list()
for(i in 1:3) {
  colnames(fimo[[i]]) = c("pattern.name", "seq.name", "start", "stop", "strand", "score", "p.val", "q.val", "matched.seq") 
  fimo_unique = unique(fimo[[i]], by = "seq.name") # drop repeated genes
  deg_hnf4a_t[[i]] = deg[[i]][deg[[i]]$entrez %in% fimo_unique$seq.name,] # deg with hnf4a motif
}

# export results 

# all the deg hnf4a
write.csv(deg_hnf4a_t[[1]], "de.hnf4a.target.all.csv") 
write.csv(deg_hnf4a_t[[2]], "pe.hnf4a.target.all.csv") 
write.csv(deg_hnf4a_t[[3]], "blc.hnf4a.target.all.csv") 

# downregulated
deg_hnf4a_t_down1 = filter(deg_hnf4a_t[[1]], logFC < 0) 
write.csv(deg_hnf4a_t_down1, "de.hnf4a.target.down.csv")
## ... and so on 

# upregulated 
deg_hnf4a_t_up1 = filter(deg_hnf4a_t[[1]], logFC > 0) 
## ... and so forth 
##########################################################

##########################################################
# Visualisation / Plots 
##########################################################

# simple venn diagram to see distribution 
# a fast (but bad plot) tool is venny: http://bioinfogp.cnb.csic.es/tools/venny/ 
# plot using VennDiagram package is better looking 

venn.pal = pal_simpsons("springfield", alpha = 0.5)(3) # palette

# ALL (not just HNF4A targets) deg across development stages 
pdf("all.deg.pdf")
all.deg <- draw.triple.venn(area1 = 33 + 12 + 33 + 11,
                            area2 = 931 + 12 + 33 + 375,
                            area3 = 1163 + 11 + 33 + 375,
                            n12 = 12 + 33,
                            n13 = 11 + 33,
                            n23 = 375 + 33,
                            n123 = 33,
                            category = c("DE", "PE", "BLC"),
                            fill = venn.pal,
                            alpha = 1,
                            fontface = "bold",
                            fontfamily = "sans",
                            cex = 3.5,
                            cat.fontface = "bold",
                            cat.fontfamily = "sans",
                            cat.cex = 4,
                            cat.default.pos = "outer",
                            lty = "blank",
                            overrideTriple = 1) 
dev.off()
grid.newpage()

# HNF4A DEG only
pdf("all.hnf4a.t.deg.pdf")
all.deg <- draw.triple.venn(area1 = 15 + 4 + 13 + 6,
                            area2 = 495 + 4 + 13 + 188,
                            area3 = 587 + 6 + 13 + 188,
                            n12 = 4 + 13,
                            n13 = 6 + 13,
                            n23 = 188 + 13,
                            n123 = 13,
                            category = c("DE", "PE", "BLC"),
                            fill = venn.pal,
                            alpha = 1,
                            fontface = "bold",
                            fontfamily = "sans",
                            cex = 3.5,
                            cat.fontface = "bold",
                            cat.fontfamily = "sans",
                            cat.cex = 4,
                            cat.default.pos = "outer",
                            lty = "blank",
                            overrideTriple = 1) #all deg across development stages (cut-off LFC = 0)
dev.off()

# HNF4A DEG downregulated
pdf("all.hnf4a.t.deg.down.pdf")
all.deg <- draw.triple.venn(area1 = 7 + 2 + 9 + 4,
                            area2 = 287 + 2 + 9 + 57,
                            area3 = 322 + 4 + 9 + 57,
                            n12 = 2 + 9,
                            n13 = 4 + 9,
                            n23 = 57 + 9,
                            n123 = 9,
                            category = c("DE", "PE", "BLC"),
                            fill = venn.pal,
                            alpha = 1,
                            fontface = "bold",
                            fontfamily = "sans",
                            cex = 3.5,
                            cat.fontface = "bold",
                            cat.fontfamily = "sans",
                            cat.cex = 4,
                            cat.default.pos = "outer",
                            lty = "blank",
                            overrideTriple = 1) #all deg across development stages (cut-off LFC = 0)
dev.off()

# HNF4A DEG upregulated 
pdf("all.hnf4a.t.deg.up.pdf")
all.deg <- draw.triple.venn(area1 = 8 + 2 + 4 + 2,
                            area2 = 238 + 2 + 4 + 101,
                            area3 = 295 + 2 + 4 + 101,
                            n12 = 2 + 4,
                            n13 = 2 + 4,
                            n23 = 101 + 4,
                            n123 = 4,
                            category = c("DE", "PE", "BLC"),
                            fill = venn.pal,
                            alpha = 1,
                            fontface = "bold",
                            fontfamily = "sans",
                            cex = 3.5,
                            cat.fontface = "bold",
                            cat.fontfamily = "sans",
                            cat.cex = 4,
                            cat.default.pos = "outer",
                            lty = "blank",
                            overrideTriple = 1) #all deg across development stages (cut-off LFC = 0)
dev.off()
##########################################################
