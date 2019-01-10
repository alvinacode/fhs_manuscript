#first try. A problem with the names of the chromosomes. The PCA plot doesn't quite make sense 

###################################################################
# All differentiation stages & multiple experiments (NEED WORK!!!)
###################################################################

#count matrix 
##marta
marta = read.delim("~/OneDrive/oxford/summer_internship/counts/marta.peak.counts.txt", 
                   header = TRUE, row.names = 1, check.names = FALSE, sep = "") # already ordered by stages and donor
marta = marta[,-c(1:5)] 
marta.donors = c("Ad2.1", "Ad3.1", "Neo1.1")
marta.stage = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "EN", "BLC") # 8 stages
colnames(marta) = paste(rep(marta.stage, each = 3), rep(marta.donors, 8), sep = "_")  # rename columns

##me
hnf4a = read.delim("~/OneDrive/oxford/summer_internship/counts/hnf4a.peak.counts.txt",
                   header = TRUE, row.names = 1, check.names = FALSE)
hnf4a = hnf4a[,-c(1:5)]
hnf4a.donors = c("SB", "EX1", "03A", "04A")  #data comes from four 'donors' (2 indiduals, 4 cell clones)
hnf4a.stage = c("DE", "PE", "BLC") #3 stages
order_by_donors.hnf4a = function(counts, donors. = hnf4a.donors) {
  nc = counts[, grepl("SB" , names(counts))]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in donors.[-1])  {
    i = counts[, grepl(s , names(counts))]
    nc = cbind(nc, i)
  }
  
  return(nc)
}
hnf4a = order_by_donors.hnf4a(hnf4a) #order
order_by_stages.hnf4a = function(counts, stage. = hnf4a.stage) {
  nc = counts[, grepl("DE" , names(counts))]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in stage.[-1])  {
    i = counts[, grepl(s , names(counts))]
    nc = cbind(nc, i)
  }
  
  return(nc)
}
hnf4a = order_by_stages.hnf4a(hnf4a) #order
##too many chromatin regions 
library(edgeR)
isexpr = rowSums(cpm(hnf4a) > 1) >= 2
hnf4a = hnf4a[isexpr,]

##merge common chromatin regions 
combined = merge(marta, hnf4a, by = 0, all = TRUE) 
rownames(combined) = combined$Row.names
combined = combined[,-1]
combined = na.omit(combined)

#design matrix
samples = c(rep(marta.donors, 8),
            rep(c("SB", "03A", "04A"), times = c(3,3,4)),
            rep(c("SB", "03A", "04A"), times = c(3,3,4)),
            rep(c("SB", "03A", "04A"), times = c(3,3,4)))
samples = as.factor(samples)

stages = c(rep(marta.stage, each = 3),
           rep(c("DE", "PE", "BLC"), each = 10))
stages = as.factor(stages)
stages = factor(stages, levels = c("iPSC", "DE", "PGT", "PFG", "EP", "PE", "EN", "BLC"))

study = c(rep("Perez-Alcantara", 24),
          rep("HNF4A", 30))

library(limma)
dge = DGEList(counts = combined)
normalised = calcNormFactors(dge)

design = model.matrix( ~ stages + study)

v = voom(normalised, design, plot = F) # voom normalize the read counts

batch_corrected = removeBatchEffect(v$E, study) # remove batch effects  

pca = prcomp(t(batch_corrected))

library(ggplot2)
library(ggfortify)
pca_plot <- autoplot(pca)
pca_plot #unlabelled to quickly check if all is ok

pca_out = as.data.frame(pca$x)
percentage = round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
percentage = paste(colnames(pca_out), "(", paste( as.character(percentage), "%", ")", sep="") )
run = colnames(normalised)

all1 <- cbind(pca_out, stages)
all2 <- cbind(all1, samples)
all3 <- cbind(all2, study)
all4 <- cbind(all3, run)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),axis.text.x=element_text(colour="black"),
               axis.text.y=element_text(colour="black"),
               axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

library(ggsci)
pdf("try1.pdf")
ggplot(all4, aes(PC1, PC2, color = stages, shape = study)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#FED439FF", "#709AE1FF", "#8A9197FF",
                                "#D2AF81FF", "#FD7446FF", "#D5E4A2FF", "#197EC0FF",
                                "#46732EFF", "#71D0F5FF")) + 
  theme +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  coord_fixed()
dev.off()
###################################################################

###################################################################
# All differentiation stages, HNF4A only (Fig 2C)
###################################################################
peaks = read.delim("~/OneDrive/oxford/summer_internship/atac/merged.peaks.txt") 
peaks = peaks[,-c(1:6)] #take out gene id 
peaks = order_by_stages(peaks)

# experiment information
cut.to.pieces = strsplit(names(peaks), split = "_")
get.first = lapply(cut.to.pieces, "[", 1)
experiments = unlist(get.first)
get.second = lapply(cut.to.pieces, "[", 2)
cell_lines = unlist(get.second)
get.third = lapply(cut.to.pieces, "[", 3)
stages = unlist(get.third)

# filter 
library(edgeR)
dge_list = DGEList(counts = as.matrix(peaks))
isexpr = rowSums(cpm(dge_list) > 1) >= 0.5 * ncol(peaks)
dge_list = dge_list[isexpr, ]

# normalize
norm = calcNormFactors(dge_list)
design = model.matrix( ~ stages + experiments)
v = voom(norm, design, plot = F) #voom normalize read counts 
batch_corrected = removeBatchEffect(v$E, experiments) #remove batch effects

# labelled PCA plot
pca = prcomp(t(batch_corrected))
pca_out = as.data.frame(pca$x)
percentage = round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
percentage = paste(colnames(pca_out), "(", paste( as.character(percentage), "%", ")", sep="") )

all = cbind(pca_out, stages)
all = cbind(all, experiments)
all = cbind(all, cell_lines)

theme = theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              strip.background=element_blank(),axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

library(ggsci)
pdf("atac_pca_cell_line_try2.pdf", width = 10)
ggplot(all, aes(PC1, PC2, color = stages, shape = cell_lines)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(17, 18, 1))+
  scale_color_simpsons() +
  theme +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  coord_fixed()
dev.off()


###################################################################

###################################################################
# For each differentiation stage (Figure 2D)
###################################################################
#tidy
peaks = read.delim("~/OneDrive/oxford/summer_internship/atac/merged.peaks.txt") 
peaks = peaks[,-c(1:6)] # take out gene id
stage = c("DE", "PE", "BLC")
order_by_stages = function(peaks, stage_ = stage) {
  nc = peaks[, grepl("DE", names(peaks))]
  
  for(s in stage_[-1]) {
    i = peaks[, grepl(s, names(peaks))]
    nc = cbind(nc, i)
  }
  
  return(nc)
  
}
peaks = order_by_stages(peaks)
# NB: 'missing' i.e. not in RNA-seq EXP23SBDE, EXP2603ADE, EXP23SBPE, EXP2603APE, EXP23SBBLC, EXP2603ABLC

# count matrices 
de_peaks = peaks[,1:10]
pe_peaks = peaks[, 11:20]
blc_peaks = peaks[, 21:30]
peak_list = list(de = de_peaks,
                 pe = pe_peaks,
                 blc = blc_peaks)

# metadata matrix 
cut.to.pieces = strsplit(names(peaks), split = "_")
get.first = lapply(cut.to.pieces, "[", 1)
experiments = unlist(get.first)
experiments = toupper(experiments)
get.second = lapply(cut.to.pieces, "[", 2)
cell_lines = unlist(get.second)
get.third = lapply(cut.to.pieces, "[", 3)
stages = unlist(get.third)

meta = list(de = data.frame(differentiation = experiments[1:10],
                            cell.line = cell_lines[1:10]),
            pe = data.frame(differentiation = experiments[11:20],
                            cell.line = cell_lines[11:20]),
            blc = data.frame(differentiation = experiments[21:30],
                             cell.line = cell_lines[21:30]))

# filter
library(edgeR)
filtered_list = list()
isexpr = list()
for(i in 1:3) {
  isexpr[[i]] = rowSums(cpm(peak_list[[i]]) > 1) >= 2
  filtered_list[[i]] = peak_list[[i]][isexpr[[i]],]
}

## check if there are any null values 
# for(i in 1:3) {
#  print(length(which((apply(filtered_list[[i]], 1, var) != 0) == "FALSE")))
#}


# normalisation
library(DESeq2)
vst_list = list()
for(i in 1:3) {
  vst_list[[i]] = vst(as.matrix(filtered_list[[i]]), blind = F) # variance stabilising transformation 
}

# remove batch effects 
library(sva)
batch_corrected = list()
for(i in 1:3) {
  batch = meta[[i]]$differentiation
  modcombat = model.matrix( ~1, data = meta[[i]])
  combat_mydata = ComBat(dat = vst_list[[i]],
                         batch = batch,
                         mod = modcombat,
                         par.prior = T,
                         prior.plots = F) 
  batch_corrected[[i]] = combat_mydata # store
}

# pca
percentage = list()
all = list()
for(i in 1:3) {
  pca = prcomp(t(batch_corrected[[i]]))
  pca_out = as.data.frame(pca$x)
  percentage[[i]] = round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2) 
  percentage[[i]] = paste(colnames(pca_out), "(", paste(as.character(percentage[[i]]), "%", ")", sep = ""))
  all[[i]] = cbind(meta[[i]], pca_out)
}

theme = theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              strip.background=element_blank(),axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

# plot
library(ggplot2)
library(ggsci)
pdf("atac.pca.all.dev.stages.pdf", width = 5, height = 5)
for(i in 1:3) {
  print(ggplot(all[[i]], aes(PC1, PC2, colour = differentiation, shape = cell.line)) +
          geom_point(size = 5) +
          scale_shape_manual(values = c(17, 18, 0)) +
          scale_color_simpsons() +
          theme +
          xlab(percentage[[i]][1]) +
          ylab(percentage[[i]][2]) + 
          coord_fixed())
}
dev.off()
###################################################################
