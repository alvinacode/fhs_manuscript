# see if downregulated/ upregulated genes are enriched with markers of specific cell types 

library(ggplot2)
library(ggsci)

# list of gene marker sources
pp = read.delim("gene_markers/TEAD_YAP.500genes.txt", header = FALSE) # pancreatic progenitors (Cebola, 2015)
pop.a = read.delim("gene_markers/Ramond2018.popA.txt", header = FALSE) # pancreatic progenitors (Ramond, 2018)
pop.b = read.delim("gene_markers/Ramond2018.popB.txt", header = FALSE) # endocrine differentiation pathway initiating NEUROG3 expression (Ramond, 2018)
pop.c = read.delim("gene_markers/Ramond2018.popC.txt", header = FALSE) # endocrine progenitors (Ramond, 2018)
pop.d = read.delim("gene_markers/Ramond2018.popD.txt", header = FALSE) # endocrine cells  (Ramond, 2018)
alpha = read.delim("gene_markers/Muraro2016.acell.txt", header = FALSE) # alpha cells (Muraro, 2016)
beta = read.delim("gene_markers/Muraro2016.bcell.txt", header = FALSE) # beta cells (Muraro, 2016)
delta = read.delim("gene_markers/Muraro2016.dcell.txt", header = FALSE) # delta cells (Muraro, 2016)
pancreatic.polypeptide = read.delim("gene_markers/Muraro2016.pancreaticpolypeptide.txt", header = FALSE) # pancreatic polypeptide cells (Muraro, 2016)
endothelial = read.delim("gene_markers/Muraro2016.endothelial.txt", header = FALSE) # endothelial (Muraro, 2016)

gene.marker = list(pp = pp,
                   pp2 = pop.a,
                   neurog3 = pop.b,
                   endo.prog = pop.c,
                   endo = pop.d,
                   alpha = alpha,
                   beta = beta,
                   delta = delta,
                   pancreatic.polypeptide = pancreatic.polypeptide,
                   endothelial = endothelial)


# differentially expressed genes (with hnf4a motif i.e. hnf4a 'targets')
de.deg = read.csv("de.hnf4a.target.all.csv", row.names = 1)
pe.deg = read.csv("pe.hnf4a.target.all.csv", row.names = 1)
blc.deg = read.csv("blc.hnf4a.target.all.csv", row.names = 1)

deg.down = list(de = subset(de.deg, de.deg$logFC < 0),
                pe = subset(pe.deg, de.deg$logFC < 0),
                blc = subset(blc.deg, blc.deg$logFC < 0))

deg.up = list(de = subset(de.deg, de.deg$logFC > 0),
              pe = subset(pe.deg, pe.deg$logFC > 0),
              blc = subset(blc.deg, blc.deg$logFC > 0))


# count no. of genes overlapping 
down.res = matrix(nrow = 3, ncol = 10)
for(i in 1:dim(down.res)[1]) { # go along row
  for(j in 1:dim(down.res)[2]) { # go along column
    down.res[i,j] = length(deg.down[[i]][deg.down[[i]]$GeneSymbol %in% gene.marker[[j]]$V1, 1])
  }
}

up.res = matrix(nrow = 3, ncol = 10)
for(i in 1:dim(up.res)[1]) { # row
  for(j in 1:dim(up.res)[2]) { # column
    up.res[i,j] = length(deg.up[[i]][deg.up[[i]]$GeneSymbol %in% gene.marker[[j]]$V1, 1])
  }
}


# hypergeometric test 
down.hypergeo = matrix(nrow = 3, ncol = 10)
for(i in 1:dim(down.res)[1]) { # row
  for(j in 1:dim(down.res)[2]) { # column
    if(i == 1) {
      a = down.res[1,j] # no. of overlapping genes
      b = 38 # no. differentially expressed genes 
      c = 15141 - 38 # no. of genes tested for deg
      d = length(gene.marker[[j]]$V1) # no. of gene markers
      down.hypergeo[1,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
    if(i == 2) {
      a = down.res[2,j]
      b = 700 
      c = 16283 - 700
      d = length(gene.marker[[j]]$V1)
      down.hypergeo[2,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
    if(i == 3) {
      a = down.res[3,j]
      b = 794
      c = 16961 - 794
      d = length(gene.marker[[j]]$V1)
      down.hypergeo[3,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
  }
}

up.hypergeo = matrix(nrow = 3, ncol = 10)
for(i in 1:dim(up.res)[1]) { # row
  for(j in 1:dim(up.res)[2]) { # column
    if(i == 1) {
      a = up.res[1,j] # no. of overlapping genes
      b = 38 # no. differentially expressed genes 
      c = 15141 - 38 # no. of genes tested for deg
      d = length(gene.marker[[j]]$V1) # no. of gene markers
      up.hypergeo[1,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
    if(i == 2) {
      a = up.res[2,j]
      b = 700 
      c = 16283 - 700
      d = length(gene.marker[[j]]$V1)
      up.hypergeo[2,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
    if(i == 3) {
      a = up.res[3,j]
      b = 794
      c = 16961 - 794
      d = length(gene.marker[[j]]$V1)
      up.hypergeo[3,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
  }
}

# visualise log2fc of gene markers 
pe.deg = transform(pe.deg, GeneSymbol = reorder(GeneSymbol, -logFC)) # reorder
plot.data = pe.deg[pe.deg$GeneSymbol %in% pop.c$V1, c(2,4)]

p = ggplot(plot.data, aes(x = GeneSymbol, y = logFC)) + 
  geom_bar(stat = "identity", aes(fill = logFC < 0)) +
  coord_flip() +
  scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values = c("#709AE1FF", "#FED439FF")) +
  theme(panel.background = element_blank(), panel.border=element_rect(fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"), plot.margin = unit(c(1,1,1,1), "line"))

pdf("neurog3+_markers_pe.pdf")
plot(p)
dev.off()

