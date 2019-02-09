# see if module are enriched with downregulated/upregulated genes at a developmental stage

# hypergeometric test 
# phyper(a, b, c, d, lower.tail = FALSE)
# a = success in sample (no. of white balls drawn from an urn containing both black and white balls)
# b = success in population/ max success (no. of white balls in the urn)
# c = failure in population/ max failure (no. of black balls in the urn)
# d = sample size (no. of balls drawn from the urn)


# list of module genes
skyblue = read.delim("skyblue.txt", header = TRUE)
darkgrey = read.delim("darkgrey.txt", header = TRUE)
bisque4 = read.delim("bisque4.txt", header = TRUE)
cyan = read.delim("cyan.txt", header = TRUE)
brown4 = read.delim("brown4.txt", header = TRUE)
turquoise = read.delim("turquoise.txt", header = TRUE)

gene.marker = list(skyblue = skyblue, 
                   darkgrey = darkgrey, 
                   bisque4 = bisque4,
                   cyan = cyan,
                   brown4 = brown4,
                   turquoise = turquoise)

# differentially expressed genes (with hnf4a motif i.e. hnf4a 'targets')
de.deg = read.csv("~/OneDrive/oxford/summer_internship/differential_gene_expression/de.hnf4a.target.all.csv", row.names = 1)
pe.deg = read.csv("~/OneDrive/oxford/summer_internship/differential_gene_expression/pe.hnf4a.target.all.csv", row.names = 1)
blc.deg = read.csv("~/OneDrive/oxford/summer_internship/differential_gene_expression/blc.hnf4a.target.all.csv", row.names = 1)

deg.down = list(de = subset(de.deg, de.deg$logFC < 0),
                pe = subset(pe.deg, de.deg$logFC < 0),
                blc = subset(blc.deg, blc.deg$logFC < 0))

deg.up = list(de = subset(de.deg, de.deg$logFC > 0),
              pe = subset(pe.deg, pe.deg$logFC > 0),
              blc = subset(blc.deg, blc.deg$logFC > 0))

# count no. of genes overlapping 
down.res = matrix(nrow = 3, ncol = 6)
for(i in 1:dim(down.res)[1]) { # go along row
  for(j in 1:dim(down.res)[2]) { # go along column
    down.res[i,j] = length(deg.down[[i]][deg.down[[i]]$GeneSymbol %in% gene.marker[[j]]$geneSymbol, 1])
  }
}

up.res = matrix(nrow = 3, ncol = 6)
for(i in 1:dim(up.res)[1]) { # row
  for(j in 1:dim(up.res)[2]) { # column
    up.res[i,j] = length(deg.up[[i]][deg.up[[i]]$GeneSymbol %in% gene.marker[[j]]$geneSymbol, 1])
  }
}

# hypergeometric test 
down.hypergeo = matrix(nrow = 3, ncol = 6)
for(i in 1:dim(down.res)[1]) { # row
  for(j in 1:dim(down.res)[2]) { # column
    if(i == 1) {
      a = down.res[1,j] # no. of overlapping genes
      b = 38 # no. differentially expressed genes 
      c = 15141 - 38 # no. of genes tested for deg
      d = length(gene.marker[[j]]$geneSymbol) # no. of genes in the module 
      down.hypergeo[1,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
    if(i == 2) {
      a = down.res[2,j]
      b = 700 
      c = 16283 - 700
      d = length(gene.marker[[j]]$geneSymbol)
      down.hypergeo[2,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
    if(i == 3) {
      a = down.res[3,j]
      b = 794
      c = 16961 - 794
      d = length(gene.marker[[j]]$geneSymbol)
      down.hypergeo[3,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
  }
}

up.hypergeo = matrix(nrow = 3, ncol = 6)
for(i in 1:dim(up.res)[1]) { # row
  for(j in 1:dim(up.res)[2]) { # column
    if(i == 1) {
      a = up.res[1,j] # no. of overlapping genes
      b = 38 # no. differentially expressed genes 
      c = 15141 - 38 # no. of genes tested for deg
      d = length(gene.marker[[j]]$geneSymbol) # no. of genes in the module 
      up.hypergeo[1,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
    if(i == 2) {
      a = up.res[2,j]
      b = 700 
      c = 16283 - 700
      d = length(gene.marker[[j]]$geneSymbol)
      up.hypergeo[2,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
    if(i == 3) {
      a = up.res[3,j]
      b = 794
      c = 16961 - 794
      d = length(gene.marker[[j]]$geneSymbol)
      up.hypergeo[3,j] = phyper(a, b, c, d, lower.tail = FALSE)
    }
  }
}
