library(pheatmap)
library(RColorBrewer)

# select gene and its log2FC
pe_deg = read.csv("~/OneDrive/oxford/summer_internship/differential_gene_expression/dream.deg.pe.07122018.csv")
pe_heat = data.frame(log2FC = pe_deg[pe_deg$GeneSymbol %in% c("HNF4A", "PTF1A", "NEUROD1", "NEUROG3", "NKX2-2", "NKX6-1",
                                                              "KCNQ1", "KCNH2", "SLC2A2", "KCNJ8", "ABCC8", "ABCC9"), 4],
                     row.names = pe_deg[pe_deg$GeneSymbol %in% c("HNF4A", "PTF1A", "NEUROD1", "NEUROG3", "NKX2-2", "NKX6-1",
                                                                 "KCNQ1", "KCNH2", "SLC2A2", "KCNJ8", "ABCC8", "ABCC9"), 3])
pe_heat = as.matrix(pe_heat)

blc_deg = read.csv("~/OneDrive/oxford/summer_internship/differential_gene_expression/dream.deg.blc.07122018.csv")
blc_heat = data.frame(log2FC = blc_deg[blc_deg$GeneSymbol %in% c("HNF4A", "PTF1A", "NEUROD1", "NEUROG3", "NKX2-2", "NKX6-1",
                                                                 "KCNQ1", "KCNH2", "SLC2A2", "KCNJ8", "ABCC8", "ABCC9"), 4],
                      row.names = blc_deg[blc_deg$GeneSymbol %in% c("HNF4A", "PTF1A", "NEUROD1", "NEUROG3", "NKX2-2", "NKX6-1",
                                                                    "KCNQ1", "KCNH2", "SLC2A2", "KCNJ8", "ABCC8", "ABCC9"), 3])
blc_heat = as.matrix(blc_heat)

# combine
combined_heat = merge(pe_heat, blc_heat, by = 0, all = TRUE)
rownames(combined_heat) = combined_heat[,1]
colnames(combined_heat) = c("delete", "PE", "BLC")
combined_heat = combined_heat[,-1]
target = c("HNF4A", "PTF1A", "NEUROD1", "NEUROG3", "NKX2-2", "NKX6-1",
           "KCNQ1", "KCNH2", "SLC2A2", "KCNJ8", "ABCC8", "ABCC9")
combined_heat = combined_heat[match(target, rownames(combined_heat)),]
combined_mat = as.matrix(combined_heat)

# plot matrix 
breaksList = seq(-3, 3, by = 1)
pdf("l2fc.heatmap.12022019.pdf", width = 2.5)
pheatmap(combined_mat,
         cluster_rows = FALSE, # retain order of genes as given
         cluster_cols = FALSE, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,
         display_numbers = TRUE)
dev.off()