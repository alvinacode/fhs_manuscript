# hypergeometric test 
# phyper(a, b, c, d, lower.tail = FALSE)
# a = success in sample (no. of white balls drawn from an urn containing both black and white balls)
# b = success in population/ max success (no. of white balls in the urn)
# c = failure in population/ max failure (no. of black balls in the urn)
# d = sample size (no. of balls drawn from the urn)


# testing for enrichment of PE DEG (not necessarily HNF4A motif-bearing genes) in brown4

a = 99 # 99 genes overlap
b = 1351 # 1351 DEG
c = 16283 - 1351 # 16283 genes tested for DEG 
d = 321 # 321 genes in brown4 module
phyper(a, b, c, d, lower.tail = FALSE) # 7.302569e-33
