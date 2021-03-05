library(dplyr)
library(Seurat)
library(patchwork)

##########################################################################
ex = read.csv('expression-for-correlation.csv', header = TRUE, row.names = 'X')
ext = t(ex)

exp <- CreateSeuratObject(counts = ext, project = "exp")
exp



exp <- NormalizeData(exp, normalization.method = "LogNormalize", scale.factor = 10000)
no = exp[["RNA"]]@data

non <- as.matrix(no[-1,])
non <- as.data.frame(non)
write.csv(non, 'normalised_expression.csv', quote = FALSE)

################################################################

