library(ggpubr)

#######READING EXPRESSION FILE###########

exp <- read.csv('wild_exp_zscore.csv')

e <- as.matrix(exp[,-1])
row.names(e) <- exp[,1]
e <- as.data.frame(e)

et <- t(e)

############# CALCULATING PEARSON CORRELATION###############
x = cor(et, method = 'pearson')

write.csv(x, 'pearson-correlation-matrix_zscore.csv')

#############READING UMAP COORDINATES FILE###############

library(philentropy)


u <- read.csv('umap-coord.csv')
uu <- as.matrix(u[,-1])
row.names(uu) <- u[,1]
uu <- as.data.frame(uu)

###################EUCLIDEAN DISTANCE###################

ud <- distance(uu, method = 'euclidean', use.row.names = TRUE)

udd <- as.data.frame(ud)
write.csv(udd, 'euclidean.csv', quote = FALSE)

