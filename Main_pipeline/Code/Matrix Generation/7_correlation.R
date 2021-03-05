library(ggpubr)


exp <- read.csv('wild_exp_zscore.csv')

e <- as.matrix(exp[,-1])
row.names(e) <- exp[,1]
e <- as.data.frame(e)

et <- t(e)

#############
x = cor(et, method = 'pearson')

write.csv(x, 'pearson-correlation-matrix_zscore.csv')



##################################correlation of spliced and unspliced#################



sp <- read.csv('spliced-mature.csv')

spl <- as.matrix(sp[,-1])
row.names(spl) <- sp[,1]
spl <- as.data.frame(spl)



usp <- read.csv('unspliced-immature.csv')
uspl <- as.matrix(usp[,-1])
row.names(uspl) <- usp[,1]
uspl <- as.data.frame(uspl)

spl = t(spl)
uspl = t(uspl)

k <- cor(spl, uspl, method = 'pearson')
write.csv(k, 'spl-unsp-correlation.csv')


l <- read.csv('spl-unsp-correlation.csv')

##################################EUCLIDEAN DISTANCE###################
library(philentropy)


u <- read.csv('umap-coord.csv')
uu <- as.matrix(u[,-1])
row.names(uu) <- u[,1]
uu <- as.data.frame(uu)

ud <- distance(uu, method = 'euclidean', use.row.names = TRUE)

udd <- as.data.frame(ud)
write.csv(udd, 'euclidean.csv', quote = FALSE)

