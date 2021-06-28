setwd('/home/sidrah19220/nmd_review/reviewer2/GSVA/close')
tdf5<- read.csv('ttOSNs.csv')
t_exp<-setDT(tdf5, keep.rownames = "Cell_names")[]
exp<-t_exp[,-1]
colnames(exp)[1] <- "Cell_names"
write.table(exp,file="tOSN_input_for_doublet.csv",sep=",",quote=FALSE,row.names = FALSE,col.names = TRUE)

