rm(list=ls())
cat("\014")

library(readxl)
library(readr)

#setwd('/Users/chenx/Documents/Exp/A-TEAMS/lxy')

setwd('C:/AÅÌ-ÏîÄ¿/project-4-ESPCA/Issues_20180510')

source('ESPCA/fun_ESPCA.R')

gse <- data.frame(read_excel("db/breast_cancer/GSE3494.xlsx"))
gse_data <- as.matrix(gse[,-1]) # remove label

ppi <- data.frame(read_delim("db/breast_cancer/PPI.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
edges<-ppi
ppi$weight <- 1

genes <- unique(c(as.vector(t(ppi$X1)),as.vector(t(ppi$X2))))

ppi.adj<- edgelist2adjmatrix(genes = genes,edgelist = ppi,directed = FALSE)
m <- as.matrix(ppi.adj)

edges = list()
edges.number <- dim(ppi)[1]
gse_data_columns <- colnames(gse_data)
for (i in 1:edges.number){
  nodea<-ppi$X1[i];nodea.idx<-which(gse_data_columns == nodea)
  nodeb<-ppi$X2[i];nodeb.idx<-which(gse_data_columns == nodeb)
  if(length(nodea.idx) >0 && nodea.idx >0 && length(nodeb.idx) >0 && nodeb.idx > 0){
    edges[[length(edges)+1]] <- c(nodea.idx,nodeb.idx)
  }
}

colnames(gse_data) <- NULL
out =  ESPCA(gse_data, k = 2, edges, k.group=10)

sessionInfo()