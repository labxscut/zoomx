#!/usr/bin/env Rscript
library(plyr)
library(Matrix)
X_i=rep(1:61928,630)
X_j=sample(1:554715,length(X_i),replace=TRUE)
X_data=data.frame(i=X_i,j=X_j)
X_data=X_data[complete.cases(X_data),]
X_data=count(X_data, vars=c("i","j"))
X=sparseMatrix(i=X_data$i, j=X_data$j, x=X_data$freq)
Y_i=X_j
Y_j=X_i
Y_data=data.frame(i=Y_i,j=Y_j)
Y_data=Y_data[complete.cases(Y_data),]
Y_data=count(Y_data, vars=c("i","j"))
Y_data$freq=(Y_data$freq>=1)
Y=sparseMatrix(i=Y_data$i, j=Y_data$j, x=Y_data$freq)
M=X %*% Y
save.image(file="MatrixMultiplyOptimize.RData")
