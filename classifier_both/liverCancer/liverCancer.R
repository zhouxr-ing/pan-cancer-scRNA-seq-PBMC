library(caret)
library(glmnet)

load("../train/fit_lasso.Rdata")
load("../train/best_lambda_lasso.Rdata")
#scaleMinMax<-function(x){return((x-min(x))/(max(x)-min(x))*10)}

#######smart
##geneMat
load("../../classifier_geneExpr/liverCancer_test/testSet_liverCancer_smart.Rdata")
load("../../classifier_geneExpr/train/fit_lasso.Rdata")

names(testSet_liverCancer_smart)<-rownames(fit_lasso$beta)
for(i in 1:(dim(testSet_liverCancer_smart)[1])){
	testSet_liverCancer_smart[i,which(is.na(testSet_liverCancer_smart[i,]))]<-0
}
#testSet_CRC_smart<-log2(testSet_CRC_smart+0.001)
#testSet_CRC_smart<-t(apply(testSet_CRC_smart,1,scaleMinMax))

load("liver_testSet1.Rdata")
rownames(testSet1)<-sapply(rownames(testSet1),function(x){strsplit(x,"-")[[1]][1]})

if(all(rownames(testSet1)==rownames(testSet_liverCancer_smart))){
	testSet_liverCancer_smart<-cbind(testSet_liverCancer_smart,testSet1)
}

load("../train/fit_lasso.Rdata")
load("../train/best_lambda_lasso.Rdata")
testSet_liverCancer_smart<-testSet_liverCancer_smart[,rownames(fit_lasso$beta)]

pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx = as.matrix(testSet_liverCancer_smart),
			family=binomial(link = "logit"),type="response")
pred_lasso


#######10X
##geneMat
load("../../classifier_geneExpr/liverCancer_test/testSet_liverCancer_10X.Rdata")
load("../../classifier_geneExpr/train/fit_lasso.Rdata")
names(testSet_liverCancer_10X)<-rownames(fit_lasso$beta)
for(i in 1:(dim(testSet_liverCancer_10X)[1])){
	testSet_liverCancer_10X[i,which(is.na(testSet_liverCancer_10X[i,]))]<-0
}
#testSet_CRC_10X<-log2(testSet_CRC_10X+0.001)
#testSet_CRC_10X<-t(apply(testSet_CRC_10X,1,scaleMinMax))

load("liver_testSet2.Rdata")
if(all(rownames(testSet2)==rownames(testSet_liverCancer_10X))){
	testSet_liverCancer_10X<-cbind(testSet_liverCancer_10X,testSet2)
}

load("../train/fit_lasso.Rdata")
load("../train/best_lambda_lasso.Rdata")
testSet_liverCancer_10X<-testSet_liverCancer_10X[,rownames(fit_lasso$beta)]

pred_lasso <- predict(fit_lasso,s = best_lambda_lasso, 
			newx = as.matrix(testSet_liverCancer_10X),
			family=binomial(link = "logit"),type="response")
pred_lasso

