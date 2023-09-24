library(caret)
library(glmnet)

load("../train/fit_lasso.Rdata")
load("../train/best_lambda_lasso.Rdata")

#######smart
##geneMat
load("../../classifier_geneExpr/CRC_test/testSet_CRC_smart.Rdata")
load("../../classifier_geneExpr/train/fit_lasso.Rdata")
testSet_CRC_smart<-testSet_CRC_10X
names(testSet_CRC_smart)<-rownames(fit_lasso$beta)
for(i in 1:(dim(testSet_CRC_smart)[1])){
	testSet_CRC_smart[i,which(is.na(testSet_CRC_smart[i,]))]<-0
}
#testSet_CRC_smart<-log2(testSet_CRC_smart+0.001)
#testSet_CRC_smart<-t(apply(testSet_CRC_smart,1,scaleMinMax))

load("../../classifier_cellType/testSet1.Rdata")
if(all(rownames(testSet1)==rownames(testSet_CRC_smart))){
	testSet_CRC_smart<-cbind(testSet_CRC_smart,testSet1)
}

load("../train/fit_lasso.Rdata")
load("../train/best_lambda_lasso.Rdata")
testSet_CRC_smart<-testSet_CRC_smart[,rownames(fit_lasso$beta)]

pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx = as.matrix(testSet_CRC_smart),
			family=binomial(link = "logit"),type="response")
pred_lasso


#######10X
##geneMat
load("../../classifier_geneExpr/CRC_test/testSet_CRC_10X.Rdata")
load("../../classifier_geneExpr/train/fit_lasso.Rdata")
names(testSet_CRC_10X)<-rownames(fit_lasso$beta)
for(i in 1:(dim(testSet_CRC_10X)[1])){
	testSet_CRC_10X[i,which(is.na(testSet_CRC_10X[i,]))]<-0
}
#testSet_CRC_10X<-log2(testSet_CRC_10X+0.001)
#testSet_CRC_10X<-t(apply(testSet_CRC_10X,1,scaleMinMax))

load("../../classifier_cellType/testSet2.Rdata")
if(all(rownames(testSet1)==rownames(testSet_CRC_smart))){
	testSet_CRC_smart<-cbind(testSet_CRC_smart,testSet2)
}

load("../train/fit_lasso.Rdata")
load("../train/best_lambda_lasso.Rdata")

pred_lasso <- predict(fit_lasso,s = best_lambda_lasso, 
			newx = as.matrix(testSet_CRC_10X),
			family=binomial(link = "logit"),type="response")
pred_lasso

