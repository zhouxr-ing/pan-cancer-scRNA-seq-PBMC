library(caret)
library(glmnet)

load("../train/fit_lasso.Rdata")
load("../train/best_lambda_lasso.Rdata")
scaleMinMax<-function(x){return((x-min(x))/(max(x)-min(x))*10)}

#######smart
##geneMat
load("testSet_CRC_smart.Rdata")
names(testSet_CRC_10X)<-rownames(fit_lasso$beta)
for(i in 1:(dim(testSet_CRC_10X)[1])){
	testSet_CRC_10X[i,which(is.na(testSet_CRC_10X[i,]))]<-0
}
testSet_CRC_10X<-log2(testSet_CRC_10X+0.001)
testSet_CRC_10X<-t(apply(testSet_CRC_10X,1,scaleMinMax))


pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx = as.matrix(testSet_CRC_10X),
			family=binomial(link = "logit"),type="response")
pred_lasso


#######10X
##geneMat
load("testSet_CRC_10X.Rdata")
names(testSet_CRC_10X)<-rownames(fit_lasso$beta)
for(i in 1:(dim(testSet_CRC_10X)[1])){
	testSet_CRC_10X[i,which(is.na(testSet_CRC_10X[i,]))]<-0
}
testSet_CRC_10X<-log2(testSet_CRC_10X+0.001)
testSet_CRC_10X<-t(apply(testSet_CRC_10X,1,scaleMinMax))

# 模型评估
pred_lasso <- predict(fit_lasso,s = best_lambda_lasso, 
			newx = as.matrix(testSet_CRC_10X),
			family=binomial(link = "logit"),type="response")
pred_lasso

