library(caret)
library(glmnet)

load("../train/fit_lasso.Rdata")
load("../train/best_lambda_lasso.Rdata")
scaleMinMax<-function(x){return((x-min(x))/(max(x)-min(x))*10)}

#######smart
##geneMat
load("testSet_liverCancer_smart.Rdata")
names(testSet_liverCancer_smart)<-rownames(fit_lasso$beta)
for(i in 1:(dim(testSet_liverCancer_smart)[1])){
	testSet_liverCancer_smart[i,which(is.na(testSet_liverCancer_smart[i,]))]<-0
}

testSet_liverCancer_smart<-t(apply(testSet_liverCancer_smart,1,scaleMinMax))


pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx = as.matrix(testSet_liverCancer_smart),
			family=binomial(link = "logit"),type="response")
pred_lasso


#######10X
##geneMat
load("testSet_liverCancer_10X.Rdata")
names(testSet_liverCancer_10X)<-rownames(fit_lasso$beta)
for(i in 1:(dim(testSet_liverCancer_10X)[1])){
	testSet_liverCancer_10X[i,which(is.na(testSet_liverCancer_10X[i,]))]<-0
}

testSet_liverCancer_10X<-t(apply(testSet_liverCancer_10X,1,scaleMinMax))

pred_lasso <- predict(fit_lasso,s = best_lambda_lasso, 
			newx = as.matrix(testSet_liverCancer_10X),
			family=binomial(link = "logit"),type="response")
pred_lasso

