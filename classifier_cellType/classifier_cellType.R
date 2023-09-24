library(caret)
library(glmnet)

##geneMat
load("trainSet.Rdata")

#scaleMinMax<-function(x){return((x-min(x))/(max(x)-min(x))*10)}
#ov2<-data.frame(apply(trainSet[,1:5],2,scaleMinMax),label=trainSet[,6])
ov2<-data.frame(apply(trainSet[,1:5],2,function(x){scale(x,scale=F)}),label=trainSet[,6])

# 构建训练集和测试集
index <- sample(1:nrow(ov2), size = 0.5*nrow(ov2))
train <- as.data.frame(ov2[index,])
test <- as.data.frame(ov2[-index,])


######lasso
fit_lasso <-  glmnet(x = as.matrix(train[,-which(names(train)=="label")]), y = train[,"label"], alpha = 1)
pdf("fit_lasso.pdf",4.5,4.5)
plot(fit_lasso,xvar = 'lambda',label = T)
dev.off()

fit_lasso_cv <-  cv.glmnet(x = as.matrix(train[,-which(names(train)=="label")]),
				y = train[,"label"], alpha = 1,nfold=5,
				family = "binomial",type.measure="auc")
pdf("fit_lasso2.pdf",4.5,4.5)
plot(fit_lasso_cv,xvar = 'lambda',label = T)
dev.off()
best_lambda_lasso <- fit_lasso_cv$lambda.min
log(best_lambda_lasso)
save(best_lambda_lasso,file="best_lambda_lasso.Rdata")

#best_lambda_lasso=0.2
model_lasso_min<-glmnet(x = as.matrix(train[,-which(names(train)=="label")]), y = train[,"label"],alpha=1,lambda=best_lambda_lasso)
model_lasso_min$beta

# 根据最佳lambda构建模型
fit_lasso <- glmnet(x = as.matrix(train[,-which(names(train)=="label")]), 
			y = train[,"label"], family=binomial(link = "logit"),
			alpha = 1)
save(fit_lasso,file="fit_lasso.Rdata")
coeff_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			family=binomial(link = "logit"),type = 'coefficients')
save(coeff_lasso,file="coeff_lasso.Rdata")

# 模型评估
pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx = as.matrix(train[,-which(names(train)=="label")]),
			family=binomial(link = "logit"),type="response")
###roc
library(pROC)
roc<-roc(train[,"label"],pred_lasso)
pdf("train_cellType.pdf",4,4)
plot(roc,col="blue",print.thres=T,print.auc=T,grid=c(0.2,0.2))
dev.off()
table(data.frame(train$label,pred_lasso>0.5))

write.csv(coeff_lasso[which(coeff_lasso[,1]!=0),],
	file="coeff.csv",quote=F)

##测试集

pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx = as.matrix(test[,-which(names(test)=="label")]),
			family=binomial(link = "logit"),type="response")
roc<-roc(test[,"label"],pred_lasso)
pdf("test_cellType.pdf",4,4)
plot(roc,col="blue",print.auc=T,grid=c(0.2,0.2))
dev.off()
table(data.frame(test$label,pred_lasso>0.5))




##CRC
load("testSet1.Rdata")
scaleMinMax2<-function(x){
	minValue<-apply(trainSet[,-6],2,min)
	maxValue<-apply(trainSet[,-6],2,max)
	x2<-x
	for(i in 1:length(minValue)){
		x2[,i]<-(x2[,i]-minValue[i])/(maxValue[i]-minValue[i])*10
	}
	x2
}
testSet1<-data.frame(scaleMinMax2(testSet1))

testSet1<-testSet1[,rownames(fit_lasso$beta)]
pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx =as.matrix(testSet1),
			family=binomial(link = "logit"),type="response")
pred_lasso

load("testSet2.Rdata")
testSet2<-testSet2[,rownames(fit_lasso$beta)]
testSet2<-data.frame(scaleMinMax2(testSet2))
pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx =as.matrix(testSet2),
			family=binomial(link = "logit"),type="response")
pred_lasso


##CRC
load("testSet1.Rdata")
center2<-function(x){
	meanValue<-apply(trainSet[,-6],2,mean)
	sdValue<-apply(trainSet[,-6],2,sd)
	x2<-x
	for(i in 1:length(meanValue)){
		x2[,i]<-(x2[,i]-meanValue[i])/sdValue[i]
	}
	x2
}
testSet1<-data.frame(center2(testSet1))

testSet1<-testSet1[,rownames(fit_lasso$beta)]
pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx =as.matrix(testSet1),
			family=binomial(link = "logit"),type="response")
pred_lasso

load("testSet2.Rdata")
testSet2<-testSet2[,rownames(fit_lasso$beta)]
testSet2<-data.frame(center2(testSet2))
pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx =as.matrix(testSet2),
			family=binomial(link = "logit"),type="response")
pred_lasso

