library(caret)
library(glmnet)

##geneMat
load("../../classifier_geneExpr/train/trainSet2.Rdata")
ov_1<-trainSet
ov_1<-ov_1
load("../../classifier_cellType/trainSet.Rdata")
ov_2<-trainSet
ov_2<-ov_2[rownames(ov_1),]
ov2<-cbind(ov_2[,1:5],ov_1)

load("../../classifier_geneExpr/train/coeff_lasso.Rdata")
coeff_lasso_geneExpr<-coeff_lasso[-1,]
load("../../classifier_cellType/coeff_lasso.Rdata")
coeff_lasso_cellType<-coeff_lasso[-1,]

ov2<-ov2[,c(names(coeff_lasso_geneExpr)[which(coeff_lasso_geneExpr!=0)],
		names(coeff_lasso_cellType)[which(coeff_lasso_cellType!=0)],
		"label")]

index <- sample(1:nrow(ov2), size = 0.5*nrow(ov2))
train <- as.data.frame(ov2[index,])
test <- as.data.frame(ov2[-index,])

#mean_train<-apply(train,2,mean)
#sd_train<-apply(train,2,sd)
#for(i in 1:18){
#	train[,i]<-scale(train[,i],scale=F)
#}
#for(i in 1:18){
#	test[,i]<-(test[,1]-mean_train[i])
#}


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




