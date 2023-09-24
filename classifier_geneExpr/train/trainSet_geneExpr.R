library(caret)
library(glmnet)

##geneMat
#load("trainSet.Rdata")
load("trainSet2.Rdata")
ov2<-trainSet

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

model_lasso_min<-glmnet(x = as.matrix(train[,-which(names(train)=="label")]), y = train[,"label"],alpha=1,lambda=best_lambda_lasso)
model_lasso_min$beta

# 根据最佳lambda构建模型
fit_lasso <- glmnet(x = as.matrix(train[,-which(names(train)=="label")]), 
			y = train[,"label"], family=binomial(link = "logit"),
			alpha = 1)
save(fit_lasso,file="fit_lasso.Rdata")
coeff_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			family=binomial(link = "logit"),type = 'coefficients')
coeff_lasso

# 模型评估
pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx = as.matrix(train[,-which(names(train)=="label")]),
			family=binomial(link = "logit"),type="response")
###roc
library(pROC)
roc<-roc(train[,"label"],pred_lasso)
pdf("roc_train.pdf",4,4)
plot(roc,col="blue",print.thres=T,print.auc=T,grid=c(0.2,0.2))
dev.off()
table(data.frame(train$label,pred_lasso>0.65))

write.csv(coeff_lasso[which(coeff_lasso[,1]!=0),],
	file="coeff.csv",quote=F)

##测试集

pred_lasso <- predict(fit_lasso, s = best_lambda_lasso, 
			newx = as.matrix(test[,-which(names(test)=="label")]),
			family=binomial(link = "logit"),type="response")
roc<-roc(test[,"label"],pred_lasso)
pdf("roc_test.pdf",4,4)
plot(roc,col="blue",print.thres=F,print.auc=T,grid=c(0.2,0.2))
dev.off()
table(data.frame(test$label,pred_lasso>0.62))

