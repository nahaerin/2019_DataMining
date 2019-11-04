####################################
# Statistical Data Mining Example2 #
####################################

# Q1 ------------------------------------------------

library(MASS)

# Dataset
Boston[1:5,]

dim(Boston)

?Boston

# Q2 ------------------------------------------------

n = nrow(Boston)
p = 13
X = Boston[,1:13]
y = Boston[,14]

best.m = vector('list', length=p)
aic = numeric(p)

for (k in 1:p)
{
	index = combn(1:p,k)
	m = ncol(index)
	
	RSS = numeric(m)
	for (i in 1:m)
	{
		subset = index[,i]
		dat = data.frame(X[,subset],y)
		model = lm(y~.,data=dat)
		RSS[i] = sum(model$residuals^2)
	}
	loc = which.min(RSS)
	best.m[[k]] = index[,loc]

	sigma2.hat = RSS[loc]/n
	aic[k] = n*(log(2*pi)+log(RSS[loc])-log(n)) + 2*(n+k+1)
}

loc1 = which.min(aic)
bestmodel = best.m[[loc1]]
colnames(X)[bestmodel]

# Q3 ------------------------------------------------

install.packages('leaps')
library(leaps)

?regsubsets

# Best subset selection ========================
fit.full = regsubsets(medv~., data=Boston, nvmax=13, nbest=3)
summary(fit.full)

fit.full = regsubsets(medv~., data=Boston, nvmax=13)
summary(fit.full)

bs = summary(fit.full)
names(bs)

bs$rss

par(mfrow=c(2,2))

plot(bs$rss,xlab="Number of Variables",ylab="RSS",type="l")

plot(bs$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="l")
l = which.max(bs$adjr2)
points(l,bs$adjr2[l], col="red",cex=2,pch=20)

plot(bs$cp,xlab="Number of Variables",ylab="Cp",type='l')
l = which.min(bs$cp)
points(l,bs$cp[l],col="red",cex=2,pch=20)

plot(bs$bic,xlab="Number of Variables",ylab="BIC",type='l')
l = which.min(bs$bic)
points(l,bs$bic[l],col="red",cex=2,pch=20)

plot(fit.full,scale="r2")
plot(fit.full,scale="adjr2")
plot(fit.full,scale="Cp")
plot(fit.full,scale="bic")

# The model with 11 predictors was selected. 
coef(fit.full,11)

# Forward Stepwise Selection =========================

fit.fwd = regsubsets(medv~.,data=Boston,nvmax=13,method="forward")
summary(fit.fwd)

fs = summary(fit.fwd)
plot(fs$cp,xlab="Number of Variables",ylab="Cp",type='l')
l = which.min(fs$cp)
points(l,fs$cp[l],col="red",cex=2,pch=20)

# Backward Stepwise Selection =========================

fit.bwd=regsubsets(medv~.,data=Boston,nvmax=13,method="backward")
summary(fit.bwd)

fs = summary(fit.bwd)
plot(fs$cp,xlab="Number of Variables",ylab="Cp",type='l')
l = which.min(fs$cp)
points(l,fs$cp[l],col="red",cex=2,pch=20)

# Q4 ------------------------------------------------


X = as.matrix(X)
y = as.vector(y)

l = sample(1:nrow(X),200)
tr.X = X[l,]
te.X = X[-l,]

tr.y = y[l]
te.y = y[-l]

# Standardized inputs
tr.Xs = scale(tr.X)

apply(tr.Xs,2,mean)
apply(tr.Xs,2,sd)

tr.Xbar = apply(tr.X,2,mean)
tr.sdx = apply(tr.X,2,sd)

te.Xs = NULL
for (j in 1:p)
{
	temp = (te.X[,j] - tr.Xbar[j])/tr.sdx[j]
	te.Xs = cbind(te.Xs,temp)
}

# Centered output

tr.ys = scale(tr.y,center=T,scale=F)

tr.ybar = mean(tr.y)


te.ys = te.y - tr.ybar



# Q5 ------------------------------------------------


b.ridge = function(lambda, Xs, ys) 
{
	p = ncol(Xs)	
	beta = solve(t(Xs) %*% Xs + lambda * diag(rep(1,p))) %*% t(Xs) %*% ys
	return(beta)
}


br = b.ridge(lambda=0.1, tr.Xs, tr.ys)
te.yhat = te.Xs %*% br
testMSE = sum((te.ys - te.yhat)^2)/length(te.ys)
testMSE

br = b.ridge(lambda=1, tr.Xs, tr.ys)
te.yhat = te.Xs %*% br
testMSE = sum((te.ys - te.yhat)^2)/length(te.ys)
testMSE


br = b.ridge(lambda=10, tr.Xs, tr.ys)
te.yhat = te.Xs %*% br
testMSE = sum((te.ys - te.yhat)^2)/length(te.ys)
testMSE


# Q6 ------------------------------------------------

install.packages('glmnet')
library(glmnet)

?glmnet
# if alpha=0 -> ridge; if alpha=1 -> lasso.
# Default: standardize = TRUE


# ridge regression
ridge = glmnet(tr.X,tr.y,alpha=0,lambda=c(0.1,1,10))
coef(ridge,s=c(0.1,1,10))

te.yhat = predict(ridge,newx = te.X, s=c(0.1,1,10))

sres = (te.yhat - te.y)^2
apply(sres,2,mean)


# lasso regression

lasso=glmnet(tr.X,tr.y,alpha=1,lambda=c(0.1,1,10))
coef(lasso,s=c(0.1,1,10))

te.yhat = predict(lasso,newx = te.X, s=c(0.1,1,10))

sres = (te.yhat - te.y)^2
apply(sres,2,mean)


# Q7 ------------------------------------------------

svd.Xs=svd(tr.Xs)

svd.Xs

svd.Xs$d

var.pc = svd.Xs$d^2 / nrow(tr.Xs)

explained = var.pc / sum(var.pc)
explained

V=svd.Xs$v

Z = tr.Xs %*% V
Z

Z1 = Z[,1]
Z2 = Z[,2]

theta1.pcr = (t(Z1) %*% tr.ys)/(t(Z1)%*%Z1)
theta2.pcr = (t(Z2) %*% tr.ys)/(t(Z2)%*%Z2)

b.pcr = theta1.pcr %*% V[,1] + theta2.pcr %*% V[,2]
b.pcr = as.vector(b.pcr)
b.pcr

te.yshat.pcr = te.Xs %*% b.pcr
testMSE = mean((te.ys - te.yshat.pcr)^2)
testMSE

# Q8 ------------------------------------------------

install.packages('pls')
library(pls)

train = data.frame(tr.X,tr.y)

# Principal component regression

pcr.fit = pcr(tr.y~.,data=train,scale=TRUE)
summary(pcr.fit)

mse = numeric(p)
for (npc in 1:p)
{
	te.yhat = predict(pcr.fit,te.X,ncomp=npc)
	mse[npc] = mean((te.y - te.yhat)^2)
}

plot(1:p,mse,type='b',col='red')

# Partial least squares

pls.fit=plsr(tr.y~., data=train,scale=TRUE)
summary(pls.fit)

mse = numeric(p)
for (npc in 1:p)
{
	te.yhat = predict(pls.fit,te.Xs,ncomp=npc,scale=F)
	mse[npc] = mean((te.y - te.yhat)^2)
}

plot(1:p,mse,type='b',col='red')







