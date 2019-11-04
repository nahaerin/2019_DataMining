#####################################
# Statistical Data Mining Example 3 #
#####################################

# Part (1) ==============================================

hd.data = read.table("D:/Class/SKKU/2015/Datamining_Graduate/Example3/heartdisease.txt",sep="\t",header=TRUE)
training.data = hd.data[hd.data$train==TRUE,]
test.data = hd.data[hd.data$train==FALSE,]

tr.age = training.data$age
tr.fh = training.data$famhist=="Present"
tr.y = training.data$chd

te.age = test.data$age
te.fh = test.data$famhist=="Present"
te.y = test.data$chd

X = cbind(1, tr.fh, tr.age)
Y = cbind(1-tr.y, tr.y)

te.X = cbind(1,te.fh,te.age)

Ntr = nrow(training.data)
Nte = nrow(test.data)


# Part (2) ==============================================

B = solve(t(X)%*%X)%*%t(X)%*%Y
te.Y.hat = te.X%*%B
G.hat = te.Y.hat[,2]>.5
test.error = 1-sum(G.hat==te.y)/Nte
test.error            # Misclassification rate

# Part (3) ==============================================

# estimation of pi
pi.0 = sum(tr.y==0)/Ntr
pi.1 = sum(tr.y==1)/Ntr

pi.0
pi.1

# estimation of mu
mu.0 = c(mean(tr.fh[tr.y==0]),mean(tr.age[tr.y==0]))
mu.0

mu.1 = c(mean(tr.fh[tr.y==1]),mean(tr.age[tr.y==1]))
mu.1


Sigma = 0
for (i in 1:Ntr){
  xi = c(tr.fh[i],tr.age[i])
  if (tr.y[i]==0)
  Sigma = Sigma+(xi-mu.0)%*%t(xi-mu.0)
  else
  Sigma = Sigma+(xi-mu.1)%*%t(xi-mu.1)
 }
Sigma = Sigma/(Ntr-2)
Sigma

install.packages('mvtnorm')
library(mvtnorm)
?dmvnorm


LDA.y.hat = numeric(Nte)
for (i in 1:Nte){
  xi = c(te.fh[i],te.age[i])
  LDA.y.hat[i] = pi.1*dmvnorm(xi,mu.1,Sigma)/(pi.1*dmvnorm(xi,mu.1,Sigma)+pi.0*dmvnorm(xi,mu.0,Sigma))
 }
LDA.G.hat = LDA.y.hat>.5
LDA.test.error = 1-sum(te.y==LDA.G.hat)/Nte
LDA.test.error


# Built-in function

library(MASS)
?lda

tr.dat = data.frame(tr.y,tr.fh,tr.age)
colnames(tr.dat) = c('Y','FH','AGE')

fit = lda(Y~FH+AGE,data=tr.dat)


te.dat = data.frame(te.y,te.fh,te.age)
colnames(te.dat) = c('Y','FH','AGE')

te.yhat = predict(fit, newdata=te.dat)
te.Ghat = te.yhat$class
LDA.test.error = 1-sum(te.y==te.Ghat)/Nte
LDA.test.error

# QDA 

fit = qda(Y~FH+AGE,data=tr.dat)
te.yhat = predict(fit, newdata=te.dat)
te.Ghat = te.yhat$class
QDA.test.error = 1-sum(te.y==te.Ghat)/Nte
QDA.test.error


# Part (4) ==============================================

y = tr.y

# Iteratively reweighted least squares
beta.new = c(0,0,0)
beta.old = c(-1,-1,-1)
count = 1
while (max(abs(beta.new-beta.old)) > 1e-4)
{
	beta.old = beta.new
	p = exp(X%*%beta.old)/(1+exp(X%*%beta.old))
	W = diag(c(p*(1-p)))
	beta.new = beta.old+solve(t(X)%*%W%*%X)%*%t(X)%*%(y-p)
	cat("Iteration=",count," beta=",beta.new,"\n")
	count = count+1
}

LR.y.hat = exp(te.X%*%beta.new)/(1+exp(te.X%*%beta.new))

LR.G.hat=LR.y.hat > .5

LR.test.error = 1-sum(LR.G.hat==te.y)/Nte

LR.test.error
LDA.test.error

# Built-in function

fit = glm(Y~FH+AGE,data=tr.dat,family=binomial)
te.yhat = predict(fit,data=te.dat,type='response')
te.Ghat = te.yhat > .5
LR.test.error = 1-sum(LR.G.hat==te.y)/Nte
LR.test.error







