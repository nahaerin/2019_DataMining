#####################################
# Statistical Data Mining Example 4 #
#####################################

install.packages('np')
library(np)

data(cps71)

attach(cps71)

# (1) Polynomial Regression ====================================

fit=lm(logwage~poly(age,3),data=cps71)
coef(summary(fit))

# poly(age,3) generates orthogonal polynomials with degree 3.
# i.e., each column from poly(age,3) is a linear combination of (age, age^2, age^3).
# Also, columns are orthogonal.

agelims=range(age)
age.grid=seq(from=agelims[1],to=agelims[2])
preds=predict(fit,newdata=list(age=age.grid),se=TRUE)
se.bands=cbind(preds$fit+2*preds$se.fit,preds$fit-2*preds$se.fit)
plot(age,logwage,xlim=agelims,cex=.5,col="darkgrey")
lines(age.grid,preds$fit,lwd=2,col="blue")
matlines(age.grid,se.bands,lwd=1,col="blue",lty=3)

# (2) Regression Splines ======================================

library(splines)

kn = quantile(age,prob=c(0.25,0.5,0.75))

fit = lm(logwage~bs(age,knots=kn), data=cps71)

pred = predict(fit, newdata=list(age=age.grid), se=T)
plot(age, logwage, col="gray")
lines(age.grid, pred$fit, lwd=2)
lines(age.grid, pred$fit+2*pred$se, lty="dashed")
lines(age.grid, pred$fit-2*pred$se, lty="dashed")

dim(bs(age,knots=kn))	
# =>  4 regions * 4 parameters - 3 knots * 3 constraints 
#	= 7 (intercept + 6 basis functions)

dim(bs(age,df=6))
attr(bs(age,df=6),"knots")

# (3) Natrual Splines =======================================

fit2 = lm(logwage~ns(age,df=4),data=cps71)
# => 2 regions * 4 parameters + 2 regions * 2 parameters
#    - 1 knot * 3 constraints - 2 knots * 2 constraints 
#    = 5 (intercept + 4 basis functions)

attr(ns(age,df=4),"knots")

pred2 = predict(fit2,newdata=list(age=age.grid),se=T)
lines(age.grid, pred2$fit,col="red",lwd=2)


# (4) Smoothing Splines =====================================

fit = smooth.spline(age, logwage, df=16)
fit2 = smooth.spline(age, logwage, cv=TRUE)
fit2$df

plot(age,logwage,xlim=agelims,cex=.5,col="darkgrey")
lines(fit,col="red",lwd=2)
lines(fit2,col="blue",lwd=2)
legend("topright",legend=c("16 DF","8.51 DF"),col=c("red","blue"),lty=1,lwd=2,cex=.8)

# (5) Nonparametric logistic regression =====================

Iwage = as.numeric(logwage > 13.5)

fit = glm(Iwage ~ ns(age,df=4), family=binomial)

pred3 = predict(fit,newdata=list(age=age.grid),se=T,type='response')





