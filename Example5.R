#####################################
# Statistical Data Mining Example 5 #
#####################################

# Q1 ==============================================

?Nile

Time=1871:1970
plot(Time,Nile)
x=seq(1870,1971,by=.1)

# Q2 ==============================================

# Use Epanechnikov kernel smoother with lambda=5 
# red=Epanechnikov

D = function(t) ifelse(abs(t)<=1,3/4*(1-t^2),0)

K=function(x0,x,h) D(abs(x-x0)/h)

lambda=5
L=length(x)
f.hat=numeric(L)
n=length(Nile)

K0i=numeric(n)
for (j in 1:L)
{
 for (i in 1:n) K0i[i]=K(x[j],Time[i],lambda)
 f.hat[j]=sum(K0i*Nile)/sum(K0i)
}
lines(x,f.hat,col="red")

# Use tri-cube kernel smoother with lambda=5 
# black=tri-cube

D=function(t) ifelse(abs(t)<=1,(1-abs(t)^3)^3,0)

for (j in 1:L)
{
  for (i in 1:n) K0i[i]=K(x[j],Time[i],lambda)
  f.hat[j]=sum(K0i*Nile)/sum(K0i)
}
lines(x,f.hat,col="black")


# Use Gaussian kernel smoother with lambda=5 
# green=gaussian
D=function(t) 1/sqrt(2*pi)*exp(-t^2/2)

for (j in 1:L)
{
 for (i in 1:n) K0i[i]=K(x[j],Time[i],lambda)
 f.hat[j]=sum(K0i*Nile)/sum(K0i)
}
lines(x,f.hat,col="green")


# Q3 ==============================================

D = function(t) 1/sqrt(2*pi)*exp(-t^2/2)

plot(Time,Nile)

col1 = c('blue','red','green','black')
lambda1 = c(1,3,5,10)

for (mm in 1:4)
{
	f.hat=numeric(L)
	K0i=numeric(n)
	for (j in 1:L)
	{
		 for (i in 1:n) K0i[i]=K(x[j],Time[i],lambda1[mm])
		 f.hat[j]=sum(K0i*Nile)/sum(K0i)
	}
	lines(x,f.hat,col=col1[mm])
}

# Q4 ==============================================

# Use Gaussian kernel smoother with lambda=5,

D = function(t) 1/sqrt(2*pi)*exp(-t^2/2)

plot(Time,Nile)
lambda=5

for (j in 1:L)
{
	for (i in 1:n) K0i[i]=K(x[j],Time[i],lambda)
	f.hat[j]=sum(K0i*Nile)/sum(K0i)
}
lines(x,f.hat,col="green")

# Use Gaussian local linear regression with lambda=5

B=cbind(1,Time)

for (j in 1:L)
{
  for (i in 1:n) K0i[i]=K(x[j],Time[i],lambda)
  W=diag(K0i)
  newx=c(1,x[j])
  f.hat[j]=t(newx)%*%solve(t(B)%*%W%*%B)%*%t(B)%*%W%*%Nile
}
lines(x,f.hat,,lty=2,col='blue')


# Built-in function =========================================

ks = ksmooth(Time, Nile, kernel='normal', bandwidth=5)

plot(Time, Nile)
lines(ks,col='blue')

install.packages('KernSmooth')
library(KernSmooth)

?locpoly

lp = locpoly(Time,Nile,bandwidth=5)
lines(lp,col='red')




