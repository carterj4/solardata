#
library(nlme)

solar <- read.csv("SolarSavings.csv",as.is=TRUE)
plot(solar$PowerBill[solar$Solar=='N'],xlim=c(0,52),ylim=c(0,240),type="l",xaxt="n", ylab="Power Bill ($)", xlab="Date", main= "Solar Savings")
axis(1,at=1:nrow(solar),labels=solar$Date)

lines((1:nrow(solar))[solar$Solar=='Y'],solar$PowerBill[solar$Solar=='Y'],col="blue")
legend("topright",col=c("black","blue"),legend=c("Regular", "Solar"),lty=1)

# create time indicator
solar$ztime <- 1:length(solar$Date) - 1
solar$rtime <- c(1:sum(solar$Solar == "N") - 1, 1:sum(solar$Solar == "Y") - 1)

# Create Month column
solar$month.n <- as.numeric(sapply(strsplit(solar$Date, "/"), "[",1)) - 1
solar$month.n[solar$month.n == 0] <- 12
solar$month <-  as.factor(sprintf("%02d",solar$month.n))

# seasonality by month
scatter.smooth(solar$month.n[solar$Solar == "N"],solar$PowerBill[solar$Solar == "N"],degree=2)

scatter.smooth(solar$month.n[solar$Solar == "Y"],solar$PowerBill[solar$Solar == "Y"],degree=2)

# make plot for seasonality by month
# par(mfrow=c(2,1))
pdf('smoothregular.pdf')
scatter.smooth(solar$month.n[solar$Solar == "N"],
               solar$PowerBill[solar$Solar == "N"],degree=2,
               ylab='Payment (Regular)',xlab='Month')#,xaxt='n')
dev.off()
pdf('smoothsolar.pdf')
scatter.smooth(solar$month.n[solar$Solar == "Y"],
               solar$PowerBill[solar$Solar == "Y"],degree=2,
               ylab='Payment (Solar)',xlab='Month')
dev.off()


# Different seasonality for Solar and Non Solar
sn <- c(0,0,0,0,0,1,1,1,1,0,0,0) 
sy <- c(1,1,0,0,0,0,0,0,0,0,1,1) 


solar$season <- sapply(1:length(solar$Solar), function(i) ifelse(solar$Solar[i] == "N",  sn[solar$month.n[i]], sy[solar$month.n[i]])) + 1


# Explore Correlation
acf(solar$PowerBill[solar$Solar == "N"])
acf(solar$PowerBill[solar$Solar == "Y"])



#########################
# Ztime models by month # 
#########################

# Fit AR(1) model
fitarZ.m <- gls(PowerBill ~ -1 + month:Solar, correlation = corAR1(form=~ztime),data=solar,method="ML")
summary(fitarZ.m)
plot(predict(fitarZ.m),type="l")

# Fit MA(1)
fitmaZ.m <- gls(PowerBill ~ -1 + month:Solar, correlation = corARMA(form=~ztime,p=0,q=1),data=solar,method="ML")
summary(fitmaZ.m)
plot(predict(fitmaZ.m),type="l")

# Fit ARMA(1,1)
fitarmaZ.m <- gls(PowerBill ~ -1 + month:Solar, correlation = corARMA(form=~ztime,p=1,q=1),data=solar,method="ML")
summary(fitarmaZ.m)
plot(predict(fitarmaZ.m),type="l")

######################
# Ztime Models By PS #
######################
# effect for 'peak' season

# fit AR(1)
fitarZ.ps <- gls(PowerBill ~ -1 + as.factor(season):Solar, correlation = corAR1(form=~ztime),data=solar,method="ML")
summary(fitarZ.ps)
plot(predict(fitarZ.ps),type="l")


# Fit MA(1)
fitmaZ.ps <- gls(PowerBill ~ -1 + as.factor(season):Solar, correlation = corARMA(form=~ztime|Solar,p=0,q=1),data=solar,method="ML")
summary(fitmaZ.ps)
plot(predict(fitmaZ.ps),type="l")

# Fit ARMA(1,1)
fitarmaZ.ps <- gls(PowerBill ~ -1 + as.factor(season):Solar, correlation = corARMA(form=~ztime|Solar,p=1,q=1),data=solar,method="ML")
summary(fitarmaZ.ps)
plot(predict(fitarmaZ.ps),type="l")


############################
# Time by Solar, Non Solar #
# 'rtime'                  #
############################
#effect for month
fitar.m <- gls(PowerBill ~ -1 + as.factor(month.n):Solar, correlation = corARMA(form=~rtime|Solar,p=1,q=0),data=solar,method="ML")
summary(fitar.m)
plot(predict(fitar.m),type="l")

fitma.m <- gls(PowerBill ~ -1 + as.factor(month.n):Solar, correlation = corARMA(form=~rtime|Solar,p=0,q=1),data=solar,method="ML")
summary(fitma.ps)
plot(predict(fitma.m),type="l")

fitarma.m <- gls(PowerBill ~ -1 + as.factor(month.n):Solar, correlation = corARMA(form=~rtime|Solar,p=1,q=1),data=solar,method="ML")
summary(fitarma.m)
plot(predict(fitma.m),type="l")


# effect for peak season
fitar.ps <- gls(PowerBill ~ -1 + as.factor(season):Solar, correlation = corARMA(form=~rtime|Solar,p=1,q=0),data=solar,method="ML")
summary(fitar.ps)
plot(predict(fitar.ps),type="l")

fitma.ps <- gls(PowerBill ~ -1 + as.factor(season):Solar, correlation = corARMA(form=~rtime|Solar,p=0,q=1),data=solar,method="ML")
summary(fitma.ps)
plot(predict(fitma.ps),type="l")

fitarma.ps <- gls(PowerBill ~ -1 + as.factor(season):Solar, correlation = corARMA(form=~rtime|Solar,p=1,q=1),data=solar,method="ML")
summary(fitarma.ps)
plot(predict(fitma.ps),type="l")

AIC(fitarma.ps,fitarma.m,fitar.ps,fitma.ps,fitar.m,fitma.m)

# Prediction into the future
sigma2 <- fitar.ps$sigma^2
phi <- coef(fitar.ps$modelStruct$corStruct,unconstrained=FALSE)

# AR(1)
intervals(fitar.ps)

n.sn <- sum(solar$Solar == "N")
n.sy<- sum(solar$Solar == "Y")


# calculating R squared for the fitar.ps
# ss.res <- sum((solar$PowerBill - predict(fitar.ps))^2)
# ss.total <- sum((solar$PowerBill - mean(solar$PowerBill))^2)
# Rsquared <- 1-(ss.res/ss.total)
# Rsquared

#make a plot of data vs fitar.ps model
plot(1:nrow(solar),solar$PowerBill,type='l',lwd=1,
     xaxt = 'n',xlab='',ylab='Power Bill',main='AR(1) Model Fit')
axis(1, at = seq(1, 51, by = 5),labels=solar$Date[seq(1, 51, by = 5)], las=2)
colors <- c('blue','red')
my.pch <- c(1,16)
abline(v=which(solar$Solar=='Y')[1],col='red',lty=2)
lines(predict(fitar.ps),type="l",col='blue',lwd = 2)
legend("topright",col=c("black","blue"),legend=c("Data", "Model"),lty=1,lwd=c(1,2))

#make a table of the model summary
library(xtable)
xtable(t(t(fitar.ps$coefficients)))

intervals(fitar.ps) #gets confidence intervals for model values

plot(fitar.ps$fitted,fitar.ps$residuals)

### De-correlate model to check assumptions
my.X<- model.matrix(fitar.ps,data=solar)
dim(my.X)

#de-correlate model to be able to check assumptions
# L.inv * Y ~ N(L.inv * X * beta, L.inv * L * t(L) * t(L.inv))
#make R
w <- "N"
n <- n.sn
k <- n.sy 
R <- diag(k+n)
R <- sigma2*phi^(abs(row(R) - col(R)))

L<- t(chol(R))
Ystar <- solve(L) %*% solar$PowerBill

XX.star <- solve(L) %*%  my.X
decorr.mod <- lm(Ystar ~ XX.star) #fit de-correlated model
#check assumptions
#residuals
pdf('residplot.pdf')
plot(decorr.mod$fitted.values,decorr.mod$residuals,
     main='Fitted Values vs. Residuals',xlab='Fitted Values',ylab='Residuals')
abline(h=0,col='red')
dev.off()
pdf('residhist.pdf')
hist(decorr.mod$residuals,prob=T,
     main='Histogram of Residuals',
     xlab='Residuals')
curve(dnorm(x,0,sd(decorr.mod$residuals)),col='red',add=T)
dev.off()



#linearity
library(car)
avPlots(decorr.mod)

##




# 1) predict how much has been saved already
w <- "N"
n <- n.sn
k <- n.sy 
R <- diag(k+n)
R <- sigma2*phi^(abs(row(R) - col(R)))
X.all <- model.matrix(~-1 + as.factor(sn[solar$month.n]+1))
Xstar <- X.all[n+1:k,]
X <- X.all[1:n,]
bhat <- coef(fitar.ps)
bw <- 1:2

pred.savings <- Xstar %*% bhat[bw] + R[n+1:k,1:n] %*% solve(R[1:n,1:n]) %*% (solar[solar$Solar == w,"PowerBill"][1:n] - X %*% bhat[bw])
pred.var <- R[n+1:k,n+1:k] - R[n+1:k,1:n] %*% solve(R[1:n,1:n]) %*% R[1:n,n+1:k]
  list(pred.mn,pred.var)
  
A <- matrix(1,ncol=k)
# Expected amount saved
m.sav <- A %*% (pred.savings - solar[solar$Solar == "Y","PowerBill"])
m.sav <- as.numeric(m.sav)

# Confidence interval on amount saved
sd.sav <- as.numeric(sqrt(A %*% pred.var %*% t(A)))
m.sav + c(-1,1) * qnorm(0.975) * sd.sav


prediction <- function(n,k,R,X.all,bw,w) {
  Xstar <- X.all[n+1:k,]
  X <- X.all[1:n,]
  pred.mn <- Xstar %*% bhat[bw] + R[n+1:k,1:n] %*% solve(R[1:n,1:n]) %*% (solar[solar$Solar == w,"PowerBill"][1:n] - X %*% bhat[bw])
  pred.var <- R[n+1:k,n+1:k] - R[n+1:k,1:n] %*% solve(R[1:n,1:n]) %*% R[1:n,n+1:k]
  list(pred.mn,pred.var)
}


# General prediction for k months into the future into the future
k <- 12

predict.sav <- function(k) {
	months <- c(solar$month.n,rep(c(2:12,1),times=ceiling(k/12))[1:k])
	
	R.n <- diag(k+n.sn+n.sy)
	R.n <- sigma2*phi^(abs(row(R.n) - col(R.n)))
	nn <- n.sn
	kn <- k+n.sy
	xn <- as.factor(sn[months]+1)
	Xn.all <- model.matrix(~-1 + xn)
	
	R.s <- diag(k+n.sy)
	R.s <- sigma2*phi^(abs(row(R.s) - col(R.s)))
	ns <- n.sy
	ks <- k
	xs <- as.factor(sy[months[(n.sn+1):length(months)]]+1)
	Xs.all <- model.matrix(~-1 + xs)
	
	if (k != 0) {
	  p.s <- prediction(ns,ks,R.s,Xs.all,3:4,"Y")
	  As <- matrix(1,ncol=k)
	  v.s <- As %*% p.s[[2]] %*% t(As) 
	} else {
	  p.s <- list(0,0)
	  v.s <- 0
	}
	p.n <- prediction(nn,kn,R.n,Xn.all,1:2,"N")
	
	m.sav <- sum(p.n[[1]]) - sum(p.s[[1]]) - sum(solar[solar$Solar == "Y","PowerBill"])
	
	An <- matrix(1,ncol=n.sy+k)
	sd.sav <- as.numeric(sqrt(v.s + An %*% p.n[[2]] %*% t(An)))
	c(est=m.sav, m.sav + c(lower=-1,upper=1) * qnorm(0.975) * sd.sav)
}

# show savings for years
years <- 5
m.out <- seq(0,years*12,by=12)

out <- sapply(m.out,predict.sav)
plot(0:years,out[1,],type="l",ylim=c(min(out),max(out)))
lines(0:years,out[2,],lty=2)
lines(0:years,out[3,],lty=2)

# show savings by months
m <- 12*8
out <- sapply(0:m,predict.sav)
plot(0:m,out[1,],type="l",ylim=c(min(out),max(out)))
lines(0:m,out[2,],lty=2)
lines(0:m,out[3,],lty=2)
abline(h=8000,col="red")
which(out[2,] > 8000)[1]
which(out[3,] > 8000)[1]

# predict number of months to recoop $8000 dollars
single.draw <- function(k,mu,R) {
  vv <- t(chol(R))
  vv %*% rnorm(k) + mu
}

predict.sav.sim <- function(m=5000,k) {
	months <- c(solar$month.n,rep(c(2:12,1),times=ceiling(k/12))[1:k])
	
	R.n <- diag(k+n.sn+n.sy)
	R.n <- sigma2*phi^(abs(row(R.n) - col(R.n)))
	nn <- n.sn
	kn <- k+n.sy
	xn <- as.factor(sn[months]+1)
	Xn.all <- model.matrix(~-1 + xn)
	
	R.s <- diag(k+n.sy)
	R.s <- sigma2*phi^(abs(row(R.s) - col(R.s)))
	ns <- n.sy
	ks <- k
	xs <- as.factor(sy[months[(n.sn+1):length(months)]]+1)
	Xs.all <- model.matrix(~-1 + xs)

	p.s <- prediction(ns,ks,R.s,Xs.all,3:4,"Y")
	p.n <- prediction(nn,kn,R.n,Xn.all,1:2,"N")
	
	m.sav <- sum(p.n[[1]]) - sum(p.s[[1]]) - sum(solar[solar$Solar == "Y","PowerBill"])
	
	single.draw <- function(k,mu,R) {
     vv <- t(chol(R))
     p <- vv %*% rnorm(k) + mu
     p[p<0] <- 0
     p
    }
    pred.s <- sapply(1:m, function(i) {
      pb.s <- c(solar[solar$Solar == "Y","PowerBill"], single.draw(ks,p.s[[1]],p.s[[2]]))
      pb.n <- single.draw(kn,p.n[[1]],p.n[[2]])
      cumsum(pb.n - pb.s)
    } )
    
   rbind(est=apply(pred.s,1,mean),apply(pred.s,1,quantile,p=c(0.025,0.975)))
}

predict.days.sim <- function(m=5000) {
	k <- 12*8
	months <- c(solar$month.n,rep(c(2:12,1),times=ceiling(k/12))[1:k])
	
	R.n <- diag(k+n.sn+n.sy)
	R.n <- sigma2*phi^(abs(row(R.n) - col(R.n)))
	nn <- n.sn
	kn <- k+n.sy
	xn <- as.factor(sn[months]+1)
	Xn.all <- model.matrix(~-1 + xn)
	
	R.s <- diag(k+n.sy)
	R.s <- sigma2*phi^(abs(row(R.s) - col(R.s)))
	ns <- n.sy
	ks <- k
	xs <- as.factor(sy[months[(n.sn+1):length(months)]]+1)
	Xs.all <- model.matrix(~-1 + xs)

	p.s <- prediction(ns,ks,R.s,Xs.all,3:4,"Y")
	p.n <- prediction(nn,kn,R.n,Xn.all,1:2,"N")
	
	m.sav <- sum(p.n[[1]]) - sum(p.s[[1]]) - sum(solar[solar$Solar == "Y","PowerBill"])
	
	single.draw <- function(k,mu,R) {
     vv <- t(chol(R))
     p <- vv %*% rnorm(k) + mu
     p[p<0] <- 0
     p
    }
    pred.s <- sapply(1:m, function(i) {
      pb.s <- c(solar[solar$Solar == "Y","PowerBill"], single.draw(ks,p.s[[1]],p.s[[2]]))
      pb.n <- single.draw(kn,p.n[[1]],p.n[[2]])
      which(cumsum(pb.n - pb.s) > 8000)[1] - n.sy
    } )
    
   c(est=mean(pred.s),quantile(pred.s,p=c(0.025,0.975)))
}

predict.days.sim()







# Simulate data to get confidence interval on Amount Saved (validate CI)

m <- 5000
preds <- sapply(1:m, function(i) single.draw(k,pred.savings,pred.var))
ci.pred.s <- t(apply(preds,1,quantile,p=c(0.025,0.975)))


# via simulation
preds.sav <- apply(preds,2,function(x) sum(x - solar[solar$Solar == "Y","PowerBill"]))

mean(preds.sav)
sd(preds.sav)
quantile(preds.sav,c(0.025,0.975))

# Confidence Interval On Prediction
ci.pred <- cbind(pred.savings,pred.savings) + qnorm(0.975)*sqrt(diag(pred.var)) %*% cbind(-1,1)

#
plot(solar$PowerBill[solar$Solar=='N'],xlim=c(0,52),ylim=c(0,240),type="l",xaxt="n", ylab="Power Bill ($)", xlab="Date", main= "Solar Savings")
axis(1,at=1:nrow(solar),labels=solar$Date)

lines((1:nrow(solar))[solar$Solar=='Y'],solar$PowerBill[solar$Solar=='Y'],col="blue")
legend("topright",col=c("black","blue"),legend=c("Regular", "Solar"),lty=1)

lines(n+1:k,pred.savings)
lines(n+1:k,ci.pred[,1],lty=2)
lines(n+1:k,ci.pred[,2],lty=2)



# Test with no seasonal effect
###
fitar <-  gls(PowerBill ~ -1 + Solar, correlation = corARMA(form=~rtime|Solar,p=1,q=0),data=solar,method="ML")
summary(fitar)

sigma2 <- fitar$sigma^2
phi <- 0.5185

n <- n.sn
k <- n.sy 
R <- diag(k+n)
R <- (sigma2/(1-phi^2))*phi^(abs(row(R) - col(R)))
#X.all <- model.matrix(~ 1, data=solar)
Xstar <- matrix(1,nrow=k)
X <- matrix(1,nrow=n)
bw <- c(1)
bhat <- coef(fitar)

pred.savings <- Xstar %*% bhat[bw] + R[n+1:k,1:n] %*% solve(R[1:n,1:n]) %*% (solar[solar$Solar == w,"PowerBill"][1:n] - X %*% bhat[bw])
pred.var <- R[n+1:k,n+1:k] - R[n+1:k,1:n] %*% solve(R[1:n,1:n]) %*% R[1:n,n+1:k]

m <- 5000
preds <- sapply(1:m, function(i) single.draw(k,pred.savings,pred.var))
ci.pred <- t(apply(preds,1,quantile,p=c(0.025,0.975)))

ci.pred <- cbind(pred.savings,pred.savings) + qnorm(0.975)*sqrt(diag(pred.var)) %*% cbind(-1,1)


plot(solar$PowerBill[solar$Solar=='N'],xlim=c(0,52),ylim=c(0,240),type="l",xaxt="n", ylab="Power Bill ($)", xlab="Date", main= "Solar Savings")
axis(1,at=1:nrow(solar),labels=solar$Date)

lines((1:nrow(solar))[solar$Solar=='Y'],solar$PowerBill[solar$Solar=='Y'],col="blue")
legend("topright",col=c("black","blue"),legend=c("Regular", "Solar"),lty=1)

lines(n+1:k,pred.savings)
lines(n+1:k,ci.pred[,1],lty=2)
lines(n+1:k,ci.pred[,2],lty=2)

# Months correlated
mn <- sapply(table(solar$month.n[solar$Solar == "N"]),function(x) 1:x)
my <- sapply(table(solar$month.n[solar$Solar == "Y"]),function(x) 1:x)
for (i in 1:12) {
  solar$mc[solar$month.n == i & solar$Solar=="N"] <- mn[[i]]
  solar$mc[solar$month.n == i & solar$Solar=="Y"] <- my[[i]]
}

fitar.m <- gls(PowerBill ~ -1 + as.factor(season):Solar, correlation = corAR1(form=~mc|Solar/month), data=solar)

# TO DO:
# Check assumptions
# - decorrelate and check

# Make predictions 
# - errors on
# - simulation study?

# Predictive accuracy




