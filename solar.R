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
fitar.ps$modelStruct$corStruct
sigma2 <- fitar.ps$sigma^2
phi <- -0.01488195

# AR(1)
intervals(fitar.ps)

n.sn <- sum(solar$Solar == "N")
n.sy<- sum(solar$Solar == "Y")

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


# Simulate data to get confidence interval on Amount Saved (validate CI)
single.draw <- function(k,mu,R) {
  vv <- t(chol(R))
  vv %*% rnorm(k) + mu
}
m <- 5000
preds <- sapply(1:m, function(i) single.draw(k,pred.savings,pred.var))
ci.pred.s <- t(apply(preds,1,quantile,p=c(0.025,0.975)))

A <- matrix(1,ncol=k)
# Expected amount saved
m.sav <- A %*% (pred.savings - solar[solar$Solar == "Y","PowerBill"])
m.sav <- as.numeric(m.sav)

# Confidence interval on amount saved
sd.sav <- as.numeric(sqrt(A %*% pred.var %*% t(A)))
m.sav + c(-1,1) * qnorm(0.975) * sd.sav

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




