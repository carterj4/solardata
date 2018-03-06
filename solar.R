
#
library(nlme)
library(astsa)

solar <- read.csv("SolarSavings.csv",as.is=TRUE)
plot(solar$PowerBill[solar$Solar=='N'],xlim=c(0,52),ylim=c(0,240),type="l",xaxt="n", ylab="Power Bill ($)", xlab="Date", main= "Solar Savings")
axis(1,at=1:nrow(solar),labels=solar$Date)


lines((1:nrow(solar))[solar$Solar=='Y'],solar$PowerBill[solar$Solar=='Y'],col="blue")
legend("topright",col=c("black","blue"),legend=c("Regular", "Solar"),lty=1)

# create time indicator
solar$ztime <- 1:length(solar$Date) - 1

# add columns for season
s4 <- c("W","W","Sp","Sp","Sp","Su","Su","Su","F","F","F","W")
s4.adj <- c("W","W","W","Sp","Sp","Sp","Su","Su","Su","F","F","F")
s2 <- ifelse(s4 == "Sp" | s4 == "Su", "Su", "W" )
s2.adj <- c("W","W","Sp","Sp","Sp","Su","Su","Su","F","F","F","W")

# Different seasonality for Solar and Non Solar
sn <- c(0,0,0,0,0,1,1,1,1,0,0,0) 
sy <- c(1,1,0,0,0,0,0,0,0,0,1,1) 

# 
solar$month.n <- as.numeric(sapply(strsplit(solar$Date, "/"), "[",1)) - 1
solar$month.n[solar$month.n == 0] <- 12
solar$month <-  as.factor(sprintf("%02d",solar$month.n))

solar$season4 <- s4[solar$month.n]
solar$season4ad <- s4.adj[solar$month.n]
solar$season2 <- s2[solar$month.n]

solar$season <- sapply(1:length(solar$Solar), function(i) ifelse(solar$Solar[i] == "N",  sn[solar$month.n[i]], sy[solar$month.n[i]])) + 1


# Explore Correlation
acf(solar$PowerBill[solar$Solar == "N"])
acf(solar$PowerBill[solar$Solar == "Y"])

# seasonality by month
scatter.smooth(solar$month.n[solar$Solar == "N"],solar$PowerBill[solar$Solar == "N"],degree=2)

scatter.smooth(solar$month.n[solar$Solar == "Y"],solar$PowerBill[solar$Solar == "Y"],degree=2)

boxplot(PowerBill ~ season4, data=solar[solar$Solar == "N",])

sarmaN <- sarima(solar$PowerBill[solar$Solar == "N"],1,1,1,1,1,1,12)
sarmaY <- sarima(solar$PowerBill[solar$Solar == "Y"],1,1,1,1,1,1,12)


# Fit AR 1 model
# effect for month
fitar.m <- gls(PowerBill ~ -1 + month:Solar, correlation = corAR1(form=~1|ztime),data=solar,method="ML")
summary(fitar.m)
plot(predict(fitar.m),type="l")

# effect for 'peak' season
fitar.ps <- gls(PowerBill ~ -1 + as.factor(season):Solar, correlation = corAR1(form=~1|ztime),data=solar,method="ML")
summary(fitar.ps)
plot(predict(fitar.ps),type="l")

# s4ad
fitar.s4a <- gls(PowerBill ~ -1 + as.factor(season4ad):Solar, correlation = corAR1(form=~1|ztime),data=solar,method="ML")
plot(predict(fitar.s4a),type="l")

# Fit MA(1)
# effect for month
fitma.m <- gls(PowerBill ~ -1 + month:Solar + ztime, correlation = corARMA(form=~1|ztime,p=0,q=1),data=solar,method="ML")
summary(fitma.m)
plot(predict(fitma.m),type="l")

# effect for 'peak' season
fitma.ps <- gls(PowerBill ~ -1 + as.factor(season):Solar +ztime, correlation = corARMA(form=~1|ztime,p=0,q=1),data=solar,method="ML")
summary(fitma.ps)
plot(predict(fitma.ps),type="l")

fitma.ps <- gls(PowerBill ~ -1 + as.factor(season):Solar +ztime, correlation = corARMA(form=~1|ztime,p=0,q=1),data=solar,method="ML")
summary(fitma.ps)
plot(predict(fitma.ps),type="l")

# s4ad
fitma.s4a <- gls(PowerBill ~ -1 + as.factor(season4ad):Solar, correlation = corARMA(form=~1|ztime,p=0,q=1),data=solar,method="ML")
pred <- predict(fitma.s4a,interval="prediction")
plot(pred,type="l")

AIC(fitar.m)
AIC(fitar.s4a)
AIC(fitar.ps)

AIC(fitma.m)
AIC(fitma.s4a)
AIC(fitma.ps)








fit1 <- gls(PowerBill ~ -1 + month:Solar, correlation = corARMA(form=~1|ztime,p=1,q=1),data=solar)
summary(fit1)

plot(predict(fit1),type="l")
