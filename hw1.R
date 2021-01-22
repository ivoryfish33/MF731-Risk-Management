library(tseries)
library(fGarch)
library(fBasics)
library(timeDate)
library(zoo)
library(roll)
library(caTools)


df <- read.csv("/Users/ivy/Documents/hw1/Nasdaq_Data.csv")

price <- df[,2]
df['LogReturn'] <- NaN
df[2:length(df$Date),3] <- log(price[2:length(df$Date)] / price[1:length(df$Date)-1]) 

df$Date <- as.Date(df$Date,format="%Y/%m/%d")
rets <- df$LogReturn[df$Date<"2014-01-01"]


#################### MA EWMA ##################################################
windowSize <- 100
lambda1 <- 0.94
lambda2 <- 0.97

logReturn <- as.matrix(df$LogReturn[2:length(df$LogReturn)])

# ma <- runsd(df$LogReturn,windowSize,endrule='NA')

df["MA"] <- NaN
df['EWMA(0.97)'] <- NaN
df['EWMA(0.94)'] <- NaN

df[2:length(df$Date),4] <- roll_var(logReturn,100)

#################### EWMA ####################################################
df$`EWMA(0.97)`[windowSize+1] <- df$MA[windowSize+1]
df$`EWMA(0.94)`[windowSize+1] <- df$MA[windowSize+1]

for (t in 1:length(df$Date)) {
  if(t > windowSize + 1){
    df$`EWMA(0.97)`[t] <- lambda1 * (df$`EWMA(0.97)`[t-1]) + (1-lambda1) * (df$LogReturn[t]^2)
    df$`EWMA(0.94)`[t] <- lambda2 * (df$`EWMA(0.94)`[t-1]) + (1-lambda2) * (df$LogReturn[t]^2)
  }
}

df[,4:6] <- df[,4:6] * 252

plot(df$Date[windowSize+1:length(df$Date)],sqrt(df$MA[windowSize+1:length(df$Date)]),'l',xlab='Date',ylab='Volatility',ylim=c(0,0.8),main='MA and EWMA')
lines(df$Date[windowSize+1:length(df$Date)],sqrt(df$`EWMA(0.97)`[windowSize+1:length(df$Date)]),'l',col='Blue')
lines(df$Date[windowSize+1:length(df$Date)],sqrt(df$`EWMA(0.94)`[windowSize+1:length(df$Date)]),'l',col='Red')
legend('topright',legend=c("MA", "EWMA(0.97)","EWMA(0.94)"), col=c("black", "blue","red"),lty=1,box.lty=0)

#################### Compute Garch Model #####################################
garch11 <- garch(rets[2:length(rets)])
plot(df$Date[2:length(rets)], garch11$fitted.values[,1]*sqrt(252),'l',xlab='Date',ylab='Volatility',main='Garch(1,1)')

a0 <- garch11$coef[1]
a1 <- garch11$coef[2]
b1 <- garch11$coef[3]

Vl <- a0/(1-a1-b1) * 252


################### Simulation ###############################################
nsim <- 50
nDay <- length(df$Date[df$Date>'2013-12-31'])
sim <- matrix(0,nDay+1,nsim)
x <- matrix(0,nrow = nDay+1,ncol=nsim)
forcastDate <- df$Date[df$Date>'2013-12-30']

# 2013-12-31 volatility
sim[1,] <- garch11$fitted.values[2264]

for(t in 1:nDay+1){
  x[t,] <- sim[t-1,] * rnorm(nsim,0,1)
  sim[t,] <- sqrt(a0 + a1* x[t,]^2 + b1 *sim[t-1,]^2)
}

for (i in 1:50){
  if(i==1){
    plot(forcastDate,sim[,i]*sqrt(252),type='l',col=1,ylim=c(0.05,0.6),main='Garch Simulation Plot',xlab='Date',ylab='Volatility')
  }
  else
    lines(forcastDate,sim[,i]*sqrt(252))
}



################ simple variance  ###############################
sd(df$LogReturn[df$Date>'2013-12-31'])^2 * 252

###########################################

