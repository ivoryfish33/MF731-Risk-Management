
setwd("/Users/ivy/Documents/学习/MF731/HW/hw3")
# Part1
# Exercise2
library(janitor)
# Price Data
Price <- read.csv("Prices.csv") # Have changed the date format 
Price$Date <- as.Date(Price$Date, format="%m/%d/%y")
N <- dim(Price)[1]
tail(Price)
# Log Return Data
# LRet <- diff(log(as.numeric(Price[,2:5])))
# head(log(Price[,2:5]))

LRet <- Price[2:N,]
LRet$Boeing <- log(Price$Boeing[2:N]/Price$Boeing[1:(N-1)])
LRet$McDonalds <- log(Price$McDonalds[2:N]/Price$McDonalds[1:(N-1)])
LRet$Nike <- log(Price$Nike[2:N]/Price$Nike[1:(N-1)])
LRet$Walmart <- log(Price$Walmart[2:N]/Price$Walmart[1:(N-1)])
head(LRet)
ret <- as.matrix(LRet[,2:5])
head(ret)
# EWMA
# Initialization
M <- 30 

ewma_mu <- as.matrix(apply(ret[1:30,],2,mean))
ewma_cov <- as.matrix(var(ret[1:30,]))

N <- dim(LRet)[1]
lambda <- 0.94
theta <- 0.97 
for(i in (M+1):(N)){
  ewma_cov <- theta * ewma_cov + (1 - theta) * (as.matrix(ret[i,]) - ewma_mu) %*% t(as.matrix(ret[i,]) - ewma_mu)
  ewma_mu <- lambda * ewma_mu + (1 - lambda) * as.matrix(ret[i,])
}
# estimated mean and variance for t + delta
ewma_mu
ewma_cov

require(MASS)

ewma_x1 <- matrix(0,50000,4)
ewma_x2 <- matrix(0,50000,4)
ewma_x3 <- matrix(0,50000,4)
ewma_x4 <- matrix(0,50000,4)
ewma_x5 <- matrix(0,50000,4)
ewma_x6 <- matrix(0,50000,4)
ewma_x7 <- matrix(0,50000,4)
ewma_x8 <- matrix(0,50000,4)
ewma_x9 <- matrix(0,50000,4)
ewma_x10 <- matrix(0,50000,4)
ewma_x1 <- mvrnorm(50000, ewma_mu, ewma_cov)
dim(ewma_x1)
for(i in 1:50000){
	cov <- ewma_cov
	mu <- ewma_mu
	cov <- theta * cov + (1 - theta) * (as.matrix(ewma_x1[i,]) - mu) %*% t(as.matrix(ewma_x1[i,]) - mu)
	mu <- lambda * mu + (1 - lambda) * as.matrix(ewma_x1[i,])
	ewma_x2[i,] <- mvrnorm(1, mu, cov)
	cov <- theta * cov + (1 - theta) * (as.matrix(ewma_x2[i,]) - mu) %*% t(as.matrix(ewma_x2[i,]) - mu)
	mu <- lambda * mu + (1 - lambda) * as.matrix(ewma_x2[i,])
	ewma_x3[i,] <- mvrnorm(1, mu, cov)
	cov <- theta * cov + (1 - theta) * (as.matrix(ewma_x3[i,]) - mu) %*% t(as.matrix(ewma_x3[i,]) - mu)
	mu <- lambda * mu + (1 - lambda) * as.matrix(ewma_x3[i,])
	ewma_x4[i,] <- mvrnorm(1, mu, cov)
	cov <- theta * cov + (1 - theta) * (as.matrix(ewma_x4[i,]) - mu) %*% t(as.matrix(ewma_x4[i,]) - mu)
	mu <- lambda * mu + (1 - lambda) * as.matrix(ewma_x4[i,])
	ewma_x5[i,] <- mvrnorm(1, mu, cov)
	cov <- theta * cov + (1 - theta) * (as.matrix(ewma_x5[i,]) - mu) %*% t(as.matrix(ewma_x5[i,]) - mu)
	mu <- lambda * mu + (1 - lambda) * as.matrix(ewma_x5[i,])
	ewma_x6[i,] <- mvrnorm(1, mu, cov)
	cov <- theta * cov + (1 - theta) * (as.matrix(ewma_x6[i,]) - mu) %*% t(as.matrix(ewma_x6[i,]) - mu)
	mu <- lambda * mu + (1 - lambda) * as.matrix(ewma_x6[i,])
	ewma_x7[i,] <- mvrnorm(1, mu, cov)
	cov <- theta * cov + (1 - theta) * (as.matrix(ewma_x7[i,]) - mu) %*% t(as.matrix(ewma_x7[i,]) - mu)
	mu <- lambda * mu + (1 - lambda) * as.matrix(ewma_x7[i,])
	ewma_x8[i,] <- mvrnorm(1, mu, cov)
	cov <- theta * cov + (1 - theta) * (as.matrix(ewma_x8[i,]) - mu) %*% t(as.matrix(ewma_x8[i,]) - mu)
	mu <- lambda * mu + (1 - lambda) * as.matrix(ewma_x8[i,])
	ewma_x9[i,] <- mvrnorm(1, mu, cov)
	cov <- theta * cov + (1 - theta) * (as.matrix(ewma_x9[i,]) - mu) %*% t(as.matrix(ewma_x9[i,]) - mu)
	mu <- lambda * mu + (1 - lambda) * as.matrix(ewma_x9[i,])
	ewma_x10[i,] <- mvrnorm(1, mu, cov)
}

weights <- as.matrix(c(196.94,125.86,131.57,282.87)/sum(c(196.94,125.86,131.57,282.87)))

dim(weights)
dim(ewma_x1)

# K = 1
Loss_1 <- (-1000000)* (exp(ewma_x1)%*%weights -1)
L_1 <- sort(Loss_1)
VaR_1 <- L_1[50000*0.95] # 14552.71
sqrt(10)*VaR_1
ES_1 <- sum(L_1[(50000*0.95):50000])/(50000*(1-0.95)) # 18580.4
sqrt(10)*ES_1
Spectral_1 <- 0
n <- 50000
for (i in 1:n){
	Spectral_1 <- Spectral_1 + L_1[i]*(exp(i/n*30)-exp((i-1)/n*30))/(exp(30)-1)
}
Spectral_1 # 18283.53
sqrt(10)*Spectral_1

# K = 10 
Loss <-(-1000000)* ((exp(ewma_x1)%*%weights )*(exp(ewma_x2)%*%weights)*(exp(ewma_x3)%*%weights)*(exp(ewma_x4)%*%weights)*(exp(ewma_x5)%*%weights)*(exp(ewma_x6)%*%weights)*(exp(ewma_x7)%*%weights)*(exp(ewma_x8)%*%weights)*(exp(ewma_x9)%*%weights)*(exp(ewma_x10)%*%weights)-1)

quantile(Loss,0.95)
L <- sort(Loss)
VaR <- L[50000*0.95] # 50423.23
ES <- sum(L[(50000*0.95):50000])/(50000*(1-0.95)) # 66710.53

Spectral_10 <- 0
n <- 50000
for (i in 1:n){
	Spectral_10 <- Spectral_10 + L[i]*(exp(i/n*30)-exp((i-1)/n*30))/(exp(30)-1)
}
Spectral_10 # 65608.52


# Part2
FS_Price <- read.csv("/Users/ivy/Documents/学习/MF731/HW/hw3/Five_Stock_Prices.csv")
head(FS_Price)
N <- dim(FS_Price)[1]
M <- 50
Pos <- as.matrix(rep(3000000,5))
dim(Pos)
FS_Ret <- FS_Price[2:N,2:6]
FS_Ret <- log(FS_Price[2:N,2:6]/FS_Price[1:(N-1),2:6])
head(FS_Ret)
N <- dim(FS_Ret)[1]
# Parameters 
alpha <- 0.99
lambda <- 0.94
theta <- 0.96
# Mean and Variance Estimation
ewma_mu <- as.matrix(apply(FS_Ret[1:M,],2,mean))
ewma_cov <- as.matrix(var(FS_Ret[1:M,]))
CVar <- matrix(0,(N-M),5) # tp store percent contribution to the loss variance
RCVaR <- matrix(0,(N-M),5) # to store percent component value at risk
RCES <- matrix(0,(N-M),5)

dim(ewma_mu)
# 
VaR_Z <- qnorm(alpha)
ES_Z <- 1/(1-alpha)*(1/sqrt(2*pi)*exp(-1/2*qnorm(alpha)^2))

for(i in (M+1):(N)){
	VaR <- -t(Pos)%*%ewma_mu + sqrt(t(Pos)%*%ewma_cov%*%Pos)*VaR_Z
	tempt1 <- - t(ewma_mu) + t(ewma_cov%*%Pos)*VaR_Z/as.numeric(sqrt(t(Pos)%*%ewma_cov%*%Pos))
	RCVaR[(i-M),] <- 100*tempt1*(3*10^6)/as.numeric(VaR)
	ES <- -t(Pos)%*%ewma_mu + sqrt(t(Pos)%*%ewma_cov%*%Pos)*ES_Z
	tempt2 <- - t(ewma_mu) + t(ewma_cov%*%Pos)*ES_Z/as.numeric(sqrt(t(Pos)%*%ewma_cov%*%Pos))
	RCES[(i-M),] <- 100*tempt2*(3*10^6)/as.numeric(ES)
	CVar[(i-M),] <- t(ewma_cov%*%(Pos^2))/as.numeric(t(Pos)%*%ewma_cov%*%Pos)*100
	ewma_cov <- theta * ewma_cov + (1 - theta) * (t(as.matrix(FS_Ret[i,])) - ewma_mu) %*% t(t(as.matrix(FS_Ret[i,])) - ewma_mu)
	ewma_mu <- lambda * ewma_mu + (1 - lambda) * t(as.matrix(FS_Ret[i,]))
}
dim(RCES)
# Expected shortfall 
ymax <- max(RCES)
ymin <- min(RCES)
plot(c(1:703),RCES[,1],ylim =c(ymin,ymax),type='l',lty=1,xlab='703 data points',ylab = 'Component ES %', main='Percentage component ES')
lines(c(1:703),RCES[,2],ylim =c(ymin,ymax),type='l',lty=2,col='red')
lines(c(1:703),RCES[,3],ylim =c(ymin,ymax),type='l',lty=3,col='blue')
lines(c(1:703),RCES[,4],ylim =c(ymin,ymax),type='l',lty=4,col='green')
lines(c(1:703),RCES[,5],ylim =c(ymin,ymax),type='l',lty=5,col='brown')
legend("topright",legend=c("Walmart","Target","Costco","Citigroup","JP.Morgan"),col=c('black','red','blue','green','brown'),lty=c(1,2,3,4,5))
# VaR
ymax <- max(RCVaR)
ymin <- min(RCVaR)
plot(c(1:703),RCVaR[,1],ylim =c(ymin,ymax),type='l',lty=1,xlab='703 data points',ylab = 'Component VaR %',main='Percentage component VaR')
lines(c(1:703),RCVaR[,2],ylim =c(ymin,ymax),type='l',lty=2,col='red')
lines(c(1:703),RCVaR[,3],ylim =c(ymin,ymax),type='l',lty=3,col='blue')
lines(c(1:703),RCVaR[,4],ylim =c(ymin,ymax),type='l',lty=4,col='green')
lines(c(1:703),RCVaR[,5],ylim =c(ymin,ymax),type='l',lty=5,col='brown')
legend("topright",legend=c("Walmart","Target","Costco","Citigroup","JP.Morgan"),col=c('black','red','blue','green','brown'),lty=c(1,2,3,4,5))
# Variance 
ymax <- max(CVar)
ymin <- min(CVar)
plot(c(1:703),CVar[,1],ylim =c(ymin,ymax),type='l',lty=1,xlab='703 data points',ylab = 'Component Variance %',main='Percentage variance contribution')
lines(c(1:703),CVar[,2],ylim =c(ymin,ymax),type='l',lty=2,col='red')
lines(c(1:703),CVar[,3],ylim =c(ymin,ymax),type='l',lty=3,col='blue')
lines(c(1:703),CVar[,4],ylim =c(ymin,ymax),type='l',lty=4,col='green')
lines(c(1:703),CVar[,5],ylim =c(ymin,ymax),type='l',lty=5,col='brown')
legend("topright",legend=c("Walmart","Target","Costco","Citigroup","JP.Morgan"),col=c('black','red','blue','green','brown'),lty=c(1,2,3,4,5))

