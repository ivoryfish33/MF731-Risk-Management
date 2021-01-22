# # HW4 
# part 1-2
library(MASS)

setwd("/Users/ivy/Documents/11学习/MF731/HW/hw4/")
df <- read.csv("MSFT_AAPL_Log_Returns.csv",header = FALSE)
names(df)[1] <- paste("Date")
names(df)[2] <- paste("MSFT")
names(df)[3] <- paste("AAPL")

Mkt_Cps <- c(448.77,575.11)
W_Mkt_Cps <- as.matrix(Mkt_Cps/sum(Mkt_Cps))
W_Mkt_Cps
dim(W_Mkt_Cps)

# (1) Estimating mean vector and covariance matrix 
M <- 100
mu0 <- apply(df[1:M,2:3],2,mean)
cov0 <- var(df[1:M,2:3])
N <- dim(df)[1]
mu <- mu0
cov <- cov0 
theta <- 0.97
lambda <- 0.97
for(i in (M+1):N){
  cov <- theta *cov + (1-theta)*t(as.matrix(df[i,2:3]-mu))%*%as.matrix(df[i,2:3]-mu)
  mu <- lambda*mu + (1-lambda)*df[i,2:3]
}
mu <- as.matrix(mu)
cov <- as.matrix(cov)
mu;cov


# (2)
M <- 10000 # number of simulation
V <- 1000000
V_pos <- V*W_Mkt_Cps
VaR_loss <- -t(V_pos) %*% t(mu) + sqrt(t(V_pos)%*% cov %*%V_pos)*qnorm(0.95) # 14383.85
sqrt(10)*VaR_loss # 45485.72
3*sqrt(10)*VaR_loss # 136457.1


# (3)
# 
mu; cov 
# assume negative return for Apple
X_2 <- mu[2]-5*sqrt(cov[2,2])
X_2 # -.0.0605

rho <- cov[1,2]/(sqrt(cov[1,1]*cov[2,2]))
mu_1 <- mu[1]+rho*sqrt(cov[1,1])/sqrt(cov[2,2])*(X_2-mu[2])
var_1 <- cov[1,1]*(1-rho^2)
dim(mu_KT)
dim(mu)
dim(est_X)
K <- 10
M <- 50000 # the number of simulation
loss_KT <- c() 
V_pos <- V*W_Mkt_Cps

for( j in 1:M){
  cov_KT <- cov
  mu_KT <- mu
  est_X_1 <- rnorm(1,mu_1,var_1)
  est_X <- as.matrix(c(est_X_1,X_2))
  est_X_KT <- est_X
  for(i in 1:(K-1)){
    cov_KT <- theta*cov_KT + (1-theta)*(est_X-t(mu_KT))%*%(t(est_X)-mu_KT)
    mu_KT <- lambda*mu_KT + (1-lambda)*t(est_X)
    est_X <- as.matrix(mvrnorm(1,mu_KT,cov_KT))
    est_X_KT <- est_X_KT + est_X
  }
  loss_KT <- c(loss_KT,-t(V_pos)%*%est_X_KT)
}

mean(loss_KT)
loss_KT_sorted <- sort(loss_KT)
loss_KT_sorted[ceiling(50000*0.95)]
sqrt(10)*VaR_loss
sum(ifelse( loss_KT_sorted > rep(sqrt(10)*VaR_loss,500000),1,0)) # 43270
freq1 <- sum(ifelse( loss_KT_sorted > rep(sqrt(10)*VaR_loss,500000),1,0))/500000
freq1 # 0.44504
freq2 <- sum(ifelse( loss_KT_sorted > rep(3*sqrt(10)*VaR_loss,500000),1,0))/500000
freq2 # 0.00744