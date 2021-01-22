library(OptionPricing)
# ------------------------- Initialize Parameters -------------------------------
r <- 0.0132
mu <- .15475
sigma <- .2214
kappa <- 170
T <- .25
delta <- 1/252
S0 <- 158.12
K <- 10
alpha <- .95

# ---------------------European Call Option Price using BS -----------------------
require(OptionPricing)
BS <- BS_EC(TT,kappa,r,sigma,S0)

bspx <- function(t,x){
  return(BS_EC(TT-t,kappa,r,sigma,x)[1])
}

# ---------------------European Call Option Delta using BS -----------------------
bsdelta <- function(t,x){
  return(BS_EC(TT-t,kappa,r,sigma,x)[2])
}

# -------------------------------- Portfolio Loss --------------------------------
# S_t+delta = S_t * exp(X), S_t+delta - S_t = S_t * ( exp(X)-1 )
Loss <- function(S0,TT,kappa,r,mu,sigma,t,delta){
  BS <- BS_EC(TT,kappa,r,sigma,S0)
  S_delta <- S0 * exp((mu - 0.5*sigma^2)*delta + sigma*sqrt(delta)*rnorm(1))
  callPrice <- BS[1]
  BS2 <- BS_EC(TT-delta,kappa,r,sigma,S_delta)
  Call_delta <- BS2[1]
  h <- BS[2]
  lss <- h*(S0 - S_delta) - (callPrice - Call_delta)
}

Lss <- c()
for(i in 1:10000){
  Lss <- c(Lss,100*Loss(158.12,0.25,170,0.0132,0.15475,0.2214,0,1/252))
}
# -------------------------------- VaR --------------------------------
VaR_10_sqrt <- sqrt(10) * quantile(Lss,0.95)  # 40.98


# --------------------------------  10-day VaR ------------------------

S_est <- function(S,mu,sigma,t,delta){
  S_delta <- S * exp((mu - 1/2*sigma^2)*delta + sigma * sqrt(delta)* rnorm(1))
}

S <- 158.12; k <- 170; r <- 0.0132; mu <- 0.15475; sigma <- 0.2214;delta <- 1/252; TT <- 0.25; t <- 0 ; K <- 10

result <- c()
for(ii in 1:10000){
  S_K <- c(S)
  for(i in 1:K){
    tmpt <- S_est(S_K[i], mu, sigma, t, delta)
    S_K <- c(S_K,tmpt)
  }
  
  Call_K <- c()
  for(i in 1:(K+1)){
    BS <- BS_EC((TT-(i-1)*delta),k,r,sigma,S_K[i])
    tmpt <- BS[1]
    Call_K <- c(Call_K, tmpt) 
  }
  
  h_K <- c()
  for(i in 1:(K+1)){
    BS <- BS_EC((TT-(i-1)*delta),k,r,sigma,S_K[i])
    tmpt <- BS[2]
    h_K <- c(h_K,tmpt)
  }
  
  Residual_K <- c(0)
  Value_1 <- h_K[1]*S - Call_K[1]
  Value_K <- c(Value_1)
  Loss_K <- c()
  for(i in 1:K){
    tmpt1 <- h_K[i]*S_K[i+1] - Call_K[i+1] + Residual_K[i]*exp(r*delta)# for value process
    
    Value_K <- c(Value_K,tmpt1)
    tmpt2 <- Value_K[i+1] - (h_K[i+1]*S_K[i+1] - Call_K[i+1]) # for residual 
    Residual_K <- c(Residual_K, tmpt2)
    tmpt3 <- Value_K[i]-Value_K[i+1]
    Loss_K <- c(Loss_K,tmpt3)
  }
  result <- c(result,sum(Loss_K))
}
100*quantile(result,0.95) # 39.00
