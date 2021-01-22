# HW4 
# part 1-1

library(OptionPricing)
K <- 3
S <- 158.12
L <- 158.12
sigma <- 0.2214
r <- 0.0132 
TT <- 0.25

EC_t <- BS_EC(TT,(L+K),r,sigma,S)
EP_t <- BS_EP(TT,(L-K),r,sigma,S)

# compute value
V_t <- as.numeric((EC_t['delta']+EP_t['delta'])*S - (EC_t['price']+EP_t['price'])) 

Delta <- EC_t['delta']+EP_t['delta']
X_T <- c(0.6,-0.6,0.4,-0.4,0.2,-0.2)
sigma_T <- c(0.5,0.75,1.25,1.5,1.75,2)

W_X <- c(0.5,0.5,0.75,0.75,1,1)
W_sigma <- c(0.5,0.75,1,1,0.75,0.5)
V_T <- c()
L_T <- c()
W <- c()


for(i in 1:6){
  tempt_X <- X_T[i]
  S_T <- S*exp(tempt_X)
  W_X_tmpt <- W_X[i]
  for(j in 1:6){
    tempt_sigma <- sigma_T[j]*sigma
    EC_T <- BS_EC((TT-5/252),(L+K),r,tempt_sigma,S_T)
    EP_T <- BS_EP((TT-5/252),(L-K),r,tempt_sigma,S_T)
    V_T <- c(V_T,as.numeric(Delta)*S_T - (EC_T['price']+EP_T['price']))
    W <- c(W,W_X_tmpt*W_sigma[j])
  }
}

L_T <- -(V_T-V_t)
V_T
max(L_T)
sum(W)
length(W)
max(W*L_T) # 53
L_TT <- W*L_T 

# find the index for the minimum value 
loc <- which(L_TT == min(L_TT)) # the 31th value, (6,1) 
loc
L_T[loc]
X_T[6];sigma_T[1]





