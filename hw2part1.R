library(MASS)

df_file <- "/Users/ivy/Documents/学习/MF731/HW/hw2/CAT_TSLA_Data.csv"
df <- read.csv(df_file)
df$Date <- as.Date(df$Date) 
N <- length(df[,1])
New_df <- df[1:N-1, ]
New_df$CAT.Log.Return <- log(df$CAT[1:(N-1)]/df$CAT[2:N])
New_df$TSLA.Log.Return <- log(df$TSLA[1:(N-1)]/df$TSLA[2:N])
logReturn <- New_df[,4:5]
# ------------------------- EWMA -------------------------------
head(New_df)
# ------------------------- Initial Value ----------------------
mu0 <- apply(New_df[,4:5],2,mean)
var0 <- var(New_df[,4:5])
# ------------------------- estimate mean vector and covariance matrix
lambda <- 0.97
theta <- 0.97
var_est <- matrix(c(0,0,0,0),2,2)
mu_est <- matrix(c(0,0),1,2)
N <- dim(New_df)[1]
for(i in 1:N){
  var_est <- theta * var_est + (1- theta)*t(as.matrix(New_df[i,4:5]-mu_est))%*%as.matrix(New_df[i,4:5]-mu_est)
  mu_est <- lambda * mu_est + (1-lambda)*New_df[i,4:5]
}

# -------------------------The Market Cap of CAT, TSLA -------------------------
MktCap <- c(82.52,51.46)
weights <- MktCap/sum(MktCap)
portfolioValue <- 1000000 

################################ (a)  ######################################
# ------------------------- Compute Empirical portfolio loss -------------------
loss_Emp <- -portfolioValue * (exp(as.matrix(logReturn))-1) %*% as.matrix(weights)

################################ (b)  ######################################
# number of samples
K <- 25000
# Simulate Normal Distribution 
normal_EWMA <- mvrnorm(K,as.matrix(mu_est),var_est)

# ------------------------- Compute Simulated portfolio loss -------------------
loss_sim <- -portfolioValue * (exp(as.matrix(normal_EWMA))-1) %*% as.matrix(weights)

# ------------------------- Compute Var ----------------------------------------
alpha <- 0.95
sort_loss_Emp <- as.matrix(loss_Emp[order(loss_Emp),])
sort_loss_Sim <- as.matrix(loss_sim[order(loss_sim),])
VaR_Emp <- sort_loss_Emp[ceiling(N*alpha)]
VaR_Sim <- sort_loss_Sim[ceiling(K*alpha)]

# ------------------------- linearized loss operator --------------------------- 
VaR_Linear <- portfolioValue * (-as.matrix(mu_est) %*% as.matrix(weights) + sqrt(t(as.matrix(weights)) %*% as.matrix(var_est) %*% as.matrix(weights)) * qnorm(alpha))

