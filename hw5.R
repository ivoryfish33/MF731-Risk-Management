# 731HW5 Huishu Xue
# Part 1 
# initiate parameters
lambda <- 100
mu <- 0.1
sig <- 0.4
# mean and variance of N
mu_N <- lambda
var_N <- lambda
# mean and variance of X
mu_X <- exp(mu+1/2*sig^2)
var_X <- exp(2*mu+sig^2)*(exp(sig^2)-1)
# mean, variance and skewness of SN
mu_SN <- mu_N*mu_X
var_SN <- mu_N*var_X+mu_X^2*var_N
skew_SN <- exp(3*mu+1/2*9*sig^2)/sqrt(lambda*(exp(2*mu+1/2*4*sig^2))^3)
mu_SN; var_SN; skew_SN
# translated gamma approximation
alpha <- (2/skew_SN)^2 # 247.5134
beta <- sqrt(alpha/var_SN) # 1.213061
k <- mu_SN - alpha/beta # -84.31853
alpha;beta;k
# empirical distribution
M <- 100000
SN_tmpt_sim <- c()
for(i in 1:M){
  N <- rpois(1,lambda)
  SN_tmpt_sim <- c(SN_tmpt_sim,sum(exp(rnorm(N,mu,sig))))
}
sorted_SN_tmpt_sim <- sort(SN_tmpt_sim)
low_val <- sorted_SN_tmpt_sim[ceiling(0.95*M)]
high_val <- sorted_SN_tmpt_sim[ceiling(0.99999*M)]

empiriclal_cdf <- ecdf(sorted_SN_tmpt_sim)

# approximation distribution
N <- 1000
step <- (high_val-low_val)/N
x_seq <- seq(low_val,high_val,by=step)
normal_approx <- 1 - pnorm(x_seq,mean=mu_SN,sd=sqrt(var_SN))
beta_approx <- 1 - pgamma(x_seq-k,shape=alpha,rate=beta)
emprical_approx <- 1 - empiriclal_cdf(x_seq)

plot(x_seq,normal_approx,type='l',log='xy',xlab='X',ylab = '1-F(SN)',title='log-log plot',col='blue',lty='twodash')
points(x_seq,beta_approx,type='l',log='xy',col='red',lty='solid')
points(x_seq,emprical_approx,type='l',log='xy',col='black',lty='dotdash')
legend('topright',legend = c("Normal","Gamma","Empirical"),col=c('blue','red','black'),lty=c('twodash','solid','dotdash'))
# parameters summary
mu_SN;var_SN;skew_SN;
alpha;beta;k

# Part 1-3
M <- 100000
alpha <- 0.99
k <- 3
zeta1 <- 0.08/100
nu1 <- 0.2/100
S0 <- 59
Delta <- 100
mu1 <- 0
sig1 <- 0.4/sqrt(252)
Lss_sim <- c()
LLss_sim <- c()
LC_sim <- c()
for(i in 1:M){
  X1 <- rnorm(1,mean=mu1,sd=sig1)
  Lss_tmpt <- -Delta*S0*(exp(X1)-1)
  Lss_sim <- c(Lss_sim,Lss_tmpt)
  LC_tmpt <- Delta*S0*exp(X1)*rnorm(1,mean=nu1,sd=k*zeta1)/2
  LLss_tmpt <- Lss_tmpt + LC_tmpt
  LLss_sim <- c(LLss_sim,LLss_tmpt)
}
sorted_LLss_sim <- sort(LLss_sim)
LVaR_sim <- sorted_LLss_sim[ceiling(alpha*M)]
sorted_Lss_sim <- sort(Lss_sim)
VaR_sim <- sorted_Lss_sim[ceiling(alpha*M)]
# Theoretical VaR
VaR_theoretical <- -Delta*S0*(exp(mu1+sig1*qnorm(1-alpha))-1)
LC_sim <- sorted_LLss_sim[ceiling(alpha*M)]-VaR_theoretical#sorted_Lss_sim[ceiling(alpha*M)]
LVaR_ind <- VaR_theoretical + 1/2*Delta*S0*(nu1+k*zeta1)
LC_ind <- 1/2*Delta*S0*(nu1+k*zeta1)
# output
cat(" Result for Part1 (3)","\n","Confidence: ",alpha, "\n","Empirical liquidity adjusted VaR:",LVaR_sim,"\n","Assuming mean is known","\n","Theoretical VaR based on sample mean and variance:",VaR_theoretical,"\n","Estimated liquidity cost: ",LC_sim,"\n", "Estimated percentage increase in the risk measure: ", 100*(LVaR_sim/VaR_theoretical-1),"\n", "The industry approximation liquidity adjusted VaR: ", LVaR_ind,"\n", "The industry liquidity cost: ", LC_ind,"\n", "The industry percentage increase in the risk measure: ", 100*(LVaR_ind/VaR_theoretical-1))




# Part 2
install.packages("janitor",dependencies = TRUE)
library(janitor)

setwd("/Users/ivy/Documents/11学习/MF731/HW/hw5/")
df <- read.csv('AAPL_Data.csv',header=FALSE)
names(df)[1] <- 'Date'
names(df)[2] <- 'Close'
N <- dim(df)[1]
Date <- excel_numeric_to_date(as.numeric(df[1:N,1]))
df['Date'] <- Date 
logRet <- as.matrix(log(df$Close[2:N]/ df$Close[1:(N-1)]))

mean_log_return  <- mean(logRet)
var_log_return <- var(logRet)
port_value <- 1000000
port_loss <- -port_value*(exp(logRet)-1)
sort_port_loss <- sort(port_loss,decreasing = FALSE)

alpha <- .97
emp_var <- sort_port_loss[ceiling(length(logRet)*alpha)]
sample_var <- port_value*(1-exp(mean_log_return+(var_log_return)^(.5)*qnorm(1-alpha,0,1)))


beta <- .02
low_chi_val <- qchisq(beta/2,length(logRet)-1)
high_chi_val <- qchisq(1-beta/2,length(logRet)-1)

low_sigma_val <- ((length(logRet)-1)*var_log_return/high_chi_val)^(.5)
high_sigma_val <- ((length(logRet)-1)*var_log_return/low_chi_val)^(.5)

low_VaR_val <- port_value*(1-exp(mean_log_return+low_sigma_val*qnorm(1-alpha,0,1)))
high_VaR_val <- port_value*(1-exp(mean_log_return + high_sigma_val*qnorm(1-alpha,0,1)))

M <- 125000
VaR_est <- matrix(0,M,1)

for(m in 1:M){
  chi_square_value <- rchisq(1,length(logRet)-1)
  sigma_est <- ((length(logRet)-1)*var_log_return/chi_square_value)^(.5)
  mu_est <- mean_log_return
  VaR_est[m] <- port_value*(1-exp(mu_est+sigma_est*qnorm(1-alpha,0,1)))
}

VaR_est_sim_mean <- mean(VaR_est)
VaR_est_sim_var <- var(VaR_est)

sort_VaR_est <- sort(VaR_est,decreasing = FALSE)
low_VaR_val_sim <- sort_VaR_est[ceiling(M*beta/2)]
high_VaR_val_sim <- sort_VaR_est[ceiling(M*(1-beta/2))]