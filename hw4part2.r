# # HW4 
# part 2-3

setwd("/Users/ivy/Documents/11学习/MF731/HW/hw4/")
df <- read.csv('SP500_Log_Returns.csv',header=FALSE)
names(df)[1] <- 'Date'
names(df)[2] <- 'Log_Ret'
N <- dim(df)[1]

alpha_seq <- seq(0.99,0.9999,by=0.000099)
# (1) empirical distribution 
V <- 1000000
loss_emp <- -V*df[,2]
loss_emp_sorted <- sort(loss_emp)
VaR_emp_seq <- c()
for(alpha in alpha_seq){
  VaR_emp_seq <- c(VaR_emp_seq,loss_emp_sorted[ceiling(N*alpha)])
}
# (2) EWMA
M <- 500
mu0 <- mean(df[1:(M-1),2])
var0 <- var(df[1:(M-1),2])
mu <- mu0
var <- var0
lambda <- 0.97
theta <- 0.97
for(i in M:N){
  var <- theta * var + (1 - theta)*(df[i,2]-mu)^2
  mu <- lambda * mu + (1 - lambda)*df[i,2]
}
VaR_ewma_seq <- c()
for(alpha in alpha_seq){
  tmpt <- -V*mu + sqrt(var)*V*qnorm(alpha)
  VaR_ewma_seq <- c(VaR_ewma_seq,tmpt)
}

# (3) GEV 
library(ismev)
library(fExtremes)
V <- 1000000
loss_emp <- -V*df[,2]
N <- length(loss_emp)
n <- 125
m <- N/n
maximum <- c()
for( i in 1:m){
  le <- (i-1)*n +1
  re <- i*n 
  maximum <- c(maximum,max(loss_emp[le:re]))
}
gevfit2 <- gev.fit(maximum)
xi <- gevfit2$mle[3]
mu <- gevfit2$mle[1]
beta <- gevfit2$mle[2]
VaR_gev_seq <- c()
for(alpha in alpha_seq){
  if(xi>0){
    tmpt <- mu-beta/xi*(1-(-n*log(alpha))^(-xi))
    VaR_gev_seq <- c(VaR_gev_seq,tmpt)
  }
}
# (4) GP
N <- 6875
mu <- loss_emp_sorted[ceiling(N*0.95)]
cdf_u_value <- ceiling(N*0.95)/N
excess_loss <- loss_emp[loss_emp > mu]-mu
N <- length(loss_emp[loss_emp > mu])
gpdfit1 <- gpd.fit(loss_emp, threshold = mu) 
beta <- gpdfit1$mle[1]
xi <- gpdfit1$mle[2]
VaR_gp_seq <- c()
alpha_seq <- seq(0.99,0.9999,by=0.000099)
for( alpha in alpha_seq){
  if(xi < 0.0001){
    tmpt <- mu + beta*log((1-cdf_u_value)/(1-alpha))
    VaR_gp_seq <- c(VaR_gp_seq,tmpt)
  } else {
    tmpt <- mu + beta/xi*(((1-cdf_u_value)/(1-alpha))^(xi)-1)
    VaR_gp_seq <- c(VaR_gp_seq,tmpt)
  }	
}
VaR_emp_seq[length(VaR_emp_seq)] # 69088.98
VaR_ewma_seq[length(VaR_ewma_seq)] # 59361.88
VaR_gev_seq[length(VaR_gev_seq)] # 63151.54
VaR_gp_seq[length(VaR_gp_seq)] # 54593.68

# plot
yl <- min(c(VaR_emp_seq,VaR_ewma_seq,VaR_gev_seq,VaR_gp_seq))
yu <- max(c(VaR_emp_seq,VaR_ewma_seq,VaR_gev_seq,VaR_gp_seq))
plot(alpha_seq,VaR_emp_seq,ylim=c(yl,yu),type='l',col='red',xlab = expression(paste(alpha)), ylab = expression(paste('Expected Loss ', L)),yaxt='n') 
mtext(expression(10^4),adj=0,padj=-1,outer=FALSE)
y <- seq(0,70000,by=10000)
axis(side=2,at=y,labels=round(y/10000,1))
points(alpha_seq,VaR_ewma_seq,type='l',col='blue',lty='twodash')
points(alpha_seq,VaR_gev_seq,type='l',col='black',lty='dotdash')
points(alpha_seq,VaR_gp_seq,type='l',col='brown',lty='dotted')
legend('topleft',legend = c("EWMA VaR","EMPIRICAL VaR","GEV VaR","GP VaR"),col=c('blue','red','black','brown'),lty=c('twodash','solid','dotdash','dotted'))
