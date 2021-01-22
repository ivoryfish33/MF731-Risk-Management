rm(list = ls())

# Problem 4
# (a)
T <- seq(0, 10, 0.05)
w1 <- 1
w1_bar <- 0
w2 <- 0
w2_bar <- 1
alpha1 <- 0.1
alpha2 <- 0.1
beta1 <- 1
beta2 <- 0.9
sigma1 <- 0.3
sigma2 <- 0.2

gamma1 <- sqrt(beta1^2 + 2 * (w1 + w1_bar) * sigma1^2)
gamma2 <- sqrt(beta2^2 + 2 * (w2 + w2_bar) * sigma2^2)

h11 <- ((2 * gamma1 * exp(T / 2 * (gamma1 + beta1))) / ((gamma1 + beta1) * (exp(gamma1 * T) - 1) + 2 * gamma1))^(2 * alpha1 / sigma1^2)
h12 <- ((2 * gamma2 * exp(T / 2 * (gamma2 + beta2))) / ((gamma2 + beta2) * (exp(gamma2 * T) - 1) + 2 * gamma2))^(2 * alpha2 / sigma2^2)
h21 <- 2 * (exp(gamma1 * T) - 1) / ((gamma1 + beta1) * (exp(gamma1 * T) - 1) + 2 * gamma1)
h22 <- 2 * (exp(gamma2 * T) - 1) / ((gamma2 + beta2) * (exp(gamma2 * T) - 1) + 2 * gamma2)

# Assume the default-free short rate r(t) and the default intensity lambda(t) are constants
def <- h11 * exp(-h21) * h12 * exp(-h22)
plot(T, def, type = 'l', main = 'T vs price of the defaultable zero-coupon bond',
     xlab = 'T', ylab = 'price of the defaultable zero-coupon bond')

# (b)
T <- 1
mu <- seq(0, 1, 0.001)
w1 <- 1
w1_bar <- mu
w2 <- 0
w2_bar <- 1 - mu
alpha1 <- 0.1
alpha2 <- 0.1
beta1 <- 1
beta2 <- 0.9
sigma1 <- 0.3
sigma2 <- 0.2

gamma1 <- sqrt(beta1^2 + 2 * (w1 + w1_bar) * sigma1^2)
gamma2 <- sqrt(beta2^2 + 2 * (w2 + w2_bar) * sigma2^2)

h11 <- ((2 * gamma1 * exp(T / 2 * (gamma1 + beta1))) / ((gamma1 + beta1) * (exp(gamma1 * T) - 1) + 2 * gamma1))^(2 * alpha1 / sigma1^2)
h12 <- ((2 * gamma2 * exp(T / 2 * (gamma2 + beta2))) / ((gamma2 + beta2) * (exp(gamma2 * T) - 1) + 2 * gamma2))^(2 * alpha2 / sigma2^2)
h21 <- 2 * (exp(gamma1 * T) - 1) / ((gamma1 + beta1) * (exp(gamma1 * T) - 1) + 2 * gamma1)
h22 <- 2 * (exp(gamma2 * T) - 1) / ((gamma2 + beta2) * (exp(gamma2 * T) - 1) + 2 * gamma2)

# Assume the default-free short rate r(t) and the default intensity lambda(t) are constants
def <- h11 * exp(-h21) * h12 * exp(-h22)
plot(mu, def, type = 'l', main = 'u vs price of the defaultable zero-coupon bond',
     xlab = 'u', ylab = 'price of the defaultable zero-coupon bond')