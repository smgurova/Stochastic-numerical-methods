# Black-Scholes model

#Black-Scholes assumes that log(S_T) are normally distributed

#log(S_T)???N(log(S_0)+(???????^2/2)T,??*sqrt(T))
#dS/S???N(????t,??sqrt(??T)),
#where S0 is the current stock price

#ex.1
#Consider astock with S_0=30, expected return=16% p.a.and volatility=20% p.a.
#The probability distribution S_T in 6 mothns will be

S0 <- 40
mu <- 16/100 
sd <- 20/100
T <- 0.5
(mean <- log(S0) + (mu-sd^2/2)*T) #[1] 3.758879
(var   <- sd^2*T)  #[1] 0.02
(stdev <- sqrt(var)) #[1] 0.1414214
# CI of ln(S_T) 
(CI <- c(mean-1.96*stdev, mean+1.96*stdev))   #[1] 3.481694 4.036065
# Plot distribution
x <- seq(3,5,by = .01)
dens <- dnorm(x,mean = mean, sd = stdev)
plot(x     , dens , type='l', main= "distribution of ln(S_t)")
i <- x >= 3.48 & x <= 4.04
polygon(c(3.48,x[i],4.04), c(0,dens[i],0), col="lightblue") 
# This implies CI of S_T = 
(CI2 <- exp(CI)) #[1] 32.51474 56.60319
plot(exp(x), dens , type='l', main= "distribution of S_t")
i <- exp(x) >= 32.5 & exp(x) <= 56.6
polygon(c(32.5,exp(x[i]),56.6), c(0,dens[i],0), col="lightblue") 

#ex.2 Pricing function
#call option
Call <- function(S, K, r, T, sigma) {
  d1  <-  (log(S/K) + (r + sigma^2/2)*T) / (sigma*sqrt(T))
  d2  <-  d1 - sigma*sqrt(T)
  S * pnorm(d1)  - K*exp(-r*T)*pnorm(d2)
}
#put option
Put <- function(S, K, r, T, sigma) {
  d1  <-  (log(S/K) + (r + sigma^2/2)*T) / (sigma*sqrt(T))
  d2  <-  d1 - sigma*sqrt(T)
  -S * pnorm(-d1) + K*exp(-r*T)*pnorm(-d2)
}
  # Parameters
  S     <- 100
  K     <- 70
  r     <- 0.05
  T     <- 1       # seq(10/52,1/52,-1/52)
  sigma <- c(.16,.21,.3, .4,.5,.6,.75,1,1.06) 
  
  d1  <-  (log(S/K) + (r + sigma^2/2)*T) / (sigma*sqrt(T))
  d2  <-  d1 - sigma*sqrt(T)
  p   <-  cbind(S, K, T, sigma,  p = -S * pnorm(-d1) + K*exp(-r*T)*pnorm(-d2))  
  c   <-  cbind(S, K, T, sigma,  c =  S * pnorm(d1)  - K*exp(-r*T)*pnorm(d2)) 
 # Probability of exercise:
  pnorm(d2)  
  #[1] 0.9930863 0.9664906 0.8860109 0.7929464 0.7134017 0.6472073 0.5664067 0.4628227
  #[9] 0.4418248
  #Expected value of payoff at maturity:
  K*pnorm(d2) 
#[1] 69.51604 67.65434 62.02076 55.50625 49.93812 45.30451 39.64847 32.39759 30.92773
  
  #?? is the time value of an option. Here we are looking at
  
#evolution of stocks and calls
  cprct <- c[,c(1,5)] / matrix(rep(c[1,c(1,5)],each=length(c[,5]) ),ncol = 2) * 100
ts.plot(cprct, col = c(1,4), main = "Evolution of Stocks and Calls in %")
legend('topleft', c(paste("Stock +", tail(cprct,1)[1]-100, "%") , 
                      paste("Long Call +", round(tail(cprct,1)[2])-100, "%")), 
         col=c(1,4), lty=1)

#evolustion of stocks and puts
pprct <- p[,c(1,5)] / matrix(rep(p[1,c(1,5)],each=length(c[,5]) ),ncol = 2) * 100
ts.plot(pprct, col=1:2, main = "Evolution of Stocks and Puts in %")

# Put-Call Parity
c[,5] + K*exp(-r*T)
#[1] 100.0229 100.1714 100.9814 102.6007 104.7023 107.0822 110.9202 117.5396 119.1218
p[,5] + S 
#[1] 100.0229 100.1714 100.9814 102.6007 104.7023 107.0822 110.9202 117.5396 119.1218

#Find Implicit Volatility
S <- 40
C <- 1.50
VperSigma <- Call(S,K,r,T,seq(0.01,0.50,0.005))
IV <- which.min(abs(VperSigma-C))/2
print(paste("Implicit Volatility is around", IV, "%)"))
#[1] "Implicit Volatility is around 45 %)

#parameters
mu <- 0
sd <- 0.2
x <- rlnorm(10000, mu,sd)
plot(density(x), xlim=c(0,4)) ; abline(v=c(1, 1.5), lty=2)

dlnorm(1.2, mu,sd,log = T)   # pdf [1] 0.09266345
plnorm(1.54, mu,sd)          # CDF [1] 0.9845715
qlnorm(0.95, mu,sd)          # Quantile Function [1] 1.389537

par(mfrow = c(1, 2))
curve(dnorm, from = -3, to = 3, n = 1000, main = "Normal Probability Density Function")
text(-2, 0.3, expression(f(x) == paste(frac(1, sqrt(2 * pi * sigma^2)),
                                       " ", e^{
                                         frac(-(x - mu)^2, 2 * sigma^2)
                                       })), cex = 1.2)
x <- dnorm(seq(-3, 3, 0.001))
plot(seq(-3, 3, 0.001), cumsum(x)/sum(x), type = "l", col = "blue",
     xlab = "x", main = "Normal Cumulative Distribution Function")
text(-1.5, 0.7, expression(phi(x) == paste(frac(1, sqrt(2 * pi)),
                                           " ", integral(e^(-t^2/2) * dt, -infinity, x))), cex = 1.2)
