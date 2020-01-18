#Vasicek Interest Rate Model

#https://en.wikipedia.org/wiki/Vasicek_model

#dr_t=K(theta-r_t)dt+betadW_t,
#W_t is a Wiener process under the risk neutral framework modelling the random market risk factor;
#theta: "long term mean level". All future trajectories of r will evolve around a mean level theta in the long run;
#K: "speed of reversion".K characterizes the velocity at which such trajectories will regroup around theta in time;
#beta : "instantaneous volatility", measures instant by instant the amplitude of randomness entering the system. 
#Higher beta implies more randomness.
#r_0- initial data
#
#basic formula: r_t=r_0*e^(-kt)+theta(1-e^(-kt))+beta*e^(-kt)int_{0}^{t} e^(ks)ds
#r_T ~ N ( theta+(r_0-theta)e^{-kT}, {beta^2}/{2k}(1-e^{-2kT})).

## Simulate Sample Paths ##

## define model parameters
r0 <- 0.03
theta <- 0.10
k <- 0.3
beta <- 0.03

## simulate short rate paths
n <- 10    # MC simulation trials
T <- 10    # total time
m <- 200   # subintervals
dt <- T/m  # difference in time each subinterval

r <- matrix(0,m+1,n)  # matrix to hold short rate paths
r[1,] <- r0

for(j in 1:n){
  for(i in 2:(m+1)){
    dr <- k*(theta-r[i-1,j])*dt + beta*sqrt(dt)*rnorm(1,0,1)
    r[i,j] <- r[i-1,j] + dr
  }
} 

## plot paths
t <- seq(0, T, dt)
rT.expected <- theta + (r0-theta)*exp(-k*t)
rT.stdev <- sqrt( beta^2/(2*k)*(1-exp(-2*k*t)))
matplot(t, r[,1:10], type="l", lty=1, main="Short Rate Paths", ylab="rt") 
abline(h=theta, col="red", lty=2)
lines(t, rT.expected, lty=2) 
lines(t, rT.expected + 2*rT.stdev, lty=2) 
lines(t, rT.expected - 2*rT.stdev, lty=2) 
points(0,r0)

## function to find ZCB price using Vasicek model
VasicekZCBprice <- 
  function(r0, k, theta, beta, T){
    b.vas <- (1/k)*(1-exp(-T*k)) 
    a.vas <- (theta-beta^2/(2*k^2))*(T-b.vas)+(beta^2)/(4*k)*b.vas^2
    return(exp(-a.vas-b.vas*r0))
  }

## define model parameters for plotting yield curves
theta <- 0.10
k <- 0.5
beta <- 0.03

r0 <- seq(0.00, 0.20, 0.05)
n <- length(r0)
yield <- matrix(0, 10, n)
for(i in 1:n){
  for(T in 1:10){
    yield[T,i] <- -log(VasicekZCBprice(r0[i], k, theta, beta, T))/T
  }
}

maturity <- seq(1, 10, 1)
matplot(maturity, yield, type="l", col="black", lty=1, main="Yield Curve Shapes")
abline(h=theta, col="red", lty=2)

#Formula vs Monte Carlo

## define model parameters
r0 <- 0.03
theta <- 0.10
k <- 0.3
beta <- 0.03

## simulate short rate paths
n <- 1000  # MC simulation trials
T <- 1     # total time
m <- 200   # subintervals
dt <- T/m  # difference in time each subinterval

r <- matrix(0,m+1,n)  # matrix to hold short rate paths
r[1,] <- r0
for(j in 1:n){
  for(i in 2:(m+1)){
    dr <- k*(theta-r[i-1,j])*dt + beta*sqrt(dt)*rnorm(1,0,1)
    r[i,j] <- r[i-1,j] + dr
  }
} 

## calculate Monte Carlo bond price and compare to Exact Vasicek solution
ss <- colSums(r[2:(m+1),]*dt)  # integral estimate
c <- exp(-ss)
estimate <- mean(c)
stdErr <- sd(c)/sqrt(n)
exact <- VasicekZCBprice(r0, k, theta, beta, T)

cat('Exact Vasicek Price:', round(exact,4), 'n')
cat('MC Price:', round(estimate,4), 'n')
cat('MC Standard Error:', round(stdErr,4), 'n')

#Exact Vasicek Price: 0.9614 n> cat('MC Price:', round(estimate,4), 'n')
#MC Price: 0.9617 n> cat('MC Standard Error:', round(stdErr,4), 'n')
#MC Standard Error: 5e-04 n