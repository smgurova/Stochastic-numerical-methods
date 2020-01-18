#Sampling from discrete distributions
#  ex. 1
#suppose that X takes values in S = {1, 2, 3} with probability mass function defined by
#the following table:
# x  1   2   3
# p  p1  p2 p3
#
#To generate from this distribution we partition (0, 1) into the three sub-intervals (0, p1),
#(p1, p1 +p2), and (p1 +p2, p1 +p2 +p3), generate a Uniform(0, 1), and check which interval
#the variable falls into. The following R code does this, and checks the results for p1 = .4,
#p2 = .25, and p3 = .35.

# n is the sample size
# p is a 3-length vector containing the corresponding probabilities
rX <- function(n, p)
{
  # generate the underlying uniform(0,1) variables
  U <- runif(n)
  # storage for the random variables, X
  X <- rep(0,n)
  # which U are in the first interval
  w1 <- which(U <= p[1])
  X[w1] <- 1
  # which U are in the second interval
  w2 <- which( (U > p[1]) & (U < sum(p[1:2])) )
  X[w2] <- 2
  # which U are in the third interval
  w3 <- which( U > sum(p[1:2]) )
  X[w3] <- 3
  return(X)
}
X = rX(10000, c(.4, .25, .35))
mean(X == 1)
#[1] 0.407
mean(X == 2)
#[1] 0.2504
mean(X == 3)
#[1] 0.3426


#ex. 2
#The probability mass function for the poisson with parameter ?? has the form
#p(x) = e^(-lambda)* (lambda^i)/(i!)
#whose sample space is all non-negative integers. The following R program generates from
#this distribution and compared the empirical mass function with the true mass function
#for ?? = 4.

# function to calculate poission CDF at an integer x
# when lambda = L. This does the same thing as ppois()
pois.cdf <- function(x, L)
{
  # grid of values at which to calculate the mass function
  v <- seq(0, x, by=1)
  # return CDF value
  return( sum( exp(-L)*(L^v)/factorial(v) ) )
}
# n: sample size, L: lambda
r.pois <- function(n, L)
{
  U <- runif(n)
  X <- rep(0,n)
  # loop through each uniform
  for(i in 1:n)
  {
    # first check if you are in the first interval
    if(U[i] < pois.cdf(0,L))
    {
      X[i] <- 0
    } else
    {
      # while loop to determine which subinterval,I, you are in
      # terminated when B = TRUE
      B = FALSE
      I = 0
      while(B == FALSE)
      {
        # the interval to check
        int <- c( pois.cdf(I, L), pois.cdf(I+1,L) )
        # see if the uniform is in that interval
        if( (U[i] > int[1]) & (U[i] < int[2]) )
        {
          # if so, quit the while loop and store the value
          X[i] <- I+1
          B = TRUE
        } else
        {
          # If not, continue the while loop and increase I by 1
          I=I+1
        }
      }
    }
  }
  return(X)
}

# generate 1000 Pois(4) random variables
V = r.pois(1000, 4)
# empirical mass function
Mf <- c(0:15)
for(i in 0:15) Mf[i+1] <- mean(V==i)
# plot the observed mass function with the
# the true mass function overlaying it
b <- c(0:15)
plot(b, Mf, xlab="x", ylab="p(x)", main="Empirical mass function", col=2)
lines(b, Mf, col=2)
# overlay with true mass function
points(b, dpois(b,4), col=4)
lines(b, dpois(b,4), col=4)

#ex.3
#suppose you want to generate from a distribution with CDF
#F(x) = x/(1+x) for x ??? (0,???).

# n is the sample size
r.harmonic <- function(n)
{
  # generate uniforms
  U <- runif(n)
  # return F^-1(U)
  return( U/(1-U) )
}
X <- r.harmonic(1000)
# empirical CDF
v <- seq(0, 20, length=1000)
emp.cdf <- rep(0, 1000)
for(i in 1:1000) emp.cdf[i] <- mean(X <= v[i])
# true CDF
true.cdf <- v/(1+v)
plot(v, emp.cdf, xlab="X", ylab="F(X)", main="Empirical vs True CDF", col=2,type="l")
lines(v,true.cdf,col=4)

#ex. 4
#suppose X has CDF F(x)=(1/(1+e^(-x))^(theta)
#where ?? > 0 is a parameter. This distribution is known as the skew logistic distribution,
#which is symmetric when ?? = 1, and skewed otherwise.
#
#
# generate n samples from the skew logistic
# distribution with parameter Theta

r.skewlog <- function(n, Theta)
{
  U <- runif(n)
  return( -log( (1/U)^(1/Theta) - 1) )
}
X <- r.skewlog(1000, 4)
# empirical CDF
v <- seq(-10, 10, length=100)
emp.cdf <- rep(0, 1000)
for(i in 1:1000) emp.cdf[i] <- mean(X <= v[i])
# true CDF
true.cdf <- (1 + exp(-v))^(-4)
plot(v, emp.cdf, xlab="X", ylab="F(X)", main="Empirical vs True CDF", col=2,type="l")
lines(v,true.cdf,col=4)
         