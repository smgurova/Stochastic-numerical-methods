### Inverse Transformation for continuous random variables
#The method needs the inverse function F^(???1)_{X}(u) to work.
#For each value required the method follows:
#1.Generate a random u from U(0,1)
#2.Return F???1X(u)
#Here’s the algorihm:

# return n samples based on the inverse of the target cdf
inv.transform <- function(inv.f, n) {
  
  # Non-vectorized version (for explanation purposes only)  
  #
  #   result.sample <- rep(NA,n)
  #   
  #   for (i in 1:n) {
  #     u <- runif(1,0,1)              # step 1
  #     result.sample[i] <- inv.f(u)   # step 2
  #   }
  #   
  #   result.sample
  
  # Vectorized version
  inv.f(runif(n,0,1))
}

#ex.1 
#Used this method to simulate a random variable which pdf is f_{X}(x)=3x2,0<x<1.
# The cdf is F_{X}(x)=int_{0}^{x} 3*t^2 dt=x^3.
# The inverse is F^(-1)(u)=u^(1/3).

inv.f <- function(u) u^(1/3)

vals <- inv.transform(inv.f, 5e4)

# Checking if it went well
hist(vals, breaks=50, freq=FALSE, main=expression("Sample vs true Density [ f(x)=" ~3*x^2~"]"))
curve(3*x^2, 0, 1, col="red", lwd=2, add=T)

#ex.2
# Let pdf f(x)=(1/3)*x^2 for -1<x<2
# The cdf is F(x)=(1+x^3)/9 and the inverse F^(-1)(u)=(9*u-1)^(1/3).

# The next R code is like this due to the fact that cubic squares of negative numbers also have complex roots
# But basically this is 
#   inv.f <- function(u) (9*u-1)^(1/3)
inv.f <- function(u) ifelse((9*u-1)>=0,  (9*u-1)^(1/3),  -(1-9*u)^(1/3))

vals <- inv.transform(inv.f, 5e4)

# Checking if it went well
hist(vals, breaks=50, freq=FALSE, main="Sample vs true Density")
curve((1/3)*x^2, -1, 2, col="red", lwd=2, add=T)

