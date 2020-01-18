#Markov Chain Monte Carlo (MCMC) Integration

# MC integration can still be applied if the generated (dependent) observations have a 
#joint density roughly the same of the joint density of a random, iid sample. 
#To achieve this it is used Markov Chains, which provides the sampler that generates 
#the dependent observations from the target distribution. 

#The Metropolis-Hastings algorithm is a MCMC method that tries to achieve this task. 
#The main idea is to generate a Markov Chain Xt|t=0,1,… such that its stationary distribution is the target distribution. 
#The algorithm, given Xt, knows how to compute Xt+1. 
#To do that it must be able to generate a candidate point Y from a proposal distribution g(???|Xt)
#which (probably) depends on the previous state. 
#This point Y may or may not be accepted. If it is accepted, then Xt+1=Y, 
#otherwise the chain remains in the same place Xt+1=Xt.The proposal distribution g must be 
#chosen so that the generated chain will converge to a stationary distribution
#in this case the target distribution f. The proposal distribution is the way we generate possible 
#good points for the target distribution. If the proposal is not well chosen, the algorithm will 
#produce lots of rejections and the time to converge to the target distribution might take more time than we have available.
#If the generation of candidates does not depend on the current region of the chain,
#the proposal distribution can be independent of xt, and the algorithm 
#will accept the new candidate y if f(y)/g(y)???f(xt)/g(xt)

#https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm

#Here’s the R code (algorithm details can be found at Rizzo’s book, pag.248):

metropolis.hastings <- function(f,  # the target distribution
                                g,  # the proposal distribution
                                rg, # a sample from the proposal distribution
                                x0, # initial value for chain, in R it is x[1]
                                chain.size=1e5,  # chain size
                                burn.perc=0.1) { # burn in percentage
  
  x <- c(x0, rep(NA,chain.size-1))  # initialize chain
  
  for(i in 2:chain.size)   {
    y <- rg(x[i-1])                 # generate Y from g(.|xt) using sampler rg
    alpha <- min(1, f(y)*g(x[i-1],y)/(f(x[i-1])*g(y,x[i-1])))
    x[i] <- x[i-1] + (y-x[i-1])*(runif(1)<alpha)  # update step
  }
  
  # remove initial part of the chain before output result
  x[(burn.perc*chain.size) : chain.size] 
}
#ex. 1
#This first eg samples from an uniform distribution (the proposal distribution) to generate a sample from a Beta(2.7, 6.3) distribution:

a<-2.7; b<-6.3; size<-1e4

f  <- function(x)   dbeta(x,a,b)
rg <- function(x)   runif(1,0,1)
g  <- function(x,y) 1 # i.e., dunif(x,0,1)

X <- metropolis.hastings(f,g,rg,x0=runif(1,0,1),chain.size=size)

par(mfrow=c(1,2),mar=c(2,2,1,1))
hist(X,breaks=50,col="blue",main="Metropolis-Hastings",freq=FALSE)
curve(dbeta(x,a,b),col="sienna",lwd=2,add=TRUE)
hist(rbeta(size,a,b),breaks=50,col="grey",main="Direct Sampling",freq=FALSE)
curve(dbeta(x,a,b),col="sienna",lwd=2,add=TRUE)

#ex.2
#Compute the expected value of 
#f(x)=c*{exp(1/2* ((x-2)^2-2x))+5exp(-1/2*(x+2)^4)}

# first plot it
#  c is 1/30.8636 necessary to make it a density, but we didn't need to know it
f <- function(x) (exp(-((x-2)^4-2*x)/2) + 5*exp(-(x+2)^4/2)) / 30.8636 

xs <- seq(-5,5,len=100)
plot(xs,f(xs),type="l",col="red",lwd=2)

#To find Ef[X] we’ll use the Metropolis-Hastings algorithm. 
#The candidate function after xt will be q(y|xt)???N(xt,??2). 
#In our first test we choose ??=0.1:
g  <- function(x, y) dnorm(x,y,0.1)
rg <- function(x)    rnorm(1,x,0.1)

set.seed(101)
X <- metropolis.hastings(f,g,rg,x0=1,chain.size=5e4)
mean(X) # the answer?
# [1] 2.482991

#This value seems wrong. Let’s compare the histogram of chain X with the true density:
hist(X,breaks=50,col="blue",xlim=c(-5,5),main="Metropolis-Hastings",freq=FALSE)
#
#What happened? Since the candidate function is a normal with a very short ?? the potential
#candidates that it produces are very close to the last xt which means the algorithm
#is unable to cross 0 to the left side. We can check what happens if we start at a negative x0.

set.seed(101)
X1 <- metropolis.hastings(f,g,rg,x0=-2,chain.size=5e4)
mean(X1)  ## [1] -1.928099

hist(X1,breaks=50,col="blue",xlim=c(-5,5),main="Metropolis-Hastings",freq=FALSE)
curve(f(x),col="red",lwd=2,add=TRUE)

#Precisely what was expected, now the chain is unable to cross to the right side.
#Let’s visualize a bit of both previous markov chains and see how they are unable to jump to the other side of the bimodal density:
par(mfrow=c(2,1),mar=c(2,2,1,1))
plot(1:3000,X[1:3000], lwd=2,type="l",ylim=c(-4,4))
plot(1:3000,X1[1:3000], lwd=2,type="l",ylim=c(-4,4))

#So let’s try a higher sigma, say ??=1
g  <- function(x, y) dnorm(x,y,1)
rg <- function(x)    rnorm(1,x,1)

set.seed(101)
X2 <- metropolis.hastings(f,g,rg,x0=runif(1,-4,4),chain.size=5e4)
mean(X2) # the answer [1] 0.797984

par(mfrow=c(2,1),mar=c(2,2,1,1))
hist(X2,breaks=50,col="blue",xlim=c(-5,5),main="Metropolis-Hastings",freq=FALSE)
curve(f(x),col="red",lwd=2,add=TRUE)
plot(1:3000,X2[1:3000], lwd=2,type="l",ylim=c(-4,4))

#Now the candidate function is able to make longer jumps, and both parts of 
#the density are visited, providing a good estimate of the true density.
#An exagerated value of sigma will have another type of disadvantage: most candidates will be
#so far the interesting area that they will be simply rejected which will result in a poor estimate.

g  <- function(x, y) dnorm(x,y,100) # sigma = 100 (!)
rg <- function(x)    rnorm(1,x,100)

set.seed(101)
X3 <- metropolis.hastings(f,g,rg,x0=runif(1,-4,4),chain.size=5e4)
mean(X3) #[1] 0.9172019

par(mfrow=c(2,1),mar=c(2,2,1,1))
hist(X3,breaks=50,col="blue",xlim=c(-5,5),main="Metropolis-Hastings",freq=FALSE)
curve(f(x),col="red",lwd=2,add=TRUE)
plot(1:3000,X3[1:3000], lwd=2,type="l",ylim=c(-4,4))

#The plateau’s above show repeated rejections, making the chain stay at the last change of xt.
