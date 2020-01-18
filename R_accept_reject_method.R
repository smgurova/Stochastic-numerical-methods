# Acceptance-Rejection Method
#Problem: Generate X???f from an arbitray pdf f (especially when it’s hard to sample from f).
#We must find Y???g under the only restriction that ???f(x)>0:f(x)<c.g(x),c>1.
#Instead of sampling from f(x) which might be difficult, we use c.g(x) to sample instead.
#For each value required the method follows:
#1.Generate a random y from g
#2.Generate a random u from U(0,1)
#3.If u<f(y)/(c.g(y)) then return y else reject y and goto 1.
#The algorihtm in R:
# generate n samples from f using rejection sampling with g (rg samples from g)
accept.reject <- function(f, c, g, rg, n) { 
  n.accepts     <- 0
  result.sample <- rep(NA, n)
  
  while (n.accepts < n) {
    y <- rg(1)               # step 1 Generate a random y from  g
    u <- runif(1,0,1)        # step 2 Generate a random u from U(0,1)
    if (u < f(y)/(c*g(y))) { # step 3  (accept)
      n.accepts <- n.accepts+1
      result.sample[n.accepts] = y
    }
  }
  
  result.sample
}

#ex.1
#Generate samples from distribution Beta(2,2), where we use the uniform has g, since f(x)<2??g(x).

f  <- function(x) 6*x*(1-x)     # pdf of Beta(2,2), maximum density is 1.5
g  <- function(x) x/x           # g(x) = 1 but in vectorized version
rg <- function(n) runif(n,0,1)  # uniform, in this case
c  <- 2                         # c=2 since f(x) <= 2 g(x)

vals <- accept.reject(f, c, g, rg, 10000) 

# Checking if it went well
hist(vals, breaks=30, freq=FALSE, main="Sample vs true Density")
xs <- seq(0, 1, len=100)
lines(xs, dbeta(xs,2,2), col="red", lwd=2)

#Let’s visualize the method accepting (green dots) or rejecting (red dots) at some specified segments.
xs <- seq(0, 1, len=100)
plot(xs, dbeta(xs,2,2), ylim=c(0,c*1.3), type="l", col="red", lwd=2, ylab="densities")
lines(xs, c*g(xs), type="l", col="blue", lwd=2)
legend("topleft",c("f(x)","c*g(x)"), col=c("red","blue"), lwd=2) 

draw.segment <- function(begin.segment, end.segment) {
  segments(c(begin.segment,end.segment,end.segment,begin.segment), c(0,0,c*1.025,c*1.025), 
           c(end.segment,end.segment,begin.segment,begin.segment), c(0,c*1.025,c*1.025,0))
  n.pts <- 100
  us <- runif(n.pts, 0, 1)
  ys <- begin.segment + rg(n.pts)*(end.segment-begin.segment)
  accepted <- us < f(ys)/(c*g(ys))
  points(ys, c*us, col=ifelse(accepted,"green","red"), pch=18)  
}

draw.segment(0.10, 0.20) 
draw.segment(0.45, 0.55)
draw.segment(0.90, 1.00)

# if c=10

c <- 10

xs <- seq(0, 1, len=100)
plot(xs, dbeta(xs,2,2), ylim=c(0,c*1.25), type="l", col="red", lwd=2, ylab="densities")
lines(xs, c*g(xs), type="l", col="blue", lwd=2)
legend("topleft",c("f(x)","c*g(x)"), col=c("red","blue"), lwd=2) 

draw.segment(0.10, 0.20)
draw.segment(0.45, 0.55)
draw.segment(0.90, 1.00)

#ex.2
#Let’s try with the initial eg for the pdf fX(x)=3x2.
#The function needs the log of the pdf, in this case, log(3x2) and its first derivative:
#d(log(3x2))/dx=2x

library(ars)

log.f <- function(x) log(3*x^2) # log(f(x))
log.f.dx <- function(x) 2/x     # d/dx log(f(x))

vals <- ars(1e4,                # how many points are required
            log.f, log.f.dx,    # the needed functions
            lb=TRUE, ub=TRUE,   # there are a lower and upper bounds for the pdf
            xlb=0, xub=1,       # and which are those bounds
            x=c(.1,.5,.9))      # some initial points inside the pdf

# Checking if it went well
hist(vals, breaks=30, freq=FALSE, main="Sample vs true Density")
xs <- seq(0, 1, len=100)
curve(3*x^2, 0, 1, col="red", lwd=2, add=T)
