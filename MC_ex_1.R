# Monte Carlo MethodsList
#
# E[g(x)]=int_{I^d} g(x)f(x) dx,  x~f(x) //pdf// 
#
#Let ??=E[g(X)]. The Monte Carlo algorithm to compute an approximation to ??, 
#denoted by ??*_{MC} is as follows:
#Algorithm:
#1. Generate n samples from f(x): X1,X2,…,Xn.
#2. The MC approximation is given by, ??*_{MC}=1/n*???_{i=1}^{n}g(Xi).

#ex.1
# Compute the integral of g(x)=(cos50x+sin20x)^2 from 0 to 1.
#
#1.Analytically

g <- function(x){
  (cos(50*x) + sin(20*x))^2
}
integrate(g, 0, 1)
### 0.9652009 with absolute error < 1.9e-10

# 2. MC
n.sam <- 1000
# step 1: generate n i.i.d samples from f (in this case uniform(0,1))
x.sam <- runif(n.sam)

# compute the MC approximation
theta.mc <- sum(sapply(x.sam, g))/n.sam
#[1] 0.9434973

# increasing n=1000
n.sam <- 10000
# step 1: generate n i.i.d samples from f (in this case uniform(0,1))
x.sam <- runif(n.sam)

# compute the MC approximation
theta.mc <- sum(sapply(x.sam, g))/n.sam
print(theta.mc)

#ex.2
# Let X???N(??,??2) and let t be a real number. Approximate the probability P(X<t) using MC integration.
#For concreteness take, ??=5, ??=1.25 and t=7.1.

#1. Directly
pnorm(7.1, mean = 5, sd = 1.25) #[1] 0.9535213

#2. MC
g.ind <- function(x, t){
  if (x < t) return(1)
  else return(0)
}
n.sam <- 1000
# Step 1: Sample of size n from the density f (in this case N(5, 1.25))
x.sam <- rnorm(n.sam, mean = 5, sd = 1.25)

# Step 2: get the monte carlo estimate
theta.mc <- sum(sapply(x.sam, g.ind, t = 7.1))/n.sam
print(theta.mc) #0.94

#  try with t=9.1 // P(x>9.1)=?
g.ind2 <- function(x, t){
  if (x > t) return(1)
  else return(0)
}
n.sam <- 1000
# Step 1: Sample of size n from the density f (in this case N(5, 1.25))
x.sam <- rnorm(n.sam, mean = 5, sd = 1.25)
# Step 2: get the monte carlo estimate
theta.mc <- sum(sapply(x.sam, g.ind2, t = 9.1))/n.sam
print(theta.mc) #[1]0.002

1-pnorm(9.1, mean = 5, sd = 1.25) #[1] 0.0005190354

#ex.3
#Let A={(x,y):x2+y2???1. Compute the Area(A) using Monte Carlo Integration.

x.sam <- runif(10000, min = -1, max = 1)
y.sam <- runif(10000, min = -1, max = 1)
joint.sam <- cbind(x.sam, y.sam)

g.indA <- function(point){
  if ((point[1]^2 + point[2]^2) <= 1) return (1.0)
  else return(0)
}
theta.mc <- sum(apply(joint.sam, 1, g.indA))/nrow(joint.sam)
4*theta.mc  #[1] 3.1448

#Visualization
joint.data <- data.frame(joint.sam)
colnames(joint.data) <- c("x", "y")
# create a factor variable to see which were accepted
joint.data$Accepted <- (joint.data$x^2 + joint.data$y^2 < 1)
library(ggplot2)
# we will create a scatter plot
# "x" values will be y (i.e samples from g)
# "y" values will be ue (u * e)

# If (u*e < f(y)) then samples are accepted otherwise rejected
# Hence, in scatter plot we will color by the Accepted variable

plot_ar <- ggplot(joint.data, aes(x=x, y=y)) + 
  geom_point(shape=20, aes(color=Accepted)) +
  stat_function(fun = function(x) sqrt(1-x^2), size=1.5) + 
  stat_function(fun = function(x) -1*sqrt(1-x^2), size=1.5) +
  coord_fixed()


# beautify : Increase font size in axes 
plot_ar + theme(axis.text.x = 
                  element_text(face = "bold",size = 12),
                axis.text.y = 
                  element_text(face = "bold", size = 12),
                axis.line = element_line(colour = "black", 
                                         size = 1, linetype = "solid"))

#ex.4 Estimate int_{2}^{4} e^(-x)dx by MC.
# pre-condition: a < b
MC.simple.est <- function(g, a, b, n=1e4) {
  xi <- runif(n,a,b)      # step 1
  g.mean <- mean(g(xi))   # step 2
  (b-a)*g.mean            # step 3
}
g <- function(x) exp(-x)

MC.simple.est(g, 2, 4)  #[1] 0.1162728
