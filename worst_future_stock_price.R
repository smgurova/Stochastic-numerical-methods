#Worst-case scenario for future stock price

#the Normal distribution can be useful in constructing Monte Carlo simulations
#and it is still commonly found in applications such as calculating the Value 
#at Risk (VaR) of a portfolio, pricing options, and estimating the liabilities 
#in variable annuity contracts. 

#ex.1
#Assume we want to calculate the worst-case scenario of a future stock price.
#This problem called value at risk is heavily used in risk management.
#By “worst-case scenario” we mean the value that the stock price 
#will exceed with 99% probability 
#(i.e., there is only 1% probability that the stock price will be below). 
#The current price of our stock is 100 $. 
#We want to see the possible future prices after 20 trading days.

#Assumptions:

#1. Let the drift over the 20 trading days is 10% (i.e. in average the price
#will go up to 110) and the volatility is 20%.
#2. Generate some random inputs, which in our case means to model the future price given the current price,
#drift and the volatility. For our purpose we do not need anything fancy so we will use 
#the standard stock price model called Geometric Brownian Motion:
  
  #S_{t+1} = S_{t} * (1 + ????t + ????srqt(??t)),
#Where S_{t} is price in time t,
#?? is our drift, 
#??t means one period (in our case 1 trading day, i.e. one twentieth of the period),
#?? is our volatility,
#?? is a random number between 0 and 1.

#The core idea of Monte Carlo method is to generate the future price (which is random)
#high number of times to simulate what are all the situations that can occur.

#' Stock price calculation
#' Calculates stock price after n periods using standard stock price model
#' @param stock_price original stock price
#' @param n number of periods
#' @param stock_mu expected percentual stock drift over n periods
#' @param stock_sigma expecter percentual stock volatility
#' @return stock price after n periods
f_stock_return <- function(stock_price, n, stock_mu, stock_sigma){
  delta_t <- 1/n # one period
  for (i in seq(n)){
    epsilon <- runif(n=1, min=0, max=1) # random generated number
    # calculate stock price (using quantile function of normal distribution)
    stock_price <- stock_price * (1 + qnorm(epsilon, 
                                            stock_mu * delta_t, 
                                            stock_sigma* sqrt(delta_t)))
  }
  return(stock_price)
}

# parameters
simulations <- 10000 # number of MC simulations
n <- 20 # trading days
stock_price <- 100
stock_mu <- .1 # drift 10%
stock_sigma <- .2 # volatility 20%

# Monte Carlo simulations
set.seed(42) # for reproducibility
stock_prices <- c()
for (i in seq(simulations)){
  stock_prices <- c(stock_prices,
                    f_stock_return(stock_price=stock_price, 
                                   n=n, 
                                   stock_mu=stock_mu, 
                                   stock_sigma=stock_sigma))
}
quantile(stock_prices, c(.01, .05))

#    1%       5% 
#  67.46501 77.69902
#We can see that in 9 900 (99%) scenarios the price is larger than  67 $;
#95% of them were higher than  77 $.
#Which gives us needful assessment of how bad our investment could go.


#ex.2 (Single Security Example)
#suppose we are interested in constructing a distribution of quarterly returns, 
#where ?? = 10% and ??= 15%.  In order to get a reasonable approximation of the distribution, 
#we will generate n = 10,000 returns.

n <- 10000
# Fixing the seed gives us a consistent set of simulated returns
set.seed(106)
z <- rnorm(n)     #epsilon   # mean = 0 and sd = 1 are defaults
mu <- 0.10
sd <- 0.15  #sigma
delta_t <- 0.25
# apply to expression (*) above
qtr_returns <- mu*delta_t + sd*z*sqrt(delta_t) 
hist(qtr_returns, breaks = 100, col = "green")

#The symmetric bell shape of the histogram is consistent with the Normal assumption.  Checking the annualized mean and variance of the simulated returns,

stats <- c(mean(qtr_returns) * 4, sd(qtr_returns) * 2)   # sqrt(4)
names(stats) <- c("mean", "volatility")
stats
  #   mean volatility 
#0.09901252 0.14975805 
#which is very close to our original parameter settings of ?? = 10% and ??= 15%.