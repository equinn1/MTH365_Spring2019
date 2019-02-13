## rethinking package exercises

# Estimating the probability of success from a binomial sample

N = 100              #sample size

n = rpois(N,10)      #number of trials

p = 0.35             #p value

y=rbinom(N,n,p)      #generate the vector of binomial random variables

df = data.frame(y,n)

# Use ulam to run the MCMC simulation
m1 <- ulam(
  alist(
    y ~ dbinom(n,p),
    p ~ beta(1,1)
  ) ,
  data=df, chains=2, cores=1 , sample=TRUE )

m1@stanfit

pd = extract(m1@stanfit)
str(pd)

quantile(pd$p,c(.005,.025,.05,.25,.5,.75,.95,.975,.995))

mean(pd$p < 0.37)      #estimated probability that p<.37

#Estimating the parameters of a beta distribution

N = 100

a=3.0
b=2.0

y=rbeta(N,a,b)

df = data.frame(y)

m <- ulam(
  alist(
    y ~ dbeta(a,b),
    c(a,b) ~ dexp(0.3)
  ) ,
  data=df, chains=2, cores=1 , sample=TRUE )

m@stanfit

pd = extract(m@stanfit)
str(pd)

quantile(pd$a,c(.005,.025,.05,.25,.5,.75,.95,.975,.995))
quantile(pd$b,c(.005,.025,.05,.25,.5,.75,.95,.975,.995))

# estimated 90th percentile of y

yp = qbeta(0.9,pd$a,pd$b)
str(yp)

quantile(yp,c(.005,.025,.05,.25,.5,.75,.95,.975,.995))

#Simple linear regression  y = a+bx+e

N = 80

a=-8
b=1.0
sigma_e = 3

x=seq(1:N)/20

y = a + b*x + rnorm(N,0,sigma_e)

df = data.frame(y,x)

m <- ulam(
  alist(
    y ~ dnorm(mu,sigma),
    mu <- A + B*x,
    c(A,B) ~ dnorm(0,10),
    sigma ~ dexp(1)
  ) ,
  data=df, chains=2, cores=1 , sample=TRUE )

m@stanfit

pd = extract(m@stanfit)
str(pd)

quantile(pd$A,c(.005,.025,.05,.25,.5,.75,.95,.975,.995))
quantile(pd$B,c(.005,.025,.05,.25,.5,.75,.95,.975,.995))

#confidence interval for fitted value at x=1.5

y1 = pd$A + 1.0*pd$B
quantile(y1,c(.005,.025,.05,.25,.5,.75,.95,.975,.995))
