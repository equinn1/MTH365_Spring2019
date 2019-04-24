## R code 15.1
# simulate a pancake and return randomly ordered sides
sim_pancake <- function() {
  pancake <- sample(1:3,1)
  sides <- matrix(c(1,1,1,0,0,0),2,3)[,pancake]
  sample(sides)
}

# sim 10,000 pancakes
pancakes <- replicate( 1e4 , sim_pancake() )
up <- pancakes[1,]
down <- pancakes[2,]

# compute proportion 1/1 (BB) out of all 1/1 and 1/0
num_11_10 <- sum( up==1 )
num_11 <- sum( up==1 & down==1 )
num_11/num_11_10

## R code 15.2
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

# points
plot( d$Divorce ~ d$MedianAgeMarriage , ylim=c(4,15) ,
      xlab="Median age marriage" , ylab="Divorce rate" )

# standard errors
for ( i in 1:nrow(d) ) {
  ci <- d$Divorce[i] + c(-1,1)*d$Divorce.SE[i]
  x <- d$MedianAgeMarriage[i]
  lines( c(x,x) , ci )
}

## R code 15.3
dlist <- list(
  D_obs = standardize( d$Divorce ),
  D_sd = d$Divorce.SE / sd( d$Divorce ),
  M = standardize( d$Marriage ),
  A = standardize( d$MedianAgeMarriage ),
  N = nrow(d)
)

m15.1 <- ulam(
  alist(
    D_obs ~ dnorm( D_true , D_sd ),
    vector[N]:D_true ~ dnorm( mu , sigma ),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ) , data=dlist , chains=4 , cores=4 )

## R code 15.4
precis( m15.1 , depth=2 )

## R code 15.5
dlist <- list(
  D_obs = standardize( d$Divorce ),
  D_sd = d$Divorce.SE / sd( d$Divorce ),
  M_obs = standardize( d$Marriage ),
  M_sd = d$Marriage.SE / sd( d$Marriage ),
  A = standardize( d$MedianAgeMarriage ),
  N = nrow(d)
)

m15.2 <- ulam(
  alist(
    D_obs ~ dnorm( D_est , D_sd ),
    vector[N]:D_est ~ dnorm( mu , sigma ),
    mu <- a + bA*A + bM*M_est[i],
    M_obs ~ dnorm( M_est , M_sd ),
    vector[N]:M_est ~ dnorm( 0 , 1 ),
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp( 1 )
  ) , data=dlist , chains=4 , cores=4 )

## R code 15.6
post <- extract.samples( m15.2 )
D_est <- apply( post$D_est , 2 , mean )
M_est <- apply( post$M_est , 2 , mean )
plot( dlist$M_obs , dlist$D_obs , pch=16 , col=rangi2 ,
      xlab="marriage rate (std)" , ylab="divorce rate (std)" )
points( M_est , D_est )
for ( i in 1:nrow(d) )
  lines( c( dlist$M_obs[i] , M_est[i] ) , c( dlist$D_obs[i] , D_est[i] ) )

## R code 15.7
library(rethinking)
data(milk)
d <- milk
d$neocortex.prop <- d$neocortex.perc / 100
d$logmass <- log(d$mass)

## R code 15.8
dat_list <- list(
  K = standardize( d$kcal.per.g ),
  B = standardize( d$neocortex.prop ),
  M = standardize( d$logmass )
)

m15.3 <- ulam(
  alist(
    K ~ dnorm( mu , sigma ),
    mu <- a + bB*B + bM*M,
    B ~ dnorm( nu , sigma_B ),
    c(a,nu) ~ dnorm( 0 , 0.5 ),
    c(bB,bM) ~ dnorm( 0, 0.5 ),
    sigma_B ~ dexp( 1 ),
    sigma ~ dexp( 1 )
  ) , data=dat_list , chains=4 , cores=4 )

## R code 15.9
precis( m15.3 , depth=2 )

## R code 15.10
obs_idx <- which( !is.na(d$neocortex.prop) )
dat_list_obs <- list(
  K = dat_list$K[obs_idx],
  B = dat_list$B[obs_idx],
  M = dat_list$M[obs_idx]
)
m15.4 <- ulam(
  alist(
    K ~ dnorm( mu , sigma ),
    mu <- a + bB*B + bM*M,
    B ~ dnorm( nu , sigma_B ),
    c(a,nu) ~ dnorm( 0 , 0.5 ),
    c(bB,bM) ~ dnorm( 0, 0.5 ),
    sigma_B ~ dexp( 1 ),
    sigma ~ dexp( 1 )
  ) , data=dat_list_obs , chains=4 , cores=4 )
precis( m15.4 )

## R code 15.11
plot( coeftab(m15.3,m15.4) , pars=c("bB","bM") )

## R code 15.12
post <- extract.samples( m15.3 )
B_impute_mu <- apply( post$B_impute , 2 , mean )
B_impute_ci <- apply( post$B_impute , 2 , PI )

# B vs K
plot( dat_list$B , dat_list$K , pch=16 , col=rangi2 ,
      xlab="neocortex percent (std)" , ylab="kcal milk (std)" )
miss_idx <- which( is.na(dat_list$B) )
Ki <- dat_list$K[miss_idx]
points( B_impute_mu , Ki )
for ( i in 1:12 ) lines( B_impute_ci[,i] , rep(Ki[i],2) )

# M vs B
plot( dat_list$M , dat_list$B , pch=16 , col=rangi2 ,
      ylab="neocortex percent (std)" , xlab="log body mass (std)" )
Mi <- dat_list$M[miss_idx]
points( Mi , B_impute_mu )
for ( i in 1:12 ) lines( rep(Mi[i],2) , B_impute_ci[,i] )

## R code 15.13
m15.5 <- ulam(
  alist(
    # K as function of B and M
    K ~ dnorm( mu , sigma ),
    mu <- a + bB*B_merge + bM*M,
    
    # M and B correlation
    MB ~ multi_normal( c(muM,muB) , Rho_BM , Sigma_BM ),
    matrix[29,2]:MB <<- append_col( M , B_merge ),
    
    # define B_merge as mix of observed and imputed values
    vector[29]:B_merge <- merge_missing( B , B_impute ),
    
    # priors
    c(a,muB,muM) ~ dnorm( 0 , 0.5 ),
    c(bB,bM) ~ dnorm( 0, 0.5 ),
    sigma ~ dexp( 1 ),
    Rho_BM ~ lkj_corr(2),
    Sigma_BM ~ exponential(1)
  ) , data=dat_list , chains=4 , cores=4 )
precis( m15.5 , depth=3 , pars=c("bM","bB","Rho_BM" ) )

## R code 15.14
set.seed(100)
x <- c( rnorm(10) , NA )
y <- c( rnorm(10,x) , 100 )
d <- list(x=x,y=y)
