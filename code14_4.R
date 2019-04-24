## R code 14.46
library(rethinking)
data(Primates301)
data(Primates301_nex)

# plot it
library(ape)
plot( ladderize(Primates301_nex) , type="fan" , font=1 , no.margin=TRUE , label.offset=1 , cex=0.5 )

## R code 14.47
d <- Primates301
d$name <- as.character(d$name)
dstan <- d[ complete.cases( d$group_size , d$body , d$brain ) , ]
spp_obs <- dstan$name

## R code 14.48
dat_list <- list(
  N_spp = nrow(dstan),
  M = standardize(log(dstan$body)),
  B = standardize(log(dstan$brain)),
  G = standardize(log(dstan$group_size)),
  Imat = diag( nrow(dstan) )
)

m14.8 <- ulam(
  alist(
    G ~ multi_normal( mu , SIGMA ),
    mu <- a + bM*M + bB*B,
    matrix[N_spp,N_spp]: SIGMA <- Imat * sigma_sq,
    a ~ normal( 0 , 1 ),
    c(bM,bB) ~ normal( 0 , 0.5 ),
    sigma_sq ~ exponential( 1 )
  ), data=dat_list , chains=4 , cores=4 )
precis( m14.8 )

## R code 14.49
tree_trimmed <- keep.tip(Primates301_nex, spp_obs )
Rbm <- corBrownian( phy=tree_trimmed )
V <- vcv(Rbm)
Dmat <- cophenetic( tree_trimmed )
plot( Dmat , V , xlab="phylogenetic distance" , ylab="covariance" )

## R code 14.50
# put species in right order
dat_list$V <- V[ spp_obs , spp_obs ]
# convert to correlation matrix
dat_list$R <- dat_list$V / max(V)

# Brownian motion model
m14.9 <- ulam(
  alist(
    G ~ multi_normal( mu , SIGMA ),
    mu <- a + bM*M + bB*B,
    matrix[N_spp,N_spp]: SIGMA <- R * sigma_sq,
    a ~ normal( 0 , 1 ),
    c(bM,bB) ~ normal( 0 , 0.5 ),
    sigma_sq ~ exponential( 1 )
  ), data=dat_list , chains=4 , cores=4 )
precis( m14.9 )

## R code 14.51
# add scaled and reordered distance matrix
dat_list$Dmat <- Dmat[ spp_obs , spp_obs ] / max(Dmat)

m14.10 <- ulam(
  alist(
    G ~ multi_normal( mu , SIGMA ),
    mu <- a + bM*M + bB*B,
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat , etasq , rhosq , 0.01 ),
    a ~ normal(0,1),
    c(bM,bB) ~ normal(0,0.5),
    etasq ~ exponential(1),
    rhosq ~ exponential(1)
  ), data=dat_list , chains=4 , cores=4 )
precis( m14.10 )

## R code 14.52
post <- extract.samples(m14.10)
plot( NULL , xlim=c(0,max(dat_list$Dmat)) , ylim=c(0,5) ,
      xlab="phylogenetic distance" , ylab="covariance" )
for ( i in 1:50 ) curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE , col=grau() )

## R code 14.53
S <- matrix( c( sa^2 , sa*sb*rho , sa*sb*rho , sb^2 ) , nrow=2 )

