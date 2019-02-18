library(rethinking)
rm(list=ls())

data(Howell1)

d<-(Howell1)

d2<-d[d$age >= 18,]

flist <- alist(
  height ~ dnorm(mu,sigma) ,
  mu ~ dnorm(178,20) ,
  sigma ~ dunif(0,50)
)

m4.1 <- map(flist, data=d2)

precis(m4.1)

m4.2 <- map2stan(flist, data=d2)

precis(m4.2)

m4.3 <- ulam(flist, data=d2)

precis(m4.3)

start <- list(
  mu=mean(d2$height),
  sigma=sd(d2$height)
)

m4.4 <- map(
  alist(
    height ~ dnorm(mu,sigma),
    mu ~ dnorm(178,0.1),
    sigma ~ dunif(0,50)
  ),
  data=d2 )
precis(m4.4)

vcov(m4.1)

diag(vcov(m4.1))

cov2cor(vcov(m4.1))

post <- extract.samples(m4.1, n=1e4)
head(post)

precis(post)