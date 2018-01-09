#' packages
install.packages("MASS")
install.packages("abind")
install.packages("mnormt")
install.packages("LaplacesDemon")
install.packages("coda")
library(MASS)
library(abind)
library(mnormt)
library(LaplacesDemon)
library(coda)


#' 1. preliminaries: 
#' simulate a gaussian mixture 

d=2
K <- 3
N <- 500
set.seed(1)
p <- c(3/10, 2/10, 5/10) 
NN <- rmultinom(n = 1, size = N, prob = p)  # dimension 3 multinomial
Mu <- rbind(c(-1,-1), c(0,1), c(1, -1))  # the u
X <- matrix(ncol = d, nrow = 0)
Sigma <- array(dim = c(2,2,K))
for(j in 1:K){
    Sigma[,,j] <- rwishart(nu = 5, S = 0.05*diag(d)) 
    # Wishart distribution
}  
for(j in 1:K){
    X <- rbind(X, mvrnorm(n=NN[j], mu = Mu[j,], Sigma=Sigma[,,j])) 
    # multivariate normal distribution
}

#' labs: vector of labels
labs <- rep(0,N)
count=1
for(j in 1:K)
{
    labs[count:(count+NN[j]-1)] <- j
    count=count + NN[j]
}

#' Plot the labelled data 
plot(X[,1], X[,2], col=labs)

#########################'
#' 2. EM  
#########################' 

#' Run EM
Kfit <- 4 ## try with Kfit= 2,3,4,5 ...
outputem <- emalgo(x=X,k=Kfit, tol=1e-6)

#' inspect the objective function (stopping criterion)
length(outputem$objective)

#' Plot the (labelled) data 
plot(X[,1], X[,2], col = labs, pch = 19)

#' Add the starting points (from kmeans) to the plot
Init <-  initem(X,Kfit)
points(Init$Mu[,1],Init$Mu[,2], col="orange",pch=18,cex=5)
#' Add the centers from EM
points(outputem$last$Mu[,1],outputem$last$Mu[,2], col="blue",pch=18,cex=5)
#' Add the true centers
points(Mu[,1],Mu[,2], col="black",pch=3,cex=5) 

#' Draw 1.64 sd level sets 
for(j in 1:Kfit){
    ellips <- draw_sd(outputem$last$Mu[j,], outputem$last$Sigma[,,j])
    lines(ellips[1,], ellips[2,])
}

for(j in 1:K){
    ellips <- draw_sd(Mu[j,], Sigma[,,j])
    lines(ellips[1,], ellips[2,], col='black')
}


#########################'
#' 3. VB
#########################' 
#' Bayesian model: 
#' p ~ dirichlet(alpha);  alpha = (alpha0, ... , alpha0)
#' [ xi | p ] ~ Multinomial(p)
#' [ mu_j | Lambda_j ] ~ Normal(m0, beta0 Lambda_j^(-1))
#' Lambda_j ~ Wishart(W0, nu0)
#' [ X| xi=j, mu, Lambda ] ~ Normal (mu_j, Lambda_j^(-1))


#' hyper-parameters : to be varied 
alpha0 <- 0.1
m0 <- rep(0,2)
beta0 <- 0.1
W0 <- 1*diag(2)
nu0 <- 10
#' Run VB 
#'
seed <- 10
set.seed(seed)
outputvb <- vbalgo(x=X,k=Kfit, alpha0 = alpha0, W0inv = solve(W0),
                 nu0 = nu0, m0 = m0, beta0=beta0, tol=1e-6)

#' plot the lowerbound over iterations 
plot(outputvb$lowerbound)

##' show a summary of VB's output
T <- ncol(outputvb$alphamat)
outputvb$alphamat[,T]
outputvb$Marray[,,T]


##' show a summary of VB's output
T <- ncol(outputvb$alphamat)
outputvb$alphamat[,T]
outputvb$Marray[,,T]

#' Visual summary of VB's output :
#' posterior expectancy of each parameter
p_vb <- outputvb$alphamat[,T] / sum(outputvb$alphamat[,T]) ## complete the code
    ## (variational posterior expectancy of mixture weights)
Mu_vb <- outputvb$Marray[,,T] ## complete the code
    ## (variational posterior expectancy of mixture centers)
Sigma_vb <- array(dim=c(d,d,Kfit))
for(j in 1:Kfit){
    Sigma_vb[,,j] <- 1/outputvb$Numat[,T][j]* (outputvb$Winvarray[,,j,T]) ## complete the code
    ## (variational posterior expectancy of mixture covariances)
}

## show the data, true centers and initial positions from K-means
graphics.off()
plot(X[,1], X[,2], col=labs)
points(Mu[,1],Mu[,2], col="black",pch=8,cex=10*p) 
set.seed(seed)
Init <-  initem(X,Kfit)
points(Init$Mu[,1],Init$Mu[,2], col="orange",pch=18,cex = 10*Init$p)
## Add a  summary of the VB solution
nonneg <- which(p_vb>0.001)
for(j in nonneg){
    points(Mu_vb[j,1], Mu_vb[j,2], col="blue",
           pch=18,cex= 10 * p_vb[j])
    ellips <- draw_sd(mu = Mu_vb[j,], 
                      sigma = Sigma_vb[,,j])
    lines(ellips[1,], ellips[2,], col='blue')
}


####################################################'
####' 4 Metropolis-Hastings
####################################################'
#' Basic testing for the MH sampler
Kmc <- Kfit ## try with different values
init <- initem(x=X, k=Kmc)

hpar <- list( alpha0=rep(alpha0, Kmc),
           m0 = rep(0, d), beta0 = beta0, 
           W0 = W0, nu0 = nu0)

ppar <- list(var_Mu = 0.001,
            nu_Sigma = 500,
            alpha_p = 500) 


set.seed(1)
pct <- proc.time()
outputmh <- MHsample(x=X, k=Kmc, nsample= 3000,
                    init=init, hpar=hpar, ppar=ppar)
newpct <- proc.time()
elapsed <- newpct - pct
elapsed
outputmh$naccept ## should not be ridiculously low. 

################'
#' Heidelberger and Welch's convergence diagnostic
#' 
#'      Stationarity start     p-value
#################' 
y_0 <- mcmc(data=cdfTrace(c(0,0), sample = outputmh))
heidel.diag(y_0)

y_1 <- mcmc(data=cdfTrace(x[100], sample = outputmh))
heidel.diag(y_1)

y_2 <- mcmc(data=cdfTrace(c(1,1), sample = outputmh))
heidel.diag(y_2)

y_3 <- mcmc(data=cdfTrace(c(1,1), sample = outputmh))
heidel.diag(y_3)

##############' 
#' 
#' Gelman and Rubin’s diagnostic
#' 
##############' 
#' the outputmh2
#' 
Kmc2 <- 2 ## try with different values
init2 <- initem(x=X, k=Kmc2)
hpar2 <- list( alpha0=rep(alpha0, Kmc2),
              m0 = rep(0, d), beta0 = beta0, 
              W0 = W0, nu0 = nu0)
set.seed(1)
pct <- proc.time()
outputmh2 <- MHsample(x=X, k=Kmc2, nsample= 3000,
                     init=init2, hpar=hpar2, ppar=ppar)
newpct <- proc.time()
elapsed2 <- newpct - pct
elapsed2
outputmh2$naccept 


#' the outputmh5
#' 
Kmc5 <- 5 ## try with different values
init5 <- initem(x=X, k=Kmc5)
hpar5 <- list( alpha0=rep(alpha0, Kmc5),
              m0 = rep(0, d), beta0 = beta0, 
              W0 = W0, nu0 = nu0)
set.seed(1)
pct <- proc.time()
outputmh5 <- MHsample(x=X, k=Kmc5, nsample= 3000,
                      init=init5, hpar=hpar5, ppar=ppar)
newpct <- proc.time()
elapsed5 <- newpct - pct
elapsed5
outputmh5$naccept

#' the outputmh3, nsample = 3000
#' 
Kmc3 <- 3 ## try with different values
init3 <- initem(x=X, k=Kmc3)
hpar3 <- list( alpha0=rep(alpha0, Kmc3),
               m0 = rep(0, d), beta0 = beta0, 
               W0 = W0, nu0 = nu0)
set.seed(1)
pct <- proc.time()
outputmh3 <- MHsample(x=X, k=Kmc3, nsample= 3000,
                      init=init3, hpar=hpar3, ppar=ppar)
newpct <- proc.time()
elapsed3 <- newpct - pct
elapsed3
outputmh3$naccept
y_4 <- mcmc(data=cdfTrace(c(1,1), sample = outputmh3))
heidel.diag(y_4)

#' the outputmh3, nsample = 10000
#' 
pct <- proc.time()
outputmh3_1 <- MHsample(x=X, k=Kmc3, nsample= 10000,
                      init=init3, hpar=hpar3, ppar=ppar)
newpct <- proc.time()
elapsed3 <- newpct - pct
elapsed3
outputmh3$naccept
y_5 <- mcmc(data=cdfTrace(c(1,1), sample = outputmh3_1))
heidel.diag(y_5)

#####' different test for init value of kfit
mcObject_0 <- mcmc(data=cdfTrace(c(1,1), sample = outputmh5))
mcObject_1 <- mcmc(data=cdfTrace(c(0,1), sample = outputmh5))
mcObject_2 <- mcmc(data=cdfTrace(c(0,0), sample = outputmh5))
mcList <- mcmc.list(mcObject_0, mcObject_1, mcObject_2)
gelman.diag(mcList)
gelman.plot(mcList)

mcObject_0 <- mcmc(data=cdfTrace(c(1,1), sample = outputmh))
mcObject_1 <- mcmc(data=cdfTrace(c(0,1), sample = outputmh))
mcObject_2 <- mcmc(data=cdfTrace(c(0,0), sample = outputmh))
mcList <- mcmc.list(mcObject_0, mcObject_1, mcObject_2)
gelman.diag(mcList)
gelman.plot(mcList)

mcObject_0 <- mcmc(data=cdfTrace(c(1,1), sample = outputmh3))
mcObject_1 <- mcmc(data=cdfTrace(c(0,1), sample = outputmh3))
mcObject_2 <- mcmc(data=cdfTrace(c(0,0), sample = outputmh3))
mcList <- mcmc.list(mcObject_0, mcObject_1, mcObject_2)
gelman.diag(mcList)
gelman.plot(mcList)

mcObject_0 <- mcmc(data=cdfTrace(c(1,1), sample = outputmh3_1))
mcObject_1 <- mcmc(data=cdfTrace(c(0,1), sample = outputmh3_1))
mcObject_2 <- mcmc(data=cdfTrace(c(0,0), sample = outputmh3_1))
mcList <- mcmc.list(mcObject_0, mcObject_1, mcObject_2)
gelman.diag(mcList)
gelman.plot(mcList)

#' Predictive density
xx <- seq(-2,2,length.out=20)
yy <- xx
dtrue <- outer(X= xx, Y=yy,
            FUN = function(x,y){
                wrapper(x=x, y=y,
                        FUN=function(u,v){
                          exp(gmllk(x = c(u,v), Mu = Mu,
                          Sigma = Sigma, p = p))
                        })
            })

dpredmh <-  outer(X= xx, Y=yy,
            FUN = function(x,y){
                wrapper(x = x, y = y,
                        FUN =function(u,v){
                          ## complete the code 
                          MHpredictive(c(u, v), outputmh3, burnin = 2950, thin = 1)
                          })
            })

breaks <- c(seq(0.01,0.09, length.out=5),seq(0.1,0.3,length.out=5))
nbreaks <- length(breaks)
contour(xx,yy, z = dtrue, nlevels=nbreaks, levels = breaks)
contour(xx,yy, z = dpredmh,  nlevels=nbreaks, levels = breaks,
          add=TRUE, col='red')

#' comments for the results 
#' 


#########################'
##' predictive Cdf's
#########################'


Pexcess <- rep(0,10)
Pexcess_em <- Pexcess; Pexcess_vb <- Pexcess; Pexcess_mh <- Pexcess
thres_vect <-  seq(-3, 3, length.out=30)
#seq(1, 5, length.out=30)
Tem <- length(outputem$objective)
T <- ncol(outputvb$alphamat)

for(i in seq_along(thres_vect)){
threshold <- rep(thres_vect[i], 2)
Pexcess[i] <- 1 - gmcdf(x = threshold, Mu = Mu, Sigma=Sigma, p=p)
Pexcess_em[i] <- 1 - gmcdf(x = threshold, Mu = outputem$Muarray[,,Tem],
                           Sigma = outputem$Sigmaarray[,,,Tem], p = outputem$parray[,-1] ) 
                  ## complete the code:
                 ##maximum likelihood estimator using EM output

Pexcess_vb[i] <- 1 - vbPredictiveCdf(x = threshold, 
                                     alpha = outputvb$alphamat[,T], 
                                     Beta = outputvb$Betamat[, T],
                                     M = outputvb$Marray[,,T], 
                                     Winv = outputvb$Winvarray[,,,T],
                                     Nu = outputvb$Numat[,T] )  ## complete the code:
    ## posterior predictive  estimator using VB output: 
    ## use vbPredictiveCdf
    
Pexcess_mh[i] <- 1 - MHpredictiveCdf(x = threshold, sample = outputmh, burnin = 2900, thin = 1) 
    ## complete the code:
    ## posterior predictive  estimator using MH output:
    ## use MHpredictiveCdf. 
}
ylim <- range(Pexcess, Pexcess_em,Pexcess_vb)
plot(thres_vect,Pexcess, ylim = ylim)
lines(thres_vect, Pexcess_vb, col='red')
lines(thres_vect, Pexcess_em, col='blue')
lines(thres_vect, Pexcess_mh, col='green')

#' 
#' tails 
#' 
Pexcess <- rep(0,10)
Pexcess_em <- Pexcess; Pexcess_vb <- Pexcess; Pexcess_mh <- Pexcess
thres_vect <-  seq(1, 5, length.out=30)
Tem <- length(outputem$objective)
T <- ncol(outputvb$alphamat)

for(i in seq_along(thres_vect)){
  threshold <- rep(thres_vect[i], 2)
  Pexcess[i] <- 1 - gmcdf(x = threshold, Mu = Mu, Sigma=Sigma, p=p)
  Pexcess_em[i] <- 1 - gmcdf(x = threshold, Mu = outputem$Muarray[,,Tem],
                             Sigma = outputem$Sigmaarray[,,,Tem], p = outputem$parray[,-1] ) 
  ## complete the code:
  ##maximum likelihood estimator using EM output
  
  Pexcess_vb[i] <- 1 - vbPredictiveCdf(x = threshold, 
                                       alpha = outputvb$alphamat[,T], 
                                       Beta = outputvb$Betamat[, T],
                                       M = outputvb$Marray[,,T], 
                                       Winv = outputvb$Winvarray[,,,T],
                                       Nu = outputvb$Numat[,T] )  ## complete the code:
  ## posterior predictive  estimator using VB output: 
  ## use vbPredictiveCdf
  
  Pexcess_mh[i] <- 1 - MHpredictiveCdf(x = threshold, sample = outputmh, burnin = 2900, thin = 1) 
  ## complete the code:
  ## posterior predictive  estimator using MH output:
  ## use MHpredictiveCdf. 
}
ylim <- range(log(Pexcess), log(Pexcess_em),log(Pexcess_vb))
plot(thres_vect,log(Pexcess), ylim = ylim)
lines(thres_vect, log(Pexcess_vb), col='red')
lines(thres_vect, log(Pexcess_em), col='blue')
lines(thres_vect, log(Pexcess_mh), col='green')
