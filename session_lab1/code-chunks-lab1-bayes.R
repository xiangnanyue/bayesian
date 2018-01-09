############################################
## part 2: bayesian linear regression
############################################

## 2.1 posterior distribution

glinear_fit <- function(Alpha, Beta, data, feature_map, target)
    ## Alpha: prior precision on theta
    ## Beta: noise precision
    ## data: the input variables (x): a matrix with n rows
#### where n is the sample size
    ##feature_map: the basis function, returning a vector of
#### size p equal to the dimension of theta
    ## target: the observed values y: a vector of size n 
{
  Phi <- ## complete:
  ## ? apply
  ## ? t 
    p = ncol(Phi)
    posterior_variance <-  ## complete 
    posterior_mean <-  ## complete 
    return(list(mean=posterior_mean, cov=posterior_variance))
}


## dummy plotting example
xx <- 1:100
yy <- sin(xx/10) * xx/10
pp <- yy + rnorm(n=100)
interv <-  sqrt(abs(yy))
lsup <- yy+ 1.96*interv
linf <- yy- 1.96*interv
plot(xx, yy, type='l', lwd=3, ylim=range(lsup,linf))
lines(xx, lsup, col='red')
lines(xx, linf, col='red')
points(xx, pp, pch=19,col='blue')
legend('top', legend=c('estimate', 'credible levels', 'data'),
       col=c('black', 'red', 'blue'),
       pch= c(NA, NA, NA, 19),
       lwd=c(3,1,1,NA)
       )    


## 2.2 predictive distribution
  glinear_pred <- function(Alpha, Beta, fitted, data, feature_map)
    ## Alpha: prior pecision for theta
    ## Beta: noise variance
    ## fitted: the output of glinear_fit: the posterior mean and
#### variance of the parameter theta.
    ## data: new input data where the predictive distribution 
#### of Y must be computed  
    ## feature map: the vector of basis functions 
{
    Phi_transpose <-  ## complete 
    pred_mean <- ## complete 
    pred_variance <- ## complete 
    return(list(mean = pred_mean, variance = pred_variance))
    
}


## 2.3 Empirical bayes

logevidence <- function(Alpha, Beta, data ,feature_map, target)
    ## Alpha: prior precision for theta
    ## Beta: noise precision
    ## data: the input points x_{1:n}
    ## feature_map: the vector of basis functions
    ## target: the observed values y: a vector of size n.  
{
    Phi_transpose <- ## complete  the code
    if(is.vector(Phi_transpose)){
        Phi_transpose = matrix(Phi_transpose,nrow=1)
      }
      ## avoids undesired matrix-> vector conversions for one
      ## dimensional feature maps
    Phi <- t(Phi_transpose)
    N <- nrow(Phi)
    p <- ncol(Phi)
    
    ### complete the code 
   
    return(res)    
}



F7 <- function(x){c(1, x, x^2, x^3, x^4,x^5, x^6, x^7)}
F6 <- function(x){c(1, x, x^2, x^3, x^4,x^5, x^6)}
F5 <- function(x){c(1, x, x^2, x^3, x^4,x^5)}
F4 <- function(x){c(1, x, x^2, x^3, x^4)}
F3 <- function(x){c(1, x, x^2, x^3)}
F2 <- function(x){c(1, x, x^2)}
F1 <- function(x){c(1, x)}
F0 <- function(x){1}
listF=list(F0,F1,F2,F3,F4,F5,F6,F7)

############################################
## part 3: Naive Bayes
############################################


## 3.1 data pre-processing
 install.packages('RCurl')
library(RCurl)
urlfile <- 'https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data'

 # download the file
downloaded <- getURL(urlfile, ssl.verifypeer=FALSE)
# treat the text data as a steam so we can read from it
connection <- textConnection(downloaded)
# parse the downloaded data as CSV
dataset <- read.csv(connection, header=FALSE, as.is=FALSE)
adult <- dataset 
indx <- sapply(adult, is.factor)

N = dim(adult)[1]
p = dim(adult)[2] - 1

Ntrain = floor(2*N/3)
training = adult[1:Ntrain,]
test = adult[(Ntrain+1):N,]

Class = training[,p+1]
Features = training[,1:p]


## 3.2 defining a prior distribution for each class using the dataset

Nclass = length(levels(Class))
classPrior=rep(0,Nclass)
for (i in 1:Nclass){
    classPrior[i] = ##complete 
}
classPrior


## 3 training the model


  train_naivebayes <- function(Class, predictors, laplace = 1)
  ## Class: vector of target classes for the classification problem
  ###### (a factor) 
  ## predictors: the features (x) of the training data
  ## laplace: an integer(0 or 1): 1 for laplace regularization,
  ###### i.e. for assigning an additional observation to each class
  ###### (avoids zeros)
{
    Nclass = length(levels(Class))
    p = dim(predictors)[2]
    indClass=list()
    ## select indices corresponding to each class k
    for(k in 1:Nclass){
        indClass[[k]]= which(Class==levels(Class)[k])
    }
    fitted=list()
    for (j in 1:p){
        if(is.factor(predictors[,j])){
            ## fit a multinomial model:
            ## store the estimated probability for a class given 
            ## a  predictor in a  matrix (nrows: number of levels,
            ## ncol: number of classes)
            nlev_j = length(levels(predictors[,j]))
            fitted[[j]]=matrix( nrow= nlev_j, ncol=Nclass)
        }
        else{ ##fit a normal model: store estimated mean and variance
            fitted[[j]]=matrix(0,nrow = 2, ncol = Nclass)
        }
        for(k in 1:Nclass){
            X = predictors[indClass[[k]], j]
            if( is.factor(X)){
             
                 fitted[[j]][,k] = ##  complete 
            }
            else{
                fitted[[j]][,k] = ## complete 
            }
        }
    }
    return(fitted)
}

## 3.4 class prediction 
  predict_naivebayes <- function(model, xnew, classPrior = NULL,
                                   labels = NULL)
{
    ## model: the output of function train_naivebayes
    ## xnew: the new input points (a matrix with as any rows
#### as test data)
    ## classPrior: the prior probabilitie of each class.
#### if NULL, a uniform prior will be imposed.
    ## labels: the target class labels. If NULL, the class labels
#### contained in the 'model' agument will be used. 
                                     
    p <- ncol(xnew)
    if(is.null(labels)){
        nclass <- ncol(model[[1]])
        labels <- as.character(1:nclass)
    }
    else{nclass<- length(labels)}
    if(is.null(classPrior)){
        classPrior <- rep(1/nclass, nclass)
      }
      
    ntest <- nrow(xnew)
    posteriorProb<- matrix(0, nrow = ntest, ncol = nclass,)
    for( k in 1:nclass){
        marglikelihoods= matrix(1,nrow= ntest, ncol=p)
        for(j in 1:p){
            if(is.factor(xnew[,j])){
                inds <- sapply(xnew[,j],
                              function(x){which(x == levels(x))})
                marglikelihoods[,j] <- ## complete 
            }
            else{
              marglikelihoods[,j] <-
              ## complete 
            }
            
        }
        posteriorProb[,k] <- classPrior[k] *
            exp(apply(log(marglikelihoods), 1, sum ))
        
    }
    normalize = apply(posteriorProb,1,sum)
    posteriorProb = posteriorProb/normalize
    ## M/v: each row M[i,] is divided by v[i] 
    return(posteriorProb)
}  


## 3.5 testing your model

fitted_naive <- train_naivebayes(Class=Class,
                                 predictors = training[,-15],
                                 laplace=1)
Pclass <-  predict_naivebayes(fitted_naive,xnew = test[,-15],
                              classPrior = classPrior,
                              labels = levels(Class))
mapClass = apply(Pclass, 1, which.max)
mapClass = factor(mapClass, labels=c(" <=50K", " >50K"))
table(mapClass,testClass)


## Compare with the output of the e1071} package:

library(e1071)
m <- naiveBayes(x = training[, 1:14], y=training[ , 15],laplace =1)
predm <- predict(m, test, type ="class")
probm <- predict(m, test, type ="raw")
table(predm, testClass)
