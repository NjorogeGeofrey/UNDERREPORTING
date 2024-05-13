

####################### SIMULATION PART ######################




# SIMULATED DATA
lambda <- 10
P <- 0.7
yp <- rpois(47, lambda * P)

## LIKELIHOOD 
L=function(data,P,lambda){
  n=length(data)
  PP1=((P*lambda)^sum(data))*exp(-(n*P*lambda))/(prod(factorial(data)))
  return(PP1)
}


# FULLP FUNCTION; computing the full probability based on data, P, lambda, a, and b
FULLP <- function(data, P, lambda, a, b) {
  n <- length(data)
  L1 <- (P^(sum(data) + a - 1)) * (1 - P)^(b - 1)  * exp(-n * lambda * P)
  return(L1)
}

# METROPOLIS FUNCTION;  calculating the acceptance ratio for a Metropolis-Hastings algorithm
MHA <- function(data, PP, P, lambda, a, b){
  R1 <- FULLP(data, PP, lambda, a, b) / FULLP(data, P, lambda, a, b)
  return(R1)
}

# TRIAL FUNCTION; implementing a Markov chain Monte Carlo (MCMC) 
# algorithm using the Metropolis-Hastings method
TRIAL <- function(data, alpha, Beta, P1, lambda1, a, b,L) {
  n <- length(data)
  lambda_chain <- numeric(L)
  P_chain <- numeric(L)
  P_chain[1]=P1
  lambda_chain[1]=lambda1
  for (i in 2:L) {
    lambda_chain[i] <- rgamma(1, (sum(data) + alpha), (Beta + (n * P_chain[i -1])))
    PP <- runif(1, 0.6,.9)
    U <- runif(1, 0, 1)
    R <- MHA(data, PP, P_chain[i - 1], lambda_chain[i], a, b)
    if (R < U | is.nan(R)==T) {
      P_chain[i] <- PP
    } else {
      P_chain[i] <- P_chain[i - 1]
    }
  }
  return(list(lambda = lambda_chain, P = P_chain))
}



# Call the TRIAL function; Calls the TRIAL function with 
# the simulated data yp and other parameters, running the MCMC algorithm for 10000 iterations
result <- TRIAL(yp, 0.0001, 0.00001, .1, 100, 0.00001, 0.0001,10000)

# Access lambda and P chains; Extracting the chains of sampled values of lambda and 
# P from the result of the MCMC algorithm
lambda_chain <- result$lambda
P_chain <- result$P
length(lambda_chain)

#Plot; Plotting the chains of lambda and P against their iteration numbers
plot(lambda_chain,type='l')
plot(P_chain,type = 'l')

# Creating sequences of values for P and lambda.
p1=seq(0.2,0.9,.01)
l1=seq(6,13,.1)
length(p1)
length(l1)


I=length(p1)
J=length(l1)

## Initializing a matrix LLL filled with NA values, 
## with dimensions determined by the lengths of p1 and l1
LLL=matrix(rep(NA,I*J),ncol =J)

# Nested loops iterating over each combination of P and lambda values.
for(i in 1:I ){
  for( j in 1:J){
    LLL[i,j]=log(L(yp,p1[j],l1[i]),base=exp(1)) ###Calculating the log likelihood for each 
    ## combination of P and lambda and storing it in the matrix LLL
  }
}


# Identifying the largest values of lambda and P; Finding the indices of the maximum 
# value in the likelihood matrix LLL
max_index <- which(LLL == max(LLL), arr.ind = TRUE)





mean(lambda_chain)
mean(P_chain)
unique(max_index[,1])
sort(unique(max_index[,1]))
p1[sort(unique(max_index[,1]))]

#Plot; Plotting the chains of lambda and P against their iteration numbers
plot(lambda_chain,type='l')
plot(P_chain,type = 'l')
mean(P_chain[-c(1:100)])

quantile(P_chain[-c(1:100)],c(.025,.975))
mean(lambda_chain[-c(1:100)])

acf(P_chain[-c(1:100)])
acf(lambda_chain[-c(1:100)])

unique(p1[max_index[,1]])
sort(unique(p1[max_index[,1]]),decreasing = T)









################ WORKING WITH REAL DATA ####################









#importing the real data
library(readr)
DATA <- read.csv("DATA.csv")
head(DATA)
dim(DATA)
CLEAN_DATA <- DATA[-48, ]
dim(CLEAN_DATA)
#DIVIDING THE DATA IN TO YEARS
DATASET14 <- subset(CLEAN_DATA, Year == "2014")
DATA14 <- DATASET14$People_Living_with_HIV

DATASET15 <- subset(CLEAN_DATA, Year == "2015")
DATA15 <- DATASET15$People_Living_with_HIV

DATASET17 <- subset(CLEAN_DATA, Year == "2017")
DATA17 <- DATASET17$People_Living_with_HIV



## LIKELIHOOD 
L=function(data,P,lambda){
  n=length(data)
  PP1=((P*lambda)^sum(data))*exp(-(n*P*lambda))/(prod(factorial(data)))
  return(PP1)
}



# FULLP FUNCTION; computing the full probability based on data, P, lambda, a, and b
FULLP <- function(data, P, lambda, a, b) {
  n <- length(data)
  L1 <- (P^(sum(data) + a - 1)) * (1 - P)^(b - 1)  * exp(-n * lambda * P)
  return(L1)
}

# METROPOLIS FUNCTION;  calculating the acceptance ratio for a Metropolis-Hastings algorithm
MHA <- function(data, PP, P, lambda, a, b){
  R1 <- FULLP(data, PP, lambda, a, b) / FULLP(data, P, lambda, a, b)
  return(R1)
}

##fINDING THE LIMITS OF P
L1 <- function(params, dta) {
  P <- params[1]
  lambda <- params[2]
  ##Scaling the data
  dta = (DATA14/5000)
  j = -log(L(dta, P, lambda))
  return(j)
}


optim(par = c(0.7, 11), fn = L1, method = "L-BFGS-B", 
      lower = c(0.35, 6), upper = c(0.9, 20), hessian = TRUE)

# TRIAL FUNCTION; implementing a Markov chain Monte Carlo (MCMC) 
# algorithm using the Metropolis-Hastings method
TRIAL <- function(data, alpha, Beta, P1, lambda1, a, b,L) {
  n <- length(data)
  lambda_chain <- numeric(L)
  P_chain <- numeric(L)
  P_chain[1]=P1
  lambda_chain[1]=lambda1
  for (i in 2:L) {
    lambda_chain[i] <- rgamma(1, (sum(data) + alpha), (Beta + (n * P_chain[i -1])))
    PP <- runif(1, 0.5,0.7)
    U <- runif(1, 0, 1)
    R <- MHA(data, PP, P_chain[i - 1], lambda_chain[i], a, b)
    if (R < U | is.nan(R)==T) {
      P_chain[i] <- PP
    } else {
      P_chain[i] <- P_chain[i - 1]
    }
  }
  return(list(lambda = lambda_chain, P = P_chain))
}



# Call the TRIAL function; Calls the TRIAL function with 
# the simulated data yp and other parameters, running the MCMC algorithm for 10000 iterations
result <- TRIAL(DATA14, 0.0001, 0.00001, 0.2, 20, 0.00001, 0.0001,10000)

# Access lambda and P chains; Extracting the chains of sampled values of lambda and 
# P from the result of the MCMC algorithm
lambda_chain <- result$lambda
P_chain <- result$P
length(lambda_chain)

#Plot; Plotting the chains of lambda and P against their iteration numbers
plot(lambda_chain[-c(1:100)], type = 'l', ylim = c(2000, 80000))
plot(P_chain[-c(1:100)],type = 'l',ylim = c(0, .9))

hist(P_chain[-c(1:100)])
hist(lambda_chain[-c(1:100)])


quantile(lambda_chain[-c(1:100)],c(.025,.975))
mean(lambda_chain[-c(1:100)])
median(lambda_chain[-c(1:100)])
quantile(P_chain[-c(1:100)],c(.025,.975))
mean(P_chain[-c(1:100)])
median(P_chain[-c(1:100)])

acf(P_chain[-c(1:100)])
acf(lambda_chain[-c(1:100)])

### Estimation
median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)])
median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)])-median(DATA14)

##Visual 
plot(P_chain[-c(1:100)]*lambda_chain[-c(1:100)], type = "l")
abline(h=median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)]), col='red')

## Proportion of underreported counts
1-median(DATA14)/median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)])


median(DATA14)
hist(DATA14)


######################### 2015
##fINDING THE LIMITS OF P
##Scaling the data
dta2=DATA15/5000
L2 <- function(params, dta2) {
  P <- params[1]
  lambda <- params[2]
  dta2 = (DATA15/5000)
  j = -log(L(dta2, P, lambda))
  return(j)
}
dta2
mean(dta2)/.5

optim(par = c(0.7, 13), fn = L2, method = "L-BFGS-B", 
      lower = c(0.4, 8), upper = c(0.9, 18), hessian = TRUE)

# TRIAL FUNCTION; implementing a Markov chain Monte Carlo (MCMC) 
# algorithm using the Metropolis-Hastings method
TRIAL <- function(data, alpha, Beta, P1, lambda1, a, b,L) {
  n <- length(data)
  lambda_chain <- numeric(L)
  P_chain <- numeric(L)
  P_chain[1]=P1
  lambda_chain[1]=lambda1
  for (i in 2:L) {
    lambda_chain[i] <- rgamma(1, (sum(data) + alpha), (Beta + (n * P_chain[i -1])))
    PP <- runif(1, 0.5,0.7)
    U <- runif(1, 0, 1)
    R <- MHA(data, PP, P_chain[i - 1], lambda_chain[i], a, b)
    if (R < U | is.nan(R)==T) {
      P_chain[i] <- PP
    } else {
      P_chain[i] <- P_chain[i - 1]
    }
  }
  return(list(lambda = lambda_chain, P = P_chain))
}



# Call the TRIAL function; Calls the TRIAL function with 
# the simulated data yp and other parameters, running the MCMC algorithm for 10000 iterations
result <- TRIAL(DATA14, 0.0001, 0.00001, 0.2, 20, 0.00001, 0.0001,10000)

# Access lambda and P chains; Extracting the chains of sampled values of lambda and 
# P from the result of the MCMC algorithm
lambda_chain <- result$lambda
P_chain <- result$P
length(lambda_chain)

#Plot; Plotting the chains of lambda and P against their iteration numbers
plot(lambda_chain[-c(1:100)], type = 'l', ylim = c(2000, 80000))
plot(P_chain[-c(1:100)],type = 'l',ylim = c(0, .9))

hist(P_chain[-c(1:100)])
hist(lambda_chain[-c(1:100)])


quantile(lambda_chain[-c(1:100)],c(.025,.975))
mean(lambda_chain[-c(1:100)])
median(lambda_chain[-c(1:100)])
quantile(P_chain[-c(1:100)],c(.025,.975))
mean(P_chain[-c(1:100)])
median(P_chain[-c(1:100)])

acf(P_chain[-c(1:100)])
acf(lambda_chain[-c(1:100)])

### Estimation
median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)])
median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)])-median(DATA15)

##Visual 
plot(P_chain[-c(1:100)]*lambda_chain[-c(1:100)], type = "l")
abline(h=median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)]), col='red')

## Proportion of underreported counts
1-median(DATA14)/median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)])


median(DATA14)
hist(DATA14)

###################DATA 3
######################### 2017
##fINDING THE LIMITS OF P
##Scaling the data
dta3 <- DATA17/5000
L3 <- function(params, dta3) {
  P <- params[1]
  lambda <- params[2]
  dta3 = (DATA17/5000)
  j = -log(L(dta3, P, lambda))
  return(j)
}

mean(dta3)/.5

optim(par = c(0.5, 13), fn = L3, method = "L-BFGS-B", 
      lower = c(0.4, 8), upper = c(0.9, 18), hessian = TRUE)

# TRIAL FUNCTION; implementing a Markov chain Monte Carlo (MCMC) 
# algorithm using the Metropolis-Hastings method
TRIAL <- function(data, alpha, Beta, P1, lambda1, a, b,L) {
  n <- length(data)
  lambda_chain <- numeric(L)
  P_chain <- numeric(L)
  P_chain[1]=P1
  lambda_chain[1]=lambda1
  for (i in 2:L) {
    lambda_chain[i] <- rgamma(1, (sum(data) + alpha), (Beta + (n * P_chain[i -1])))
    PP <- runif(1, 0.4,0.6)
    U <- runif(1, 0, 1)
    R <- MHA(data, PP, P_chain[i - 1], lambda_chain[i], a, b)
    if (R < U | is.nan(R)==T) {
      P_chain[i] <- PP
    } else {
      P_chain[i] <- P_chain[i - 1]
    }
  }
  return(list(lambda = lambda_chain, P = P_chain))
}



# Call the TRIAL function; Calls the TRIAL function with 
# the simulated data yp and other parameters, running the MCMC algorithm for 10000 iterations
result <- TRIAL(DATA17, 0.0001, 0.00001, 0.2, 20, 0.00001, 0.0001,10000)

# Access lambda and P chains; Extracting the chains of sampled values of lambda and 
# P from the result of the MCMC algorithm
lambda_chain <- result$lambda
P_chain <- result$P
length(lambda_chain)

#Plot; Plotting the chains of lambda and P against their iteration numbers
plot(lambda_chain[-c(1:100)], type = 'l', ylim = c(2000, 80000))
plot(P_chain[-c(1:100)],type = 'l',ylim = c(0, .9))

hist(P_chain[-c(1:100)])
hist(lambda_chain[-c(1:100)])


quantile(lambda_chain[-c(1:100)],c(.025,.975))
mean(lambda_chain[-c(1:100)])
median(lambda_chain[-c(1:100)])
quantile(P_chain[-c(1:100)],c(.025,.975))
mean(P_chain[-c(1:100)])
median(P_chain[-c(1:100)])

acf(P_chain[-c(1:100)])
acf(lambda_chain[-c(1:100)])

### Estimation
median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)])
median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)])-median(DATA17)

##Visual 
plot(P_chain[-c(1:100)]*lambda_chain[-c(1:100)], type = "l")
abline(h=median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)]), col='red')

## Proportion of under reported counts
1-median(DATA17)/median(P_chain[-c(1:100)]*lambda_chain[-c(1:100)])


median(DATA17)
hist(DATA17)







