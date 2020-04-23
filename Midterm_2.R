library(reshape2)
library(ggplot2)
library(dplyr)
library(rjags)

#read data
df_men <- read.csv("C:\\Users\\ashwi\\Downloads\\ST 540\\men_sampled_shoe.csv")

#omits na values
newdata_men <- na.omit(df_men)

#assign id's for each runner
newdata_men <- transform(newdata_men,nameid=as.numeric(factor(match_name)))
newdata_men <- transform(newdata_men,marathonid=as.numeric(factor(marathon)))
#add binary vaporfly yes or no
vaporbinary <- newdata_men$vaporfly*1
newdata_men <- cbind(newdata_men, vaporbinary)

#MODEL 1 (Varying intercepts and slopes)
#start mcmc process
Y <- newdata_men$time_minutes
v <- newdata_men$vaporbinary
m_id<-newdata_men$marathonid
n_id<-newdata_men$nameid
X <- cbind(v,m_id,n_id)
M <- 22
J <- 308

zero.u<-numeric(2)

mean(Y)
#model Varying intercepts and slopes
data     <- list(Y=Y,n=840, v=v, J=J, M=M,
                 m_id=m_id, n_id=n_id, zero.u=zero.u)
burn     <- 10000
n.iter   <- 30000
thin     <- 20
n.chains <- 2
model_string <- textConnection("model{
   # Likelihood
    for(i in 1:n){
      mu[i] <- (beta[1]+w[m_id[i],1]+u[n_id[i],1])
      +(beta[2]+w[m_id[i],2]+u[n_id[i],2])*v[i]
  
      Y[i] ~ dnorm(mu[i], tau.e)
    }

   # Random effects
    for(k in 1:M){
      w[k,1:2] ~ dmnorm(zero.u,invSig.w)
    }
    for(j in 1:J){
      u[j,1:2] ~ dmnorm(zero.u,invSig.u)
    }
    
   # priors
   beta[1]   ~ dnorm(140,0.001)
   beta[2]   ~ dnorm(0,0.001)
   
   #error variance
   Sigma_e ~ dgamma(20,1)
   tau.e <- pow(Sigma_e,-2)
   
   #for individual random effects
   tau.u1    ~ dgamma(1,.01)
   tau.u2    ~ dgamma(1,.01)
   sigma.u1 <- pow(tau.u1,-1/2)
   sigma.u2 <- pow(tau.u2,-1/2)
   
   R.u[1,1] <- pow(sigma.u1,2)
   R.u[2,2] <- pow(sigma.u2,2)
   R.u[1,2] <- rho.u*sigma.u1*sigma.u2
   R.u[2,1] <- rho.u*sigma.u1*sigma.u2
   invSig.u ~ dwish(R.u,2.1)
   Sigma.u <- inverse(invSig.u)
   
   rho.u     ~ dnorm(mu_rho.u,tau_rho.u)T(-1,1)
   mu_rho.u  ~ dunif(-1,1)
   tau_rho.u ~ dgamma(.1,.0001)
   
   #for marathon location random effects
   tau.w1    ~ dgamma(.1,.1)
   tau.w2    ~ dgamma(.1,.1)
   sigma.w1 <- pow(tau.w1,-1/2)
   sigma.w2 <- pow(tau.w2,-1/2)
   
   R.w[1,1] <- pow(sigma.w1,2)
   R.w[2,2] <- pow(sigma.w2,2)
   R.w[1,2] <- rho.w*sigma.w1*sigma.u2
   R.w[2,1] <- rho.w*sigma.w1*sigma.u2
   invSig.w ~ dwish(R.w,2.1)
   Sigma.w <- inverse(invSig.w)
   
   rho.w ~ dnorm(mu_rho.w,tau_rho.w)T(-1,1)
   mu_rho.w ~ dunif(-1,1)
   tau_rho.w ~ dgamma(1,.0001)
     }")
params  <- c("beta","Sigma.u","Sigma.w","Sigma_e", "rho.w", "rho.u")
model   <- jags.model(model_string,data = data, n.chains=n.chains,quiet=TRUE)
update(model, burn, progress.bar="none")
samples1 <- coda.samples(model, variable.names=params,
                        n.iter=n.iter, thin=thin, progress.bar="none")
plot(samples1)
summary(samples1)

samples1<-rbind(samples1[[1]],samples1[[2]])

DIC1 <- dic.samples(model, 
                    variable.names=c("beta"), 
                    n.iter=30000, progress.bar="none")
DIC1

jags_out <- rjags::coda.samples(model,
                                variable.names=c("beta","Sigma.u","Sigma.w","Sigma_e", "rho.w", "rho.u"),
                                n.iter=30000, progress.bar="none")
install.packages("MCMCvis")
library(MCMCvis)
MCMCsummary(jags_out, round=2)
gelman.diag(samples1)
#Model 2 (treat all subjects as having a fixed mean)

data     <- list(Y=Y,n=840, v=v, J=J, M=M,
                 m_id=m_id, n_id=n_id)
burn     <- 10000
n.iter   <- 30000
thin     <- 20
n.chains <- 2
model_string <- textConnection("model{
   # Likelihood
    for(i in 1:n){
      mu[i] <- beta[1]+(beta[2]
               +w[m_id[i]]+u[n_id[i]])*v[i]
  
      Y[i] ~ dnorm(mu[i], tau.e)
    }

   # Random effects
    for(k in 1:M){
      w[k] ~ dnorm(0,tau.w)
    }
    for(j in 1:J){
      u[j] ~ dnorm(0,tau.u)
    }
    
   # priors
   beta[1]   ~ dnorm(140,0.001)
   beta[2]   ~ dnorm(0,0.001)
   
   #error variance
   Sigma_e ~ dgamma(20,1)
   tau.e <- pow(Sigma_e,-2)
   
   #for individual random effect
   tau.u    ~ dgamma(1,.001)
   sigma.u <- pow(tau.u,-1/2)
   
  
   #for marathon location random effect
   tau.w    ~ dgamma(1,.001)
   sigma.w <- pow(tau.w,-1/2)
   
     }")
params  <- c("beta","sigma.u","sigma.w", "Sigma_e")
model   <- jags.model(model_string,data = data, n.chains=n.chains,quiet=TRUE)
update(model, burn, progress.bar="none")
samples2 <- coda.samples(model, variable.names=params,
                        n.iter=n.iter, thin=thin, progress.bar="none")
plot(samples2)
summary(samples2)

DIC2 <- dic.samples(model, 
                    variable.names=c("beta"), 
                    n.iter=30000, progress.bar="none")
DIC2


#Model 3 Fixed Effects


data     <- list(Y=Y,n=840, v=v)
burn     <- 10000
n.iter   <- 30000
thin     <- 20
n.chains <- 2
model_string <- textConnection("model{
   # Likelihood
    for(i in 1:n){
      mu[i] <- beta[1]+(beta[2]*v[i])
  
      Y[i] ~ dnorm(mu[i], tau.e)
    }

   # priors
   beta[1]   ~ dnorm(140,0.001)
   beta[2]   ~ dnorm(0,0.001)
   
   #error variance
   Sigma_e ~ dgamma(20,1)
   tau.e <- pow(Sigma_e,-2)
     }")
params  <- c("beta","Sigma_e")
model   <- jags.model(model_string,data = data, n.chains=n.chains,quiet=TRUE)
update(model, burn, progress.bar="none")
samples3 <- coda.samples(model, variable.names=params,
                         n.iter=n.iter, thin=thin, progress.bar="none")
plot(samples3)
summary(samples3)

DIC3 <- dic.samples(model, 
                    variable.names=c("beta"), 
                    n.iter=30000, progress.bar="none")
DIC3
gelman.diag(samples3)

df_women <- read.csv("C:\\Users\\ashwi\\Downloads\\ST 540\\women_sampled_shoe.csv")


DIC1
DIC2
DIC3

