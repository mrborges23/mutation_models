# ESTIMATORS OF POPULATION SIZE WITH BOUNDARY AND RECURRENT MUTATION MODELS
# Author: Rui Borges

# Sampled frequency spectrum of the recurrent mutation model
#  mu : mutation rate
#  N  : effective population size
#  M  : number of sampled individuals

sampled_recurrent_sfs <- function(mu,N,M){
  
  # creates an empty vector
  ssfs <- rep(NA,M+1)
  
  # calculates the frequencies of the monomorphic counts
  # this step guarantees that we do not encounter underflow errors
  if (abs(N*mu-1) < 0.00001){
    ssfs[1]   <- (N-M)/(M+1)
    ssfs[M+1] <- ssfs[1]
  } else {
    ssfs[1]   <- exp(lgamma(1-N-N*mu)+lgamma(1-2*N*mu-M)-lgamma(1-N-2*N*mu)-lgamma(1-N*mu-M))
    ssfs[M+1] <- ssfs[1]
  }
  
  # calculates the frequencies of the polymorphic counts
  m <- 1:(M-1)
  ssfs[m+1] <- exp(log(N) + log(mu) + lchoose(M,m)+lgamma(N*mu+m)+lgamma(2*N*mu+N)+lgamma(N*mu-m+M)-lgamma(N*mu+1)-lgamma(2*N*mu+M)-lgamma(N*mu+N))

  # normalizes and returns the sampled SFS
  return(ssfs/sum(ssfs))
  
}


# Sampled frequency spectrum of the boundary mutation model
#  mu : mutation rate
#  N  : effective population size
#  M  : number of sampled individuals

sampled_boundary_sfs <- function(mu,N,M){
  
  # creates an empty vector
  ssfs <- rep(NA,M+1)
  
  # calculates some harmonic numbers
  hn <- digamma(N)-digamma(1)
  hm <- digamma(M)-digamma(1)
    
  # calculates the frequencies of the monomorphic counts
  ssfs[1]   <- 1+mu*N*(hn-hm)
  ssfs[M+1] <- 1+mu*N*(hn-hm)
  
  # calculates the frequencies of the polymorphic counts
  m <- 1:(M-1)
  ssfs[m+1] <- mu*N*M/(m*(M-m))
  
  # normalizes and returns the sampled SFS
  return(ssfs/sum(ssfs))
}


# Likelihood of the boundary mutation model
#  counts : the empirical sampled site frequency spectrum
#  lN     : log of the effective population size
#  mu     : mutation rate

llik_sampled_boundary_sfs <- function(counts,lN,mu){
  
  # calculates the sampled frequency spectrum of the boundary mutation model assuming that
  # a sample of M individuals was taken from a population of size N
  M <- length(counts)-1
  N <- round(exp(lN))
  expected_counts <- sampled_boundary_sfs(mu,N,M)
  
  # calculates the log ikelihood
  llik <- sum(counts*log(expected_counts))
  
  # returns the likelihood
  return(llik)
  
}

# Likelihood of the recurrent mutation model
#  counts : the empirical sampled site frequency spectrum
#  lN     : log of the effective population size
#  mu     : mutation rate

llik_sampled_recurrent_sfs <- function(counts,lN,mu){
  
  # calculates the sampled frequency spectrum of the recurrent mutation model assuming that
  # a sample of M individuals was taken from a population of size N
  M <- length(counts)-1
  N <- round(exp(lN))
  expected_counts <- sampled_recurrent_sfs(mu,N,M)
  
  # calculates the log ikelihood
  llik <- sum(counts*log(expected_counts))
  
  # returns the likelihood
  return(llik)
  
}



# Metropolis-Hastings on the effective population size under the boundary mutation model
#  mcmc : list of objects (this list is automatically created by the mcmc functions)

mh_boundary_N <- function(mcmc){
  
  # samples from a normal proposal
  # while avoiding too small or big population sizes
  lN1 <- rnorm(1,mcmc$lN,mcmc$tunning)
  if (lN1<2 || lN1>120){
    return(mcmc)
  }
  
  # calculates the likelihood under the boundary mutation model
  lk1 <- llik_sampled_boundary_sfs(mcmc$counts,lN1,mcmc$mu)
  
  # Accept-reject step
  # in case of acceptance the newly proposed population size and the likelihood are updated
  alpha <- exp(lk1-mcmc$lk)
  if (alpha > runif(1,0,1)){
    mcmc$lN <- lN1
    mcmc$lk <- lk1
    mcmc$ap <- mcmc$ap +1 
  }
  
  # returns the updated mcmc list
  return(mcmc)
  
}



# Metropolis-Hastings on the effective population size under the boundary mutation model
#  mcmc : list of objects (this list is automatically created by the mcmc functions)

mh_recurrent_N <- function(mcmc){
  
  # samples from a normal proposal
  # while avoiding too small or big population sizes
  lN1 <- rnorm(1,mcmc$lN,mcmc$tunning)
  if (lN1<2 || lN1>100){
    return(mcmc)
  }
  
  # calculates the likelihood under the recurrent mutation model
  lk1 <- llik_sampled_recurrent_sfs(mcmc$counts,lN1,mcmc$mu)
  
  # Accept-reject step
  # in case of acceptance the newly proposed population size and the likelihood are updated
  alpha <- exp(lk1-mcmc$lk)
  if (alpha > runif(1,0,1)){
    mcmc$lN <- lN1
    mcmc$lk <- lk1
    mcmc$ap <- mcmc$ap +1 
  }
  
  # returns the updated mcmc list
  return(mcmc)
  
}


# Markov Chain Monte Carlo simulations to estimate the population size under the recurrent mutation model
#  I      : number of MCMC iterates (should be a multiple of 100 and higher than 1000; only every 100th iterate are kept!)
#  counts : the empirical sampled site frequency spectrum
#  mu     : mutation rate

mcmc_recurrent <- function(I,counts,mu){
  
  # sets the acceptance probability to 0 and calculates the likelihood for lN
  lN     <- -log(mu)
  ap     <- 0
  lk <- llik_sampled_recurrent_sfs(counts,lN,mu)
  
  # creates a list that includes the mutation rate, the current population size, likelihood 
  # proposal tunning and acceptance probability
  mcmc <- list(lN=lN,mu=mu,counts=counts,lk=lk,tunning=1,ap=ap)
  
  # vector with the sampled values of N
  mcmclN <- rep(NA,(I-1000)/100)
  
  # runs the MCMC
  p <- 1
  for (i in 1:I){
    
    # Metropolis-Hastings step with the recurrent mutation model
    mcmc <- mh_recurrent_N(mcmc)
    
    # saves each 100th iterate
    # forces a burning phase of 1000 iterates
    if ((i%%100 ==0) && i>1000 ){
      
      # saves the effective population size
      mcmclN[p] <- mcmc$lN
      p         <- p+1
      
      # calculates the acceptance probability
      aprob <- mcmc$ap/i
      
      # retunes the proposal variance according to the acceptance probability so that
      # the acceptance probability revolves around 0.345
      mcmc$tunning <- mcmc$tunning*(-0.8312691*aprob*aprob+2.331269*aprob+0.5)

    }
  }
  
  # returns the vector of sampled population sizes
  return(mcmclN/log(10))
  
}


# Markov Chain Monte Carlo simulations to estimate the population size under the boundary mutation model
#  I      : number of MCMC iterates (should be a multiple of 100 and higher than 1000; only every 100th iterate are kept!)
#  counts : the empirical sampled site frequency spectrum
#  mu     : mutation rate

mcmc_boundary <- function(I,counts,mu){

  # sets the acceptance probability to 0 and calculates the likelihood for lN
  lN     <- -log(mu)
  ap     <- 0
  lk <- llik_sampled_boundary_sfs(counts,lN,mu)
  
  # creates a list that includes the mutation rate, the current population size, likelihood 
  # proposal tunning and acceptance probability
  mcmc <- list(lN=lN,mu=mu,counts=counts,lk=lk,tunning=1,ap=ap)
  
  # vector with the sampled values of N
  mcmclN <- rep(NA,(I-1000)/100)
  
  # runs the MCMC
  p <- 1
  for (i in 1:I){
    
    # Metropolis-Hastings step with the boundary mutation model
    mcmc <- mh_boundary_N(mcmc)
      
    # saves each 100th iterate
    # forces a burning phase of 1000 iterates
    if ((i%%100 ==0) && i>1000 ){
        
      # saves the effective population size
      mcmclN[p] <- mcmc$lN
      p <- p+1
        
      # calculates the acceptance probability
      aprob <- mcmc$ap/i
        
      # returns the proposal variance according to the acceptance probability so that
      # the acceptance probability revolves around 0.345
      mcmc$tunning <- mcmc$tunning*(-0.8312691*aprob*aprob+2.331269*aprob+0.5)

    }
  }
  
  # returns the vector of sampled population sizes
  return(mcmclN/log(10))
  
}


# Markov Chain Monte Carlo simulations to estimate the population size under the recurrent mutation model and uncertain mu
#  I      : number of MCMC iterates (should be a multiple of 100 and higher than 1000; only every 100th iterate are kept!)
#  counts : the empirical sampled site frequency spectrum
#  mu     : mutation rate
#  sd     : standard deviation of the mutation rate: samples are taken from a Normal(mu,sd)

mcmc_recurrent_umu <- function(I,counts,mu,sd){
  
  # sets the acceptance probability to 0 and calculates the likelihood for lN
  lN     <- -log(mu)
  ap     <- 0
  lk <- llik_sampled_recurrent_sfs(counts,lN,mu)
  
  # creates a list that includes the mutation rate, the current population size, likelihood 
  # proposal tunning and acceptance probability
  mcmc <- list(lN=lN,mu=mu,sd=sd,counts=counts,lk=lk,tunning=0.1,ap=ap)
  
  # vector with the sampled values of N
  mcmclN <- rep(NA,(I-1000)/100)
  
  # runs the MCMC
  p <- 1
  for (i in 1:I){
    
    # Metropolis-Hastings step with the recurrent mutation model
    mcmc <- mh_recurrent_N(mcmc)

    # saves each 100th iterate
    # forces a burning phase of 1000 iterates
    if ((i%%100 ==0) && i>1000 ){
      
      # saves the effective population size
      mcmclN[p] <- mcmc$lN
      p <- p+1
      
      # calculates the acceptance probability
      aprob <- mcmc$ap/i
      
      # returns the proposal variance according to the acceptance probability so that
      # the acceptance probability revolves around 0.345
      mcmc$tunning <- mcmc$tunning*(-0.8312691*aprob*aprob+2.331269*aprob+0.5)

      
      # samples a new mutation rate
      # updates the mutation rate and the likelihood
      mcmc$mu <- rnorm(1,mu,sd)
      mcmc$lk <- llik_sampled_recurrent_sfs(mcmc$counts,mcmc$lN,mcmc$mu)
      mcmc    <- mh_recurrent_N(mcmc)
      
    }
  }
  
  # returns the vector of sampled population sizes
  return(mcmclN/log(10))
  
}




# Markov Chain Monte Carlo simulations to estimate the population size under the boundary mutation model and uncertain mu
#  I      : number of MCMC iterates (should be a multiple of 100 and higher than 1000; only every 100th iterate are kept!)
#  counts : the empirical sampled site frequency spectrum
#  mu     : mutation rate
#  sd     : standard deviation of the mutation rate: samples are taken from a Normal(mu,sd)

mcmc_boundary_umu <- function(I,counts,mu,sd){
  
  # sets the acceptance probability to 0 and calculates the likelihood for lN
  lN     <- -log(mu)
  ap     <- 0
  lk <- llik_sampled_boundary_sfs(counts,lN,mu)
  
  # creates a list that includes the mutation rate, the current population size, likelihood 
  # proposal tuning and acceptance probability
  mcmc <- list(lN=lN,mu=mu,sd=sd,counts=counts,lk=lk,tunning=0.1,ap=ap)
  
  # vector with the sampled values of N
  mcmclN <- rep(NA,(I-1000)/100)
  
  # runs the MCMC
  p <- 1
  for (i in 1:I){
    
    # Metropolis-Hastings step with the boundary mutation model
    mcmc <- mh_boundary_N(mcmc)
    
    # saves each 100th iterate
    # forces a burning phase of 1000 iterates
    if ((i%%100 ==0) && i>1000 ){
      
      # saves the effective population size
      mcmclN[p] <- mcmc$lN
      p <- p+1
      
      # calculates the acceptance probability
      aprob <- mcmc$ap/i
      
      # returns the proposal variance according to the acceptance probability so that
      # the acceptance probability revolves around 0.345
      mcmc$tunning <- mcmc$tunning*(-0.8312691*aprob*aprob+2.331269*aprob+0.5)
      
      
      # samples a new mutation rate
      # updates the mutation rate and the likelihood
      mcmc$mu <- rnorm(1,mu,sd)
      mcmc$lk <- llik_sampled_boundary_sfs(mcmc$counts,mcmc$lN,mcmc$mu)
      mcmc    <- mh_boundary_N(mcmc)
      
    }
  }
  
  # returns the vector of sampled population sizes
  return(mcmclN/log(10))
  
}
