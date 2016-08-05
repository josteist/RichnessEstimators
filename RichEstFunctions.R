# Functions:

# Copy of Functions_DinoBino_toweb.R taken 030816. 
# From C:\Users\josteist\Documents\R\DinoBinoGit\DinoBino



# This library must be a bit cleaned up before publication.

# get_ltt <- function(Out,tmax){
# OLD
#   # JOS: 170915
#   # JOS: 051015 Updated to include extinction events also (out[1:2])
#   # Function to tally lineages through time from the output of a simulateBDF.
#   # Out is an array of size (no_lineages) by (3)
#   # where [,1] is start of lineage
#   # [,2] end of lineages
#   # [,3] number of fossils left (not used for this function)
#   spectimes=sort(unique(c(tmax,as.vector(Out[,1:2]))))
#   N_tt = array(NA,c(length(spectimes),1));
#   for (ii in 1:length(spectimes)){
#     # N_tt[ii] = sum(Out[,2]>=spectimes[ii] & Out[,1]<=spectimes[ii])
#     N_tt[ii] = sum(Out[,2]>spectimes[ii] & Out[,1]<=spectimes[ii])
#   }
#   
#   return(list(spectimes,N_tt))
# }

get_ltt <- function(Out){
  # JOS: 170915
  # Update 250915: Now includes the true ltt, not the cum only.
  # Function to tally lineages through time from the output of a simulateBDF.
  # Out is an array of size (no_lineages) by (3)
  # where [,1] is start of lineage
  # [,2] end of lineages
  # [,3] number of fossils left (not used for this function)
  spectimes=sort(unique(as.vector(Out[,1:2])))
  N_tt = array(NA,c(length(spectimes)));
  # Need to think here; spectimes is when events happen. So we want N_tt[ii] to be
  # the number of lineages between spectimes[ii] and spectimes[ii+1]. I.e. if a lineage
  # started at spectimes[ii] we includ it. If it dies out at spectimes[ii] we do not.
  for (ii in 1:(length(spectimes)-1)){
    N_tt[ii] = sum(Out[,2]>spectimes[ii] & Out[,1]<=spectimes[ii])
  }
  N_tt[ii+1]=N_tt[ii];
  
  return(list(spectimes,N_tt))
  # list(c(spectimes,max(Out[,2])),N_tt))
}


doTRiPS <- function(Out){
  # doing the whole TRiPS estimation with the Out array from a simulation.
  # 
  # Occs = Out[Out[,3]>0,3];
  # dTs = max(Out[,2]); # assuming this is the duration.
  p_lam = estimatePoiss(rep(max(Out[,2]),sum(Out[,3]>0)),Out[Out[,3]>0,3])
  p_bino = 1-exp(-max(Out[,2])*p_lam);
  N_true = c(estimatetrue(sum(Out[,3]>0),p_bino[1])[1],
             min(estimatetrue(sum(Out[,3]>0),p_bino[3])),
             max(estimatetrue(sum(Out[,3]>0),p_bino[2])))
  
  
}

doTRiPS_abs <- function(abs,t=1){
  # Performs TRiPS on the count of observations in abs. t=1 (default) is assumed duration
  p_lam = estimatePoiss(rep(t,length(abs)),abs);
  p_bino = 1-exp(-p_lam*t);
  N_true = c(estimatetrue(length(abs),p_bino[1])[1],
             min(estimatetrue(length(abs),p_bino[3])),
             max(estimatetrue(length(abs),p_bino[2])))
  out = array(NA,c(3,3))
  
  out[1,]=p_lam
  out[2,]=p_bino
  out[3,]=N_true
  rownames(out)<-c("Sampling rate","Sampling probability","Estimated richness")
  colnames(out)<-c("MLE","lower CI","upper CI")
  return(out)
}

SimulateBDF <- function(spec,mu,lambda,lambdavar,tmax,n_init){
  ## Simulating a birth-death-fossilize model.
  # Written by Jostein Starrfelt (jostein.starrfelt[AT]ibv.uio.no)
  # Last checked 02.09.15
  #   
  #   spec = 0.01;
  #   mu = 0;
  #   lambda = 0.4;
  #   lambdavar = 0;
  #   tmax = 10;
  #   n_init = 100;
  # Assuming fixed speciation and extinction rates for whole interval.
  rho <- function(t,spec,mu){mu*t - spec*t}
  rho2 <- function(t){mu*t - spec*t}
  tmpfun <- function(t){exp(-rho2(t))*spec}
  nosp <- function(t,spec,mu,n_init){ n_init*(1 + integrate(tmpfun,0,t)$value)}
  # Rough estimate of total species richness expected by Kendall
  minspec <- round(nosp(tmax,spec,mu,n_init))
  Out = array(NA,c(minspec*10,3)); # preallocating the array for species
  done = 0;
  Out[1:n_init,1]=0;
  tix = 1;
  ntix = n_init+1; #index into first row in Out with no entries yet.
  while (done==0) {
    # When does this lineage go extinct
    if (mu>0){ # if extinction rate is nonzero
      Out[tix,2] = min(Out[tix,1]+rexp(1,mu),tmax)
    } else {
      Out[tix,2] = tmax;
    }
    # Drawing number of fossils for this lineage
    # draw one number from Poisson distribution with mean drawn from a normal distribution with mean lambda
    # and st.dev lambdavar times the duration of the taxon.
    Out[tix,3] = rpois(1,max(0,rnorm(1,lambda,lambdavar))*(Out[tix,2]-Out[tix,1]))
    
    # Drawing waiting times to possible speciation events.
    if (spec>0){
      tmp <- rexp(1e1,spec)
      tmptix = 1; 
      # The above draw might be too small;
      while (sum(tmp)<(Out[tix,2]-Out[tix,1])){
        tmp = rexp(10^(tmptix),spec);
        tmptix = tmptix+1; # increase number of draws until end of
        # lineage is traversed
      }
      spectimes = cumsum(tmp)<(Out[tix,2]-Out[tix,1]);
      if (sum(cumsum(tmp)<(Out[tix,2]-Out[tix,1]))>0){
        
        if ((ntix+sum(spectimes)-1)>nrow(Out)){
          ## need to enlarge OUT
          oldout= Out;
          # Generating new array ~50 % bigger.
          Out = array(NA,c(round(nrow(Out)*1.5),3));
          Out[1:(ntix-1),]=oldout[1:(ntix-1),];
        }
        # if some of these 'speciation times' are inside the actual duration of taxon tix
        Out[ntix:(ntix+sum(spectimes)-1),1]=Out[tix,1]+cumsum(tmp)[spectimes]
        ntix = ntix+sum(spectimes)
        
      }
    } else {
      # tmp = 0;
      # Do nothing if speciation rate is 0
    }
    
    tix = tix+1;
    if (tix==(ntix)){
      done=1;
    }
    
  }
  Out = Out[!is.na(Out[,1]),]; # removin excess NA's
  return(Out)
}
# Drawing a number/index from an empirical probability distribution.
Emprand <- function(x) {
  j <- runif(1) ## drawing random uniform number
  which(cumsum(x)/sum(x)>j)[1] #returning first entry in the cumsum/sum [empirical density function]
}

## Likelihood function of the poisson intensity conditioned on more than one occurrence.
likepoissint <- function(lambda,dt,nobs) {
  (((lambda*dt)^(nobs))/(factorial(nobs))*(exp(-lambda*dt)))/(1- exp(-lambda*dt))
}
# Can we utilize the much simpler log likelihood for Poisson, even when we want to condition on >0?
# These are the same, but perhaps with some numerical stability
# JOS;: 030915
# Calculating and storing the log(factorial(n!)) needed for Poisson estimation, 
tab_lfactorials = array(NA,c(1000,1));
tab_lfactorials[1] = log(1);
for (ii in 1:1000){
  tab_lfactorials[ii+1] = tab_lfactorials[ii] + log(ii+1)
}
likepoissint2 <-function(lambda,dt,nobs){
  nobs*log(lambda*dt) - length(nobs)*(lambda*dt) - log(1-exp(-lambda*dt)) - tab_lfactorials[nobs];
}


## Negative Log likelihood for many observations
# Make more general, i.e. date multiple nobs, but one dt value (then repeat it)
# 
loglikepoissint <- function(lambda,dtin,nobs) {
  if (length(nobs)!=length(dtin)){
    dt = rep(dt,length(nobs))
  }else{
    dt = dtin
  } # if not equal lengths, assume dt has 1 value and repeat
  
  ll <- 0
  for (ii in (1:length(dt))){
    # ll <- ll - log(likepoissint(lambda,dt[ii],nobs[ii]))
    # testing for likepoissint2 (using log likelihod for stability in calcs.)
    ll <- ll - likepoissint2(lambda,dt[ii],nobs[ii]);
  }
  return(ll)
}

# Function to return phat and pci
estimatePoiss <- function(dTs,Occs){
  # Negative log likelihood of the observations.
  nllnow <- function(lambda){ loglikepoissint(lambda,dTs,Occs)}
  # Getting maximum likelihood estimate of the poisson rate.
  fit1 <- mle(nllnow,start = list(lambda=1),nobs=length(Occs),method="Brent",lower = 1e-8,upper=30)
  # Defining the likelihood ratio function to estimate confidence intervals
  mycis2 <- function(l_x) ((2*(-nllnow(coef(fit1))+nllnow(l_x)))-qchisq(0.95,1))^2
  
  # lower CI
  ci1<-optim(0.9*coef(fit1),mycis2,method="Brent",lower=1e-6,upper=coef(fit1))
  # upper CI
  ci2<-optim(1.1*coef(fit1),mycis2,method="Brent",lower=coef(fit1),upper=40)
  return(c(coef(fit1),ci1$par,ci2$par))
}

## Function to collect occurrances and durations from PBDB output
# This is custom made for dinosaur data, it will only extract occurrences from
# pbdb intervals 112:138
createDataArrs <- function(dinos){
  
  ## Counting occurrences inside each bin for each uniqe species
  # A uniqe indexlist of occurrences matched to species rank.
  uniqspec <- unique(dinos$mid[dinos$mra==3])
  # Data has [species by interval] with number of occurrences per species in each interval.
  Data = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Times are the durations in a matrix of same size for ease of computation.
  Times = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Times[,ii] <- Bins[ii,3]
    
  }
  ## species by interval matrix
  tix = 1 # loop counter
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # For each occurrence of this species
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # if the last interval of this occurrence in this species is before Maastricthian-1
        # OR after Ladinian (WE include Ladinian here)
        
        if (j1[jj]>j2[jj]) {
          countdoubles=countdoubles+1
          bix = seq(j1[jj],j2[jj])
          x = Bins[bix-(Bins[1,4]-1),3]			
          binow <- bix[Emprand(x)]
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Data[tix,binow-(Bins[1,4]-1)] <- Data[tix,binow-(Bins[1,4]-1)]+1
        }
      }
    }
    tix <-tix+1
  }
  results <- list(Data = Data, Times=Times)
  return(results)
}


createDataArrs_v2 <- function(dinos,removedoubles=FALSE){
  
  ## Counting occurrences inside each bin for each uniqe species
  # A uniqe indexlist of occurrences matched to species rank.
  # This function adds a third list [Isbird] which is true/false if the taxon is a bird.
  # This is to 'unify' the datamatrices with dinos and theropods with and without birds.
  # if removedoubles is set to TRUE, each occurrence that spans more than one interval is ignored.
  uniqspec <- unique(dinos$mid[dinos$mra==3])
  # Data has [species by interval] with number of occurrences per species in each interval.
  Data = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  Isbird = array(NA,c(length(uniqspec),1))
  # Times are the durations in a matrix of same size for ease of computation.
  Times = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Times[,ii] <- Bins[ii,3]
    
  }
  ## species by interval matrix
  tix = 1 # loop counter this is for each unique species
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # TESTING TO INCLUDE LADINIAN_ 040615
        if (j1[jj]>j2[jj]) {
          if (removedoubles==TRUE){
            binow <- -100 #make it not count.
          }else{
            countdoubles=countdoubles+1
            bix = seq(j1[jj],j2[jj])
            x = Bins[bix-(Bins[1,4]-1),3]  		
            binow <- bix[Emprand(x)]
          }
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Data[tix,binow-(Bins[1,4]-1)] <- Data[tix,binow-(Bins[1,4]-1)]+1
          
          
        }
      }
    }    
    Isbird[tix] = any(dinos[dinos$mid==ii,]$cln==36616)
    Isbird[tix]
    tix <-tix+1
  }
  
  results <- list(Data = Data, Times=Times,Isbird=Isbird)
  return(results)
}


createDataArrs_v3 <- function(dinos){
  
  ## Counting occurrences inside each bin for each uniqe species
  # A uniqe indexlist of occurrences matched to species rank.
  uniqspec <- unique(dinos$mid[dinos$mra==3])
  # Dats has [species by interval] with number of occurrences per species in each interval.
  Dats = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  uniqnam<-unique(dinos$mna[dinos$mra==3])
  rownames(Dats)<-uniqnam;
  colnames(Dats)<-interval.names;
  # Trying to not defined Data, but do it iteratively to accurately get the ones included. earlier (bf 070815)
  # this code included rows with no occurrences used.
  # The prealloc is needed for computation, how about removing the zeros
  # Times are the durations in a matrix of same size for ease of computation.
  Tims = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Tims[,ii] <- Bins[ii,3]
    
  }
  ## species by interval matrix
  tix = 1 # loop counter
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # For each occurrence of this species
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # if the last interval of this occurrence in this species is before Maastricthian-1
        if (j1[jj]>j2[jj]) {
          countdoubles=countdoubles+1
          bix = seq(j1[jj],j2[jj])
          x = Bins[bix-(Bins[1,4]-1),3]			
          binow <- bix[Emprand(x)]
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Dats[tix,binow-(Bins[1,4]-1)] <- Dats[tix,binow-(Bins[1,4]-1)]+1
        }
      }
    }
    tix <-tix+1
  }
  tmp<-rowSums(Dats)>0;
  Data = Dats[tmp,]
  Times = Tims[tmp,]
  results <- list(Data = Data, Times=Times)
  return(results)
}


createDataArrs_v4 <- function(dinos){
  
  ## Counting occurrences inside each bin for each uniqe GENUS
  # A uniqe indexlist of occurrences matched to species rank.
  uniqspec <- unique(dinos$gnn)
  uniqspec = uniqspec[!is.na(uniqspec)] # if there's a NA
  
  # Dats has [species by interval] with number of occurrences per species in each interval.
  Dats = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  dinos[which(dinos$gnn==uniqspec[2]),]$gnl
  uniqnam = array(NA,dim=c(length(uniqspec),1))
  for (ii in 1:length(uniqspec)){
    uniqnam[ii] = unique(dinos[which(dinos$gnn==uniqspec[ii]),]$gnl)
  }
  # uniqnam<-unique(dinos$gnl)
  rownames(Dats)<-uniqnam;
  colnames(Dats)<-interval.names;
  # Trying to not defined Data, but do it iteratively to accurately get the ones included. earlier (bf 070815)
  # this code included rows with no occurrences used.
  # The prealloc is needed for computation, how about removing the zeros
  # Times are the durations in a matrix of same size for ease of computation.
  Tims = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Tims[,ii] <- Bins[ii,3]
    
  }
  ## species by interval matrix
  tix = 1 # loop counter
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[which(dinos$gnn==ii)]#dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[which(dinos$gnn==ii)]#dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # For each occurrence of this species
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # if the last interval of this occurrence in this species is before Maastricthian-1
        if (j1[jj]>j2[jj]) {
          countdoubles=countdoubles+1
          bix = seq(j1[jj],j2[jj])
          x = Bins[bix-(Bins[1,4]-1),3]			
          binow <- bix[Emprand(x)]
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Dats[tix,binow-(Bins[1,4]-1)] <- Dats[tix,binow-(Bins[1,4]-1)]+1
        }
      }
    }
    tix <-tix+1
  }
  tmp<-rowSums(Dats)>0;
  Data = Dats[tmp,]
  Times = Tims[tmp,]
  results <- list(Data = Data, Times=Times)
  return(results)
}

# To estimate maximum likelihood of the true number of species given observed number and 
# binomial sampling probability.
# pdf = 
estimatetrue <- function(nobs,binomprob) {
  if (!is.na(binomprob)){
    n <- seq(0,nobs/binomprob * 4+10)
    liks <- log(dbinom(nobs,size=n,prob=binomprob))
    tmp <- n[which((2*(max(liks)-liks))<qchisq(0.95,1))]
    
    return(c(n[which.max(liks)],min(tmp),max(tmp)))
  } else {
    return(c(NA,NA,NA))
  }
}

