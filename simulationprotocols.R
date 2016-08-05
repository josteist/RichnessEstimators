# Simulation protocol


## How to simulate fossil data.
# Can disentangle a couple of different simulation regimes:

# All components below can be simulated using identical 'sampling rates' for all lineages/individuals
# or different. If identical sampling for individuals (i.e. p(fossil per individual per 'time')), then 
# a 'rate' for the species can be used. TRiPS is estimating the rate per species. Lets assume that
# 'abundance', while possibly changing should be interpreted as a population at carrying capacity with
# level 'abundance', and any sampling from a standing abundance can be seen as a binomial probability
# with almost infitesimally small duration (i.e. ecological time is ignored)

# Still we assume that sampling at the lowest level of variability is Poisson.

dt = 1; # duration
samprate_spec = 0.18; # Sampling rate (ind. fossil finds per species per time)
samprate_ind  = 1e-4; # sampling rate (prob fossiliziation per individual per time). We interpret 'individual'
# here rather peculiar; we assume that a species has a carrying capacity and at all points in time
# has an abundance equal to this CC. We assume ecological turnover to be very rapid, i.e. fossilization
# is simulated by applying a poisson sampling for each 'individual' though these individuals are not proped
# individuals, but rather 'slots' for possible individuals.

# 1. Fixed species pool & fixed abundances (i.e. no "time")
#    1-a fixed number of species and equal abundances for all species
nsp = round(runif(1,min=10,max=1e3))
# nspf <- function() {round(runif(1,min=10,max=1e3))}
ocs = rpois(nsp,dt*samprate_spec)
doTRiPS_abs(ocs[ocs>0],dt)

#    1-b fixed number of species and different abundances for all species
nsp = round(runif(1,min=10,max=1e4));
ocs = array(NA,nsp)
# Liow&Hannisdal uses meanlog=1.75 and sdlog=2, but this is for 
# Quote: " distribution of populations among species in effect lognormal, mean c. 1.75 and sd 2"
abunds = round(rlnorm(nsp,meanlog=1.75,sdlog=2));
# But these are very low. Mean of lognormal is approx exp(mu + sd^2/2)
abunds = round(rlnorm(nsp,meanlog=10,sdlog=1.2));
for (ii in 1:nsp){
  ocs[ii]=sum(rpois(abunds[ii],samprate_ind*dt))
}
doTRiPS_abs(ocs[ocs>0],dt)

#    & With identical vs different sampling rate for each species/individual.
# 2. Fixed species pool & variable abundances (i.e. abundances of species can vary over "time")
#    2-a - This might be moved to variable species numbers, to implement a 'hat' to a lineage?
#    

# 3. Variable species pool
#    3-a All species have same abundance at any point in time, but different durations.
#    3-b All species have differenc abundances, but non varying in time.
#    3-c All species have temporally varying abundances.

# 4. Other 'imperfections'
#    4-a Some occurrences are misidentified
#         - wastebasket taxa, i.e. a couple of species are identified to be the same
#         - split taxa, some species are misidentified to be different, while being the same.
#    4-b impact of 'Lagerstätten'. At some points in time a disproportionate number of species are sampled

# Also potentials;
#   spatial dimension? Some species are somewhere but others are not?
