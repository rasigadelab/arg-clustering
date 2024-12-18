# VERSION 2019-08-17
# Collection of routines for the parametric SOM
# Parameters reflect text, where
# a = alpha, the weight of the Beta component in the arc distance distribution
# b = beta, the parameter of the Beta component in the arc distance distribution
# g = gamma, the parameter of the Beta component in the interval size distribution

#' Larger g means smaller DNA length mobilized
#' Larger a means less uniform-distributed arc distances
#' Larger b means smaller arc distances
#' 

# Arc distance PDF (Probability Density Function)
darc <- function(x, a, b) { # x = arc distance (interval [0,0.5]), a & b = parameters
  2 * a * b * (1 - 2*x)^(b-1) + 2*(1-a)
}
# plot(darc(seq(0,0.5,0.01),1,0.5)) # afar genes
# plot(darc(seq(0,0.5,0.01),1,1)) # a=1, b=1 (uniform distrib)
# plot(darc(seq(0,0.5,0.01),1,10)) # relatively close genes
# plot(darc(seq(0,0.5,0.01),1,100)) # close genes
# plot(darc(seq(0,0.5,0.01),1,1000)) # very close genes
# plot(darc(seq(0,0.5,0.01),2,10)) # idem form but higher weight


# Arc distance CDF (Cumulative Distribution Function)
parc <- function(x, a, b) {
  a * (1 - (1 - 2*x)^b) + (1 - a)*2*x
}
# plot(parc(seq(0,0.5,0.01),0.01,0.01)) # a=0.01, b=0.01
# plot(parc(seq(0,0.5,0.01),0.5,0.5)) # a=0.5, b=0.5
# plot(parc(seq(0,0.5,0.01),1,1)) # a=1, b=1

# Arc distance mean and variance
muvar <- function(a, b) {
  mu <- .5 * (a / (b+1) + (1 - a) / 2)
  s2 <- 1/4 * ( 2 * a / ((b+1) * (b+2)) + (1-a)/3) - mu^2
  return(c(mu = mu, var = s2))
}

# Retrieve arc distance parameters from mean (m) and variance (s)
absolve <- function(m, s) { # m = mean, s = variance

  alpha <- ( -24 * m * (m^2 + s) + 10*m*(3*m - 1) + 6*s + 1 ) / (18*(m^2+s) - 10*m +1)
  beta  <- (-12*(m^2 + s) + 8*m - 1) / (6 * (m^2 + s) - 2*m)
  
  alpha[!is.finite(alpha)] <- 0
  beta[!is.finite(beta)] <- 1
  
  return(c(alpha = pmax(alpha, 0), beta = pmax(beta, 1))) # return a=alpha & b=beta of arc distance beta-distribution
}

# Conditional probability of cotransfer for fixed s, P(A|B,s)
comob <- function(s, a, b) { # s = size transfered DNA
  a*(1 - (1 - (1 - 2*s)^(b+1)) / ( 2 * s * (b + 1))) + (1 - a) * s
}

# Unconditional probability of cotransfer for fixed s, p(AB|s)
comobu <- function(s, a, b) { # s = size transfered DNA
  a*(s - (1 - (1 - 2*s)^(b+1)) / ( 2 * (b + 1))) + (1 - a) * s^2
}

# Unconditional probability of cotransfer, Beta-distributed s
comobg <- function(g, a, b) { # g = distribution parameter gamma, a = alpha, b = beta
  .5 * g / (b+1) * (
    a / (b+g+1) +
      (a * (g+1) * (b-g-1) + b + 1) / (g*(g+1)*(g+2))
  )
}
#plot(sapply(seq(1,500,10), function(x) {comobg(x,1,1E4)}))
#plot(sapply(seq(1,500,10), function(x) {comobg(x,1,1)}))


# Conditional mean and variance, Beta-distributed s
muvarg <- function(g, a, b) { # g = s distribution parameter gamma, a = alpha & b = beta beta-distribution parameters
  pab <- comobg(g, a, b) # Unconditional probability of cotransfer, Beta-distributed s
  m <- 1 / (4*(g+1)) * (
    (a*b)/prod(g + b + 1:2) + (1-a) / prod(g + 2:3)
  ) / pab
  s2 <- 1 / (4*(g+1)) * (
    (a*b)/prod(g + b + 1:3) + (1-a) / prod(g + 2:4)
  ) / pab - m^2
  return(c(mu = m, var = s2))
}

# Periodic antibiotic treatment
# freq <- function(t, period) { # periodic antibiotic treatment (both atb)
#   as.integer((ceiling(t) - 1) %% (2*period) >= period) # 0 (no antibiotic) or 1 (during treatment), start at 0
# }

freq <- function(t, period1, period2=period1) { # different period length w/ & w/o atb
  as.integer((((t-1)/(period1+period2)) %% 1) < (period1/(period1+period2)))
}
#plot(sapply(1:1000, function(x){freq_bis(x, 100, 200)}))
#plot(sapply(1:1000, function(x){freq_bis(x, 100)}))

# deathRateFunc <- function(t) {
#   return(c(
#     d00 = c(0.2, 0.3)[1 + freq(t, deathPeriod)],
#     d01 = c(0.2, 0.3)[1 + freq(t, deathPeriod)],
#     d10 = c(0.2, 0.3)[1 + freq(t, deathPeriod)],
#     d11 = 0.2
#   ))
# }

deathRateFunc <- function(t, deathBasal=0.2, deathIncrease=0.1, deathPeriod=100) {
  return(c(
    d00 = c(deathBasal, deathBasal+deathIncrease)[1 + freq(t, deathPeriod)],
    d01 = c(deathBasal, deathBasal+deathIncrease)[1 + freq(t, deathPeriod)],
    d10 = c(deathBasal, deathBasal+deathIncrease)[1 + freq(t, deathPeriod)],
    d11 = deathBasal
  ))
}

deathRateFunc_bis <- function(t, deathBasal, deathIncrease, deathPeriod, treatment=1, lifePeriod=deathPeriod) {
  if(treatment==1){ # synchronized antibiotics
    return(c(
      d00 = c(deathBasal, deathBasal+deathIncrease)[1 + freq(t, deathPeriod)],
      d01 = c(deathBasal, deathBasal+deathIncrease)[1 + freq(t, deathPeriod)],
      d10 = c(deathBasal, deathBasal+deathIncrease)[1 + freq(t, deathPeriod)],
      d11 = deathBasal
    ))
  }
  if(treatment==2){ # alterned antibiotics
    return(c(
      d00 = c(deathBasal, deathBasal+deathIncrease)[2],
      d01 = c(deathBasal, deathBasal+deathIncrease)[1 + freq(t, deathPeriod)],
      d10 = c(deathBasal, deathBasal+deathIncrease)[2 - freq(t, deathPeriod)],
      d11 = deathBasal
    ))
  }
  if(treatment==3){ # synchronized antibiotics, uneven periods atb/no atb
    return(c(
      d00 = c(deathBasal, deathBasal+deathIncrease)[1 + freq(t, deathPeriod, lifePeriod)],
      d01 = c(deathBasal, deathBasal+deathIncrease)[1 + freq(t, deathPeriod, lifePeriod)],
      d10 = c(deathBasal, deathBasal+deathIncrease)[1 + freq(t, deathPeriod, lifePeriod)],
      d11 = deathBasal
    ))
  }
} # include different treatment types

deathRateFunc_tris <- function(t, db, di, dp, trt=3, lp=dp) {
  if(trt==1){ # no antibiotic
    return(c(
      d00 = db,
      d01 = db,
      d10 = db,
      d11 = db
    ))
  }
  if(trt==2){ # constant antibiotic
    return(c(
      d00 = db+di,
      d01 = db+di,
      d10 = db+di,
      d11 = db
    ))
  }
  if(trt==3){ # synchronized antibiotics
    return(c(
      d00 = c(db, db+di)[1 + freq(t, dp)],
      d01 = c(db, db+di)[1 + freq(t, dp)],
      d10 = c(db, db+di)[1 + freq(t, dp)],
      d11 = db
    ))
  }
  if(trt==4){ # alterned antibiotics
    return(c(
      d00 = c(db, db+di)[2],
      d01 = c(db, db+di)[1 + freq(t, dp)],
      d10 = c(db, db+di)[2 - freq(t, dp)],
      d11 = db
    ))
  }
  if(trt==5){ # synchronized antibiotics, uneven periods atb/no atb
    return(c(
      d00 = c(db, db+di)[1 + freq(t, dp, lp)],
      d01 = c(db, db+di)[1 + freq(t, dp, lp)],
      d10 = c(db, db+di)[1 + freq(t, dp, lp)],
      d11 = db
    ))
  }
} # identical to deathRateFunc_bis but with different parameter names

#plot(sapply(1:500, function(x){deathRateFunc(x)[2]}))
#plot(sapply(1:1000, function(x){deathRateFunc_bis(x, deathBasal=0.2, deathIncrease=0.1, deathPeriod=100, treatment=3, lifePeriod=300)[3]})) # x
#plot(sapply(1:1000, function(x){deathRateFunc_bis(x, deathBasal=0.1, deathIncrease=0.2, deathPeriod=50, treatment=3, lifePeriod=200)[3]})) # x
#plot(sapply(1:5000, function(x){deathRateFunc_bis(x, deathBasal, deathIncrease, deathPeriod, treatment=3, lifePeriod)[3]}))
#plot(sapply(1:5000, function(x){deathRateFunc_tris(x, deathBasal, deathIncrease, deathPeriod, treatmentType, lifePeriod)[3]}))
