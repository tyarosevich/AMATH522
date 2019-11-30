###########################################################################
#                                                                         #
#   ESTIMATING THE REPRODUCTION NUMBER OF EBOLA VIRUS (EBOV)              #
#   DURING THE 2014 OUTBREAK IN WEST AFRICA                               #
#                                                                         #
#   REFERENCE: ALTHAUS CL. PLOS CURR. 2014;6.                             #
#   DOI: 10.1371/currents.outbreaks.91afb5e0f279e7f29e7056095255b288      #
#                                                                         #
#   Copyright (C) 2014 by Christian L. Althaus                            #
#   (christian.althaus@alumni.ethz.ch)                                    #
#                                                                         #
#   This script is free software; you can redistribute it and/or modify   #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation; either version 2 of the License, or     #
#   (at your option) any later version.                                   #
#                                                                         #
#   This script is distributed in the hope that it will be useful,        #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#   GNU General Public License for more details.                          #
#                                                                         #
###########################################################################

# Necessary libraries
library(deSolve)
library(chron)
library(bbmle)

# Read the data file
ebola <- read.csv("2014-2016 epidemic case counts_inverted.csv")
ebola$Date <- chron(as.character(ebola$Date),format=c(dates = "day-mon-year"))

# Definition of the SEIR model
SEIR <- function(t, x, parms) {
  with(as.list(c(parms,x)),{
    if(t < tau1) beta <- beta0
    else beta <- beta0*exp(-k*(t-tau1))
    N <- S + E + I + R
    dS <- - beta*S*I/N
    dE <- beta*S*I/N - sigma*E
    dI <- sigma*E - gamma*I
    dR <- (1-f)*gamma*I
    dD <- f*gamma*I
    dC <- sigma*E
    der <- c(dS,dE,dI,dR,dD,dC)
    list(der)
  })
}



# Negative log-likelihood
nll <- function(beta0,k,f,tau0,tau1,sigma,gamma) {
  pars <- c(beta0=beta0,k=k,f=f,tau0=tau0,tau1=tau1,sigma=sigma,gamma=gamma)
  pars <- trans(pars)
  times <- c(0,data$times+pars["tau0"])
  # Time steps the SEIR model using current parameters.
  simulation <- as.data.frame(ode(init,times,SEIR,parms=pars))
  simulation <- simulation[-1,]
  # This assumes a Poisson distribution to the data, and then I THINK passes the 
  # actual case information with the simulation result passed as the lambda value, or
  # the mean. So it's going to return some kind of difference value. This is then returned 
  # by the nll function with the negative sign, which is a convention of this optimization 
  # approach - you maximize the negative log likelihood, thus minimizing the error.
  ll <- sum(dpois(data$cases,simulation$C,log=TRUE)) + sum(dpois(data$deaths,simulation$D,log=TRUE))
  return(-ll)
}

# Parameter transformation
trans <- function(pars) {
  pars["beta0"] <- exp(pars["beta0"])
  pars["k"] <- exp(pars["k"])
  pars["f"] <- plogis(pars["f"])
  pars["tau0"] <- exp(pars["tau0"])
  pars["tau1"] <- exp(pars["tau1"])
  return(pars)
}

# 2014-16 Guinea Cases: Prepare the data, set the initial values and fit the model
data <- na.omit(ebola[c("Date","guinea_cases","guinea_deaths")])
names(data) <- c("times","cases","deaths")
begin <- min(data$times)
data$times <- data$times - data$times[1]
N <- 1e6
init <- c(S = N - 1, E = 0, I = 1, R = 0, D = 0, C = 1)
fixed <- c(tau1 = -Inf, sigma = 1/5.3, gamma = 1/5.61)
free <- c(beta0 = log(0.2), k = log(0.001), f = 0, tau0 = log(60))
fit <- mle2(nll,start=as.list(free),fixed=as.list(fixed),method="Nelder-Mead",control=list(maxit=1e3))
sim_pars <- trans(coef(fit))

# Case and Death numerical solution
#pars <- c(beta0=sim_pars[1],k=sim_pars[2],f=sim_pars[3],tau0=sim_pars[4],tau1=sim_pars[5],sigma=sim_pars[6],gamma=sim_pars[7])

times <- c(data$times+sim_pars[[4]])
num_soln <- as.data.frame(ode(init,times,SEIR,parms=sim_pars))



# Plotting 
graph_deaths <- function(data, num_soln) {
  plot(data$times, data$deaths,col="red")
  par(new = TRUE)
  plot(data$times, num_soln$D,axes = FALSE, xlab = "", ylab = "")
}

graph_cases <- function(data, num_soln) {
  plot(data$times, data$cases, col="red")
  par(new = TRUE)
  plot(data$times, num_soln$C, axes = FALSE, xlab = "", ylab = "")
  
}
graph_deaths(data = data, num_soln = num_soln)