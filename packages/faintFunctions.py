# -*- coding: iso-8859-1 -*-
#
# Functions used to correct the faint end of observations 
#
from scipy.stats import norm
from scipy.special import erf
import numpy as np

def calculate_mean(psfFlux, psfFluxErr):
    xObs = psfFlux / psfFluxErr
    
    # Initialize as zeros.  if xObs is so small ( eg < -50 ) that 
    # we're getting norm.sf(-xObs) = 0, then we can't divide by that 
    # number --> only divide if sf is nonzero
    xMean = np.zeros_like(xObs)  
    
    # calculate the survival function 
    sf = norm.sf(-xObs )
    ind = np.where(sf !=0 )

    # only divide by sf where the divisor is nonzero ... 
    xMean[ind] = (1/ (sf[ind]*np.sqrt(2*np.pi))) * np.exp(-(xObs[ind]**2.0) / 2.0) + xObs[ind]
    F_mean = xMean * psfFluxErr
    return F_mean

def calculate_median(psfFlux, psfFluxErr):
    xObs = psfFlux / psfFluxErr
    return psfFluxErr*norm.ppf((1+norm.cdf(-xObs))/2.0) + psfFlux


def calculate_2sigma(psfFlux, psfFluxErr):
    xObs = psfFlux / psfFluxErr
    return psfFlux + psfFluxErr * norm.isf(0.05 * norm.sf(-xObs))

def calculate_rms(psfFlux, psfFluxErr):
    xObs  = psfFlux / psfFluxErr
    xMean = calculate_mean(psfFlux, psfFluxErr) / psfFluxErr
    delX  = xObs - xMean
    I1    = norm.sf(-xObs)
    I0bysig2 = 0.5*erf(xObs/np.sqrt(2)) + (1.0/np.sqrt(2*np.pi))*np.exp(-(xObs**2.0) / 2.0)*(2*delX - xObs) + 0.5 + delX*delX*norm.sf(-xObs)

    xRMS = np.ones_like(xObs)
    xRMS[I0bysig2 > 0] = np.sqrt(I0bysig2[I0bysig2 > 0] / I1[I0bysig2 > 0])
    return  xRMS * psfFluxErr