# -*- coding: iso-8859-1 -*-
#
# Functions used to correct the faint end of observations 
#
from scipy.stats import norm
from scipy.special import erf
import numpy as np

def calculate_mean(psfFlux, psfFluxErr):
    xObs = psfFlux / psfFluxErr
    # Calculate my results
    xMean = (1/ (norm.sf(-xObs )*np.sqrt(2*np.pi))) * np.exp(-(xObs**2.0) / 2.0) + xObs
    F_mean = xMean * psfFluxErr
    return F_mean

def calculate_median(psfFlux, psfFluxErr):
    x0= -psfFlux / psfFluxErr
    return psfFluxErr*norm.ppf((1+norm.cdf(x0))/2.0) + psfFlux


def calculate_2sigma(psfFlux, psfFluxErr):
    x0= -psfFlux / psfFluxErr
    return psfFlux + psfFluxErr * norm.isf(0.05 * norm.sf(x0))

def calculate_rms(psfFlux, psfFluxErr):
    xObs = psfFlux / psfFluxErr
    xMean = (1/ (norm.sf(-xObs )*np.sqrt(2*np.pi))) * np.exp(-(xObs**2.0) / 2.0) + xObs
    delX = xObs - xMean
    I1 = norm.sf(-xObs)
    I0bysig2 = 0.5*erf(xObs/np.sqrt(2)) + (1.0/np.sqrt(2*np.pi))*np.exp(-(xObs**2.0) / 2.0)*(2*delX - xObs) + 0.5 + delX*delX*norm.sf(-xObs)
    xRMS = np.sqrt(I0bysig2 / I1)
    return  xRMS * psfFluxErr