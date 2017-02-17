# -*- coding: iso-8859-1 -*-
#
# Functions used to calculate variability parameters 
#
# needed by LC_var_stats.py 


from astroML.stats import median_sigmaG
import pandas as pd
import numpy as np
import sys 

def gaussgauss_logL(xi, ei, mu, sigma):
    """Equation 5.63: gaussian likelihood with gaussian errors
    Source code https://github.com/astroML/astroML/blob/master/astroML/
    """
    ndim = len(np.broadcast(sigma, mu).shape)

    xi = xi.reshape(xi.shape + tuple(ndim * [1]))
    ei = ei.reshape(ei.shape + tuple(ndim * [1]))

    s2_e2 = sigma ** 2 + ei ** 2
    return -0.5 * np.sum(np.log(s2_e2) + (xi - mu) ** 2 / s2_e2,
                         -1 - ndim)


def approximate_mu_sigma(xi, ei, axis=None):
    """Estimates of mu0 and sigma0 via equations 5.67 - 5.68"""

    if axis is not None:
        xi = np.rollaxis(xi, axis)
        ei = np.rollaxis(ei, axis)
        axis = 0

    mu_approx, sigmaG = median_sigmaG(xi, axis=axis)
    e50 = np.median(ei, axis=axis)
    var_twiddle = (sigmaG ** 2 + ei ** 2 - e50 ** 2)
    sigma_twiddle = np.sqrt(np.maximum(0, var_twiddle))

    med = np.median(sigma_twiddle, axis=axis)
    mu = np.mean(sigma_twiddle, axis=axis)

    zeta = np.ones_like(mu)
    zeta[mu != 0] = med[mu != 0] / mu[mu != 0]

    var_approx = zeta ** 2 * sigmaG ** 2 - e50 ** 2
    sigma_approx = np.sqrt(np.maximum(0, var_approx))

    return mu_approx, sigma_approx


def get_pdf_info(sigma=None, p_sigma=None):
    '''
    A convenience function to calculate the properties of 
    the p_sigma : it's mean, median, Gaussianity, standard 
    deviation... 
    Used directly in get_mu_sigma()

    '''
    # initialize dictionary to store p_sigma info : 
    stats = {}

    # calculate mean (p_sigma is normalized )
    E = np.sum(p_sigma*sigma) / np.sum(p_sigma)   

    # calculate standard deviation 
    stdev = np.sqrt(np.sum(sigma*sigma*p_sigma)/ np.sum(p_sigma)  - E*E)  

    # calculate weighted median 
    median = calcWeightedPercentile(sigma,p_sigma,50)

    # weighted interquartile sigmaG 
    q75 = calcWeightedPercentile(sigma, p_sigma, 75)
    q25 = calcWeightedPercentile(sigma, p_sigma, 25) 
    sigmaG_weighted  = 0.7413 * (q75-q25)

    # probability of being withing +/- 2 * stdev 
    dsigma = sigma[1]-sigma[0]
    mask = (sigma < E+2*stdev) * (sigma > E-2*stdev)
    prob = np.sum(p_sigma[mask]*dsigma)
   
    # if that probability is at least 90% that of a Gaussian : 
    # that's how I got p_Gaussian : 
    # from scipy.stats import norm
    # mean, var, skew, kurt = norm.stats(moments='mvsk')
    # p_Gaussian = norm.sf(mean-2*np.sqrt(var))-norm.sf(mean+2*np.sqrt(var))
    p_Gaussian = 0.95449973610364158
    if prob >= 0.9*p_Gaussian:
        GaussLike = True
    else:
        GaussLike = False 

    # store all info : 
    stats['mean'] = E
    stats['median'] = median 
    stats['stdev'] = stdev
    stats['sigmaG'] = sigmaG_weighted
    stats['probInTwoStdev'] = prob
    stats['GaussLike'] = GaussLike
    
    return stats 

def get_mu_sigma(xi,ei, N_boot=1000, return_plot_data = False, gridsize=70, sigma_gridsize=None, mu_gridsize=None, 
    return_sigma_pdf_info = False):
    ''' A short function
    to calculate a full mu, sigma, based 
    on N_boot bootstraps of the given sample
    with approximate method, 
    and then calculating the 2D maximum of the
    log-likelihood calculated on the space spanning
    from the minimum to the maximum value of the bootstrapped 
    result. 
    
    Input:
    -------------
    xi : array of measurement values (assumed flux, i.e. order 1e-27)
    ei ; array of measurement errors (same order of mag as xi)
    N_boot : integer, representing the number of bootstraps 
    return_plot_data : if True,  then part from mu, sig, we return also a
        dictionary of quantities sampling the log-likelihood used to find  maximum likelihood estimate of mu, sigma .
    gridsize : the length of grid on which sigma, mu likelihood is evaluated 
    return_sigma_pdf_info : if True , then we also calculate various parameters that describe the shape of the pdf

    Returns:
    ---------------
    pdf_stats, plot_data,  mu_max, sig_max   in that order.  
    --> plot_data is only returned if return_plot_data == True
    --> pdf_stats is only returned if return_sigma_pdf_info == True
    --> mu_max, sig_max are always returned, even if they are not found 
        (then np.nan is returned in case of catastrophic failure, 
         but hopefully all cases have been tested in the code)

    mu_max : mu calculated as a maximum of the 2D log-likelihood 
    sig_mag : sigma calculated as a maximum of the 2D log-likelihood 

    if return_plot_data == True , we also return a dictionary with the 
        following keys : 
          mu, sigma : the grids spanning searched mu and sigma values
          p_mu, p_sigma : pdf of mu, sigma : log-Likelihood marginalized over the other variable
          mu_boot, sigma_boot : results of bootstrapped resamples, that can be used to plot the histogram 
          logL : the actual 2D log-likelihood in mu, sigma space, that can be used
          to plot the contour plot of log-likelihood, analoguous to Fig. 5.8 AstroML
          mu_max, sigma_max : values of mu, sigma that lie at the 2D maximum of the log-likelihood  
    if return_sigma_pdf_info == True , we return a dictionary with the
    following keys : 
     :  
        mean : the mean value;   
        stdev :  the standard deviation; 
        median : the weighted median, 
        sigmaG : the robust Gaussian interquartile deviation, 
        probInTwoStdev : the integral under  p(sigma)  between mean +/- 2 stdev  
        GaussLike : the Gaussian-likeness : if the probability of being within 
            2 standard deviations from the mean is greater than 0.9 of what we 
            would expect for a pure Gaussian (0.95449973610364158), then the 
            p(sigma) has a GaussLike  shape 
    '''

    # Calculate bootstrapped approximate.... 
    #print N_boot
    indices = np.random.randint(0, len(xi), (len(xi), N_boot))

    xi_boot = xi[indices]
    ei_boot = ei[indices]

    mu_boot, sigma_boot = approximate_mu_sigma(xi_boot, ei_boot, 0)


    # Calculate marginalized likelihood sigma and mu
    # using as boundaries the minimum and maximum result of the 
    # bootstrapped approximate calculation 
 
    if sigma_gridsize is None : 
        sigma_gridsize = gridsize 
    if mu_gridsize is None : 
        mu_gridsize = gridsize 

    # max_factor = 1.0  # 
    # 
    # new addon : force sigma grid to go a bit higher... 
    # that way we avoid infinities in case all 
    # bootstrapped resamples of sigma yield 0 ... 
    if max(sigma_boot) < np.std(xi) : 
        # I choose the regular standard deviation, 
        # not the weighted  standard deviation...
        max_sigma =  np.std(xi)
        sigma = np.linspace(0 , max_sigma, sigma_gridsize)
    else:  # how it used to be ...
        sigma = np.linspace(0, max(sigma_boot), sigma_gridsize)
    mu = np.linspace(min(mu_boot), max(mu_boot), mu_gridsize)
    
    if (len(mu) == 0) | (len(sigma) == 0) : 
        # for some reason we may be completely missing it... 
        # this branch would not happen if we prevent the 
        # upper limit to be the same as lower limit 
        # (in case when bootstrapped resamples  all pile up at 
        # 0, so that min(sigma_boot) = max(sigma_boot) = 0 ) 

        return np.nan, np.nan 
    else : 
        logL = gaussgauss_logL(xi, ei, mu, sigma[:, np.newaxis])
        logL -= logL.max()
        ind = np.where(logL == np.max(logL))
        L = np.exp(logL)

        # normalize p(sigma) by the integral over all distribution
        p_sigma = L.sum(1)
        p_sigma /= (sigma[1] - sigma[0]) * p_sigma.sum() 


        # follow this branch if want to get plot quantities  
        if return_plot_data == True : 
            # L = np.exp(logL)

            # normalize p(sigma) by the integral over all distribution
            p_sigma = L.sum(1)
            p_sigma /= (sigma[1] - sigma[0]) * p_sigma.sum() 

            
            # normalize p(mu) by the integral over all distribution
            p_mu = L.sum(0)
            p_mu /= (mu[1] - mu[0]) * p_mu.sum()

            plot_data = {}
            plot_data['mu'] = mu 
            plot_data['p_mu'] = p_mu
            plot_data['sigma'] = sigma
            plot_data['p_sigma'] = p_sigma
            plot_data['mu_boot'] = mu_boot
            plot_data['sigma_boot'] = sigma_boot
            plot_data['logL'] = logL
            plot_data['mu_max'] = mu[ind[1]][0]
            plot_data['sigma_max'] =  sigma[ind[0]][0]
    	  
            if  len(ind) < 2 : 
                # for some reason, we may be unable to find the maximum...
                # but we still want the stats... 
                if return_sigma_pdf_info == True: 
                    stats = get_pdf_info(sigma, p_sigma)
                    return stats, plot_data, np.nan, np.nan

                else: 
                    # we don't need stats, only plot data 
                    return plot_data, np.nan, np.nan

            else : 
                # check if last element of array...
                #if ind[1][0] == len(mu)-1:
                #    print('mu at the last grid point')
                #if ind[0][0] ==  len(sigma)-1:
                #    print('sigma at the last grid point')
    	        # return the mu and sigma at the maximum of the likelihood
    	        # (note : I assume log-likelihood is smooth, and has only 
    	        # one maximum )
                if return_sigma_pdf_info == True  : 
                    # for some reason, we may be unable to find the maximum...
                    stats = get_pdf_info(sigma, p_sigma)
                    return stats, plot_data, mu[ind[1]][0], sigma[ind[0]][0]
                else: 
                    return plot_data, mu[ind[1]][0], sigma[ind[0]][0]
        



        # follow this branch if plot quantities not needed 
        if len(ind) < 2  : 
            # for some reason, we may be unable to find the maximum...
            if return_sigma_pdf_info == True  : 
                stats = get_pdf_info(sigma, p_sigma)
                return stats, np.nan, np.nan
            else: 
                return np.nan, np.nan
        else :  
            #check if last element of array...
            #if ind[1][0] == len(mu)-1 :
            #        print('mu at the last grid point')
            #if ind[0][0] == len(sigma)-1 :
            #        print('sigma at the last grid point')
            # return the mu and sigma at the maximum of the likelihood
            # (note : I assume log-likelihood is smooth, and has only 
            # one maximum )
            if return_sigma_pdf_info == True  : 
                stats = get_pdf_info(sigma, p_sigma)
                return stats, mu[ind[1]][0], sigma[ind[0]][0] 
            else : 
                return mu[ind[1]][0], sigma[ind[0]][0]

def calcChi2raw(y, yerr):
    """Compute simple reduced chi2  (mean-based) if more than 1 datapoints present
    chi2 = np.sum(((flux-meanFlux)**2.0) / (fluxErr ** 2.0)) / (N-1.0)
    """
    N = len(y)
    if N < 2:
        return np.nan
    else:
        chi2 = np.sum(((y-np.mean(y))/yerr)**2)
        return chi2/(N-1)
    
def calcChi2robust(y, yerr):
    """Compute simple robust reduced  chi2 (percentile-based) if more than 1 datapoints present
    """
    N = len(y)
    if N < 2:
        return np.nan
    else:
        Z = (y - np.median(y)) / yerr
        # calculating the robust interquartile range width  of the 
        # quantity Z : chi2R = calcSigmaG(Z)
        chi2R = 0.7413 * (np.percentile(Z,75)-np.percentile(Z,25))
        return chi2R

def calcMedian(y):
    ''' Calculate the median as 50-th percentile '''
    N = len(y)
    if N == 1 : 
        return float(y)    
    elif N == 0 : 
        return np.nan
    else : 
        return float(np.percentile(y,50))

def calcWeightedMean(y,yerr):
    ''' Calculate the weighted mean '''
    N = len(y)
    if N == 1 : 
        return float(y)    
    elif N == 0 : 
        return np.nan
    else: 
        # weights = 1 / (yerr ** 2.0)  
        # wMean = np.sum(weights * flux) / np.sum(weights)
        return float(np.add.reduce(y / (yerr * yerr)) / np.add.reduce((1/yerr)*(1/yerr)))

def calcWeightedMeanErr(yerr):
    ''' The error on the weighted mean '''
    N = len(yerr)
    if N == 1 : 
        return float(yerr)    
    elif N == 0 : 
        return np.nan
    else : 
        # weights = 1 / (yerr ** 2.0)
        # sigma_mean = 1 / sqrt(sum(weights))
        return np.sqrt(1. / np.add.reduce((1 / yerr)*(1/yerr)))


def calcWeightedPercentile(data, weights=None, percentile = 50):
    """Calculate the weighted percentile of an array
    Given a vector V, a q-th percentile of V is the value q/100 
    of the way from minimum to maximum of the sorted array V.  
    With weights , we find the q/100 point of the sorted weights,
    and take that value of the data. 
    """
    import numpy as np
    # if there are actually no weights, then pass it on to regular percentile 
    if weights is None:
        return np.percentile(np.array(data).flatten(), percentile)
    data, weights = np.array(data).flatten(), np.array(weights).flatten()
    if any(weights > 0):
        # 1) sort the data and weights: 
        # first : zip(data,weights) gives a list of tuples: [(data1, weight1),(data2,weight2), ()....]
        # second: sorted sorts these according to data : 
        # example: data = [  2,   5,   1,   7,    3]  ,  
        #       weights = [0.1, 0.2, 0.3, 0.1, 0.05]
        # sorted(zip(data,weights))  yields [(1, 0.3), (2, 0.1), (3, 0.05), (5, 0.2), (7, 0.1)] 
        sorted_data, sorted_weights = map(np.array, zip(*sorted(zip(data, weights))))
        
        # 2) calculate the fraction to which we need to go. For median fraction=0.5
        fraction =  percentile / 100.0
        midpoint = fraction * sum(sorted_weights)
        
        # 3) check here if there is any weight that is way bigger than the rest, i.e. 
        # given the sum of sorted_weights, is there any weight that contributes more than the 
        # fraction point q/100 ? If there is , then q/100 will be on the datapoint corresponding to 
        # that weight  
        # example : data = [  2,   5,   1,   7,    3] ,  
        #         weights = [0.1, 10, 0.3, 0.2, 0.05])
        # for 50-th percentile, midpoint = 5.32 , so we see that the datapoint 5 has such a large weight 
        # that indeed any(weights > midpoint) == True ,  so that we want to return 
        # the datapoint with that huge weight . 
        # Note : It does not happen if there are two or more points 
        # in a dataset with large weight, so this if statement is only catching the case of one datapoint
        # with weightmuch greater than other points 
        if any(weights > midpoint):
            return (data[weights == np.max(weights)])[0]
        
        # 4) if there isn't any datapoint with exceedingly large weight, 
        # then calculate the cumulative sum of weights, and find out the 
        # datapoint corresponding to the position just below the midpoint 
        cumulative_weight = np.cumsum(sorted_weights)
        below_midpoint_index = np.where(cumulative_weight <= midpoint)[0][-1]
        
        # 5) if this is almost the same, then return that value 
        if cumulative_weight[below_midpoint_index] - midpoint < sys.float_info.epsilon:
            return np.mean(sorted_data[below_midpoint_index:below_midpoint_index+2])
        return sorted_data[below_midpoint_index+1]
    


def calcWeightedStDev(y, yerr, yWmean):
    ''' Calculate the  weighted standard deviation
    Needs the weighted mean on y to be calculated 
    beforehand, using eg. calcWeightedMean()
    I'm using Bessel's correction to make it unbiased ...
    
    # calculate  weighted standard deviation corrected for intrinsic scatter 
    # using http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf
    # Yusra uses 1/N-1 instead of N/N-1.... calcWStdCorrAndMean
    # I'm pretty confused having read https://en.wikipedia.org/wiki/Bessel's_correction
    # and https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
    # after that http://stats.stackexchange.com/questions/6534/how-do-i-calculate-a-weighted-standard-deviation-in-excel 
    # I'm done. 


    # Okay, update  : using Bessel's correction ONLY makes sense with unweighted samples. 
    # Otherwise we are unnecessarily multiplying the standard deviation by 
    # sqrt(1/(N-1)).  
    # It is incorrect, because in the limit of identical errors.
    # we need to have std_w(x) --> std(x)
    # Fully explained in VariabilityReport ! 

    '''
    N = len(y)
    if N == 1:
        return float(yerr)
    elif N == 0:
        return np.nan 
    else :     
        weights=1.0 / ( yerr *yerr)
        return np.sqrt( (np.sum(weights * ((y - yWmean) ** 2.0)) / np.sum(weights)))  
        # used to be :
        # np.sqrt( (np.sum(weights * ((y - yWmean) ** 2.0)) / np.sum(weights)))  

def calcSigmaG(y):
    ''' Calculate the  interquartile sigma.
    Using an expression for a perfect Gaussian distribution
    from AstroML book, eq. (3.36) 
    '''
    N = len(y)
    if N == 1:
        return float(y)
    elif N == 0:
        return np.nan  
    else: 
        q25,q75 = np.percentile(y, (25,75))
        return 0.7413 * (q75-q25)

def computeVarMetrics(group):
    ''' Variability metrics to compute for each object on full lightcurve or 
    all points in a given season.  

    Takes a grouped object (outcome of pandas groupby : 
    grouped = fp_data.groupby('objectId')  

    Returns : 

    a pandas series with 

      'N' = number of points in  the lightcurve 
      'psfFluxMean':  Mean of Flux  ,
      'psfFluxMeanErr' : Mean of Flux Error , 
      'psfFluxMedian': Median of Flux , 
      'psfFluxMedianErr': np.sqrt(np.pi / 2)*FluxMeanErr,
      'psfFluxSkew' : group['psfFlux'].skew(),
      'psfFluxSigG' : calcSigmaG(Flux),
      'psfFluxStDev':psfFluxStDev,
      'psfFluxStDevUnw' : np.std(Flux),   # temporary : unweighted stdev of Flux 
      'chi2DOF' : calcChi2raw(Flux,FluxErr),
      'chi2R' : calcChi2robust(Flux,FluxErr),
      'sigmaFull' :sigma,
      'muFull' :mu,
      'avgMJD' : group['mjd'].mean(),
      'rangeMJD' : rangeMJD,
      'flagLtTenPts' : flagLtTenPts,
      'sigmaBootMax': sigmaBootMax ,  # temporary : maximum of the bootstrapped resamples 
      'sigmaMean' : stats['mean'],
      'sigmaStDev' : stats['stdev'],
      'sigmaMedian' : stats['median'],
      'sigmaSigmaG' : stats['sigmaG'],
      'sigmaProb2StDev' : stats['probInTwoStdev'],
      'pSigmaGaussLike' : stats['GaussLike']

    '''
    # print diagnostic for figuring out error...
    #print('objectId= %d'% group['objectId'].values[0])
    
    # even though I drop NaNs before, I do it here explicitly to save 
    # me from headaches 
    # eg. obj  216199180459189485   
    # have one row with NaN , not caught by other filters... 
    # and for some reason, can't use here  group.dropna(..., inplace=True) !   
    # group.dropna(subset=['psfFlux', 'psfFluxErr'], inplace=True)
    group = group.replace([np.inf, -np.inf, 0], np.nan)
    group = group.dropna(subset=['psfFlux', 'psfFluxErr'])
    
    
    # calculate range  of  dates in a given lightcurve    
    rangeMJD = group['mjd'].values.max() - group['mjd'].values.min() 
    
    # grab Flux and FluxErr values 
    Flux = group['psfFlux'].values
    FluxErr = group['psfFluxErr'].values
    
    # calculate Weighted  Mean
     
    FluxMean = calcWeightedMean(Flux,FluxErr)
    FluxMeanErr = calcWeightedMeanErr(FluxErr)
    psfFluxStDev = calcWeightedStDev(Flux,FluxErr, FluxMean)

    # calculate Median error : not necessary (since it's just const * meanErr)
    # medianErr = np.sqrt(np.pi / 2.0) * FluxMeanErr

    # note  : get_mu_sigma() gets digestion problems with NaN's - make sure 
    # that Flux and FluxErr are free of NaN's ! 
    # I multiply flux by a factor for stability issues    
    N = len(Flux)
    if N == 0 : 
        mu = np.nan
        sigma = np.nan
    elif N == 1  :
        mu, sigma = Flux[0], 0
    if N == 0 or N == 1 : 
        # if there are no points in the lightcurve, then 
        # initialize empty stats and return that ... 
        stats = {'mean':np.nan,'median':np.nan, 'stdev':np.nan, 'sigmaG':np.nan,
                 'probInTwoStdev':np.nan ,'GaussLike':np.nan}
        

    else : 
        # f = 1e27  # seems not needed since I solved the sigma normalization issue ... 
        #stats, mu, sigma = get_mu_sigma(Flux, FluxErr,1000, return_sigma_pdf_info = True)
        stats, mu, sigma = get_mu_sigma(Flux, FluxErr,1000, return_plot_data=False, return_sigma_pdf_info = True)

    # set the flag about length...
    if N > 10 : 
        flagLtTenPts = np.nan
    else:
        flagLtTenPts = 1 
   
    
    return pd.Series({'N':group['psfFlux'].count(),
                      'psfFluxMean': FluxMean,
                      'psfFluxMeanErr' : FluxMeanErr,
                      'psfFluxMedian': calcMedian(Flux),
                      'psfFluxMedianErr': np.sqrt(np.pi / 2)*FluxMeanErr,
                      'psfFluxSkew' : group['psfFlux'].skew(),
                      'psfFluxSigG' : calcSigmaG(Flux),
                      'psfFluxStDev':psfFluxStDev,
                      'psfFluxStDevUnw' : np.std(Flux),     # temporary : unweighted stdev of Flux 
                      'chi2DOF' : calcChi2raw(Flux,FluxErr),
                      'chi2R' : calcChi2robust(Flux,FluxErr),
                      'sigmaFull' :sigma,
                      'muFull' :mu,
                      'meanMJD' : group['mjd'].mean(),
                      'rangeMJD' : rangeMJD,
                      'flagLtTenPts' : flagLtTenPts,
                      'sigmaMean' : stats['mean'],
                      'sigmaStDev' : stats['stdev'],
                      'sigmaMedian' : stats['median'],
                      'sigmaSigmaG' : stats['sigmaG'],
                      'sigmaProb2StDev' : stats['probInTwoStdev'],
                      'pSigmaGaussLike' : stats['GaussLike']
                     })


def ComputeVarFullBinned(group): 
    
    ''' A function to calculate averages for the full lightcurve binned into seasons 
    '''    
    Flux= group['psfFluxMean'].values
    FluxErr =group['psfFluxMeanErr'].values

    N = len(Flux)
    if N == 0 : 
        mu = np.nan
        sigma = np.nan
    elif N == 1  :
        mu, sigma = Flux, 0
    else : 
        mu, sigma = get_mu_sigma(Flux*1e27, FluxErr*1e27, 1000)

    FluxMean = calcWeightedMean(Flux,FluxErr)
    FluxMeanErr = calcWeightedMeanErr(FluxErr)
    psfFluxStDev = calcWeightedStDev(Flux,FluxErr, FluxMean)

    return pd.Series({'Nseasons':group['psfFluxMean'].count(),
                      'psfFluxMeanMean': FluxMean, 
                      'psfFluxMeanMeanErr' :FluxMeanErr,
                      'chi2DOFmean' : calcChi2raw(Flux,FluxErr),
                      'chi2DOFmedian' : calcChi2raw(group['psfFluxMedian'].values,group['psfFluxMedianErr'].values),
                      'chi2Rmean' : calcChi2robust(Flux,FluxErr),
                      'chi2Rmedian' : calcChi2robust(group['psfFluxMedian'].values,group['psfFluxMedianErr'].values),
                      'sigmaFull' :sigma,
                      'muFull' :mu
                     })


# http://stackoverflow.com/questions/18603270/progress-indicator-during-pandas-operations-python


