# -*- coding: iso-8859-1 -*-
#
# Functions used  for the analysis in frequency  domain
#

import numpy as np 
 
def frequency_grid(times, method = 'AstroML'):
    ''' A function to determine the frequency grid for Lomb Scargle analysis. 
    
    I set the maximum frequency based on the average Nyquist factor. 
    I set the minimum frequency based on the size of the baseline. For 
    this reason, this is equivalent to setting in Astropy 
    minimum_frequency=None, maximum_frequency=None .  

    Note for the astropy option (taken from the astropy documentation) : 
     

        "Note that this assumes the peak width is driven by the observational
        baseline, which is generally a good assumption when the baseline is
        much larger than the oscillation period.
        If you are searching for periods longer than the baseline of your
        observations, this may not perform well.

        Even with a large baseline, be aware that the maximum frequency
        returned is based on the concept of "average Nyquist frequency", which
        may not be useful for irregularly-sampled data. The maximum frequency
        can be adjusted via the nyquist_factor argument, or through the
        maximum_frequency argument.""
        
        (source: http://docs.astropy.org/en/stable/_modules/astropy/stats/lombscargle/
                      core.html#LombScargle.autopower
         accessed: 11/27/16, 7:20 pm )

     Parameters
        ----------
        times : array (required)
             An array of observational times, used to calculate the frequency grid 

        method : which method to use when calculating the frequency grid . Possible options:
                 'AstroPy' , 'AstroML',  'EBellm'

        astropy : boolean (optional, default = True )
             If set, then use  the astropy method of finding the number of frequency 
             bins. Otherwise, use the method from Ivezic+2014, pg. 436, linking frequency 
             step in proportion to the minimum frequency : step_size * n = minimum_frequency 
        n : float (optional, default = 10 )
             If astropy = False, then when using Ivezic+2014 method , this sets the 
             frequency step size .  n = 1 / eta , if compared to Ivezic+2014, pg. 436.
             E. Bellm (iPTF Summer School 2016) used n=5 . 
             Ivezic+2014 suggests n=10  (eta = 1 / 10) (pg. 436) 
        samples_per_peak : float (optional, default=5)
            The approximate number of desired samples across the typical peak
        nyquist_factor : float (optional, default=5)
            The multiple of the average nyquist frequency used to choose the
            maximum frequency if maximum_frequency is not provided.

        Returns
        -------
        frequency : ndarray 
            The frequency array: linearly spaced from minimum_frequency, to 
            maximum_frequency, with n_bins  . 
    '''

    times = np.sort(times)
    baseline = abs(min(times) - max(times))
    
    if method == 'AstroPy' : 
        N = len(times)
        S = 5
        yF = 5
        f_min = 1.0 / (2 * S * baseline)
        f_max = f_min + ( yF * N ) / (2 * baseline)
        N_bins = 0.5 * N * S * yF 

    if method == 'EBellm'  : 

        delta_t  = np.zeros(len(times)-1) 
        for j in range(len(times)-1):
            delta_t[j] = (times[j+1]-times[j])

        f_min = 1.0 / baseline
        f_max = 0.5 / np.median(delta_t)
        N_bins = 2.5 * (baseline / np.median(delta_t)) - 5 
    
    if method == 'AstroML' : 

        delta_t  = np.zeros(len(times)-1) 
        for j in range(len(times)-1):
            delta_t[j] = (times[j+1]-times[j])

        f_min = 1.0 / baseline
        f_max = 0.5 * np.median(1.0 / delta_t)
        N_bins = 5 * np.median(1.0 / delta_t) * baseline - 10 
    
    # increase maximum frequency for small N points 
    # using the lookup table.  
    if method == 'AstroML' or method == 'EBellm':
        N = len(times)
        if N < 30 : 
            a = - 0.16
            b = 5.8
            f_max_corr= a * np.arange(30) + b 
            f_max = f_max_corr[N] * f_max

    # say all in angular frequency 
    omega_min = 2 * np.pi * f_min
    omega_max = 2 * np.pi * f_max 
    
    return np.linspace(omega_min, omega_max, N_bins)
 