# -*- coding: iso-8859-1 -*-
#
# Note: in comparison to processPatch.py : 
# 
# 1) here we use first AstroPy's supreme 
# Table to read the filter-patch file,
# since it is quicker than the pandas 
# read_csv.  This also enables quicker 
# processing. At the end we have to revert 
# to pandas to do group_by.apply(func)
# to calculate variability metrics.
# 2) At the level of processPatch()
# we calculate the metrics for 
# either full data, or just the 
# bright points (i.e. where S/N > 2 )
# In addition to processPatch.py,  
# we also calculate the meanSN per 
# lightcurve, to be able to split points
# into high SN and weak SN sample 
# 3) Finally, this version of processPatch 
# is compatible with the rewritten 
# variabilityFunctions.computeVarMetrics()
# 
# NEW kwargs:  
#- ability to specify flux_column='psfFlux' , 
#  error_column='psfFluxErr', 
#  time_column ='mjd', 
#- ability to say here whether we want to 
#  print objectIds or not with  verbose=True,
#- ability to decide whether we need 
#  intrinsic sigma AstroML 5.8 calculation
#  or not with  calc_sigma_pdf =True - adding 
#  this makes the calculation 3 times longer 
#  using %timeit  I find, assuming 71 rows per 
#  object,  490369 objects per filter-patch,
#  with 8.13 ms per lightcurve with 
#  calc_sigma_pdf =False , the  calculation
#  would take 1.1 hrs 
#  whereas with calc_sigma_pdf=True, it takes
#  25 ms / lightcurve, i.e. > 3 hrs per 
#  filter-patch. 
#
# My Creator, I am now willing that You should have all of me, 
# good and bad. I pray that You now remove from me 
# every single defect of character which stands in the way 
# of my usefulness to You and my fellows. 
# Grant me strength, as I go out from here, to do Your bidding.
# Amen
#
#


import pandas as pd
import numpy as np
import os 
import sys
import datetime
from astropy.time import Time
from astropy.table import Table
from astropy.table import Column


execution_environment = 'typhoon'  # or 'mac'
site = 'NCSA'


if execution_environment == 'mac' : 
    path_to_home = '/Users/chris/GradResearch/'

elif execution_environment == 'typhoon' :
    path_to_home = '/astro/users/suberlak/' 

sys.path.insert(0, path_to_home + 'SDSS_S82_FP_research/packages/')

pd.options.mode.chained_assignment = None

import faintFunctions as faintF 
import variabilityFunctions as varF
import imp


def flux2absigma(flux, fluxsigma):
    """Compute AB mag sigma given flux and flux sigma

    Here units of flux,  fluxsigma  don't matter 
    as long as they are consistent, since we are dividing 
    one by the other, so all the units cancel out.
    """
    FIVE_OVER_2LOG10 = 1.085736204758129569
    return FIVE_OVER_2LOG10 * fluxsigma / flux;


def flux2ab(flux, unit = 'Jy'):
    """Compute AB mag given flux. 

    Accept two unit types :  
    *  'cgs', meaning flux is in  ergs / s / Hz / cm2
    *  'Jy', meaning flux is in Jy.  1 Jy = 1E-23 * ergs/s/Hz/cm2
    """
    if unit == 'Jy':
        return -2.5 * np.log10(flux) + 8.90
    elif unit == 'cgs':
        return -2.5 * np.log10(flux) - 48.6


def process_patch_seasonally(name, DirIn, DirOut, pre='VarD_', 
    calc_sigma_pdf=False, limitNrows=None, verbose = None ):
    ''' Code to  calculate aggregate metrics on seasonally - binned 
    forced photometry light curves from Stripe 82. It is assumed 
    that the data is stored in filter-patch files, i.e. we process 
    only one filter at a time, with data for many objectIds, at the 
    same time. 

    Detailed step-by-step workflow : 
    --------------------------------
    - start with psfFlux , psfFluxErr,  and convert to Jansky, 
      remove all NaN or missing rows... 
    - aggregate by objectId, and within each aggregate,  define 
    in which season is each epochal forced photometry measurement
    - aggregate by objectId and season, and calculate seasonal 
      averages  : psfFluxSMean,  psfFluxSMeanErr,  
      psfFluxSMedian ,  psfFluxSMedianErr , using  
      variabilityFunctions.calcWeightedMean(), 
      variabilityFunctions.calcWeightedMeanErr(),
      variabilityFunctions.calcMedian() , 
      as well as chi2DOF,  chi2R, season start date, season end date,
      Number of points per season ... 
    - treating one of  {psfFluxSMean, psfFluxSMedian } as signal S, 
      and {psfFluxSMeanErr,psfFluxSMedianErr } as noise N,  
      calculate S / N ,  and if  S/N < 2  : perform faint point 
      pipeline, and then exchange S as faintMean, and N as faintRMS
    - now that we've ensured that each seasonal average is non-negative,
      calculate seasonal magnitudes : 
      psfFluxSMean --> psfSMean ,  
      psfFluxSMeanErr --> psfSMeanErr ,  etc.  
    - if for any season the psfSMeanErr < 0.003 mag, update fluxError : 
      fluxErrNew = sqrt(fluxErr^2 + (0.003 * flux) ^2 )
    - aggregating by objectId,  calculate statistics using N~4-5 points
      (seasonal averages)

    Parameters :
    ---------------------
    name : a string : name of the patch file, assumed to be of form  
           g00_21.csv.gz
    DirIn: a string : directory storing the raw S82 forced photometry 
           lightcurves
    DirOut : a string : name of output  directory where we save 
           aggregate information 
    limitNrows : integer number of rows to process from the patch file, 
          if not the  full file  . By default limitNrows=None, and 
          we calculate metrics for all lightcurves per patch file  
    verbose : boolean, if not  None, will print some extended 
          diagnostic information. 
    
    Returns:
    ----------
    None (there is nothing that the program explicitly returns - all 
          output is saved as text files in a specified directory ).


    '''

    # Note :  1.1 to 1.3 are exactly the same as in process_patch()
    # (reading  in the patch-file, converting to Janskys, removing 
    # all bad data... )


    print('\n Processing filter_patch file %s' % name)
    
    # read in the raw lightcurve... 
    # NOTE : if nrows = None , pd.read_csv() reads the entire file 
    raw_data_df = pd.read_csv(DirIn+name+'.gz', compression='gzip',  
                 usecols=['objectId', 'mjd', 'psfFlux', 'psfFluxErr'], 
                 nrows=limitNrows )
    raw_data = Table.from_pandas(raw_data_df)
 
    ##########  STEP 1 : single-epoch data ###########  
    #         1 Jy = 1.0E-26 W/m^2/Hz = 1.0E-23 erg/s/cm^2/Hz
    # 1.1  :  convert Flux from erg/cm2/sec/Hz  to Jansky
    # make new columns with data in Jy, to follow good 
    # programming practice (otherwise, changes are either
    # on a view or a copy and it's not always clear)

    raw_data['psfFluxJy'] = raw_data['psfFlux'] * 1E23 
    raw_data['psfFluxErrJy'] = raw_data['psfFluxErr'] * 1E23 


    # 1.2  : drop all rows which have NaNs in psfFlux .... 
    m1  = np.isnan(raw_data['psfFlux'].data)  # true if NaN 
     # true if not finite... 
    m2 = np.bitwise_not(np.isfinite(raw_data['psfFlux'].data))  


    # logical or : true if either condition satisfied 
    m = m1 | m2  

    if np.sum(m) > 0 :  # only apply if there is anything to drop ... 
        print('In  file %s there are : \n.... %d rows where psfFlux \
            is NaN'%(name, np.sum(m1)))
        print('.... %d rows where psfFlux is not finite  '% np.sum(m2))
        print('All such rows are dropped')
        indices = np.arange(len(raw_data))
        remove_rows= indices[m]
        raw_data.remove_rows(remove_rows)

    # 1.3 : check psfFluxErr : drop all rows which have NaN, or 0 ,
    # to avoid getting error when calculating S/N  = psfFlux / 
    # psfFluxErr 

    m1  = np.isnan(raw_data['psfFluxErr'].data)  # true if NaN 
    # true if not finite... 
    m2 =  np.bitwise_not(np.isfinite(raw_data['psfFluxErr'].data))  
    m3 =  raw_data['psfFluxErr'].data == 0 
    # logical or : true if either condition satisfied 
    m = m1 | m2  | m3

    if np.sum(m) > 0 :  # only apply if there is anything to drop ... 
        print('In  file %s there are : \n.... %d rows where psfFluxErr \
            is NaN'%(name, np.sum(m1)))
        print('.... %d rows where psfFluxErr is not finite'% 
               np.sum(m2))
        print('.... %d rows where psfFluxErr = 0 '% np.sum(m3))
        print('All such rows are dropped')
        indices = np.arange(len(raw_data))
        remove_rows= indices[m]
        raw_data.remove_rows(remove_rows)



    
    # 1.3 : check psfFluxErr : drop all rows which have NaN, or 0 ,
    # to avoid getting error when calculating S/N  = psfFlux / 
    # psfFluxErr 

    m1  = np.isnan(raw_data['psfFluxErr'].data)  # true if NaN 
    # true if not finite... 
    m2 =  np.bitwise_not(np.isfinite(raw_data['psfFluxErr'].data))  
    m3 =  raw_data['psfFluxErr'].data == 0 
    # logical or : true if either condition satisfied 
    m = m1 | m2  | m3

    if np.sum(m) > 0 :  # only apply if there is anything to drop ... 
        print('In  file %s there are : \n.... %d rows where psfFluxErr \
            is NaN'%(name, np.sum(m1)))
        print('.... %d rows where psfFluxErr is not finite'% 
               np.sum(m2))
        print('.... %d rows where psfFluxErr = 0 '% np.sum(m3))
        print('All such rows are dropped')
        indices = np.arange(len(raw_data))
        remove_rows= indices[m]
        raw_data.remove_rows(remove_rows)


    # 2.1 :  make a new column to designate seasons...
    raw_data['season'] = np.nan

    # I make a list of boundaries between season start / end : 
    # season_bounds
    # first it's a list with 2005-08-01 , 2006-08-01, etc...
    season_bounds = [str(year)+'-08-01 00:00:00.000' for year in 
                                               np.arange(2005,2009)]
    # Then I insert a very early date at the beginning of the list, 
    season_bounds.insert(0,'1990-08-01 00:00:00.000' )
    # so that all epochs between
    # 1990 and 2005 are clustered together.
    # Thus 1990-08-01 - 2005-08-01  is Season1
    # 2005-2006  : Season2 ;  2006-2007 : Season3, etc.  
    # raw photometry is averaged here WITHIN each Season. 

    # Use AstroPy Time module
    cutDates = Time(season_bounds, format='iso')
    seasons = np.arange(len(cutDates))+0 
        
    # Assign value of a season for each row...
    for i in range(len(cutDates.mjd)-1):
        mask = (raw_data['mjd'].data > cutDates[i].mjd) * \
               (raw_data['mjd'].data < cutDates[i+1].mjd)
        raw_data['season'][mask]  = seasons[i]  
    # this is better than raw_data[mask]['season'] : 
    # somehow the order  matters 

    # 2.2 Calculate seasonal  aggregates 
    # At the moment, it is impossible to apply complicated 
    # functions to groupby objects....
    # Easier at this point to just convert AstroPy table 
    # to Pandas,  and do it from here on in Pandas ...

    # But it only takes a few seconds 
    # even on mac,  on the full patch  file 
    raw_data_df = raw_data.to_pandas()

    # calculate seasonal averages 
    grouped = raw_data_df.groupby(['objectId','season'])

    # calculate averages within each season
    # since these are aggregates within each
    # season, call it by that name to make 
    # all things easier... 

    season_agg= grouped.apply(varF.computeVarMetrics, 
                                     flux_column='psfFluxJy',
                                     error_column = 'psfFluxErrJy',
                                     time_column = 'mjd', 
                                     calc_sigma_pdf = False, 
                                     verbose=False,
                                     seasonal_average = True)
    #
    # returns a pd.Series with the following columns :
    #
    # 'objectId', 'season', 'N', 'chi2DOF', 'chi2R', 'flagLtTenPts',
    # 'maxMJD', 'meanMJD', 'meanSN', 'minMJD', 'psfFluxMean', 
    # 'psfFluxMeanErr',  'psfFluxMedian', 'psfFluxMedianErr', 
    # 'psfFluxSigG',  'psfFluxSkew', 'psfFluxStDev'

    # 2.3 faint Pipeline for S/N  <2 
    # Calculate S/N : first using mean, and then median, and in 
    # each case if they
    # are below 2,  replace by faintMean, faintMedian, faintRMS 

    # read in the lookup interpolation table, 
    # which makes it very fast thanks to the 
    # smoothness of functions used  
    table_address = 'flat_prior_const_counts_lookup.csv'
    lookup = Table.read(table_address)


    # 3: calculate faint quantities for all rows  where S/N < 2 
    for avg in ['Mean', 'Median']:
        
        # temporary aliases to make it easier 
        S = season_agg['psfFlux'+avg].values
        N = season_agg['psfFlux'+avg+'Err'].values
        SN = S/N
        mask_SN = SN < 2 
        # make a new column storing that mask : 
        # True if SN < 2 , otherwise False 
        season_agg['flagFaint'+avg] = mask_SN
        print('There are %d seasons of %d that have %s-based  S/N < 2' %
            (np.sum(mask_SN),len(mask_SN), avg))

        # interpolate at the xObs = SN  locations .... 
        faint_SN = SN[mask_SN]
        xMean =  np.interp(faint_SN,  lookup['xObs'].data, 
                                     lookup['xMean'].data)
        xMed =  np.interp(faint_SN,   lookup['xObs'].data, 
                                   lookup['xMedian'].data)
        xSigma =  np.interp(faint_SN, lookup['xObs'].data, 
                                    lookup['xSigma'].data)
        xRMS = np.interp(faint_SN,    lookup['xObs'].data, 
                                      lookup['xRMS'].data)
        
        # need to multiply by error to find flux 
        faintMean = xMean * N[mask_SN]
        faintRMS = xRMS * N[mask_SN]

        # Not sure if need to save these values for anything ... 
        # can always do it later,
        # but for now, replace the 
        # Signal  by faintMean , 
        # and Noise by faintRMS 
        season_agg['psfFlux'+avg][mask_SN] = faintMean
        season_agg['psfFlux'+avg+'Err'][mask_SN] = faintRMS

    # 3.1 Calculate seasonal magnitudes 
    # add in quadrature 0.003 mag to 
    # those magErrs which are < 0.003 mag
    # and update fluxErrs for the same 
    # measurements using the fact that 
    # for small errors, (<0.2 mag or so), 
    # they are the same as fractional
    # error in flux. So, given flux and fluxErr, 
    # we redefine error to be
    # fluxErrNew  = sqrt[ fluxErr^2 + (0.003*flux)^2 ]

    for avg in ['Mean', 'Median'] : 
        fCol = 'psfFlux'+avg
        fErrCol = 'psfFlux'+avg+'Err'
        mCol = 'psf'+avg
        mErrCol = 'psf'+avg+'Err' 

        season_agg[mCol]= \
            procP.flux2ab(season_agg[fCol], unit='Jy')    
        season_agg[mErrCol] = \
            procP.flux2absigma(season_agg[fCol], season_agg[fErrCol])    

        # Wherever the error in magnitudes is smaller than the 
        # minimum value, 
        # add that minimum value
        # How many seasons have Mean Err < min_err ? 
        e_min = 0.003
        mask_err = season_agg[mErrCol].values < e_min
        print('Number of seasons with error < %.3f is %d'%
              (e_min,np.sum(mask_err)))

        # add to magnitude errors min_err for those very small ones....
        e_old = season_agg[mErrCol].values[mask_err]
        e_new = np.sqrt(e_min**2.0 + e_old**2.0 )

        # replace the originals
        season_agg[mErrCol].values[mask_err] = e_new

        # old flux error for these few seasons where 
        # magnitude error < 0.003  
        f_err_old =  season_agg[fErrCol].values[mask_err]
        f = season_agg[fCol].values[mask_err]

        # new flux error 
        f_err_new =  np.sqrt(f_err_old**2.0 + (e_min * f)**2.0)

        # replace the originals
        season_agg[fErrCol].values[mask_err] = f_err_new


    # 3.2 need to change the multindex to columns, to be able to 
    # aggregate again...
    season_agg.reset_index(inplace=True)  

    #### Save the seasonal averages
    path = DirOut+ 'S_'+name
    season_agg.to_csv(path)

    # 4 Aggregate by objectId
    grouped = season_agg.groupby('objectId')

    # calculate stats using seasonal averages... 
    # do it only using Mean... 
    # next time around can also do that using Median, and 
    # hstack the two ... 

    # the problem here is that we have 
    # mean-based and median-based seasonal averages,
    # and then we calculate averages of that 
    # so in a way we'd have mean of mean, median of mean,
    # mean of median,  median of median,
    # depending on what we choose to represent
    # seasonal fluxes below as 
    # flux_column  and error_column
    varMetricsSeasons  = grouped.apply(varF.computeVarMetrics, 
                                     flux_column='psfFluxMean',
                                     error_column = 'psfFluxMeanErr',
                                     time_column = 'meanMJD', 
                                     calc_sigma_pdf = False, 
                                     verbose=False,
                                     seasonal_average = False)


    # need to rename the columns to reflect the fact that we used 
    # seasonal Mean-based averages,
    # as opposed to Median-based averages 
    # Or just use mean-based averages for now 

    # calculate magnitudes for full light curve aggregates based 
    # on seasonal averages 
    varMetricsSeasons['psfMean'] = procP.flux2ab(
                                    varMetricsSeasons['psfFluxMean'], 
                                                           unit='Jy')
    varMetricsSeasons['psfMeanErr'] = procP.flux2absigma(
                                    varMetricsSeasons['psfFluxMean'], 
                                 varMetricsSeasons['psfFluxMeanErr'])
    #### Save the full LC Metrics based on averaged seasons .... 
    path = DirOut+ 'VarS_'+name
    varMetricsSeasons.to_csv(path)    

    return

def process_patch(name, DirIn, DirOut, pre='VarD_', 
    calc_sigma_pdf=False, limitNrows=None, calc_bright_part = False, 
    verbose = None ):
    '''  A code to perform our basic processing on the raw 
    Forced Photometry data from Stripe 82, 
    performed patch-by-patch.  One clone of data lives in 
    https://lsst-web.ncsa.illinois.edu/~yusra/S13Agg/rawDataFPSplit/
    and another in 
    /astro/store/pogo4/....
    
    Detailed step-by-step workflow : 
    --------------------------------
    - drop all  rows that have either psfFlux or psfFluxErr as NaN or inf 
    - make a column 'flagFaint' ==1 if  S/N < 2 for a given point 
    - one these 'faint' points, recalculate faintMean, faintMedian, 
      faintTwoSigma, faintRMS,  using faintFunctions.py  code 
    - replace all psfFlux  where S/N < 2  with  faintMean 
    - group by 'objectId', and calculate variability metrics using 
      variabilityFunctions.computeVarMetrics().  I call the data product:
      ( N_objects x  X_columns  )  varMetricsFull , since these metrics
      are computed on full lightcurves 
    - calculate magnitudes on lightcurve aggregate (varMetricsFull) : 
       psfMean  psfMedian  psfMeanErr  psfMedianErr 

    Parameters :
    ---------------------
    name : a string : name of the patch file, assumed to be of form  
           g00_21.csv.gz
    DirIn: a string : directory storing the raw S82 forced photometry 
           lightcurves
    DirOut : a string : name of output  directory where we save 
           aggregate information 
    limitNrows : integer number of rows to process from the patch file, 
          if not the  full file  . By default limitNrows=None, and 
          we calculate metrics for all lightcurves per patch file  
    calc_bright_part : boolean, if False, we use per light curve 
          all points, regardless of their signal-to-noise ratio. Otherwise
          we also calculate all metrics per light curve using only 
          bright points (with S/N > 2)
    verbose : boolean, if not  None, will print some extended 
          diagnostic information. 
    
    Returns:
    ----------
    None (there is nothing that the program explicitly returns - all 
          output is saved as text files in a specified directory ).
    '''

    print('\n Processing filter_patch file %s' % name)
    
    # read in the raw lightcurve... 
    # NOTE : if nrows = None , pd.read_csv() reads the entire file 
    raw_data_df = pd.read_csv(DirIn+name+'.gz', compression='gzip',  
                 usecols=['objectId', 'mjd', 'psfFlux', 'psfFluxErr'], 
                 nrows=limitNrows )
    raw_data = Table.from_pandas(raw_data_df)
 
    ##########  STEP 1 : single-epoch data ###########  
    #         1 Jy = 1.0E-26 W/m^2/Hz = 1.0E-23 erg/s/cm^2/Hz
    # 1.1  :  convert Flux from erg/cm2/sec/Hz  to Jansky
    # make new columns with data in Jy, to follow good 
    # programming practice (otherwise, changes are either
    # on a view or a copy and it's not always clear)

    raw_data['psfFluxJy'] = raw_data['psfFlux'] * 1E23 
    raw_data['psfFluxErrJy'] = raw_data['psfFluxErr'] * 1E23 


    # 1.2  : drop all rows which have NaNs in psfFlux .... 
    m1  = np.isnan(raw_data['psfFlux'].data)  # true if NaN 
     # true if not finite... 
    m2 = np.bitwise_not(np.isfinite(raw_data['psfFlux'].data))  


    # logical or : true if either condition satisfied 
    m = m1 | m2  

    if np.sum(m) > 0 :  # only apply if there is anything to drop ... 
        print('In  file %s there are : \n.... %d rows where psfFlux \
            is NaN'%(name, np.sum(m1)))
        print('.... %d rows where psfFlux is not finite  '% np.sum(m2))
        print('All such rows are dropped')
        indices = np.arange(len(raw_data))
        remove_rows= indices[m]
        raw_data.remove_rows(remove_rows)

    # 1.3 : check psfFluxErr : drop all rows which have NaN, or 0 ,
    # to avoid getting error when calculating S/N  = psfFlux / psfFluxErr 

    m1  = np.isnan(raw_data['psfFluxErr'].data)  # true if NaN 
    # true if not finite... 
    m2 =  np.bitwise_not(np.isfinite(raw_data['psfFluxErr'].data))  
    m3 =  raw_data['psfFluxErr'].data == 0 
    # logical or : true if either condition satisfied 
    m = m1 | m2  | m3

    if np.sum(m) > 0 :  # only apply if there is anything to drop ... 
        print('In  file %s there are : \n.... %d rows where psfFluxErr \
            is NaN'%(name, np.sum(m1)))
        print('.... %d rows where psfFluxErr is not finite'% np.sum(m2))
        print('.... %d rows where psfFluxErr = 0 '% np.sum(m3))
        print('All such rows are dropped')
        indices = np.arange(len(raw_data))
        remove_rows= indices[m]
        raw_data.remove_rows(remove_rows)



    # 1.4 : select points that have S/N  < 2 , flag as Faint...
    # initialize a new column with all values set to False :
    #flagFaint = Column(data = np.zeros(len(raw_data), dtype=bool), 
    # name = 'flagFaint')
    # this would also work, but doesn't allow to select good data type 
    # raw_data['flagFaint'] = True : makes a new col with int64

    # instead of two-steps, make just one : 
    SN = raw_data['psfFluxJy'].data / raw_data['psfFluxErrJy'].data
    mask_SN = SN < 2 
    raw_data['flagFaint'] = mask_SN
    print('There are %d points of %d that have S/N < 2' %(np.sum(mask_SN),\
        len(mask_SN)))


    # 1.5  calculate faint quantities for all rows  where S/N < 2 

    # make new columns ...
    for faintColname in ['faintMean', 'faintMedian','faintTwoSigma',
    'faintRMS']:
        raw_data[faintColname] = np.nan
        
    # temporary assignment 
    flux, flux_err = raw_data['psfFluxJy'][mask_SN].data, \
      raw_data['psfFluxErrJy'][mask_SN].data
    

    # new way of doing it with the lookup table - and very fast thanks 
    # to the smoothness of the employed functions ! 

    table_address = 'flat_prior_const_counts_lookup.csv'
    lookup = Table.read(table_address)

    # interpolate at the xObs = SN  locations .... 
    faint_SN = SN[mask_SN]
    xMean =  np.interp(faint_SN, lookup['xObs'].data, 
                        lookup['xMean'].data)
    xMed =  np.interp(faint_SN,lookup['xObs'].data, 
                        lookup['xMedian'].data)
    xSigma =  np.interp(faint_SN,lookup['xObs'].data,
                        lookup['xSigma'].data)
    xRMS = np.interp(faint_SN, lookup['xObs'].data, 
                        lookup['xRMS'].data)

    # add this info to the table 
    raw_data['faintMean'][mask_SN]     = xMean * flux_err
    raw_data['faintMedian'][mask_SN]   = xMed * flux_err
    raw_data['faintTwoSigma'][mask_SN] = xSigma * flux_err
    raw_data['faintRMS'][mask_SN]      = xRMS * flux_err



    # 1.6  replace psfFluxJy  by  faintMean ,  psfFluxErrJy  by  
    # faintRMS 
    # make sure we are only taking values  that are not NaN ... 
    #mask_NaNs = np.bitwise_not(np.isnan(raw_data['faintRMS']
    # [mask_SN]))

    raw_data['psfFluxJy'][mask_SN] = \
                 raw_data['faintMean'][mask_SN].data.data
    raw_data['psfFluxErrJy'][mask_SN] = \
                  raw_data['faintRMS'][mask_SN].data.data
 
    #
    ##########  STEP 2 : Derived Quantities ###########  
    #
    # Here the main choice is whether we want to introduce 
    # selection of which points to use per lightcurve based 
    # on S/N before grouping, or after grouping . 
    #

    # At the moment, it is impossible to apply complicated 
    # functions to groupby objects....
    # Easier at this point to just convert AstroPy table 
    # to Pandas,  and do it from here on in Pandas ...
    # Bummer. 

    # But thankfully it only takes a few seconds 
    # even on mac,  on the full patch  file 
    raw_data_df = raw_data.to_pandas()

    if calc_bright_part : 
        # 2.2 Calculate stats for LC using only bright points 
        print('Calculating the  LC statistics \
            using S/N > 2  points only ...')
        mask_bright = \
         np.bitwise_not(raw_data_df['flagFaint'].values.astype(bool))
        bright_grouped = raw_data_df[mask_bright].groupby('objectId')
        varMetricsFull_bright  = \
            bright_grouped.apply(varF.computeVarMetrics, 
                                    flux_column='psfFluxJy',
                                    error_column = 'psfFluxErrJy',
                                    time_column = 'mjd', 
                                    calc_sigma_pdf =calc_sigma_pdf)
        # otherwise just use all points...

    # 2.3 Calculate stats for LC using all points 
    print('Calculating the  LC statistics using all S/N  points  ...')
    all_grouped = raw_data_df.groupby('objectId')
    varMetricsFull_all  = all_grouped.apply(varF.computeVarMetrics, 
                                    flux_column='psfFluxJy',
                                    error_column = 'psfFluxErrJy',
                                    time_column = 'mjd',
                                    calc_sigma_pdf =calc_sigma_pdf) 


    # 2.4  calculate magnitudes from averaged fluxes ... 
   

    # Calculate magnitudes based on average fluxes :
    # psfMean  psfMedian  psfMeanErr  psfMedianErr 
    varMetricsFull_all['psfMean'] = \
            flux2ab(varMetricsFull_all['psfFluxMean'], unit='Jy')
    varMetricsFull_all['psfMedian'] = \
            flux2ab(varMetricsFull_all['psfFluxMedian'], unit='Jy')
    varMetricsFull_all['psfMeanErr'] = \
            flux2absigma(varMetricsFull_all['psfFluxMean'],
                         varMetricsFull_all['psfFluxMeanErr'])
    varMetricsFull_all['psfMedianErr'] = \
            flux2absigma(varMetricsFull_all['psfFluxMedian'],
                         varMetricsFull_all['psfFluxMedianErr'])

    if calc_bright_part : 
        varMetricsFull_bright['psfMean'] = \
               flux2ab(varMetricsFull_bright['psfFluxMean'], unit='Jy')
        varMetricsFull_bright['psfMedian'] = \
               flux2ab(varMetricsFull_bright['psfFluxMedian'], 
                           unit='Jy')
        varMetricsFull_bright['psfMeanErr'] = \
               flux2absigma(varMetricsFull_bright['psfFluxMean'],
                            varMetricsFull_bright['psfFluxMeanErr'])
        varMetricsFull_bright['psfMedianErr'] = \
               flux2absigma(varMetricsFull_bright['psfFluxMedian'],
                            varMetricsFull_bright['psfFluxMedianErr'])

    print('Calculating magnitudes from fluxes is finished')


    # 2.5  change colnames to reflect which subset of points per 
    # lightcurve
    # was used the easiest way to do it is to add_suffix in pandas 

    varMetricsFull_all = varMetricsFull_all.add_suffix('_all')

    if calc_bright_part : 
        varMetricsFull_bright = varMetricsFull_bright.add_suffix('_bright')
        # 2.6 combine the two ... 
        varMetricsFull_combined = pd.concat([varMetricsFull_all,
                                             varMetricsFull_bright], axis=1)
    else:
        varMetricsFull_combined  = varMetricsFull_all

    ######################### SAVING OUTPUT        ######################## 
    # 
    path = DirOut + 'VarD_'+name  
    print('Saving varMetricsFull to  %s '%path)
    varMetricsFull_combined.to_csv(path)

    print('The output dataframe header is : ')
    print(varMetricsFull_combined.head())
    return 


