# -*- coding: iso-8859-1 -*-
import pandas as pd
import numpy as np
import os 
import sys

execution_environment = 'typhoon'  # or 'mac'
site = 'NCSA'


if execution_environment == 'mac' : 
    path_to_home = '/Users/chris/GradResearch/'

elif execution_environment == 'typhoon' :
    path_to_home = '/astro/users/suberlak/' 

sys.path.insert(0, path_to_home + 'SDSS_S82_FP_research/packages/')

pd.options.mode.chained_assignment = None

# God, Relieve me of the bondage of self,
# that I may better do Thy will.
# Take away my difficulties,
# that victory over them may bear witness
# to those I would help of Thy Power,
# Thy Love, and Thy Way of life.
# May I do Thy will always!

# faint source treatment 
import faintFunctions as faintF 

# variability 
#print('I see that line...') 
import variabilityFunctions as varF
import imp
imp.reload(varF)
from astropy.time import Time

def process_patch(name, DirIn, DirOut, limitNrows=None, calc_seasonal_metrics=None, calc_seas_binned_metrics=None ,
                  verbose = None ):
    '''  A code to perform our basic processing on the raw 
    Forced Photometry data from Stripe 82, 
    performed patch-by-patch.  One clone of data lives in 
    https://lsst-web.ncsa.illinois.edu/~yusra/S13Agg/rawDataFPSplit/
    and another in 
    /astro/store/pogo4/....
    
    - drop all  rows that have either psfFlux or psfFluxErr as NaN or inf 
    - make a column 'flagFaint' ==1 if  S/N < 2 for a given point 
    - one these 'faint' points, recalculate faintMean, faintMedian, 
      faintTwoSigma, faintRMS,  using faintFunctions.py  code 
    - replace all psfFlux  where S/N < 2  with  faintMean 
    - group by 'objectId', and calculate variability metrics using 
      variabilityFunctions.computeVarMetrics().  I call the data product:
      ( N_objects x  X_columns  )  varMetricsFull , since these metrics
      are computed on full lightcurves (not seasonally binned or averaged 
      in any fashion). X_columns include [''
    - calculate magnitudes on lightcurve aggregate (varMetricsFull) : 
       psfMean  psfMedian  psfMeanErr  psfMedianErr 
    - identify variable candidates, based on lightcurve chi2dof , chi2robust, 
       sigmaFull . The criteria are  : 
       m1 = sigmaFull > 0
       m4 = chi2dof > (1 + 3.0 * np.sqrt(2 / N )
       m5 = chi2robust > (1 + 3.0 * np.sqrt(2 / N)) 
       m = m1 * (m4 | m5 )    
    - select variable candidates, and identify seasons.  
    - group the variable candidates by 'objectId', 'season'
    - calculate variability metrics using  variabilityFunctions.computeVarMetrics()
    (same as for full lightcurves),  but on seasonally-binned lightcurves .
    Return seasonal variability metrics varMetricsSeasonal

    varMetricsFull  and varMetricsSeasonal are saved as csv files, and so 
    there is nothing that this function explicitly returns.  

    '''

    print('\n Processing filter_patch file %s' % name)
    
    # read in the raw lightcurve... 
    if limitNrows is not None:
        fp_data = pd.read_csv(DirIn+name+'.gz', compression='gzip',  
                     usecols=['objectId', 'mjd', 'psfFlux', 'psfFluxErr'], nrows=limitNrows)
    else : 
        fp_data = pd.read_csv(DirIn+name+'.gz', compression='gzip',  
                     usecols=['objectId', 'mjd', 'psfFlux', 'psfFluxErr'])
    #   
    ##########  STEP 1 : single-epoch data ###########  
    #
    ### at the very beginnig : 
    # convert  Flux from erg/cm2/sec/Hz  to Jansky
    # 1 Jy = 1.0E-26 W/m^2/Hz = 1.0E-23 erg/s/cm^2/Hz

    fp_data['psfFlux'] = fp_data['psfFlux']  * 1E23 
    fp_data['psfFluxErr'] = fp_data['psfFluxErr'] * 1E23 

    ####  first drop all NaNs  in psfFlux...      
    m1  = np.isnan(fp_data['psfFlux'])  # True if NaN  
    m2 =  ~np.isfinite(fp_data['psfFlux']) #  True if not finite  
    m  = m1 | m2  # a logical or 
    if np.sum(m) > 0 :  # only apply if there is anything to drop ... 
        fp_data.drop(m.index[m], inplace=True)  # drop entire rows 
        print('Okay, we dropped %d rows where psfFlux is NaN or inf'%np.sum(m))

    #### check to make sure that there are no NaN or 0 psfFluxErr... 
    m1  = np.isnan(fp_data['psfFluxErr'])  # True if NaN  
    m2 =  ~np.isfinite(fp_data['psfFluxErr']) #  True if not finite
    m3 =   fp_data['psfFluxErr'].values == 0 # True if Err == 0  (IN2P3 problem...)
    m  = m1 | m2 | m3  # a logical or 
    if np.sum(m) > 0 :  # only apply if there is anything to drop ... 
        fp_data.drop(m.index[m], inplace=True)
        print('Okay, we dropped %d rows where psfFluxErr is NaN or inf'%np.sum(m))

    # make a new column, fill with 0's
    fp_data['flagFaint'] = 0

    # mask those rows that correspond to SNR < 2
    mask = (fp_data['psfFlux'].values / fp_data['psfFluxErr'].values) < 2

    # print info how many points are affected
    print('There are %d points of %d that have SNR<2' %(np.sum(mask),len(mask)))

    # set flag at those rows to 1
    fp_data.ix[mask, 'flagFaint'] = 1

    # make new columns for  Mean  Median  2 sigma...
    fp_data['faintMean'] = np.nan
    fp_data['faintMedian'] = np.nan
    fp_data['faintTwoSigma'] = np.nan
    fp_data['faintRMS'] = np.nan

    # calculate the faint replacement only for faint points...
    print('Faint points treatment...')
    fp_data.ix[mask, 'faintMean'] = faintF.calculate_mean(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
    fp_data.ix[mask, 'faintMedian'] = faintF.calculate_median(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
    fp_data.ix[mask, 'faintTwoSigma'] = faintF.calculate_2sigma(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
    fp_data.ix[mask, 'faintRMS'] = faintF.calculate_rms(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
   
    # 
    ######################### SAVING OUTPUT (unneccesarily takes time...)
    #

    if verbose is not None : 
        print('Saving the diagnostics...')
        path = DirOut+'Proc_'+name
        diagFile = path+'.diag'
        file = open(diagFile, "w")
        file.write('Input :  \n')
        s = '    '+ DirIn + '\n' 
        file.write(s)
        s = '    '+ name + '\n\n'
        file.write(s)
        s = 'There are '+str(np.sum(mask)) + ' points out of ' + str(len(mask)) + ' that have SNR < 2 \n '
        file.write(s)
        s = '( SNR = psfFlux / psfFluxErr ) \n \n'
        file.write(s)
        file.write('Output :  \n')
        s = '    '+ DirOut + 'Var'+name+ '\n'  
        file.write(s)
        file.close()   

    #
    ##########  STEP 2 : Derived Quantities ###########  
    #

    ####  replace all psfFlux  where SNR < 2  with  faintMean  
    rows = fp_data['flagFaint'] == 1
    fp_data.ix[rows, 'psfFlux'] = fp_data.ix[rows, 'faintMean']

    # group by objectId to calculate full LC variability characteristics 
    print('Calculating the full LC statistics...')
    grouped = fp_data.groupby('objectId')

    #
    #  Calculate low-order statistical properties of full lightcurves:
    # 
     
    varMetricsFull = grouped.apply(varF.computeVarMetrics)
    print('Calculating metrics for full lightcurves is finished')

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

    # Calculate magnitudes based on average fluxes :
    # psfMean  psfMedian  psfMeanErr  psfMedianErr 
    varMetricsFull['psfMean'] = flux2ab(varMetricsFull['psfFluxMean'], unit='Jy')
    varMetricsFull['psfMedian'] = flux2ab(varMetricsFull['psfFluxMedian'], unit='Jy')
    varMetricsFull['psfMeanErr'] = flux2absigma(varMetricsFull['psfFluxMean'],varMetricsFull['psfFluxMeanErr'])
    varMetricsFull['psfMedianErr'] = flux2absigma(varMetricsFull['psfFluxMedian'],varMetricsFull['psfFluxMedianErr'])
    print('Calculating magnitudes from fluxes is finished')
    #
    ######################### SAVING OUTPUT        ######################### 
    # 
    path = DirOut + 'Var'+name
    print('Saving varMetricsFull to  %s '%path)
    varMetricsFull.to_csv(path)
   
    #
    ##########  STEP 3 : Variable Candidates ###########  
    # 
    if calc_seasonal_metrics is None : 
        return 
    print('Identifying Variable Candidates...')
    #
    # for variable candidates we calculate  metrics for 
    # points per season, as well as bin the lightcurve 
    # into seasons, and calculate the metrics of binned 
    # lightcurve 
    # 

    # Check how many sources  are variable... 
    # The  '3 sigma' definition (more strict) ...
    m1 = (varMetricsFull['sigmaFull'] > 0)
    N = varMetricsFull['N'] #  number of pts per lightcurve 
    m4 = varMetricsFull['chi2DOF'] > (1 + 3.0 * np.sqrt(2 / N )) ##  170517
    m5 = varMetricsFull['chi2R'] > (1 + 3.0 * np.sqrt(2 / N))   ## 37464
 
    # save the diagnostics about variability to a file...
    diagFile = path+'.diag'
    file = open(diagFile, "w")
    file.write('Input :  \n')
    s = '    '+ DirIn + '\n' 
    file.write(s)
    s = '    '+ name + '\n\n'
    file.write(s)
    s = "chi2DOF is  varMetricsFull['chi2DOF'] > (1 + 3.0 * np.sqrt(2 / N )) \n" 
    file.write(s)
    s = "chi2R is  varMetricsFull['chi2R'] > (1 + 3.0 * np.sqrt(2 / N )) \n " 
    file.write(s)
    m = m1|m4|m5
    s = 'There are '+str(np.sum(m)) + ' points out of ' + str(len(m)) + \
        ' that fulfill sigma OR  chi2R  OR chi2DOF  \n '
    file.write(s)
    #m = m1*(m4|m5)
    s = 'There are '+str(np.sum(m)) + ' points out of ' + str(len(m)) + \
        ' that fulfill sigma  AND  (chi2R OR chi2DOF)  \n '
    file.write(s)
    file.write('\n')
    file.write('Output :  \n')
    s = '    '+ DirOut + '\n' 
    file.write(s)
    s = '    '+ 'Var'+name + '\n\n'  
    file.write(s)
    file.close()   
    print('Saved variable candidates info to %s' %diagFile)

    if calc_seasonal_metrics is None : 
        return 
    #
    ################ SEASONAL METRICS #################
    #
    #m = m1*(m4|m5) (already defined when saving the var stats above) 

    # Grab names of variable objects, that fulfill the criteria above... 
    varObjectIds = varMetricsFull[m].index

    # Grab only those rows that correspond to variable objects...
    rows = np.in1d(fp_data['objectId'].values, varObjectIds)
    fp_var = fp_data.ix[rows] 

    # make a new column to designate seasons...
    fp_var['season'] = np.nan

    # I insert a very early date at the beginning of the list, so that all obs between
    # 1990 and 2005 are clustered together, and 2005-2006, 2006-2007, etc are 
    # averaged seasonally (season starts August 1st)

    dates = [str(year)+'-08-01 00:00:00.000' for year in np.arange(2005,2009)]
    dates.insert(0,'1990-08-01 00:00:00.000' )
    cutDates = Time(dates, format='iso')   # use AstroPy Time module... 
    seasons = np.arange(len(cutDates))+0 

    # Assign value of a season for each row...
    for i in range(len(cutDates.mjd)-1):
        mask = (fp_var['mjd'].values > cutDates[i].mjd) * (fp_var['mjd'].values < cutDates[i+1].mjd)
        fp_var.ix[mask, 'season'] = seasons[i]  

    # Calculate seasonal metrics for objects that are variable (based on their full LC...) 
    grouped = fp_var.groupby(['objectId','season'])
    print('Calculating Seasonal metrics for %s ...'%name)
    varMetricsSeasonal = grouped.apply(varF.computeVarMetrics)
    #
    ######################### SAVING OUTPUT        ######################### 
    #
    path = DirOut + 'SeasVar'+name
    varMetricsSeasonal.to_csv(path)
    print('Saving Seasonal statistics to %s'%path)

    if calc_seas_binned_metrics is None : 
        return 
    #
    ################ SEASONALLY-BINNED LIGHTCURVE METRICS #################
    #
    # Calculate binned lightcurve metrics (binned by seasons)
    grouped  = varMetricsSeasonal.groupby(level=0)
    #grouped.get_group(grouped.groups.keys()[0])['psfFluxMean'].values
    print('Calculating Seasonally Binned LC statistics')
    varMetricsFullSeasonal = grouped.apply(varF.ComputeVarFullBinned)
    #
    ######################### SAVING OUTPUT        ######################### 
    #
    path = DirOut + 'FullSeasVar'+name
    varMetricsFullSeasonal.to_csv(path)
    print('Saving Full Seasonally binned LC statistics to %s'%path)



