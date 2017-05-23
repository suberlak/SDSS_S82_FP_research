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


def process_patch(name, DirIn, DirOut, pre='VarD_', calc_sigma_pdf=False, 
                  limitNrows=None, calc_seasonal_metrics=None, 
                  calc_seas_binned_metrics=None, verbose = None ):
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
    - identify variable candidates, based on lightcurve chi2dof , 
       chi2robust,  sigmaFull . The criteria are  : 
       m1 = sigmaFull > 0
       m4 = chi2dof > (1 + 3.0 * np.sqrt(2 / N )
       m5 = chi2robust > (1 + 3.0 * np.sqrt(2 / N)) 
       m = m1 * (m4 | m5 )    
    - select variable candidates, and identify seasons.  
    - group the variable candidates by 'objectId', 'season'
    - calculate variability metrics using  
     variabilityFunctions.computeVarMetrics()
    (same as for full lightcurves),  but on seasonally-binned lightcurves .
    Return seasonal variability metrics varMetricsSeasonal

    varMetricsFull  and varMetricsSeasonal are saved as csv files, and so 
    there is nothing that this function explicitly returns.  



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
    calc_seasonal_metrics : boolean, if not None,  the program will 
          group  all light curve points  by season, (usually few points 
          per season), and calculate statistics on these few points 
    calc_seas_binned_metrics: boolean, if not  None , the program will 
          bin the entire lightcurve  nto seasons, by averaging flux 
          etc per season. 
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
        print('In  file %s there are : \n.... %d NaN  psfFlux \
            rows'%(name, np.sum(m1)))
        print('.... %d not finite  psfFlux rows'% np.sum(m2))
        print('All such rows are dropped')
        indices = np.arange(len(test))
        rows_to_remove = indices[mask]
        raw_data.remove_rows(remove_indices)

    # 1.3 : check psfFluxErr : drop all rows which have NaN, or 0 ,
    # to avoid getting error when calculating S/N  = psfFlux / psfFluxErr 

    m1  = np.isnan(raw_data['psfFluxErr'].data)  # true if NaN 
    # true if not finite... 
    m2 =  np.bitwise_not(np.isfinite(raw_data['psfFluxErr'].data))  
    m3 = raw_data['psfFluxErr'].data == 0 
    # logical or : true if either condition satisfied 
    m = m1 | m2  | m3

    if np.sum(m) > 0 :  # only apply if there is anything to drop ... 
        print('In  file %s there are : \n.... %d NaN  psfFluxErr \
            rows'%(name, np.sum(m1)))
        print('.... %d not finite  psfFluxErr rows'% np.sum(m2))
        print('.... %d psfFluxErr = 0  rows'% np.sum(m3))
        print('All such rows are dropped')

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
    xMean =  np.interp(faint_SN, lookup['xObs'].data, lookup['xMean'].data)
    xMed =  np.interp(faint_SN, lookup['xObs'].data, lookup['xMedian'].data)
    xSigma =  np.interp(faint_SN, lookup['xObs'].data, lookup['xSigma'].data)
    xRMS = np.interp(faint_SN, lookup['xObs'].data, lookup['xRMS'].data)

    # add this info to the table 
    raw_data['faintMean'][mask_SN]     = xMean * flux_err
    raw_data['faintMedian'][mask_SN]   = xMed * flux_err
    raw_data['faintTwoSigma'][mask_SN] = xSigma * flux_err
    raw_data['faintRMS'][mask_SN]      = xRMS * flux_err



    # 1.6  replace psfFluxJy  by  faintMean ,  psfFluxErrJy  by  faintRMS 
    # make sure we are only taking values  that are not NaN ... 
    #mask_NaNs = np.bitwise_not(np.isnan(raw_data['faintRMS'][mask_SN]))

    raw_data['psfFluxJy'][mask_SN] = raw_data['faintMean'][mask_SN].data.data
    raw_data['psfFluxErrJy'][mask_SN] = raw_data['faintRMS'][mask_SN].data.data

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


    # 2.2 Calculate stats for LC using only bright points 
    print('Calculating the  LC statistics using S/N > 2  points only ...')
    mask_bright = np.bitwise_not(raw_data_df['flagFaint'].values.astype(bool))
    bright_grouped = raw_data_df[mask_bright].groupby('objectId')
    varMetricsFull_bright  = bright_grouped.apply(varF.computeVarMetrics, 
                                                  flux_column='psfFluxJy',
                                                  error_column = 'psfFluxErrJy',
                                                  time_column = 'mjd', calc_sigma_pdf =calc_sigma_pdf)

    # 2.3 Calculate stats for LC using all points 
    print('Calculating the  LC statistics using all S/N  points  ...')
    all_grouped = raw_data_df.groupby('objectId')
    varMetricsFull_all  = all_grouped.apply(varF.computeVarMetrics, 
                                                  flux_column='psfFluxJy',
                                                  error_column = 'psfFluxErrJy',
                                                  time_column = 'mjd',calc_sigma_pdf =calc_sigma_pdf) 


    # 2.4  calculate magnitudes from fluxes ... 
   

    # Calculate magnitudes based on average fluxes :
    # psfMean  psfMedian  psfMeanErr  psfMedianErr 
    varMetricsFull_all['psfMean'] = flux2ab(varMetricsFull_all['psfFluxMean'], unit='Jy')
    varMetricsFull_all['psfMedian'] = flux2ab(varMetricsFull_all['psfFluxMedian'], unit='Jy')
    varMetricsFull_all['psfMeanErr'] = flux2absigma(varMetricsFull_all['psfFluxMean'],
                                                    varMetricsFull_all['psfFluxMeanErr'])
    varMetricsFull_all['psfMedianErr'] = flux2absigma(varMetricsFull_all['psfFluxMedian'],
                                                      varMetricsFull_all['psfFluxMedianErr'])


    varMetricsFull_bright['psfMean'] = flux2ab(varMetricsFull_bright['psfFluxMean'], unit='Jy')
    varMetricsFull_bright['psfMedian'] = flux2ab(varMetricsFull_bright['psfFluxMedian'], unit='Jy')
    varMetricsFull_bright['psfMeanErr'] = flux2absigma(varMetricsFull_bright['psfFluxMean'],
                                                    varMetricsFull_bright['psfFluxMeanErr'])
    varMetricsFull_bright['psfMedianErr'] = flux2absigma(varMetricsFull_bright['psfFluxMedian'],
                                                      varMetricsFull_bright['psfFluxMedianErr'])

    print('Calculating magnitudes from fluxes is finished')


    # 2.5  change colnames to reflect which subset of points per lightcurve was used 
    # the easiest way to do it is to add_suffix in pandas 

    varMetricsFull_all = varMetricsFull_all.add_suffix('_all')
    varMetricsFull_bright = varMetricsFull_bright.add_suffix('_bright')

    # 2.6 combine the two ... 
    varMetricsFull_combined = pd.concat([varMetricsFull_all,varMetricsFull_bright], axis=1)

    ######################### SAVING OUTPUT        ######################### 
    # 
    path = DirOut + 'VarD_'+name
    print('Saving varMetricsFull to  %s '%path)
    varMetricsFull_combined.to_csv(path)

    print('The output dataframe header is : ')
    print(varMetricsFull_combined.head())
    return 


