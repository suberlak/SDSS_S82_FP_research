# place to have code that can be directly pasted into 
# ipython   interactive environment in magneto or 
# typhoon,  
# to help trace down any rare errors ...


# unfortunately, it doesn't trace changes in the original
# processPatch2.py,  it's just a temporary snapshot

import numpy as np
from itertools import product
import sys
import os 
import datetime
import argparse
import pandas as pd
from astropy.time import Time
from astropy.table import Table
from astropy.table import Column


site=  'IN2P3'
path_to_home = '/astro/users/suberlak/' 
DirIn = '/astro/store/pogo4/s13_stripe82/forced_phot_lt_23/'+site+'/'
DirOut = '/astro/store/scratch/tmp/suberlak/s13_S82_2017/'+site+'/'
sys.path.insert(0, path_to_home + 'SDSS_S82_FP_research/packages/')
pd.options.mode.chained_assignment = None
import faintFunctions as faintF 
import variabilityFunctions as varF
import processPatch2 as procP


name = 'u155_176.csv'

print('\n Processing filter_patch file %s' % name)
   
limitNrows = None
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
        rows'%(fname, np.sum(m1)))
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
        rows'%(fname, np.sum(m1)))
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
mask_SN = SN.data < 2 
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






