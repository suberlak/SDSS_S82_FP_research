# -*- coding: iso-8859-1 -*-
#
#  Chris Suberlak 04/30/2016 
#  
# OVERVIEW : 
# Code to process the S82 forced  photometry lightcurves by calculating  
# 2sigma, mean, median, based on the truncated Gaussian distribution 
# 
# INPUT : raw forced photometry lightcurves, arranged by objectId , with 
# ['objectId', 'mjd', 'psfFlux', 'psfFluxErr' ] 
#
# OUTPUT: processed forced photometry lightcurves,  with unchanged flux, etc,
# but with additional columns  :
# ['objectId', 'mjd', 'psfFlux', 'psfFluxErr', 'flag','faintMean','faintMedian','faintTwoSigma']  



# Note : to run on magneto do
# %cpaste  in jupyter notebook 
#
# Note:  
# %cpaste is sensitive to tabs : 
# make sure that in whatever editor you do tabs
# it actually puts in 4 space



# Note   2/16/17 :  many updates, most importantly, 
# made it possible to develop this code on mac, test it, 
# and then port via github to the 
# workstation (typhoon), 
# so that we can easily run the variability 
# code using magneto , but develop locally, on a mac  

################################################################################################

# make all necessary imports....
import os 
import numpy as np


# make python aware of my packages...
import sys

###
### setting these two is CRUCIAL 
###
execution_environment = 'mac'  # or 'typhoon'
site = 'NCSA'


if execution_environment == 'mac' : 
    path_to_home = '/Users/chris/GradResearch/'

elif execution_environment == 'typhoon' :
    path_to_home = '/astro/users/suberlak/' 

sys.path.insert(0, path_to_home + 'SDSS_S82_FP_research/packages/')

import processPatch as procP

#####
#####  PROCESSING FILES
#####  

if execution_environment == 'mac' : 
    DirIn = '/Users/chris/GradResearch/SDSS_S82_FP_research/raw_data/rawDataFPSplit/'
    DirOut = '/Users/chris/GradResearch/SDSS_S82_FP_research/data_products/'
elif  execution_environment == 'typhoon' : 
    DirIn = '/astro/store/pogo4/s13_stripe82/forced_phot_lt_23/'+site+'/'
    DirOut = '/astro/store/scratch/tmp/suberlak/s13_S82_2017/'+site+'/'


# Define patches which we will process in this run ... 
lToDoFilt= []
if site == 'NCSA':
    patches = ['00_21','22_43','44_65', '66_87' ,'88_109','110_131', '132_153', 
               '154_175',  '176_181', '365_387', '388_409']  # NCSA 
if site == 'IN2P3':
	patches = ['155_176', '176_197','197_218', '218_239', '239_260', '260_281', '281_302', 
               '302_323','323_344', '344_365', '365_386']  # IN2P3
# 
for patch in patches  :
    for filter in 'ugriz':
        lToDoFilt.append(filter + patch + '.csv')
    
# Run this for processing : calculation of over 27 metrics per lightcurve per band 
for name in lToDoFilt[16:]:
    procP.process_patch_minimal(name, DirIn, DirOut)





# 
# obsolete : if no use for that, remove entirely... 
# 
#lProc = []
#lProc += [each for each in os.listdir(DirOut) if each.endswith('.csv')]

#lProc = [name[5:-5] for name in lProc]
#lToDo = os.listdir(DirIn)
#lToDoComp = [name[:-3] for name in lToDo]
#if len(lProc) > 0 : 
#    maskDone = np.in1d(lToDoComp,lProc)
#    lToDoFilt =  np.array(lToDoComp)[~maskDone]
#else:
#    lToDoFilt = lToDoComp


