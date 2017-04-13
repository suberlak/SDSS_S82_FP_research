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


import sys
import datetime
###  logging to a text file...

logname = datetime.datetime.now().strftime('%Y-%m-%d')
te = open(logname+'_log.txt','w')  # File where you need to keep the logs

class Unbuffered:

   def __init__(self, stream):

       self.stream = stream

   def write(self, data):

       self.stream.write(data)
       #self.stream.flush()
       te.write(data)    # Write the data of stdout here to a text file as well
   def flush(self):
       self.stream.flush()

sys.stdout=Unbuffered(sys.stdout)


# use either console arguments to make this things easier, 
# or allow to set things in stone in the code ... 

# in the future: 
# make more user-friendly https://docs.python.org/3.3/library/argparse.html 

arg1_m = ['m', 'mac', 'macbook']
arg1_t = ['t', 'typhoon', 'workstation']
arg2_n = ['1', 'NCSA']
arg2_i = ['2', 'IN2P3']

execution_environment = None
site = None
limitNrows = None

if len(sys.argv) == 2 : 
    print('Insufficient arguments : need  to call with   arg1   arg2  *arg3')
    print('arg1 : which computer we are on' )
    print(arg1_m); print(arg1_t)
    print('arg2: which site should we process patches from ')
    print(arg2_n); print(arg2_i)
    print('arg3 : limit number of rows ? (optional)')
   
    print('\nFor example:')
    print('python LC_processing.py  m 1')
    
    sys.exit()

if len(sys.argv) > 2 :
    ex = sys.argv[1]
    if ex in  arg1_m:
        execution_environment = 'mac'
    elif ex in arg1_t:
        execution_environment = 'typhoon'
    else : 
        print('Invalid argument 1 - need to choose from ')
        print(arg1_m); print(arg1_t)
        sys.exit()

    s = sys.argv[2]
    if s in arg2_n:
        site = 'NCSA'
    elif s in arg2_i:
        site = 'IN2P3'
    else  : 
        print('Invalid argument 2 - need to choose from ')
        print(arg2_n); print(arg2_i)
        sys.exit()
    
    if site is not None and execution_environment is not None : 
        print('Using execution_environment = %s and site=%s'%(execution_environment, site))

    if len(sys.argv) > 3 :
         limitNrows = int(sys.argv[3])
         print('Limiting rows used from each patch file to  %d'%limitNrows)
    
if len(sys.argv) == 1 : 
    # Need to set things manually in the code , since no arguments are provided
    print('No arguments provided from the console ... ')
    ###
    ### setting these two is CRUCIAL 
    ###
    execution_environment = 'mac'
    site = 'NCSA'
    print('Using execution_environment = %s and site=%s'%(execution_environment, site))



# Need this to locate my packages ...
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
filter_patch_files = []
if site == 'NCSA':
    patches = ['00_21', '22_43','44_65', '66_87' ,'88_109','110_131', '132_153', 
               '154_175',  '176_181', '365_387', '388_409'] 

if site == 'IN2P3':
    patches = ['155_176', '176_197','197_218', '218_239', '239_260', '260_281', '281_302', 
               '302_323','323_344', '344_365', '365_386']  


for patch in patches  :
    for filter in 'ugriz':  
        filter_patch_files.append(filter + patch + '.csv')
    
# Run this for processing : calculation of over 27 metrics per lightcurve per band 
for name in filter_patch_files :
    if limitNrows is not None  : 
        procP.process_patch(name, DirIn, DirOut,limitNrows=limitNrows)
    else : 
        procP.process_patch(name, DirIn, DirOut)



