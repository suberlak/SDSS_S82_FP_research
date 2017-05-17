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
# ['objectId', 'mjd', 'psfFlux', 'psfFluxErr', 'flag','faintMean','faintMedian',
#  'faintTwoSigma']  



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


# Note  5/15/17 :  substantial update : introduced argparse module
# to make usage much much easier and more intuitive 
# 

##########################################################################
# make all necessary imports....
import os 
import numpy as np
from itertools import product
import sys
import datetime
import argparse

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


### Initialize the ArgumentParser   
desc  = 'Process the patch files from S82 ... '
parser = argparse.ArgumentParser(description=desc)

# Args setting the environment variables : 
# location of input and  output, etc. 

# -env :  environment : either mac , m ,   or typhoon, t
#      by default, it is mac, hence optional , 
#      but strongly recommended since it will be only 
#      correct if it is actually run on the Mac
parser.add_argument("-e", "-environment", help="set the execution \
                    environment", action="store", default = 'm', 
                    choices=['t', 'typhoon', 'm', 'mac'])


# -site  : site: NCSA ,1    or IN2P3 , 2  
parser.add_argument("-s", "-site",help="set the data processing center from \
                    which to use the data",  action ='store', default='1', 
                    choices=['NCSA', '1', 'IN2P3', '2'])

# -nlines : if want to process only n lines from each band-patch file 
parser.add_argument("-n", "-nlines", help="limit the number of rows to \
                    process in each patch-file", action="store", default=None, 
                    type=int)

#
# The following args set which patches should be processed  : 
#
# 1) a single_patch : need to provide a name 
# 2) knowing the list_of_patches per processing center, 
#    process only  list_of_patches[patch_start : patch_end]
#    by default,  we process all patches 
# 3) process only patches which do not yet have the Var.. 
#    metrics files in the output directory  : set -cd ,  and
#    provide the -pre  prefix for output filenames to 
#    check.  We assume that output of LC_processing.py 
#    follows the following naming scheme : 
#    pre + filter + patch + .csv ,   eg. VarD_i00_21.csv 

# -single_patch : if we want to process only a single patch and then stop 
parser.add_argument("-sp", "-single_patch",help="set the patch which we \
                    should process",  action ='store', default=None, 
                    type=str)



# -patch_start : which patch to start the processing with? Useful if we want to 
#    process from N to end of the alphabetic list of patches,
#    and not from 0 to end.  
parser.add_argument("-ps", "-patch_start", help='set which patch to \
                    start from, given their alphabetical ordering', 
                    action="store", default=None, type=int, 
                    choices = range(0,11))

# -patch_end : if only want to merge patches from 0 to N  ....  
# only merge N patches instead of all for which 
# aggregate metrics are available ? 
parser.add_argument("-pe", "-patch_end", help="set how many patches to merge \
                    if not all for which data is available", 
                    action='store', default = None, type=int, 
                    choices = range(0,11))

# -cd : check outDir  if yes  (default no),  it would run the check of files 
#    that startwith  -pre
parser.add_argument("-cd", "-check_dir", help='check the output directory for \
                    which patch-files have already been processed? If so, also \
                    need to set the -pre  variable indicating the prefix of \
                    the outfiles to be checked. By default it is VarD_ ', 
                    action='store_true')

# -pre : if before processing the raw lightcurve files you would like to 
#   check the outDir which processed files already exits,  you need to 
#   provide a prefix for output files, eg  VarD_ .  This is the string before 
#     g00_21.csv 
parser.add_argument("-pre", "-prefix", help = 'set the prefix for output \
                    files to be checked for which patches have already been \
                    processed', action='store', default='VarD_', type=str)



# parse all arguments : do it only once in an entire program ... 
args = parser.parse_args()

print('Environment:  %s'%args.e)
print('Site %s'%args.s)

if args.n:
    print('Nlines : %d'%args.n)

if args.ps : 
    print('Starting from the patch #%d '%args.ps)

if args.pe : 
    print('Ending from the patch #%d '%args.pe)

if args.sp : 
    print('Only processing patch %s'%args.sp)

# Need this to locate my packages ...
# Also , set the input and output
# directories depending on the 
# workspace ... 

if args.s in ['1', 'NCSA'] : 
    site = 'NCSA'
elif args.s in ['2', 'IN2P3']:
    site = 'IN2P3'

if args.e in  ['m', 'mac'] : 
    path_to_home = '/Users/chris/GradResearch/'
    DirIn = '/Users/chris/GradResearch/SDSS_S82_FP_research/raw_data/rawDataFPSplit/'
    DirOut = '/Users/chris/GradResearch/SDSS_S82_FP_research/data_products/varMetrics/'

elif args.e in ['t','typhoon'] :
    path_to_home = '/astro/users/suberlak/' 
    DirIn = '/astro/store/pogo4/s13_stripe82/forced_phot_lt_23/'+site+'/'
    DirOut = '/astro/store/scratch/tmp/suberlak/s13_S82_2017/'+site+'/'


# we import the custom-written module to process 
# each patch file   
# only once path_to_home is added as a searchable path
sys.path.insert(0, path_to_home + 'SDSS_S82_FP_research/packages/')
import processPatch2 as procP

#####
#####  PROCESSING FILES
#####  


# Define patches which we will process in this run ... 
if args.s in ['1', 'NCSA']:
    patches = ['00_21', '22_43','44_65', '66_87' ,'88_109','110_131', '132_153', 
               '154_175',  '176_181', '365_387', '388_409'] 

if args.s in ['2', 'IN2P3']:
    patches = ['155_176', '176_197','197_218', '218_239', '239_260', '260_281', 
               '281_302', '302_323','323_344', '344_365', '365_386']  

# patch selection... 
if args.ps :
    if args.pe : 
        # select patches m to n 
        patches = patches[args.ps:args.pe]
    else:
        # select from m-th patch onwards
        patches = patches[args.ps:] 
elif args.pe : 
    # select up to n-th patch 
    patches = patches[:args.pe]

if args.sp  : 
    filter_patch_files = [args.sp]
else : 
    filter_patch_files = []
    for patch, filter in product(patches,'ugriz') :
        filter_patch_files.append(filter + patch + '.csv')

#  check for already processed files, 
#  need to provide -pre argument 
#  as well if want to check for any 
#  custom prefix .... 
if args.cd : 
    print('Checking the output directory  %s'%DirOut)
    print('for files with prefix %s'%args.pre)

    lista = np.array(os.listdir(DirOut))
    mask = np.array([l.startswith(args.pre) for l in lista])

    # remove only if there are any matching files in the directory ... 
    if np.sum(mask) > 0 : 
        full_filenames_to_remove = lista[mask]
        length = len(args.pre)
        filter_patch_to_remove = [name[length:] for name in full_filenames_to_remove]
        print('Removing these files from patch-files to process...')
        print(filter_patch_to_remove)
        for name in filter_patch_to_remove : 
            filter_patch_files.remove(name)


print('Starting to process the following patch-files:')
print(filter_patch_files)
print('From the input directory %s'%DirIn)

# Check if the input files are present ... 
# assume that the files end with .gz, since 
# the raw FP lightcurves are compressed ... 
list_input = os.listdir(DirIn)
list_csv = [name[:-3] for name in list_input]
mask_missing_input = np.in1d(filter_patch_files, list_csv)
if np.sum(~mask_missing_input) > 0 : 
    print('%d of these are not in the input directory'%np.sum(~mask_missing_input))
    print('So we process only the present patch-files:')
    filter_patch_files = np.array(filter_patch_files)[mask_missing_input]
    print(filter_patch_files)

# Run this for processing : calculation of over 27 metrics per lightcurve per band 
for name in filter_patch_files :
    if args.n :
        procP.process_patch(name, DirIn, DirOut, calc_sigma_pdf=False, 
                            limitNrows=args.n)
    else : 
        procP.process_patch(name, DirIn, DirOut, calc_sigma_pdf=False)



