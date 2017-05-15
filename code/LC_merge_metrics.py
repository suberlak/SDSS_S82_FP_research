# -*- coding: iso-8859-1 -*-
# Necessary imports 
import pandas as pd
import numpy as np 
from pandas import compat

#
# A program to merge ugriz variability statistics for NCSA and IN2P3, 
# that are a result of running LC_processing.py,  stored in 
# /astro/store/scratch/tmp/suberlak/s13_stripe82/forced_phot_lt_23/
#
# It combines the photometry and variability parameters on full LCs 
# with the E(B-V), correcting for extinction, and 
# adding the RA, DEC, extendedness information from DeepSource files,
# located in  
# /astro/store/scratch/tmp/suberlak/S13Agg/repo_fls/'


# Purpose :  to provide ugriz extinction-corrected brightness measures 
# (fpsfMeanCorr, N), with a few basic variability parameters 
# ($\chi^{2}_{DOF}(f)$, $f\mu_{full}(f)$, $\sigma_{full}(f)$), 
# saved as FaintFP_NCSA_variability.csv 

# To do that we :
# - read in 6 cols per patch-filter 
# - merge ugriz data from the patch 
# - add ebv info
# - correct for extinction
# - keep only extinction-corrected mean magnitudes, N, and var params
# - save all patches as  N_objects (rows) x 6 columns * 5 filters as a single file 
 
#
# Example usage : 
# 
# -- for testing  on mac using NCSA , using only 1000 first rows of each filter-patch file... 
# python VarStat_merge_ugriz.py m 1  1000
# 
# -- for full run on typhoon with NCSA 
# python VarStat_merge_ugriz.py t 1 
# 
# -- for full run on typhoon with IN2P3 
# python VarStat_merge_ugriz.py t 2 
#
# -- note : by default , narrow_cols = True, 
# which only saves a subset of columns. 
# If need to save all columns, need to set that to narrow_cols = None 
# If need to save all patches, keep test_on_N_patches = None 
# otherwise, we will only run a subset of patches  


execution_environment = None
site = None
limitNrows = None
narrow_cols = True
test_on_N_patches = None  


import sys 
import datetime 

class Unbuffered:

   def __init__(self, stream):

       self.stream = stream

   def write(self, data):

       self.stream.write(data)
       #self.stream.flush()
       te.write(data)    # Write the data of stdout here to a text file as well
   def flush(self):
       self.stream.flush()


 # lines between here and XXX are straight from LC_processing.py
arg1_m = ['m', 'mac', 'macbook']
arg1_t = ['t', 'typhoon', 'workstation']
arg2_n = ['1', 'NCSA']
arg2_i = ['2', 'IN2P3']


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
    execution_environment = 'typhoon'
    site = 'NCSA'
    print('Using execution_environment = %s and site=%s'%(execution_environment, site))



# Need this to locate my packages ...
if execution_environment == 'mac' : 
    dir_info = '../raw_data/repo_fls/'
    dir_var  = '../data_products/varMetrics/' 
    dir_save = '../data_products/varMetricsMerged/'

elif execution_environment == 'typhoon' :
    dir_info = '/astro/users/suberlak/Desktop/deep_source/'
    dir_var  = '/astro/store/scratch/tmp/suberlak/s13_S82_2017/'+site+'/'
    dir_save = '/astro/store/scratch/tmp/suberlak/s13_S82_2017/'+site+'/varMetricsMerged/'

###  logging to a text file...
# https://docs.python.org/2/library/datetime.html#strftime-strptime-behavior
logdate = datetime.datetime.now().strftime('%Y-%m-%d-%H.%M.%S')
logname = dir_save + 'VarStat_merge_'+logdate+'_log.txt'
te = open(logname,'w')  # File where you need to keep the logs
sys.stdout=Unbuffered(sys.stdout)


print('Combining ugriz variability results for forced photometry lightcurves from %s'%site)

# Get E(B-V) : once for all objects across all patches 
ebv_file = 'ebv_'+site+'_lt235.dat'
ebv = pd.read_table(dir_info+ebv_file, delimiter=' ', usecols=[0,1])
ebv.columns = ['objectId','ebv']


def add_patch(patch='00_21', ebv = ebv, varPatchesDF = None, dir_var=dir_var, 
    ebv_file =  ebv_file, limitNrows = limitNrows, narrow = None):
    '''

    Input
    ------------
    patch - name of the patch 
    ebv - pandas table with objectId , and ebv value . It has to be passed to 
          add_patch, because it is a one big file that does not need to 
          be read each time a different patch is processed 
    varPatchesDF - a data frame storing the results of calculation. Initially there is 
         none, and each consecutive patch gets appended to that data frame, which 
         is the main result of running this function 
    dir_var -  Directory storing the results of full LC statistics...
   
    ebv_file - name of a file with E(B-V) values for each object in a given 
          data source . It's solely used to print out below information about
          how many objects in a given patch have extinction information 
    limitNrows - for testing, do we want to only take N rows from each filter-file, 
          to speed up calculation (eg, merge ugriz on a given patch for only 
          N=100 objects to test the code )
    narrow - if not None  (eg. True), we limit the output file as to how many columns
          we grab, to make the output file smaller, especially when merging across many
          patches... 

    Returns
    ------------
    varPatchesDF -  a dataframe with N_objects (rows) x ( 6 x 5 filters) (cols) 

    '''
    print('\nCombining ugriz variability results from patch %s'%patch)
    
    # only read columns that we use to speed up execution.... 
    #columns = ['objectId','N','chi2DOF','muFull','sigmaFull','psfMean']

    # Read in all filters per patch ... 
    # these files have results of variability metrics calculation 
    # from computeVarMetrics() in ../packages/variabilityFunctions.py

    varPatch = {}

    if limitNrows is not None : 
        for filter in 'ugriz':
            File = 'Var'+filter+patch+'.csv'
            varPatch[filter] = pd.read_csv(dir_var+File, nrows=limitNrows, low_memory=False)  #) , usecols = columns)
    else: 
        for filter in 'ugriz':
            File = 'Var'+filter+patch+'.csv'
            varPatch[filter] = pd.read_csv(dir_var+File,low_memory=False)  #) , usecols = columns)
    

    # Check if each patch-filter file has exactly the same number of objects... 
    for filter in 'ugriz':
        print('Number of unique objectId in %s is %d'%(filter, len(np.unique(varPatch[filter]['objectId'].values))))

    # add prefix for each filter, apart from the objectId column  
    for filter in 'ugriz':
        varPatch[filter].columns = [filter+col  if col != 'objectId' else col for col in varPatch[filter]]
        
    # merge ugriz  
    howmerge='inner' # to avoid those objects which were in one filter but not in the other...
    varPatchug = pd.merge(varPatch['u'], varPatch['g'], how=howmerge, on='objectId', copy=True, indicator=False)
    varPatchugr = pd.merge(varPatchug, varPatch['r'], how=howmerge, on='objectId', copy=True, indicator=False)
    varPatchiz = pd.merge(varPatch['i'], varPatch['z'], how=howmerge, on='objectId', copy=True, indicator=False)
    varPatchugriz = pd.merge(varPatchugr, varPatchiz , how=howmerge, on='objectId', copy=True, indicator=False)

    # check how many objects have ugriz photometry and extinction information 
    withEBV = np.sum(np.in1d(varPatchugriz['objectId'].values, ebv['objectId'].values))
    allOBJ = len(varPatchugriz['objectId'].values)
    print('Of all %d objects with ugriz info, %d have E(B-V) values from %s'%(allOBJ, withEBV, ebv_file))

    # Now this can be a left merge - I only want objects that can be extinction-corrected 
    varPatchAll =pd.merge(varPatchugriz, ebv, how='inner', on='objectId', copy=True, indicator=False)

    # Correct for extinction 
    A = [5.155, 3.793, 2.751, 2.086, 1.479]
    filters = 'ugriz'

    for i in range(len(A)):
        label = filters[i] + 'psfMean'
        varPatchAll[label+'_corr'] = varPatchAll[label] +  varPatchAll['ebv'] * A[i]

    # Drop columns with uncorrected magnitudes.... 
    # not necessary - unless one is absolutely sure 
    # that the E(B-V) correction is spotless
    # perhaps there may be some use for the uncorrected
    # magnitudes, even just to check how good (uniform)
    # is the E(B-V) correction : 
    # after all, it should not change within a very small area of the sky .... 


    if narrow is not None : 
        # instead of dropping what we don't want, I choose 
        # explicitly columns to keep ... 
        filters = 'ugriz'
        cols = ['N', 'chi2DOF', 'chi2R', 'muFull', 'psfMeanErr', 'psfMean_corr']

        cols_save = [f+c for f in filters for c in cols]
        cols_save.append('ebv')
        cols_save.append('objectId')
        varPatchSave = varPatchAll.loc[:, cols_save]
    else:
        varPatchSave = varPatchAll
    # if needed to only drop uncorrected mags... 
    #varPatchSave = varPatchAll.drop(['u'+'psfMean'], axis=1)
    #for filter in 'griz':
    #    varPatchSave = varPatchSave.drop(filter+'psfMean', axis=1)

    # add a column saying which patch we are adding data from ... 
    varPatchSave.loc[:,'patch'] = patch
    
    if varPatchesDF is not None : 
        varPatchesDF = varPatchesDF.append(varPatchSave)
        print('Now we have %d objects total'%len(varPatchesDF))
        
    else : 
        varPatchesDF = varPatchSave
        
    return varPatchesDF


if site == 'NCSA':  # NCSA patches (11)
    patches = ['00_21', '22_43', '44_65','66_87', '88_109','110_131', '132_153', 
               '154_175',  '176_181', '365_387', '388_409']

if site == 'IN2P3': # IN2P3 patches (11)
    patches = ['155_176', '176_197','197_218', '218_239', '239_260', '260_281', 
               '281_302',  '302_323','323_344', '344_365', '365_386']


if test_on_N_patches is not None :  
    patches = patches[:test_on_N_patches]
    print('Using only %d patches for testing : '%test_on_N_patches)
    print(patches)
#  
# Run the first patch to start the storage DF 
varPatchesDF=  add_patch(patch=patches[0], ebv = ebv, varPatchesDF = None, narrow=narrow_cols)

# Loop over the rest of the patches to append to that DF 
for patch in patches[1:]:
    varPatchesDF=  add_patch(patch=patch, ebv = ebv, varPatchesDF = varPatchesDF, narrow = narrow_cols)
    
# At this point, we have the following columns : 
# ( if narrow  = True )
# np.ravel(varPatchesDF.columns) = array(['uN', 'uchi2DOF', 'uchi2R', 'umuFull', 'upsfMeanErr',
#       'upsfMean_corr', 'gN', 'gchi2DOF', 'gchi2R', 'gmuFull',
#       'gpsfMeanErr', 'gpsfMean_corr', 'rN', 'rchi2DOF', 'rchi2R',
#       'rmuFull', 'rpsfMeanErr', 'rpsfMean_corr', 'iN', 'ichi2DOF',
#       'ichi2R', 'imuFull', 'ipsfMeanErr', 'ipsfMean_corr', 'zN',
#       'zchi2DOF', 'zchi2R', 'zmuFull', 'zpsfMeanErr', 'zpsfMean_corr',
#       'ebv', 'objectId', 'patch'], dtype=object)

#
# Add ra, dec , extendedness information 
#
deep_source_ext = pd.read_csv(dir_info+'DeepSource'+site+'_i_lt235_extendedness.csv.gz', compression='gzip', 
                         index_col=0)
# np.ravel(deep_source_ext.columns) = array(['deepSourceId', 'extendedness'], dtype=object)


deep_source_radec = pd.read_csv(dir_info+'DeepSource'+site+'_i_lt235_narrow.csv.gz', compression='gzip', 
                         index_col=0)
# np.ravel(deep_source_radec.columns) = array(['parentDeepSourceId', 'deepCoaddId', 'ra', 'decl', 'psfMag',
#      'psfMagSigma', 'tract', 'patch', 'detect_is_primary'], dtype=object)

# only choose these two columns 
radec = deep_source_radec[['ra','decl']]


# Add extendedness, and ra,dec information :  
ext_radec = pd.merge( deep_source_ext,radec, how='left', left_on='deepSourceId', right_index=True)

# merge this in to the varPatchesDF 
# left merge : only add the information from 
# ext_radec to the objects already present in 
# varPatchesDF 
varPatchesDF1 =  pd.merge(varPatchesDF,ext_radec, how='left', left_on = 'objectId', 
                      right_on = 'deepSourceId') 

# Select out those objects that had parents brighter than iPsfMag  (uncorrected for extinction)
# < 17 mag , because for those objects the deblender was not working properly,
# and variability may be spurious 

# Make sure we are keeping only objects without bright parents
good_sources = np.load(dir_info+site+'_source_without_bright_parent.npy')
mask_keep = np.in1d(varPatchesDF1.objectId.values , good_sources)

varPatchesDF_save = varPatchesDF1[mask_keep]
varPatchesDF_discard = varPatchesDF1[~mask_keep]


print('Out of total number of %d objects, with combined metrics across ugriz filters and  %d patches '%
      (len(varPatchesDF1), len(patches)))

print('There are %d that have a parent i<17 mag, which are discarded' % len(varPatchesDF_discard))

# http://stackoverflow.com/questions/28535067/unable-to-remove-unicode-char-from-column-names-in-pandas 
# this thing prevents a disaster before I can get a hang of 
# how to run Python 3.5 on typhoon...
if sys.version_info >= (3,0,0):
    compat.PY3 = True
else : 
    compat.PY3 = False 


if len(varPatchesDF_discard) > 0 : 
    file_discard = 'Var_ugriz_'+str(len(patches))+'_patches_'+site+'_discarded.csv.gz'
    print('\nWe save these objects separately, to  %s '%file_discard)
    varPatchesDF_discard.to_csv(dir_save+file_discard , compression='gzip')

if narrow_cols is not None : 
    file_save = 'Var_ugriz_'+str(len(patches))+'_patches_'+site+'_narrow.csv.gz'
    print('\nSaving only selected columns : ')
    print(np.ravel(varPatchesDF_save.columns))
else:
    file_save = 'Var_ugriz_'+str(len(patches))+'_patches_'+site+'.csv.gz'

print('\We save the  %d objects without bright parents to %s'%(len(varPatchesDF_save), dir_save+file_save))

# This is the main product : across filters and patches merged file... 
varPatchesDF_save.to_csv(dir_save+file_save, compression='gzip' ) 



