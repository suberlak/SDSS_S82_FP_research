# Extract  sources without bright  parents (i.e. those to keep) from either processing center 

# Note  : this code is based on Yusra's comments in README file 
# I wrote a similar code, but it was a bit slower, so I decided to use this one, 
# but their results have been tested, and they work in an identical fashion 

# Both my old code and Yusra's example code perfectly agree : 
# we have in NCSA 585920 sources that have an unreliably bright (i<17) parent. 
# That corresponds to  (585920/5474350)∗100%≈10%
# of all NCSA sources.


# I feel that in the long run this code should be streamlined, 
# and become part of VarStat_merge_ugriz.py,  
# so that there are as little steps and pieces of code 
# to execute, files to read-in, keep updating and passing on, 
# as possible . After all, running this for all 
# NCSA and IN2P3 sources takes not more than a minute 
# even on a mac... 

# In NCSA,  there are 4888430 good sources (i.e. without a bright parent ) (=  5474350 - 585920 )
# In IN2P3, there are 4602931 good sources (i.e. without a bright parent ) ()

import pandas as pd
import numpy as np

# set where we run the code, which processing center, a
# nd whose code shall we use 
site = 'IN2P3'  # NCSA 
execution_environment = 'typhoon'
use_yusra_code   = True


if execution_environment == 'mac' : 
    dir_info = '../raw_data/repo_fls/'
   
elif execution_environment == 'typhoon' :
    dir_info = '/astro/users/suberlak/Desktop/deep_source/'
   

if use_yusra_code is True : 
    print('Using code 1')
    print('Using deep source files from %s'%dir_info)
    source = pd.read_csv(dir_info+'DeepSource'+site+'_i_lt235_narrow.csv.gz', compression='gzip', 
                         index_col=0)
    parent = pd.read_csv(dir_info+'DeepSource'+site+'_i_lt235_narrow_not_primary.csv.gz', 
                         compression='gzip', index_col=0)
    parent_bright = parent.loc[parent.psfMag < 17.0]
    # we are adding information about the parent source if the parent is brighter 
    # than 17 mag. Thus all sources that do not have a bright parent will 
    # have these new columns 'null'
    parent_source_joined = pd.merge(source, parent_bright[['psfMag', 'psfMagSigma', 'ra', 'decl']],
                             how='left', suffixes=['', '_parent'],
                             left_on='parentDeepSourceId', right_index=True)
    # len(parent_source_joined.index) = 5474350  (5.47 mln)
    # these are all sources to begin with ...

    source_without_bright_parent = parent_source_joined.loc[parent_source_joined.psfMag_parent.isnull()]

    # the opposite of that are all sources that have bright parent : 

    source_with_bright_parent = parent_source_joined.loc[~parent_source_joined.psfMag_parent.isnull()]
    # len(source_without_bright_parent.index) = 4888430,
    # i.e. there are 4888430  of such sources (4.88 mln), 
    # they are all unique
    print('Number of sources without a parent brighter than 17 mag')
    print(len(np.unique(source_without_bright_parent.index.values)))

    # thus we have 5474350-4888430 = 585920  sources  (0.59 mln)  
    # that have a too bright parent in NCSA, and should not be used . 

    # Save in the quickest-to-load format ....
    print('Note : we are saving only the objectIds ')
    np.save(dir_info+site+'_source_without_bright_parent.npy', source_without_bright_parent.index.values)


if use_yusra_code is not True : 
    # Chris code ! 

    # My old code from filter_bright_parent_sources.py  , 
    # to confirm that it works  same way .....

    # Note  : I decide here to use the same DeepSource files as in 
    # /astro/store/pogo4/s13_stripe82/deep_source/  , 
    # which are cut at i < 23.5 mag 

    primary_file = 'DeepSourceNCSA_i_lt235_narrow.csv.gz'
    not_primary_file = 'DeepSourceNCSA_i_lt235_narrow_not_primary.csv.gz'

    source = pd.read_csv(dir_info+primary_file, compression='gzip')
    # primary == source 
    # NCSA has 5474350 rows (5.4 mln) , each row is a unique deepSourceId 


    parent = pd.read_csv(dir_info+not_primary_file, compression='gzip')
    # not_primary == parent 
    # NCSA has 1957486 rows here (1.9 mln) : these are deblender parents and 
    # secondary detections in overlap region 

    source_parent  =  pd.merge(source,parent, how='left', left_on = 'parentDeepSourceId', 
                      right_on = 'deepSourceId', suffixes=('_primary','_parent')) 
    # add info about parent magnitude to remove those primary sources that have a parent with
    # iPsfMag < 17

    # apart from -1 , which means that there is no parent, there may be a single parent to various 
    # primary detections.  (eg. parentDeepSourceId_primary = 1398579193184483  has three 
    # deepSourceId_primary

    # remove NaN rows (this would correspond to eg. primary sources which were not blended, i.e.
    # have no parent , and so during merge, their parentDeepSourceId = -1  would match no rows 

    mask_finite = np.isfinite(source_parent['psfMag_parent'].values)

    # How many have a parent brighter than 17 mag ? 
    mask_bright_parents = source_parent['psfMag_parent'].values[mask_finite] < 17
    print('Number of sources with a parent brighter than 17 mag')
    print(np.sum(mask_bright_parents))


    # How many deepSourceIds have too bright parents source to be reliably handled by the deblender  ? 
    # NCSA:  -  that many rows, which corresponds to 
    print('Number of unique deep source Ids that have a parent brighter than 17 mag')
    print(len(np.unique(source_parent['deepSourceId_primary'][mask_finite].values[mask_bright_parents])))
    # It's 585920  for NCSA ;  i.e. 0.58 mln out of 5.4 mln - 10 % 


    #  Save those objects (to be ignored in all future analysis as unreliable)
    #  to a file "IN2P3_lt235_primaries_with_parents_lt_17.csv" 

    merged_too_bright = source_parent[mask_finite][mask_bright_parents]

    # If needed to save them ...
    # choose only columns that are needed 
    # merged_too_bright_narrow = merged_too_bright[['deepSourceId_primary',
    # 'parentDeepSourceId_primary','psfMag_primary', 'psfMag_parent']]
    # save as a csv ...
    # merged_too_bright_narrow.to_csv(Dir+too_bright_file[:-4] + '_narrow.csv')
        
