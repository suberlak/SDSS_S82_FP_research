# -*- coding: iso-8859-1 -*-
# Necessary imports 

import numpy as np

from astropy.table import Table
from astropy.table import vstack
from astropy.table import Column
from astropy.io import ascii

#
# A program to extract lightcurves for objects in the variable sample 
# in Strip 82 forced photometry data. Initially the variable sample 
# consists of 40646 objects (sample 3a)
# 


# We tested this code on the Mac on single patch across filters : 
# DirIn,  and DirOut are here for the workstation (typhoon)
	
# Read-in the master file with info about objects in the variable sample :
fname = '/astro/users/suberlak/Desktop/variable_sample_3a_master_info.csv'
sample_data = Table.read(fname)

# find out in which patches are the objects from the variable sample ... 
patches_to_read = np.unique(sample_data['patch']).data


# Iterate over patches, extracting from each patch-filter file photometry
# for our variable objects ... 


for patch in patches_to_read : 
#patch = patches_to_read[0]   # for testing , only 00_21, otherwise we'd run over all ... 
	# choose objects that we want to extract from that patch ... 
	mask_rows= sample_data['patch'].data == patch
	objects_in_patch = sample_data['objectId'].data[mask_rows]
	print('There are %d objects of interest in patch %s ' % (len(objects_in_patch), patch))

	# for each patch, there is a different site.  On the mac I only have 00_21 which is from NCSA,
	# but in general there are two sites - for each patch it's either NCSA or IN2P3 
	site = sample_data['site'].data[mask_rows][0]


	# Initialize the dictionary to store all photometry to keep
	# from a given filter-patch file : 
	subpatch = {}

	# read all relevant rows from all filters in a given patch... 
	for filter in 'ugriz' : 
	    fname = filter+patch+'.csv.gz'
	    print('\n Reading in %s'%fname)
	    DirIn = '/astro/store/pogo4/s13_stripe82/forced_phot_lt_23/'+site+'/'
	    filter_patch = Table.read(DirIn+fname, format='csv')
	    
	    # How many rows we will extract in this filter ... 
	    n_raw = len(filter_patch)
	    print('In this filter patch there are %d raw photometric measurements '%n_raw)

	    mask_grab = np.in1d(filter_patch['objectId'].data, objects_in_patch)
	    n_raw_grab = np.sum(mask_grab)
	    print('%d of these pertain to objects of  interest, and we extract them '%n_raw_grab)

	    # select only specific columns from each patch-filter file, to speed up execution 
	    columns_of_interest = ['objectId', 'mjd', 'psfFlux', 'psfFluxErr']
	    subpatch[filter] = filter_patch[mask_grab][columns_of_interest]

	    new_col = Column(name='filter',data= np.zeros_like(subpatch[filter]['objectId'], dtype='str'))
	    new_col[:] = filter
	    subpatch[filter].add_column(new_col)

	# Save as one-file-per-object ...

	# Iterate over all objectId's ,
	# and pull photometry from  all filters
	# into one lightcurve,
	# column 'filter'  defines which 
	# filter is the measurement from 
	# data types here match precisely those in subpatch[filter],  inherited from individual patch-files ... 
	LC_stack = Table(names = ['objectId', 'mjd', 'psfFlux', 'psfFluxErr', 'filter'], 
	                            dtype=('i8', 'f8', 'f8','f8','S1'))

	cols_save = ['mjd', 'psfFlux', 'psfFluxErr', 'filter']

	for objectId in  objects_in_patch : 
	    
	    # start the stack : empty Table 
	    object_LC_stack = Table(names = ['objectId', 'mjd', 'psfFlux', 'psfFluxErr', 'filter'], 
	                            dtype=('i8', 'f8', 'f8','f8','S1'))

	    # append epochs for each filters row-by-row
	    # using very efficient vstack : columns are 
	    # exactly the same 

	    for f in 'ugriz' :
	        m =  subpatch[f]['objectId'] == objectId
	        object_LC_stack = vstack([object_LC_stack, subpatch[f][m]])
	    
	    # Save stacked LC for a given object ...
	    DirOut = '/astro/store/scratch/tmp/suberlak/s13_S82_2017/sample_LC/'
	    lc_name = str(objectId) + '.csv'
	    ascii.write(object_LC_stack[cols_save], DirOut + lc_name,  format='csv', overwrite=True)

	    # Add LC for this object to the combined stack 
	    LC_stack = vstack([LC_stack, object_LC_stack])
	    
	# Now, save the combined lightcurves too - in case it was easier to use such file...

	DirOut = '/astro/store/scratch/tmp/suberlak/s13_S82_2017/sample_LC_combined/'
	ascii.write(LC_stack,DirOut + 'Sample_LC_'+patch+'_'+str(len(len(objects_in_patch)))+'-objects.csv', 
	            format=csv,  overwrite=True )