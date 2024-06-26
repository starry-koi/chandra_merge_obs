#!/usr/bin/env python
from ciao_contrib.runtool import *
import numpy as np
import glob
import os
import regions as reg
import region #yes, this is different than the regions package. I know.
import csv
import sys
import re
import math
import collections
import warnings
import pandas as pd
import subprocess as sp
import commentjson as cjson
import poissonstats as ps #user-made file, has functions for finding the error for counts
import copy
from sh import gunzip
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy import units as u
from decimal import Decimal
from distutils.util import strtobool
from pathlib import Path

# Auther: Lilikoi Latimer
# Created November 2022
# Currently works under CIAO 4.14

'''
======================================================================================================


### Quirks and possible issues ###

- This code is meant to be run after chandra_xray_analysis, and uses several data products that result from running that script. So it
  won't work if you haven't run that beforehand.

- Need to have the folder containing this code (folder chandra_merge_obs) in the same directory as the obsids (same as for
  chandra_xray_analysis).

- The user will also need to copy-paste the params.config file from the run of xray_flux.py, as we use some of the same parameters when
  calculating the merged flux.
  
- The band_check and filter_check variables in params_merged.config will tell the code which repro folder (in the obsid folders) to look 
  in for files. So if you want to get the merged flux for e.g. filter_check = '2-7' and band_check = '2-10', you'll need to have run 
  xray_flux.py with the appropriate filter_check and band_check variable settings first.
  
- If the event files are not astrometrically aligned (which the code will check for you) then the script will exit, as merging event
  files that have astrometric differences could result in bad flux measurements (see the CIAO thread for this at
  https://cxc.cfa.harvard.edu/ciao/threads/fluxes_multiobi/ )

- The resulting merged event file unfortunately can't be used for very much (basically just wavdetect), which means that for the
  source regions we have to get clever. We make source regions (in the exact same manner as for xray_flux.py) at the position of any
  detected merged sources, but on the event files from the separate obsids, then take the resulting source radii and combine them to
  get the final source radii that we use. This works well when the objects are on the same part of the chip for all the observations
  (e.g. if all the observations are targeting the same object) but can run into issues with wildly differing source radii if for instance
  the object is at the S3 aimpoint for one observation but at the edge of the chip for another observation.

- The user should first run chandra_xray_analysis on all their obsids. After it's completed, then identify which obsids you want to merge,
  and copy-paste the obsid folders into a new directory containing only the chandra_merge_obs folder. Additionally, copy-paste over the
  params.config file from chandra_xray_analysis - the .config file you used for the run, specifically! We pull parameters out of it and
  if any of those parameters were changed between the run of xray_flux.py and merge_obs.py then you might not be getting the correct
  results. Then adjust the parameters in params_merged.config, then run merge_obs.py

### Broad Overview ###

This program merges multiple observations (2+), detects sources from the resulting merged data, and finds the x-ray fluxes for those merged sources. It was intended to be run after a run of xray_flux.py from chandra_xray_analysis, with the user copy-pasting over the params.config file and the relevant obsids.

We first detect the cleaned, astrometrically corrected evt2 files (new_evt2.fits) from the OBSID folders. We then merge the observations using merge_obs. Next, we create a new galaxy region by reading in the galaxy regions from the obsid folders and averaging them out to get the ra, dec, and radius. If the different obsids' values are too different (standard deviation above max_std) we print that to the screen. 

We then run wavdetect on the merged files from merge_obs and exclude wavdetect images from outside the galaxy region from further analysis. We then check the astrometry alignment between the different obsids' event files. We use the wavdetect sources as starting points and search around these areas for source centroids in the different evt2 files. We compare these source centroids to see if they differ by too much - if they do then we exit the code (since merging evt2 files with too big an astrometric offset can give incorrect results). If they're within the limit then we proceed with analysis.

We then create our source and background regions. For source regions we make source regions using psfsize_srcs using each obsid's evt2 file, then average them out to get our final source region data. If the different source regions differ over a user-defined limit, then we print that to the screen as well as the different source radii that were calculated for each obsid. For the background regions, we use annuli (same as for xray_flux.py) with the outer radius determined from a multiple of the source region radius (using the bg_radius_factor variable).

We finally run srcflux and save the results to a text file of the form 'results_merged_[flux band]-band_[filter band]-filter.txt'.




======================================================================================================
'''

#=====================================================================================================
#=====================================================================================================

#############################################################
######################## VARIABLES ##########################
#############################################################

# There are a few parameters that we use that will be the same as the ones used in chandra_xray_analysis. We get them from the copy-pasted
# params.config file, which the user should have copied from the run of xray_flux.py they did (without changing any of these values).
with open('params.config', 'r') as paramfile:
    
    json_obj                = cjson.load(paramfile)
    
    coord_frame             = json_obj['coord_frame']
    wavdet_analysis_scales  = json_obj['wavdet_analysis_scales']
    fluximage_psfecf        = json_obj['fluximage_psfecf']
    analysis_psfecf         = json_obj['analysis_psfecf']
    analysis_energy         = json_obj['analysis_energy']
    bg_radius_factor        = json_obj['bg_radius_factor']
    srcflux_PhoIndex        = json_obj['srcflux_PhoIndex']
    round_dec               = json_obj['round_dec']
    flux_ref                = float(json_obj['flux_ref'])

# We need additional parameters specific to merging, which we import from the params_merged.config file.
with open('params_merged.config', 'r') as paramfile:
    
    json_obj        = cjson.load(paramfile)
    
    distance        = json_obj['distance']
    filter_check    = json_obj['filter_check']
    band_check      = json_obj['band_check']
    aperture_correct= json_obj['aperture_correct']
    max_sep         = json_obj['max_sep']
    max_std         = json_obj['max_std']
    fine_correct    = json_obj['fine_correct']
    reference_obsid = str(json_obj['reference_obsid'])
    obsids_to_correct = json_obj['obsids_to_correct']
    force_continue  = bool(strtobool(json_obj['force_continue']))
    

# While the following variables should be declared in the params_merged.config file, we keep this description here so that anyone
# working on the code has a convenient explanation without having to open another file.
'''
distance            = 139.8     # Distance to the galaxy, in Mpc. Since this script was meant to be run on multiple observations of a
                    # single galaxy, this needs to be user input in params_merged.config, as pulling it from the copied
                    # params.config would make life more difficult.

filter_check        = '2-7'     # Same as the filter_check variable the user should be familiar with from xray_flux.py. Since users might
                    # have run xray_flux.py multiple times in different bands, there could be more than one repro directory
                    # in each obsid folder. This tells the code which one to focus on, as well as being used for the fluximage
                    # part of merge_obs.

band_check          = '2-10'    # Similar to above - this tells the code which repro directory to look for, and also gets plugged into
                    # srcflux later in the code. Same as the band_check variable in xray_flux.py.

max_sep             = 0.1       # In arcseconds, this is the maximum separation tolerated between source centroids in different obsids.
                    # Different observations of the same object might have slightly different astrometry, so in this code
                    # we run wavdetect on the merged event file, then use this position to search for source centroids (in a
                    # small radius around the position) in all the different obsid event files. This variable is the maximum
                    # separation we can detect between the positions of those source centroids in different obsid event files
                    # before terminating the script. The rule of thumb for this is that we need the astrometric differences
                    # to be much less than a chandra pixel (0.492 arcsec), so the default value for is is 0.1 as ~1/5th of
                    # a pixel. It probably shouldn't be changed to be much higher. If the code does detect too big an offset
                    # and stops, the user will have to do astrometry realignment on the offending obsids.

max_std             = 0.02      # Slightly complicated to explain, but essentially this is a measure of how different the calculated
                    # source radii are allowed to be before we start caring about it. In xray_flux.py, we find source radii
                    # by using psfsize_srcs to create a region that captures a certain fraction of the psf energy (given by
                    # the variable analysis_psfecf) at a certain energy (given by analysis_energy). Since psfsize_srcs 
                    # cannot be run on the merged event file, we instead run it on all the different event files from the
                    # different obsids and average the calculated radii together to arrive at a final source region radius.
                    # During this finding of the average, we also calculate the standard deviation (std) of all of the
                    # calculated radii. This can likely be set higher than the current value of 0.02, but less variation
                    # is naturally better. Note that we also use this as the limit for the std between the different galaxy 
                    # regions.

fine_correct        = 'make'    # Controls what the code does with fine astrometric corrections. Options are 'make' and 'use'. If 
                    # set to 'make': assuming fine corrections are detected as needed (e.g. the source centroids differ by a non-
                    # negligible amount) then the code will carry out fine astrometric corrections (assuming reference_obsid  has 
                    # been specified) save the corrected files, then quit. If set to 'use': when starting, instead of searching out 
                    # and using the evt2 files from the obsid directories, the code will search out and use the finely corrected
                    # images (which were saved when running with 'make'). 

reference_obsid     = 25274     # The obsid to use as 'correct' when making fine astrometric corrections between observations.
                    # Required to be set to a valid value to use fine_correct='make'. Can be specified as either an int or a string,
                    # e.g. '12345' or 12345.

obsids_to_correct   = [],       # A list of obsids to calculate and apply fine corrections to (when fine_correct='make'). If set
                    # to a list of obsids, e.g. obsids_to_correct=[12345,23456], fine corrections will only be calculated for the 
                    # listed obsids. If set to an empty list, e.g. obsids_to_correct=[], fine corrections will be calculated for
                    # all obsids (save for reference_obsid). When fine_correct='use', will be used to determine which files to use
                    # for merging, wavdetect, etc. Since there will be a mix of corrected and uncorrected obsids (for instance,
                    # the reference_obsid will remain uncorrected), this list is used to figure out whether to use the corrected
                    # or uncorrected event files for a specific obsid. For example, if obsids_to_correct=[12345,23456], then the code
                    # would only use corrected event files for obsids 12345,23456, and the original uncorrected event files for all
                    # other obsids. Basically - if you want to correct everything, leave it blank. If you want to correct not all
                    # files, set the list to the obsids you want to correct. In general, leave the list the same for both 
                    # fine_correct='make' and fine_correct='use'.

force_continue      = False     # Forces continuation of analysis, regardless of current astrometry. Should be used sparingly, as 
                    # setting to True can lead to incorrect flux results if the astrometry is bad.
'''

# Quickly checking that there isn't a mismatch between the filtering and flux bands
if filter_check == '0.5-2':
    if band_check != '0.5-2':
        print('Filter and band mismatch - please fix')
        print('Filter: '+filter_check)
        print('Band:   '+band_check)
        sys.exit()        
if filter_check == '2-7':
    if band_check != '2-10':
        print('Filter and band mismatch - please fix')
        print('Filter: '+filter_check)
        print('Band:   '+band_check)
        sys.exit()

# Using the filter_check and band_check variables to get our srcflux and fluximage (filtering) bands
if filter_check == '0.5-2':
    filt_band = '0.5:2:1.56'
elif filter_check == '0.5-7':
    filt_band = '0.5:7:2.3'
elif filter_check == '2-7':
    filt_band = '2:7:4'
else:
    print("filter_check variable not recognized")
    sys.exit()

# Quickly getting the srcflux band, which includes both the band and the specific energy.
if band_check == '0.5-2':                # The band and specific energy to use when
    srcflux_band = '0.5:2:1.56' # soft        # calculating fluxes. It's of the form
elif band_check == '0.5-7':                # (lower):(upper):(specific_energy), with
    srcflux_band = '0.5:7:2.3' # broad        # everything in keV. Since the specific
elif band_check == '0.5-8':                # energy changes depending on which band is
    srcflux_band = '0.5:8:2.3' # ~broad        # used, it seemed easier to just do if 
elif band_check == '2-10':                # statements then letting it be user-driven.
    srcflux_band = '2:10:4' # hard            # If you add a new energy band, you'll need
else:                            # to update this part as well.
    print("band_check variable not recognized")
    sys.exit()

# Directory things
code_dir = os.getcwd() # Get the current working directory, and save the paths to various files for later use
os.chdir('..') # Move up a directory out of the chandra_merge_obs folder, so we can access the obsid folders
over_dir = os.getcwd() # Store the overall working directory

# Get a list of the obsids directories only (excluding non-numeric directories)
obsids = next(os.walk('.'))[1]
obsids = [a for a in obsids if a.isdigit()]
obsids.sort()

# If no obsids were listed in obsids_to_correct, then fill it with all the obsids save for reference_obsid. Also turn all the obsids into strings.
if len(obsids_to_correct) < 1:
    obsids_to_correct = copy.deepcopy(obsids)
    obsids_to_correct.remove(reference_obsid)
obsids_to_correct = [str(obsid_indiv) for obsid_indiv in obsids_to_correct]

# Create an analysis/ folder to store all of the files and everything in
analysis_dir = over_dir + '/analysis_'+band_check+'-band__'+filter_check+'-filter'
sp.run(["mkdir", "-p", analysis_dir])

# Fill dictionaries where the keys are the obsids and values are the various directory/file names. If fine_correct has been set to
# 'use', we use the corrected event files instead.
repro_dirs  = {}
evt_files   = {}
asol_files  = {}
bpix_files  = {}
mask_files  = {}

for obsid in obsids:
    repro_dir   = over_dir + '/' + obsid + '/' + 'repro_'+band_check+'-band__'+filter_check+'-filter/'
    # If we're using the corrected files - use the new ones for all obsids in obsids_to_correct, and use the original uncorrected
    # ones otherwise
    if fine_correct  == 'use': 
        if obsid in obsids_to_correct:
            evt_file    = analysis_dir + '/' + 'finecor/' + obsid + '_new_evt2_finecor.fits'
            asol_temp   = dmkeypar(evt_file, 'ASOLFILE', echo=True)
            asol_file  = analysis_dir + '/' + 'finecor/' + asol_temp
        else:
            evt_file    = repro_dir + 'new_evt2.fits'
            asol_temp   = dmkeypar(evt_file, 'ASOLFILE', echo=True)
            asol_file   = repro_dir + asol_temp
    else: # Otherhwise, use original, uncorrected event files
        evt_file    = repro_dir + 'new_evt2.fits'
        asol_temp   = dmkeypar(evt_file, 'ASOLFILE', echo=True)
        asol_file   = repro_dir + asol_temp

    # Bad pixel and mask files will be the same regardless, so we use the ones from the original directory
    bpix_temp   = dmkeypar(evt_file, 'BPIXFILE', echo=True)
    bpix_file   = repro_dir + bpix_temp

    mask_temp   = dmkeypar(evt_file, 'MASKFILE', echo=True)
    mask_file   = repro_dir + mask_temp

    repro_dirs[obsid]   = repro_dir
    evt_files[obsid]    = evt_file
    asol_files[obsid]   = asol_file
    bpix_files[obsid]   = bpix_file
    mask_files[obsid]   = mask_file

print("Analyzing OBSIDs "+','.join(obsids))
print("Filter band: ",filter_check,"kev")
print("Flux band:   ",band_check,"kev")
print("fine_correct:",fine_correct)
print("reference_obsid:",reference_obsid)
print("Using event files:")
for obsid,evt_file in evt_files.items():
    print(f"For {obsid}: {evt_file}")
print("Using asol files:")
for obsid,asol_file in asol_files.items():
    print(f"for {obsid}: {asol_file}")
print("Using bpix files:")
for obsid,bpix_file in bpix_files.items():
    print(f"for {obsid}: {bpix_file}")
print("Using mask files:")
for obsid,mask_file in mask_files.items():
    print(f"for {obsid}: {mask_file}")

# We define a function here for getting the centroid of a source (around a point). Note that the point needs to be in the format of
# 'hh:mm:ss.ss,dd:mm:ss.ss'. We'll later use this for checking the astrometry alignment of the different observations.
def get_centroid(infile, point):
    
    #Get the total dmstat infile string
    dm_infile = infile + '[energy=500:7000,sky=circle('+point+',4")][bin sky=1]'
    
    #Run dmstat (have to run in terminal otherwise it doesn't work)
    sp.run(["dmstat",dm_infile,"centroid=yes","verb=0"])
    
    #Use pget to get coords
    proc = sp.check_output(["pget","dmstat","out_cntrd_phys"])
    output = str(proc, 'utf-8').rstrip()
    coords = output.split(',')
    xval = coords[0]
    yval = coords[1]
    
    #Run dmcoords to turn physical (pixel) into ra,dec
    #Have to run in terminal otherwise it doesn't work
    sp.run(["dmcoords",infile,"op=sky","x="+xval,"y="+yval,"celfmt=hms","v=0"])
    
    #Retrieve coords using pget again
    proc = sp.check_output(["pget","dmcoords","ra","dec"])
    output = str(proc, 'utf-8').rstrip()
    radec = output.split('\n')
    ra = radec[0]
    dec= radec[1]
    return ra,dec


#############################################################
################## MERGING and WAVDETECT ####################
#############################################################

# Move into the analysis directory/folder
os.chdir(analysis_dir)

# Set the folder for storing merged things, depending on if fine_correct has been set to 'use'
if fine_correct == 'use':
    merge_folder = 'merged_finecor/'
else:
    merge_folder = 'merged/'

# Merge the two event files
merge_obs.punlearn()
merge_obs.infiles       = ','.join(list(evt_files.values())) # gets list of 'evtfile1,evtfile2,...'
merge_obs.outroot       = merge_folder
merge_obs.bands         = filt_band
merge_obs.binsize       = 1
merge_obs.asolfiles     = ','.join(list(asol_files.values()))
merge_obs.badpixfiles   = ','.join(list(bpix_files.values()))
merge_obs.maskfiles     = ','.join(list(mask_files.values()))
merge_obs.psfecf        = fluximage_psfecf
merge_obs.clobber       = True
if merge_folder[:-1] not in os.listdir():
    print("Running merge_obs...",end=' ',flush=True)
    merge_obs()
    print("Done")
else:
    print("merge_obs already run - skipping")


# Now that merging is done, we'll need to define our galaxy region (as we use this to decide what wavdetect regions we should focus on
# later). We first retrieve the coordinates of the galaxy regions from the observations we're merging, then average these together to
# get the final coordinates. We also check whether these are the same between observations (to within a limit).
# First initialize some dictionaries
gal_ra  = {}
gal_dec = {}
gal_rad = {}

for obsid in obsids:
    
    gal_reg_file = repro_dirs[obsid] + 'gal_ds9_region.reg'
    
    with open(gal_reg_file,'r') as galfile:
        galstring = galfile.read() # Read in file
        galstring = galstring.split('\n') # Split on newlines
        galstring = galstring[2] # the 3rd line in the file is the coordinates, so focus on that
        
        # The region string looks like 'circle(ra,dec,radius)' (all in degrees) so we remove the circle() part to focus on
        # the coordinates. We then split along the commas
        galstring = galstring[ galstring.find('(')+1: galstring.rfind(')') ]
        galstring = galstring.split(',')
        
        # galstring is now a list of [ra,dec,rad] (all as strings, in degrees), which we assign to the dicts as appropriate
        gal_ra[obsid]  = float(galstring[0])
        gal_dec[obsid] = float(galstring[1])
        gal_rad[obsid] = float(galstring[2])

        
# Average all the values out to get our final values
gal_ra_avg  = np.mean(list(gal_ra.values()))
gal_dec_avg = np.mean(list(gal_dec.values()))
gal_rad_avg = np.mean(list(gal_rad.values()))

# Find the standard deviation (std) of the values
gal_ra_std  = np.std(list(gal_ra.values()))
gal_dec_std = np.std(list(gal_dec.values()))
gal_rad_std = np.std(list(gal_rad.values()))

# If std is above the max value, print that to the terminal
if gal_ra_std > max_std:
    print("Galaxy regions had differing RAs above user-specified limit set by max_std")
    for obsid in obsids:
        print("OBSID:",obsid,"RA:",gal_ra[obsid],"deg")
if gal_dec_std > max_std:
    print("Galaxy regions had differing Decs above user-specified limit set by max_std")
    for obsid in obsids:
        print("OBSID:",obsid,"Dec:",gal_dec[obsid],"deg")
if gal_rad_std > max_std:
    print("Galaxy regions had differing radii above user-specified limit set by max_std")
    for obsid in obsids:
        print("OBSID:",obsid,"radius:",gal_rad[obsid],"deg")

print("Using galaxy region with RA:",gal_ra_avg,"Dec:",gal_dec_avg,"Radius (deg):",gal_rad_avg)


# Turn the galaxy ra/dec into pixel coordinates (based on the _thresh image that we'll use for wavdetect)
# First run dmcoords to turn into pixel coords, then run pget to extract the coordinates
sp.run(["dmcoords","infile="+merge_folder+filter_check+"_thresh.img","op=cel","ra="+str(gal_ra_avg),"dec="+str(gal_dec_avg)])
proc = sp.check_output(["pget","dmcoords","x","y"])
output = str(proc, 'utf-8').rstrip() # Turn pget output into string
gal_xy = output.split('\n') # Split into x and y entries
gal_x = round(float(gal_xy[0]),4)
gal_y = round(float(gal_xy[1]),4)

# Turn the galaxy region radius into pixels. The 3600 takes it from degrees to arcsec, then / 0.492 takes it from arcsec to chandra
# pixels, since 0.492 arcsec = 1 chandra pixel.
gal_rad_pix = round(gal_rad_avg * 3600 / 0.492,4)

# Finally, create the galaxy region file
gal_file = open("final_gal_reg.reg","w")
gal_file.write("Circle("+str(gal_x)+","+str(gal_y)+","+str(gal_rad_pix)+")")
gal_file.close()


# Now we can run wavdetect on the merged event file
if fine_correct == 'use':
    wavdet_folder = 'wavdet_finecor/'
else:
    wavdet_folder = 'wavdet/'
wavdetect.punlearn()
wavdetect.infile        = merge_folder+filter_check+'_thresh.img'
wavdetect.psffile       = merge_folder+filter_check+'_thresh.psfmap'
wavdetect.expfile       = merge_folder+filter_check+'_thresh.expmap'
wavdetect.outfile       = wavdet_folder+'wav_srcs.fits'
wavdetect.scellfile     = wavdet_folder+'scell.fits'
wavdetect.imagefile     = wavdet_folder+'imgfile.fits'
wavdetect.defnbkgfile   = wavdet_folder+'nbgd.fits'
wavdetect.scales        = wavdet_analysis_scales
wavdetect.sigthresh     = 1e-06
wavdetect.clobber       = True
if wavdet_folder[:-1] not in os.listdir():
    sp.run(["mkdir", "-p", os.getcwd() + '/' + wavdet_folder])
    print("Running merged wavdetect...",end=' ',flush=True)
    wavdetect()
    print("Done")
else:
    print("Wavdetect already run - skipping")


# Now we have to analyze the wavdetect sources. We first filter the sources to keep only the ones inside the galaxy, as those are the
# ones we care about. We then extract the positions of these sources and use the positions to check how well the event files are aligned
# astrometrically. 

# We start by copy-pasting some code from xray_flux.py (also up on github). This bit filters the wavdetect regions to only keep ones
# within the galaxy. Comments have been included for clarification.
#----
# Now finding out which of the wavdetect sources is within our galaxy region. This involves
# turning the wavdetect regions and the (corrected) galaxy region into python objects via the
# region package to examine their properties. We get the ids of the wavdetect regions inside
# the galaxy region, and copy the wavdetect source region file over using dmcopy to keep
# only those regions (the [1:-1] is used to get rid of the brackets from converting the list
# to a string).
wav_regions = region.CXCRegion(wavdet_folder+'wav_srcs.fits')
gal_cor_reg = region.CXCRegion('final_gal_reg.reg')

id_list = list([])
for m in range(0,len(wav_regions)):
    # Using .shapes gives us a list of objects, and we pick out one to focus on.
    indiv_reg = wav_regions.shapes[m]
    
    xp = indiv_reg.xpoints
    yp = indiv_reg.ypoints #getting x,y coords of the center of the wavdetect source
    
    # If the wavdetect region's center is inside the galaxy region, add its id to the list
    if gal_cor_reg.is_inside(xp,yp):
        id_list.append(indiv_reg.component)

dmcopy.punlearn()
dmcopy.infile   = wavdet_folder+'wav_srcs.fits[#row='+str(id_list)[1:-1]+']' 
dmcopy.outfile  = wavdet_folder+'wav_srcs_inside.fits'
dmcopy.clobber  = True
dmcopy()
#----

# If we found no wavdetect sources, then there is no reason to continue so we exit out of the script. Otherwise, we print out how many
# wavdetect sources we found.
if len(id_list) == 0:
    print("No wavdetect sources found in galaxy region.")
    os.chdir(code_dir)
    summary_file    = open('results_merged_'+band_check+'-band__'+filter_check+'-filter.txt','w')
    summary_file.write("No sources detected!")
    summary_file.close()
    sys.exit()
else:
    print(str(len(id_list))+' wavdetect sources found in galaxy region')


# We now retrieve the wavdetect source locations (ra,dec) and turn them into a ciao-region-friendly coordinate format e.g. 
# 'hh:mm:ss.ss,dd:mm:ss.ss' and store them in a list. We'll then use these positions as starting points for checking the centroids of
# each source in each separate event file.
pos_string = []
for m in range(0,len(id_list)):
        #First we focus in on a single source, making a new file with only that source in it.
        dmcopy.punlearn()
        dmcopy.infile   = wavdet_folder+'wav_srcs_inside.fits[#row='+str(m+1)+']'
        dmcopy.outfile  = wavdet_folder+'wav_src'+str(m+1)+'_inside.fits'
        dmcopy.clobber  = True
        dmcopy()
        
        #Now we get the ra and dec of that specific source and store it
        dmlist.punlearn()
        temp_dmlist = dmlist(wavdet_folder+'wav_src'+str(m+1)+'_inside.fits[cols ra,dec]', opt='data,clean')
        temp_dmlist = temp_dmlist.split()
        wavdet_ra = float(temp_dmlist[3])
        wavdet_dec = float(temp_dmlist[4])
        print(f"Merged wavdet source at {wavdet_ra} {wavdet_dec}")
        coords = SkyCoord(wavdet_ra,wavdet_dec,frame=coord_frame,unit='deg')
        coords_str = coords.to_string('hmsdms',precision=3)
        coords_str = coords_str.replace('h',':')
        coords_str = coords_str.replace('d',':')
        coords_str = coords_str.replace('m',':')
        coords_str = coords_str.replace(' ',',')
        coords_str = coords_str.replace('s','')
        pos_string.append(coords_str)


# We now check the alignment of the astrometry between different observations. Since each observation is of the same object, we should
# be detecting the same sources. The wavdetect sources that we've found so far were found by running wavdetect on the merged image - 
# which is good and necessary for getting merged fluxes etc., but we still need to be sure that the separate observations are 
# astrometrically aligned, as if they aren't it could lead to incorrect results. We check for offsets by locating the same celestial
# source in each event file, calculating the centroid of that source, then comparing the centroids across all the different event files.
# Ideally, they should differ by much less than a chandra-sized pixel (0.492 arcsec). If the offsets are greater than that, then we quit
# out of the script as proceeding would yield bad results.

# We initialize a 3d array to store the centroid separations in. Here, s[a][b][c] gives the separation (in arcseconds) of the ath source
# between the bth and cth obsid. For example, s[0][0][1] focuses on wavdetect source 0, and gives the offset between the source's 
# centroid in obsid 0 and in obsid 1. We want these offsets to be small (much less than a chandra pixel) to be sure that the merging
# will give us good results. Technically this could have been done before merging and wavdetect, but that would have required the user
# to input the source location manually. Instead we use the merged wavdetect sources as the source locations, which lets us do it all
# automatically.
# Since finding the separation between two points is symmetric (e.g. pos_b.separation(pos_c) = pos_c.separation(pos_b)) this means that
# s[a][x][y] = s[a][y][x] i.e. the separation matrix is symmetric in the latter two indices. Additionally, s[a][y][y] = 0 since the 
# distance between a point and the same point is obviously zero. To avoid doing unnecessary work, we restrict the third for loop (for
# index c) to only run over the obsids that the second for loop (index b) hasn't touched yet.

# We also store the positions of the source centroids in each event file (for each obsid). We can then later use these to automate
# fine-corrections to the astrometry of various event files.
print("--------")
print("Using wavdetect sources to check for potential astrometric offsets between OBSIDs...")
print("Maximum separation allowed:",max_sep,"arcsec")
s = np.zeros(( len(pos_string),len(obsids),len(obsids) ))
source_centroids = {}
for a in range(0,len(pos_string)):
    for b in range(0,len(obsids)):
        
        obsid_b = obsids[b]
        ra_b,dec_b = get_centroid(evt_files[obsid_b],pos_string[a])
        pos_b = SkyCoord([ra_b],[dec_b],unit=(u.hourangle,u.deg)) # Put ra, dec in [] so later to_table() will work
        
        source_centroids[obsids[b]] = pos_b # Store centroid position
        pos_b = pos_b[0] # Return to normal SkyCoord representation once we've stored it
        
        for c in range(b+1,len(obsids)):
            
            obsid_c = obsids[c]
            ra_c,dec_c = get_centroid(evt_files[obsid_c],pos_string[a])
            pos_c = SkyCoord(ra_c,dec_c,unit=(u.hourangle,u.deg))
            
            print("OBSID " + str(obsids[b]) + " centroid at " + pos_b.to_string('decimal',precision=10))
            print("OBSID " + str(obsids[c]) + " centroid at " + pos_c.to_string('decimal',precision=10))
            
            sep_bc = pos_b.separation(pos_c)
            sep_bc = round(sep_bc.arcsecond,3)
            
            s[a][b][c] = sep_bc
            
            print("Source",a+1,"had separation of",sep_bc,"arcsec between OBSIDs "+str(obsids[b])+', '+str(obsids[c]))


# Now check to see if any separation exceeds the maximum. sep_over has the shape (3,n) where n is how many excessive offsets were found. 
# Note that even if we have more than 2 obsids (so we're merging 3+ observations) sep_over will still have the shape (3,n). This is 
# because we only compare the distance between two centroids at a time. sep_over[0] gives the list of sources with bad offsets, and
# sep_over[1] and sep_over[2] give the corresponding list of obsids that were compared. If the same source was found to have bad offsets
# between multiple different obsids (e.g. if the centroids of source 0 in obsids 0, 1, and 2 were all too far apart) then the source will
# appear multiple times on the source list, which would look like (for that example):
# sep_over[0] = np.array([0, 0, 0]) -- centroids of source 0 were bad between 3 different pairs of obsids
# sep_over[1] = np.array([0, 0, 1]) -- list of the 1st of each pair of obsids with bad offsets
# sep_over[2] = np.array([1, 2, 2]) -- list of the 2nd of each pair of obsids with bad offsets
# This tells us (read from top to bottom, one column at a time):
# Source 0 had an excessive offset between obsid 0 and obsid 1
# Source 0 had an excessive offset between obsid 0 and obsid 2
# Source 0 had an excessive offset between obsid 1 and obsid 2
sep_over        = np.where(s >= max_sep)
sep_over_vals   = s[sep_over] # Get the actual separation values (in arcsec)

# If there are no excessive offsets detected, continue. If there are, then we have two options. If the user has set fine_correct
# to 'make', then we go through and automate fine corrections, correcting the obsids to reference_obsid, so the astromety will match.
# If the user has set fine_correct to anything other than 'make', then we just quit the script, since merging between images with
# offsets gives bad results. The exception to this is if the user has set force_continue to True, in which case we continue on
# with analysis regardless.
if len(sep_over_vals) > 0:
    print("Excessive offsets found!")
    print(len(sep_over_vals),"excessive offsets found between",len(np.unique(sep_over[0])),"sources")
    for a in range(0,len(sep_over_vals)):
        src     = sep_over[0][a]
        obsid1  = sep_over[1][a]
        obsid2  = sep_over[2][a]
        
        print("wavdetect source",src,"had offset of",sep_over_vals[a],"arcsec between OBSIDs",obsids[obsid1],",",obsids[obsid2])
    print("Merged source analysis using current event files will give incorrect results")
    
    # Checking whether to do fine corrections
    if fine_correct != 'make':
        print("User will need to either adjust the fine_correct and reference_obsid variables accordingly to automate astrometric correction or manually correct the astrometry of the offending OBSIDs, then run this script again.")
        if force_continue == True:
            print("force_continue set to True - proceeding regardless")
        else:
            sys.exit()
    elif fine_correct == 'make':
        print("Attempting to make fine corrections to astrometry...")
        print("Reference OBSID: " + reference_obsid)

        # Make the directory to store all the fine astrometry correction data in
        if 'finecor' not in os.listdir():
            sp.run(["mkdir", "-p", os.getcwd() + '/finecor'])
        os.chdir('finecor')

        # Make .tsv files with centroid positions in them - we'll use these in wcs_match down below
        source_centroid_files = {} # Dictionary to fill with obsid:filename
        for obsid,centroid_pos in source_centroids.items():
            # Fill dictionary
            filename = obsid+'_centroid_pos.tsv'
            source_centroid_files[obsid] = filename

            # Write out source centroid position to .tsv file. We have to do it this weird way because the standard pandas .to_csv
            # won't let us put a '#' in the header, which ciao wants to be able to read it
            t=QTable(data=[centroid_pos])
            with open(filename,'w') as file:
                file.write('#ra dec\n')
                t.to_pandas().to_csv(file,index=False,header=False,sep='\t')
            
        # Now actually matching and updating astrometry. We match every non-reference obsid to the reference obsid, using
        # the centroid positions as our source positions. We create new event files for the fine corrections, and update the
        # astrometry on those.
        reference_source_file = source_centroid_files[reference_obsid] # Setting the ref source file

        for obsid,centroid_pos_file in source_centroid_files.items():
            
            # First, skip through correcting astrometry if it's the reference obsid (as we're assuming that is correct), and skip
            # it if the obsid isn't in obsids_to_correct
            if obsid == reference_obsid:
                continue
            elif obsid not in obsids_to_correct:
                continue
            
            print("Correcting astrometry for OBSID",obsid)

            # Next, assign variables for all the different files we'll need
            # infile          = obsid + '.dat' ##centroid_pos_file
            evt             = evt_files[obsid]            # Original event file
            asol            = asol_files[obsid]             # Original asol file
            evt_cor         = obsid + '_new_evt2_finecor.fits'  # Corrected event file
            transform       = obsid + '_out.xform'
            logmatch        = obsid + '_logmatch.txt'
            logupdate       = obsid + '_logupdate.txt'
            logupdate_asol  = obsid + '_logupdate_asol.txt'

            # Match between the observations
            wcs_match.punlearn()

            wcs_match.infile        = centroid_pos_file
            wcs_match.refsrcfile    = reference_source_file
            wcs_match.wcsfile       = evt
            wcs_match.outfile       = transform
            wcs_match.radius        = 4
            wcs_match.method        = 'trans'
            wcs_match.clobber       = True
            wcs_match.verbose       = 5
            wcs_match.logfile       = logmatch

            wcs_match()

            # Create a new copy of the original (uncorrected) event file that we can apply corrections to
            dmcopy.punlearn()
            dmcopy.infile   = evt
            dmcopy.outfile  = evt_cor
            dmcopy.clobber  = True
            dmcopy()

            # Getting the name of the asolfile for the new event file
            asol_name        = dmkeypar(evt_cor, 'ASOLFILE', echo=True)
            asol_name_cor    = asol_name[:-5] + '_cor.fits'

            # Correcting the astrometry of the new event file
            wcs_update.punlearn()
            wcs_update.infile       = evt_cor
            wcs_update.outfile      = evt_cor
            wcs_update.transformfile= transform
            wcs_update.wcsfile      = evt
            wcs_update.logfile      = logupdate
            wcs_update.clobber      = True
            wcs_update.verbose      = 5
            wcs_update()

            # Updating the asol file to be used in the new event file
            wcs_update.punlearn()
            wcs_update.infile       = asol
            wcs_update.outfile      = asol_name_cor
            wcs_update.transformfile= transform
            wcs_update.wcsfile      = evt
            wcs_update.logfile      = logupdate_asol
            wcs_update.clobber      = True
            wcs_update.verbose      = 5
            wcs_update()

            # Updating the ASOLFILE keyword in the new event file to point to the corrected asol file
            dmhedit.punlearn()
            dmhedit.infile      = evt_cor
            dmhedit.filelist    = 'none'
            dmhedit.operation   = 'add'
            dmhedit.key         = 'ASOLFILE'
            dmhedit.value       = asol_name_cor
            dmhedit()
        
        print("Fine corrections to astrometry completed - user will need to rerun this script with fine_correct='use'.")
        if force_continue == True:
            print("force_continue set to True - proceeding regardless")
        else:
            sys.exit()
else:    
    print("No excessive offsets found! Continuing with source analysis.")
    print("--------")


#############################################################
###################### SOURCE ANALYSIS ######################
#############################################################

os.chdir(analysis_dir)
# Now that we're sure that merging hasn't messed up the data, we can proceed onwards to analysis. We first have to make the source and
# background regions we'll use for srcflux.

# We start by creating a new directory to store source-related files in
sp.run(["mkdir", "-p", os.getcwd() + "/srcfiles"])

# The wavdetect sources gave us the positions of the source regions we'll use to find fluxes, and now we need to find the radii. While
# it would be nice to be able to use psfsize_srcs on the merged event file, that won't work because of merging weirdness. So instead, we
# take the list of merged wavdetect sources and run psfsize_srcs on that source list for each separate event file. This gives us a list
# of source regions (with the radii calculated based on the energy and ecf we put it) from each event file. We can then retrieve these
# radii and average them together (on a per-source basis) to get the final per-source radii to use for our source regions we'll input
# to srcflux. Additionally, since psfsize_srcs uses the RA,Dec columns of the wavdetect.fits file we don't have to worry about any 
# mismatches between physical/pixel coordinates when going between individual event files and the merged event file.
for obsid in obsids:
    psfsize_srcs.punlearn()
    psfsize_srcs.infile     = evt_files[obsid]
    psfsize_srcs.pos        = wavdet_folder+'wav_srcs_inside.fits'
    psfsize_srcs.outfile    = 'srcfiles/'+obsid+'psf_srcs.fits'
    psfsize_srcs.energy     = analysis_energy
    psfsize_srcs.ecf        = analysis_psfecf
    psfsize_srcs.clobber    = True
    psfsize_srcs()

# We initialize an array that looks like (e.g. for two obsids with three merged wavdetect sources):
#         src1    src2    src3
# obsid_1 [[       0      0      0    ]
# obsid_2  [      0      0      0    ]]
# 
# We'll then extract the calculated psfsize_srcs radii values and stick them in this array, row by row.
src_radii = np.zeros( (len(obsids),len(pos_string)) )
for n in range(0,len(obsids)):
    dmlist.punlearn()
    rad_dmlist = dmlist('srcfiles/'+obsids[n]+'psf_srcs.fits[cols r]', opt='data,clean')
    rad_dmlist = rad_dmlist.split() # Split by newline, results in list of form ['#', 'R', radius1, radius2,...]
    rad_dmlist = rad_dmlist[2:] # Removes the first 2 elements, the '#' and 'R', leaving us with a list of radii
    rad_dmlist = [float(a) for a in rad_dmlist] # Turn everything into floats
    rad_dmlist = np.array(rad_dmlist)
    src_radii[n] = rad_dmlist # Fill the first row (corresponding to obsid_1) with all the different source radii
    

# We average on a per-source basis, which results in an array where src_radii_avg[n] is the averaged radius (in pixels) for the nth
# source. We also find the standard deviation on a per-source basis. If the standard deviation is over a specified limit (which would
# mean the source radii were different above some limit, which might have implications for flux) we make note of this and store the
# index for later.
src_radii_avg = np.mean(src_radii,axis=0)
src_radii_std = np.std(src_radii,axis=0)

# Find out which (if any) the source radii differ and save this to a list that we'll print out to the results file. 
src_differ = np.where(src_radii_std >= max_std)[0]
src_differ_list = np.zeros(len(src_radii_avg))
src_differ_list[src_differ] = 1
src_differ_list = [int(a) for a in src_differ_list]

# If any source radii differ over the limit, print out to the screen which sources, and each source's radii across all the obsids.
if len(src_differ) > 0:
    
    # First print out which sources had radii over the limit
    src_differ_str = list(src_differ)
    src_differ_str = [str(a) for a in src_differ_str]
    src_differ_str = ','.join(src_differ_str)
    print("Sources "+src_differ_str+" had differing per-obsid radii above user-specified limit set by max_std")
    print("max_std set to:",max_std)
    
    # Loop through all the srcs that ticked that box, printing their calculated radii for the different obsids
    for m in range(0,len(src_differ)):
        
        # This gives us both source number and the index in src_radii that we'll need
        src_ind = src_differ[m]
        print("Source",src_ind,'radii:')
        print("   coords:",pos_string[src_ind])
        
        # Print out the radii for each obsid
        for n in range(0,len(obsids)):
            obsid = obsids[n]
            print("    OBSID -",obsid,':',src_radii[n][src_ind],'pixels')    
else:
    print("No sources had differing per-obsid radii above user-specified limit set by max_std")


# Now to make the actual source region files. We find the background annulus outer radius from the bg_radius_factor, and use the source
# radius as the inner radius. We'll also want both of these in arcsec for the region files. However, since we'll still want the region
# areas in pixels for later, we make the arcsecond conversion on a copy of the originals.
# Now to make the actual source region files. We first convert the source radii to arcseconds, and find the background annulus outer
# radius from that (the inner radius will be the source radii).
src_radii = src_radii_avg
bg_radii_out = src_radii * bg_radius_factor

src_radii_arcsec = np.round(src_radii * 0.492,3)
bg_radii_out_arcsec = np.round(bg_radii_out * 0.492,3)

# Write a region file for each source and background. Note that (since srcflux will be run on a stack of the different event files) we
# can't use pixel coordinates for these, since those will differ between the event files. Instead we have to put everything in
# celestial coordinates.
for n in range(0,len(pos_string)):
    srcreg_file = open('srcfiles/srcs.src'+str(n+1)+'.reg','w')
    srcreg_file.write('circle('+pos_string[n]+','+str(src_radii_arcsec[n])+'")\n')
    srcreg_file.close()
    
    bgreg_file = open('srcfiles/srcs.bg'+str(n+1)+'.reg','w')
    bgreg_file.write('annulus('+pos_string[n]+','+str(src_radii_arcsec[n])+'",'+str(bg_radii_out_arcsec[n])+'")\n')
    bgreg_file.close()


# We make list files of all the source and associated background region files, then run srcflux.
os.chdir('srcfiles')
src_list = glob.glob("srcs.src*.reg")
bg_list  = glob.glob("srcs.bg*.reg")

src_listfile = open('srcreg.lis','w')
bg_listfile  = open('bgreg.lis','w')

for n in range(0,len(src_list)):
    src_listfile.write(src_list[n] + '\n')
    bg_listfile.write(bg_list[n] + '\n')

src_listfile.close()
bg_listfile.close()
os.chdir('..')


# Running srcflux. Since srcflux can unfortunately not be run on the merged event file, we instead have to run it on a stack of all
# the separate event files.
srcflux.punlearn()
srcflux.infile      = ','.join(list(evt_files.values())) # Comma-separated string of event files
srcflux.pos         = wavdet_folder+'wav_srcs_inside.fits'
srcflux.outroot     = 'flux/'
srcflux.bands       = srcflux_band
srcflux.srcreg      = '@srcfiles/srcreg.lis'
srcflux.bkgreg      = '@srcfiles/bgreg.lis'
srcflux.psfmethod   = 'arfcorr'
srcflux.conf        = 0.9
srcflux.model       = 'xspowerlaw.pow1'
srcflux.paramvals   = 'pow1.PhoIndex=' + str(srcflux_PhoIndex)
srcflux.absmodel    = 'xsphabs.abs1'
srcflux.absparams   = 'abs1.nh=%GAL%' #using %GAL makes it retriev D&L maps
srcflux.mskfile     = ','.join(list(mask_files.values()))
srcflux.bpixfile    = ','.join(list(bpix_files.values()))
if 'flux' not in os.listdir():
    print('Starting merged srcflux...',end=' ',flush=True)
    srcflux()
    print('Done')
else:
    print('Srcflux already run - skipping')



# We now extract all the data from srcflux. This is a bit of a process as unfortunately the dmlist formatting needs some massaging. So
# we first output the data to a text file. This method is exactly the same as used in xray_flux.py, so the explanations here might be a 
# bit summarized.
dmlist.punlearn()
dmlist.infile   = 'flux/_'+band_check+'.flux[cols rapos,decpos,component,total_counts,total_bg_counts,num_obi,merged_net_rate_aper, umflux_cnv,merged_net_umflux_aper,merged_net_umflux_aper_lo,merged_net_umflux_aper_hi]'
dmlist.opt      = 'data,clean'
dmlist.outfile  = 'flux/data_list_init.txt'
dmlist()

# dmlist puts an unfortunate '#' at the beginning of the first line of the output file, which we have to remove or else it gums up the
# works when organizing things by columns etc.
token = open('flux/data_list_fin.txt','w')
sp.run(["tail", "-c", "+2", "flux/data_list_init.txt"], stdout = token)
token.close()

# Now we can read in the data. We do it line by line, using .split() to separate the columns.
token = open('flux/data_list_fin.txt')
linestoken = token.readlines()
resulttoken = list([])
for x in linestoken:
    resulttoken.append(x.split())
token.close()

# The previous part got us a list of lists where resulttoken[0] was all the parameter names, resulttoken[1] was all the parameter values
# for the first source, etc. We rearrange that so that data_list[0] is the first parameter name and its values for all sources, 
# data_list[1] is the second parameter name and its values for all sources, etc.
data_list = list(map(list, zip(*resulttoken)))

# Now breaking everything into a lot of smaller lists. We also get the source and background region areas (in pixels) for printing out.
param_list = resulttoken[0] # Gets the list of parameter names

ra_ind              = param_list.index('RAPOS')
dec_ind             = param_list.index('DECPOS')
comp_ind            = param_list.index('COMPONENT')
total_cnts_ind      = param_list.index('TOTAL_COUNTS')
total_bg_cnts_ind   = param_list.index('TOTAL_BG_COUNTS')
num_obi_ind         = param_list.index('NUM_OBI')
net_rate_ind        = param_list.index('MERGED_NET_RATE_APER')
umflux_cnv_ind      = param_list.index('UMFLUX_CNV')
net_umflux_ind      = param_list.index('MERGED_NET_UMFLUX_APER')
net_umflux_hi_ind   = param_list.index('MERGED_NET_UMFLUX_APER_HI')
net_umflux_low_ind  = param_list.index('MERGED_NET_UMFLUX_APER_LO')

ra              = np.asarray(data_list[ra_ind][1:],float)
dec             = np.asarray(data_list[dec_ind][1:],float)
comp            = np.asarray(data_list[comp_ind][1:],int)
total_cnts      = np.asarray(data_list[total_cnts_ind][1:],float)
total_bg_cnts   = np.asarray(data_list[total_bg_cnts_ind][1:],float)
num_obi         = np.asarray(data_list[num_obi_ind][1:],float)
net_rate        = np.asarray(data_list[net_rate_ind][1:],float)
umflux_cnv      = np.asarray(data_list[umflux_cnv_ind][1:],float)
net_umflux      = np.asarray(data_list[net_umflux_ind][1:],float)
net_umflux_hi   = np.asarray(data_list[net_umflux_hi_ind][1:],float)
net_umflux_low  = np.asarray(data_list[net_umflux_low_ind][1:],float)

src_areas       = math.pi*src_radii**2
bg_areas        = math.pi * (bg_radii_out**2 - src_radii**2)

ra_str              = data_list[ra_ind][0]
dec_str             = data_list[dec_ind][0]
comp_str            = data_list[comp_ind][0]
total_cnts_str      = data_list[total_cnts_ind][0]
total_bg_cnts_str   = data_list[total_bg_cnts_ind][0]
num_obi_str         = data_list[num_obi_ind][0]
net_rate_str        = data_list[net_rate_ind][0]
umflux_cnv_str      = data_list[umflux_cnv_ind][0]
net_umflux_str      = data_list[net_umflux_ind][0]
net_umflux_hi_str   = data_list[net_umflux_hi_ind][0]
net_umflux_low_str  = data_list[net_umflux_low_ind][0]



# Here we estimate the net counts. Since this uses the source and background areas in pixels, and assumes that the psffrac of each source
# is what we injected (analysis_psfecf), this is only a reliable estimate when the sources are near the aimpoint of the S3 chip and if
# the psffrac is indeed equal to analysis_psfecf. Usually, we can extract this information directly from the srcflux output file, but
# that only works when analyzing individual observations, instead of a stack - when analyzing a stack, srcflux doesn't output most of the
# source-specific information that we need, so we have to make do.

# Much of this next bit of code is copy-pasted from the chandra_xray_analysis repo, with a few changes to variables to account for some
# of the different variable names in this code (e.g. counts -> total_cnts). We include the comments from the original code for help
# understanding what's happening.
#------------------
#So first we'll calculate the errors on net counts. We do this by plugging in the counts
#(not background subtracted or aperture corrected) into the user-made error functions
#kraftcl and gehrelscl that automate finding errors using the methods of Kraft et al. 1991
#and Gehrels (1986). For sources with counts < 10 we use Kraft and take the background into
#account (for calculating errors). For sources with counts > 10 we use Gehrels, and assume
#the background is negligible.

#This gives us the upper and lower limits e.g. counts = 5^{+2.42}_{-1.74} (though note I just
#made those numbers up, I didn't run them through the error functions). Note that these limits
#are for the *counts*, not the net_counts, which means they still need to be background 
#subtracted and aperture corrected, so we do that next.

#For example: if we have counts = 7, bg = 0.5, aperture = 0.9. We plug into kraftcl and get
#limits of (2.97,11.87) for lower and upper (e.g. 7^{+4.87}_{-4.03}). We then background 
#subtract and aperture correct everything, counts and limits: (value - 0.5)/0.9. From this
#we get a value for the net counts of 7.22^{+5.41}_{-4.48}.


#Then, we also consider the error bars on the upper and lower intervals, and if the ratio of
#the two is less than sqrt(2), we consider them to be the same to within uncertainty, so we
#add them in quadrature and use the result for both upper and lower intervals.

#Following the previous example, the intervals are 5.41 and 4.48 for upper and lower, 
#respectively. Their ratio is 5.41/4.48 = 1.21 < sqrt(2), so we add them in quadrature to get
#4.97, which we use for both intervals, giving us a final value of 7.22^{+4.97}_{-4.97}. 

#While we're doing this we also check whether any of the source counts are below the 
#background to within the 95% confidence level (for Kraft). If for instance, the source 
#counts are below the background (to within error) we'll mark it as a bad source. We include
#it in analysis for convenience, but we'll keep track for the final stages of writing
#everything out to a text document.

bg_counts_in_src = total_bg_cnts/bg_areas*src_areas # Estimating the number of bg counts expected to be in the source aperture

net_counts_calc      = list([]) #used to store bg/aperture corrected net counts
net_counts_errors = list([]) #used to store bg/aperture corrected net count error intervals
good_src_list      = list([]) #1 means good source, 0 means bad source

temp_current_dir = os.getcwd()
os.chdir(code_dir) #need to be in code dir. for error functions to work properly

for m in range(0,np.size(ra)):
    src_counts      = total_cnts[m]
    src_bg_counts   = bg_counts_in_src[m]
    src_psffrac     = aperture_correct # Here we assume the final psffrac is equal to what was specified by the user
    
    if src_counts < 10:
        lims = np.array(ps.kraftcl(src_counts,src_bg_counts,1)) #find lims
        
        #Check whether source is above background to with 95% confidence level
        checklims = ps.kraftcl(src_counts,src_bg_counts,2)
        if checklims[0] < src_bg_counts:
            good_src_list.append(0)
        else:
            good_src_list.append(1)
    else:
        lims = np.array(ps.gehrelscl(src_counts,1))
        good_src_list.append(1) #as 95%cl test only matters for kraft
    
    #background/aperture correcting
    cor_lims        = (lims - src_bg_counts)/src_psffrac
    cor_src_counts  = (src_counts - src_bg_counts)/src_psffrac
    
    #getting upper/lower intervals
    cor_lim_up = cor_lims[1] - cor_src_counts
    cor_lim_lo = cor_src_counts - cor_lims[0]
    
    #checking whether we need to add in quadrature or not
    if max(cor_lim_up/cor_lim_lo,cor_lim_lo/cor_lim_up) < math.sqrt(2):
        single_lim = math.sqrt( (cor_lim_up**2 + cor_lim_lo**2)/2 )
        cor_lim_up = single_lim
        cor_lim_lo = single_lim
    
    #Getting final intervals and adding them (and net counts) to the lists
    fin_lims = np.array([cor_lim_lo,cor_lim_up])
    net_counts_calc.append(cor_src_counts)
    net_counts_errors.append(fin_lims)

os.chdir(temp_current_dir) #move back to normal directory after calculating errors

# Done estimating net counts, now onto saving the results.


# Get a string for the flux reference, and get the source luminosities
flxref     = '(' + str(1/flux_ref) + ')'
lum    = net_umflux*4*math.pi*(distance*3.086e24)**2
log_lum    = np.log10(lum)

# Put all this data in a dictionary so we can write it out to a results file afterwards.
od = collections.OrderedDict()
od['OBSIDS']                = ','.join(obsids)
od['SRC']                   = comp
od[ra_str]                  = ra
od[dec_str]                 = dec
od[total_cnts_str]          = total_cnts
od[total_bg_cnts_str]       = total_bg_cnts
od['EST_NET_COUNTS']        = np.round(net_counts_calc,round_dec)
od['EST_NET_COUNTS_ERR']    = [(round(x[0],round_dec),round(x[1],round_dec)) for x in net_counts_errors]
od['SRC_AREA']              = np.round(src_areas,round_dec)
od['BG_AREA']               = np.round(bg_areas,round_dec)
od['SRC_RADIUS (PIX)']      = src_radii
od['SRC_RADIUS_STD']        = src_radii_std
od['SRC_RADIUS_DIF']        = src_differ_list
od['GOOD_SOURCE']           = good_src_list
od[num_obi_str]             = num_obi
od[net_rate_str]            = net_rate
od[umflux_cnv_str]          = umflux_cnv
od[net_umflux_str+flxref]   = np.round(flux_ref*net_umflux,round_dec)
#od[net_umflux_hi_str+flxref]    = np.round(flux_ref*net_umflux_hi,round_dec)
#od[net_umflux_low_str+flxref]    = np.round(flux_ref*net_umflux_low,round_dec)
od['LOG_LUMINOSITY']        = np.round(log_lum,round_dec)

# Move back out of analysis folder to the code directory to save the results
os.chdir(code_dir)

# Write the results out to a file
summary_df      = pd.DataFrame(od)
summary_df_str  = summary_df.to_string(justify='left',index=False)
summary_file    = open('results_merged_'+band_check+'-band__'+filter_check+'-filter.txt','w')
summary_file.write(summary_df_str)
summary_file.write('\n\n\n\n'+':::OBSIDS - the obsids that were merged. \n:::SRC - source number. \n:::EST_NET_COUNTS - the estimated net counts for each source. This value is estimated using the calculated source and background areas in pixels (SRC_AREA and BG_AREA), and for the aperture correction uses the user-specified aperture_correct variable. Normally all these values (region areas and psffrac) would be extracted from srcflux, but when running srcflux on a stack of event files (as is the case when merging) these values are not listed. This estimate is only likely to be good when the sources are near the aimpoint of the S3 chips in all observations, and if the psffracs of the source regions are close to what the user specified. This can be (and should be) checked against the results from analyzing the individual observations using chandra_xray_analysis. \n:::EST_NET_COUNTS_ERR - This was calculated in the same fashion as EST_NET_COUNTS, and should have all the same caveats applied to it. It is listed here in the same manner as in the results for chandra_xray_analysis (as follows). The gist of it is we calculate the errors based off of 90% confidence levels using the methods of Kraft et al. 1991 (for counts<10) and Gehrels 1986 (for counts>10). This gives us upper and lower error intervals, e.g. if NET_COUNTS is 3.01 and NET_COUNTS_ERR is (2.61,4.05), this means the net counts value (with error) is 3.01^{+4.05}_{-2.61}.  \n:::SRC_AREA - area of the source region, in pixels. Note that this was calculated purely from the radius, not extracted from srcflux results, as srcflux does not compute region areas when finding merged flux from multiple obsids. \n:::BG_AREA - same as SRC_AREA, in pixels. \n:::SRC_RADIUS - the radius for the circular source region. This was found by computing source radii for each separate obsid, then averaging them together. \n:::SRC_RADIUS_STD - the standard deviation of the different obsid-specific calculated source radii (that were then averaged together to give SRC_RADIUS). \n:::SRC_RADIUS_DIF - whether the standard deviation of the obsid-specific calculated source radii was above the user-defined limit of max_std. A 1 here means that they fell above this limit, a 0 means they didn\'t. \n:::GOOD_SOURCE - whether or not the source counts were above the detected background counts within a 95% confidence level. A 1 here means the source was good, a 0 means the source was bad and should likely not be considered in further analysis. \n:::NUM_OBI - number of observations where the source falls in the FOV. \n:::MERGED_NET_UMFLUX_APER - the unabsorbed (merged) flux as calculated by srcflux, in erg/s/cm**2, in multiples of the number in parentheses. For example, if the number in parantheses is (1e-15), and the value is 3.62, then the actual flux would be 3.62e-15 erg/s/cm**2.')
summary_file.close()


#=====================================================================================================
#=====================================================================================================
######################################################################################################
######################################################################################################
######################################################################################################
################################################ END #################################################
######################################################################################################
######################################################################################################
######################################################################################################
