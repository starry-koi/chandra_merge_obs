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
from sh import gunzip
from astropy.io import fits
from astropy.coordinates import SkyCoord
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
	
	json_obj 		= cjson.load(paramfile)
	
	coord_frame 		= json_obj['coord_frame']
	wavdet_analysis_scales 	= json_obj['wavdet_analysis_scales']
	fluximage_psfecf 	= json_obj['fluximage_psfecf']
	analysis_psfecf 	= json_obj['analysis_psfecf']
	analysis_energy 	= json_obj['analysis_energy']
	bg_radius_factor 	= json_obj['bg_radius_factor']
	srcflux_PhoIndex	= json_obj['srcflux_PhoIndex']
	round_dec 		= json_obj['round_dec']
	flux_ref 		= float(json_obj['flux_ref'])

# We need additional parameters specific to merging, which we import from the params_merged.config file.
with open('params_merged.config', 'r') as paramfile:
	
	json_obj		= cjson.load(paramfile)
	
	distance		= json_obj['distance']
	filter_check		= json_obj['filter_check']
	band_check		= json_obj['band_check']
	max_sep			= json_obj['max_sep']
	max_std			= json_obj['max_std']
	

# While the following variables should be declared in the params_merged.config file, we keep this description here so that anyone
# working on the code has a convenient explanation without having to open another file.
'''
distance 	= 139.8 	# Distance to the galaxy, in Mpc. Since this script was meant to be run on multiple observations of a
				# single galaxy, this needs to be user input in params_merged.config, as pulling it from the copied
				# params.config would make life more difficult.

filter_check 	= '2-7'		# Same as the filter_check variable the user should be familiar with from xray_flux.py. Since users might
				# have run xray_flux.py multiple times in different bands, there could be more than one repro directory
				# in each obsid folder. This tells the code which one to focus on, as well as being used for the fluximage
				# part of merge_obs.

band_check 	= '2-10'	# Similar to above - this tells the code which repro directory to look for, and also gets plugged into
				# srcflux later in the code. Same as the band_check variable in xray_flux.py.

max_sep 	= 0.1 		# In arcseconds, this is the maximum separation tolerated between source centroids in different obsids.
				# Different observations of the same object might have slightly different astrometry, so in this code
				# we run wavdetect on the merged event file, then use this position to search for source centroids (in a
				# small radius around the position) in all the different obsid event files. This variable is the maximum
				# separation we can detect between the positions of those source centroids in different obsid event files
				# before terminating the script. The rule of thumb for this is that we need the astrometric differences
				# to be much less than a chandra pixel (0.492 arcsec), so the default value for is is 0.1 as ~1/5th of
				# a pixel. It probably shouldn't be changed to be much higher. If the code does detect too big an offset
				# and stops, the user will have to do astrometry realignment on the offending obsids.

max_std 	= 0.02		# Slightly complicated to explain, but essentially this is a measure of how different the calculated
				# source radii are allowed to be before we start caring about it. In xray_flux.py, we find source radii
				# by using psfsize_srcs to create a region that captures a certain fraction of the psf energy (given by
				# the variable analysis_psfecf) at a certain energy (given by analysis_energy). Since psfsize_srcs 
				# cannot be run on the merged event file, we instead run it on all the different event files from the
				# different obsids and average the calculated radii together to arrive at a final source region radius.
				# During this finding of the average, we also calculate the standard deviation (std) of all of the
				# calculated radii. This can likely be set higher than the current value of 0.02, but less variation
				# is naturally better. Note that we also use this as the limit for the std between the different galaxy 
				# regions.
'''

# Quickly checking that there isn't a mismatch between the filtering and flux bands
if filter_check == '0.5-2':
	if band_check != '0.5-2':
		print('Filter and band mismatch - please fix')
		print('Filter: '+filter_check)
		print('Band:   '+band_check)
		sys.exit()
if band_check == '0.5-2':
	if filter_check != '0.5-2':
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
if band_check == '0.5-2':				# The band and specific energy to use when
	srcflux_band = '0.5:2:1.56' # soft		# calculating fluxes. It's of the form
elif band_check == '0.5-7':				# (lower):(upper):(specific_energy), with
	srcflux_band = '0.5:7:2.3' # broad		# everything in keV. Since the specific
elif band_check == '0.5-8':				# energy changes depending on which band is
	srcflux_band = '0.5:8:2.3' # ~broad		# used, it seemed easier to just do if 
elif band_check == '2-10':				# statements then letting it be user-driven.
	srcflux_band = '2:10:4' # hard			# If you add a new energy band, you'll need
else:							# to update this part as well.
	print("band_check variable not recognized")
	sys.exit()

# Get the current working directory, and save the paths to various files for later use
code_dir = os.getcwd()

# Move up a directory out of the chandra_merge_obs folder, so we can access the obsid folders
os.chdir('..')

# Store the overall working directory
over_dir = os.getcwd()

# Get a list of the obsids directories only (excluding non-numeric directories)
obsids = next(os.walk('.'))[1]
obsids = [a for a in obsids if a.isdigit()]
obsids.sort()

# Fill lists with strings pointing to specific directories and files. We first get the repro directories, then save paths to the event
# files and each event file's associated asol file. We do this in a for loop so that it will work for more than combining just two
# observations. 
repro_dirs = []
evt_files  = []
asol_files = []
for n in range(0,len(obsids)):
	repro_dirs.append(over_dir + '/' + obsids[n] + '/' + 'repro_'+band_check+'-band__'+filter_check+'-filter/')
	evt_files.append(repro_dirs[n] + 'new_evt2.fits')
	asol_temp = dmkeypar(evt_files[n], 'ASOLFILE', echo=True)
	asol_files.append(repro_dirs[n] + asol_temp)

print("Analyzing OBSIDs "+','.join(obsids))
print("Filter band:",filter_check,"kev")
print("Flux band:  ",band_check,"kev")

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

# First create an analysis/ folder to store all of the files and everything in, and move into it.
analysis_dir = over_dir + '/analysis_'+band_check+'-band__'+filter_check+'-filter'
sp.run(["mkdir", "-p", analysis_dir])
os.chdir(analysis_dir)

# Merge the two event files
merge_obs.punlearn()
merge_obs.infiles 	= ','.join(evt_files) # gets list of 'evtfile1,evtfile2,...'
merge_obs.outroot 	= 'merged/'
merge_obs.bands 	= filt_band
merge_obs.binsize	= 1
merge_obs.asolfiles	= ','.join(asol_files)
merge_obs.psfecf	= fluximage_psfecf
merge_obs.clobber	= True
if 'merged' not in os.listdir():
	print("Running merge_obs...",end=' ',flush=True)
	merge_obs()
	print("Done")
else:
	print("merge_obs already run - skipping")


# Now that merging is done, we'll need to define our galaxy region (as we use this to decide what wavdetect regions we should focus on
# later). We first retrieve the coordinates of the galaxy regions from the observations we're merging, then average these together to
# get the final coordinates. We also check whether these are the same between observations (to within a limit).
# First initialize some arrays
gal_ra  = np.zeros( len(obsids) )
gal_dec = np.zeros( len(obsids) )
gal_rad = np.zeros( len(obsids) )

for n in range(0,len(obsids)):
	
	gal_reg_file = repro_dirs[n] + 'gal_ds9_region.reg'
	
	with open(gal_reg_file,'r') as galfile:
		galstring = galfile.read() # Read in file
		galstring = galstring.split('\n') # Split on newlines
		galstring = galstring[2] # the 3rd line in the file is the coordinates, so focus on that
		
		# The region string looks like 'circle(ra,dec,radius)' (all in degrees) so we remove the circle() part to focus on
		# the coordinates. We then split along the commas
		galstring = galstring[ galstring.find('(')+1: galstring.rfind(')') ]
		galstring = galstring.split(',')
		
		# galstring is now a list of [ra,dec,rad] (all as strings, in degrees), which we assign to the arrays as appropriate
		gal_ra[n]  = float(galstring[0])
		gal_dec[n] = float(galstring[1])
		gal_rad[n] = float(galstring[2])

		
# Average all the values out to get our final values
gal_ra_avg  = np.mean(gal_ra)
gal_dec_avg = np.mean(gal_dec)
gal_rad_avg = np.mean(gal_rad)

# Find the standard deviation (std) of the values
gal_ra_std  = np.std(gal_ra)
gal_dec_std = np.std(gal_dec)
gal_rad_std = np.std(gal_rad)

# If std is above the max value, print that to the terminal
if gal_ra_std > max_std:
	print("Galaxy regions had differing RAs above user-specified limit set by max_std")
	for n in range(0,len(obsids)):
		print("OBSID:",obsids[n],"RA:",gal_ra[n],"deg")
if gal_dec_std > max_std:
	print("Galaxy regions had differing Decs above user-specified limit set by max_std")
	for n in range(0,len(obsids)):
		print("OBSID:",obsids[n],"Dec:",gal_dec[n],"deg")
if gal_rad_std > max_std:
	print("Galaxy regions had differing radii above user-specified limit set by max_std")
	for n in range(0,len(obsids)):
		print("OBSID:",obsids[n],"radius:",gal_rad[n],"deg")

print("Using galaxy region with RA:",gal_ra_avg,"Dec:",gal_dec_avg,"Radius (deg):",gal_rad_avg)


# Turn the galaxy ra/dec into pixel coordinates (based on the _thresh image that we'll use for wavdetect)
# First run dmcoords to turn into pixel coords, then run pget to extract the coordinates
sp.run(["dmcoords","infile=merged/"+filter_check+"_thresh.img","op=cel","ra="+str(gal_ra_avg),"dec="+str(gal_dec_avg)])
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
wavdetect.punlearn()
wavdetect.infile	= 'merged/'+filter_check+'_thresh.img'
wavdetect.psffile	= 'merged/'+filter_check+'_thresh.psfmap'
wavdetect.expfile	= 'merged/'+filter_check+'_thresh.expmap'
wavdetect.outfile	= 'wavdet/wav_srcs.fits'
wavdetect.scellfile	= 'wavdet/scell.fits'
wavdetect.imagefile	= 'wavdet/imgfile.fits'
wavdetect.defnbkgfile	= 'wavdet/nbgd.fits'
wavdetect.scales	= wavdet_analysis_scales
wavdetect.sigthresh	= 1e-06
wavdetect.clobber	= True
if 'wavdet' not in os.listdir():
	sp.run(["mkdir", "-p", os.getcwd() + '/wavdet'])
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
wav_regions = region.CXCRegion('wavdet/wav_srcs.fits')
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
dmcopy.infile	= 'wavdet/wav_srcs.fits[#row='+str(id_list)[1:-1]+']' 
dmcopy.outfile	= 'wavdet/wav_srcs_inside.fits'
dmcopy.clobber	= True
dmcopy()
#----

# If we found no wavdetect sources, then there is no reason to continue so we exit out of the script. Otherwise, we print out how many
# wavdetect sources we found.
if len(id_list) == 0:
	print("No wavdetect sources found in galaxy region.")
	os.chdir(code_dir)
	summary_file	= open('results_merged_'+band_check+'-band__'+filter_check+'-filter.txt','w')
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
		dmcopy.infile	= 'wavdet/wav_srcs_inside.fits[#row='+str(m+1)+']'
		dmcopy.outfile	= 'wavdet/wav_src'+str(m+1)+'_inside.fits'
		dmcopy.clobber	= True
		dmcopy()
		
		#Now we get the ra and dec of that specific source and store it
		dmlist.punlearn()
		temp_dmlist = dmlist('wavdet/wav_src'+str(m+1)+'_inside.fits[cols ra,dec]', opt='data,clean')
		temp_dmlist = temp_dmlist.split()
		wavdet_ra = float(temp_dmlist[3])
		wavdet_dec = float(temp_dmlist[4])
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
print("--------")
print("Using wavdetect sources to check for potential astrometric offsets between OBSIDs...")
print("Maximum separation allowed:",max_sep,"arcsec")
s = np.zeros(( len(pos_string),len(obsids),len(obsids) ))
for a in range(0,len(pos_string)):
	for b in range(0,len(obsids)):
		for c in range(b+1,len(obsids)):
			
			ra_b,dec_b = get_centroid(evt_files[b],pos_string[a])
			ra_c,dec_c = get_centroid(evt_files[c],pos_string[a])
			
			pos_b = SkyCoord(ra_b,dec_b,unit=(u.hourangle,u.deg))
			pos_c = SkyCoord(ra_c,dec_c,unit=(u.hourangle,u.deg))
			
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
sep_over 	= np.where(s >= max_sep)
sep_over_vals 	= s[sep_over] # Get the actual separation values (in arcsec)

# If we detect excessive offsets, print them to the terminal and quit the script (since merging between images with offsets gives bad
# results). Otherwise, just continue.
if len(sep_over_vals) > 0:
	print("Excessive offsets found!")
	print(len(sep_over_vals),"excessive offsets found between",len(np.unique(sep_over[0])),"sources")
	for a in range(0,len(sep_over_vals)):
		src 	= sep_over[0][a]
		obsid1 	= sep_over[1][a]
		obsid2 	= sep_over[2][a]
		
		print("wavdetect source",src,"had offset of",sep_over_vals[a],"arcsec between OBSIDs",obsids[obsid1],",",obsids[obsid2])
	print("Merged source analysis will give incorrect results - user will need to manually correct the astrometry of the offending OBSIDs and run this script again.")
	sys.exit()
else:
	
	print("No excessive offsets found! Continuing with source analysis.")
	print("--------")


#############################################################
###################### SOURCE ANALYSIS ######################
#############################################################


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
for n in range(0,len(obsids)):
	psfsize_srcs.punlearn()
	psfsize_srcs.infile	= evt_files[n]
	psfsize_srcs.pos	= 'wavdet/wav_srcs_inside.fits'
	psfsize_srcs.outfile	= 'srcfiles/'+obsids[n]+'psf_srcs.fits'
	psfsize_srcs.energy	= analysis_energy
	psfsize_srcs.ecf	= analysis_psfecf
	psfsize_srcs.clobber	= True
	psfsize_srcs()

# We initialize an array that looks like (e.g. for two obsids with three merged wavdetect sources):
# 		src1	src2	src3
# obsid_1 [[ 	  0	  0	  0	]
# obsid_2  [	  0	  0	  0	]]
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
srcflux.infile		= ','.join(evt_files) # Comma-separated string of event files
srcflux.pos		= 'wavdet/wav_srcs_inside.fits'
srcflux.outroot		= 'flux/'
srcflux.bands		= srcflux_band
srcflux.srcreg		= '@srcfiles/srcreg.lis'
srcflux.bkgreg		= '@srcfiles/bgreg.lis'
srcflux.psfmethod	= 'arfcorr'
srcflux.conf		= 0.9
srcflux.model		= 'xspowerlaw.pow1'
srcflux.paramvals	= 'pow1.PhoIndex=' + str(srcflux_PhoIndex)
srcflux.absmodel	= 'xsphabs.abs1'
srcflux.absparams	= 'abs1.nh=%GAL%' #using %GAL makes it retriev D&L maps
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
dmlist.infile	= 'flux/_'+band_check+'.flux[cols rapos,decpos,component,total_counts,total_bg_counts,num_obi,merged_net_rate_aper, umflux_cnv,merged_net_umflux_aper,merged_net_umflux_aper_lo,merged_net_umflux_aper_hi]'
dmlist.opt	= 'data,clean'
dmlist.outfile	= 'flux/data_list_init.txt'
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

ra_ind 		  = param_list.index('RAPOS')
dec_ind 	  = param_list.index('DECPOS')
comp_ind	  = param_list.index('COMPONENT')
total_cnts_ind    = param_list.index('TOTAL_COUNTS')
total_bg_cnts_ind = param_list.index('TOTAL_BG_COUNTS')
num_obi_ind	  = param_list.index('NUM_OBI')
net_rate_ind	  = param_list.index('MERGED_NET_RATE_APER')
umflux_cnv_ind 	  = param_list.index('UMFLUX_CNV')
net_umflux_ind    = param_list.index('MERGED_NET_UMFLUX_APER')
net_umflux_hi_ind = param_list.index('MERGED_NET_UMFLUX_APER_HI')
net_umflux_low_ind= param_list.index('MERGED_NET_UMFLUX_APER_LO')

ra 		  = np.asarray(data_list[ra_ind][1:],float)
dec 		  = np.asarray(data_list[dec_ind][1:],float)
comp		  = np.asarray(data_list[comp_ind][1:],int)
total_cnts	  = np.asarray(data_list[total_cnts_ind][1:],float)
total_bg_cnts	  = np.asarray(data_list[total_bg_cnts_ind][1:],float)
num_obi		  = np.asarray(data_list[num_obi_ind][1:],float)
net_rate	  = np.asarray(data_list[net_rate_ind][1:],float)
umflux_cnv	  = np.asarray(data_list[umflux_cnv_ind][1:],float)
net_umflux 	  = np.asarray(data_list[net_umflux_ind][1:],float)
net_umflux_hi 	  = np.asarray(data_list[net_umflux_hi_ind][1:],float)
net_umflux_low 	  = np.asarray(data_list[net_umflux_low_ind][1:],float)

src_areas	  = math.pi*src_radii**2
bg_areas	  = math.pi * (bg_radii_out**2 - src_radii**2)

ra_str 		  = data_list[ra_ind][0]
dec_str 	  = data_list[dec_ind][0]
comp_str	  = data_list[comp_ind][0]
total_cnts_str	  = data_list[total_cnts_ind][0]
total_bg_cnts_str = data_list[total_bg_cnts_ind][0]
num_obi_str	  = data_list[num_obi_ind][0]
net_rate_str	  = data_list[net_rate_ind][0]
umflux_cnv_str	  = data_list[umflux_cnv_ind][0]
net_umflux_str 	  = data_list[net_umflux_ind][0]
net_umflux_hi_str = data_list[net_umflux_hi_ind][0]
net_umflux_low_str= data_list[net_umflux_low_ind][0]


# Get a string for the flux reference, and get the source luminosities
flxref 	= '(' + str(1/flux_ref) + ')'
lum	= net_umflux*4*math.pi*(distance*3.086e24)**2
log_lum	= np.log10(lum)

# Put all this data in a dictionary so we can write it out to a results file afterwards.
od = collections.OrderedDict()
od['OBSIDS']			= ','.join(obsids)
od['SRC']			= comp
od[ra_str]			= ra
od[dec_str]			= dec
od[total_cnts_str]		= total_cnts
od[total_bg_cnts_str]		= total_bg_cnts
od['SRC_AREA']			= np.round(src_areas,round_dec)
od['BG_AREA']			= np.round(bg_areas,round_dec)
od['SRC_RADIUS (PIX)']		= src_radii
od['SRC_RADIUS_STD']		= src_radii_std
od['SRC_RADIUS_DIF']		= src_differ_list
od[num_obi_str]			= num_obi
od[net_rate_str]		= net_rate
od[umflux_cnv_str]		= umflux_cnv
od[net_umflux_str+flxref]	= np.round(flux_ref*net_umflux,round_dec)
#od[net_umflux_hi_str+flxref]	= np.round(flux_ref*net_umflux_hi,round_dec)
#od[net_umflux_low_str+flxref]	= np.round(flux_ref*net_umflux_low,round_dec)
od['LOG_LUMINOSITY']		= np.round(log_lum,round_dec)

# Move back out of analysis folder to the code directory to save the results
os.chdir(code_dir)

# Write the results out to a file
summary_df 	= pd.DataFrame(od)
summary_df_str 	= summary_df.to_string(justify='left',index=False)
summary_file	= open('results_merged_'+band_check+'-band__'+filter_check+'-filter.txt','w')
summary_file.write(summary_df_str)
summary_file.write('\n\n\n\n'+':::OBSIDS - the obsids that were merged. \n:::SRC - source number. \n:::SRC_AREA - area of the source region, in pixels. Note that this was calculated purely from the radius, not extracted from srcflux results, as srcflux does not compute region areas when finding merged flux from multiple obsids. \n:::BG_AREA - same as SRC_AREA, in pixels. \n:::SRC_RADIUS - the radius for the circular source region. This was found by computing source radii for each separate obsid, then averaging them together. \n:::SRC_RADIUS_STD - the standard deviation of the different obsid-specific calculated source radii (that were then averaged together to give SRC_RADIUS). \n:::SRC_RADIUS_DIF - whether the standard deviation of the obsid-specific calculated source radii was above the user-defined limit of max_std. A 1 here means that they fell above this limit, a 0 means they didn\'t. \n:::NUM_OBI - number of observations where the source falls in the FOV. \n:::MERGED_NET_UMFLUX_APER - the unabsorbed (merged) flux as calculated by srcflux, in erg/s/cm**2, in multiples of the number in parentheses. For example, if the number in parantheses is (1e-15), and the value is 3.62, then the actual flux would be 3.62e-15 erg/s/cm**2.')
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
