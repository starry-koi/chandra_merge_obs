Merging and Analysis of Chandra Observations
===
This code merges multiple Chandra observations together, and detects and analyzes sources based on the merged data. It also calculates and applies fine astrometric corrections, if needed. This code is meant to be run after the individual observations have already been analyzed using the code in the `chandra_xray_analysis` repo.


Setup
---

You will need to download and install the current version of [CIAO](https://cxc.cfa.harvard.edu/ciao/). Additionally, this script was meant to be run after the user has already run `xray_flux.py` from the `chandra_xray_analysis` repo on all the individual observations. To this end the user should first run `chandra_xray_analysis` on all observations (without worrying about merging), then copy-paste the relevant OBSID folders to the same directory as the chandra_merge_obs/ folder. The user will also need to copy-paste the `params.config` file they used in the run of `chandra_xray_analysis` into the chandra_merge_obs/ folder.



Brief Rundown of the Code
---

* We start by retrieving the evt2, asol, bpix, and mask files associated with the observations. The locations of these files will be printed to the screen. If the `fine_correct` variable is set to 'make', the files will all the the original files from the individual analysis repro directories. If `fine_correct` is set to 'use', some of these files (event and asol) will be the ones that have had fine astrometric corrections applied to them.
* We then merge these event files (either the originals or the finecor ones) using the CIAO function `merge_obs`. All output and interim files created by the code will be placed in the created analysis/ folder, which is labelled according to the `filter_check` and `band_check` variables in the same manner as for the `chandra_xray_analysis` code.
* We define our 'galaxy region' using the already-made galaxy region files located on the OBSID repro directories. For each galaxy region file, we retrieve the RA, Dec, and radius, then average them together to get the final coordinates and radius we'll use for our galaxy region. If we find that the individual regions' coordinates/radii differ above a user-specified limit (set by the `max_std` variable), we print this to the screen.
* We run the CIAO function `wavdetect` on the merged data products to obtain a list of detected sources. We then focus in on only the sources that lie within the galaxy region.
* We check the astrometry alignment between the different observations. Using the positions of the (merged) wavdetect sources as a starting point, we search for each source's centroid in each separate `new_evt2.fits` file and compare them. If the separations exceed a user-specified limit (set by the `max_sep` variable) we print the offsets of each source between each pair of OBSIDs, and if `fine_correct` is set to 'make', make and apply the fine astrometry corrections (see this [CIAO thread](https://cxc.cfa.harvard.edu/ciao/threads/fluxes_multiobi/ for how that's done). We then quit the script - the user will need to rerun with `fine_correct` set to 'use' to use the corrected images.
* If the separations don't exceed the user-specified limit, or if `force_continue` has been set to "True", we proceed on with the analysis.
* We next create our source regions using a similar method as in `xray_flux.py`. Unfortunately this method will not work if used directly on the merged images, so instead we create one source region per source using each OBSIDs' `new_evt2.fits` file (at the coordinates of the merged wavdetect source), then combine them to find the final source region radii. If the individual radii differ above a user-specified limit (again set by the `max_std` variable) we print this to the screen, along with all the individual radii for that specific source.
* Finally, we get fluxes using `srcflux` and print out the results to a text file specific to the filtering band and flux band (`filter_check` and `band_check` variables respectively, in the same manner as in `chandra_xray_analysis`). This text file is located in the code folder.



Variables the User Will Have to Care About
---

Merging-specific parameters, as detailed in `params_merged.config`:
 
* `distance`
  * The distance to the galaxy, in Mpc. You'll have to fill this in again, instead of just pulling the distance from the copy-pasted `params.config` file.
* `band_check`
  * The same as it is for the `chandra_xray_analysis` repo, this is the energy band (in keV) you want to find source fluxes in, e.g. '0.5-2', '2-10', etc. You're required to put this in separately for the merging-specific code since there could be multiple repro/ folders in each OBSID from prior `chandra_xray_analysis` runs, so the code needs you to tell it where to look for files and the like.
* `aperture_correct`
  * The fractional aperture correction to apply (0.9 corresponds to 90%, etc.). Best way to determine this is look at what PSFFRAC is in the detected sources in the unmerged observations (in src_all.txt) and use that value here.
* `filter_check`
  * Similar as to `band_check`, except the filtering band instead of the flux band.
* `max_sep`
  * In arcseconds, this determines the maximum separation between detected source centroids from individual OBSIDs allowed before stopping the code due to astrometric alignment issues. Generally, this should be set to much less than a Chandra pixel (0.492 arcsec) and the default of 0.1 arcsec seems to do alright.
* `max_std`
  * A measure of how different averaged region radii can be before we start taking notice. When we average the different calculated radii, we also take their standard deviation (std), and if that exceeds this variable then we start printing information to the screen.
* `fine_correct`
  * Options of 'use' or 'make'. Used to automate running fine astrometric corrections between to-be-merged observations. Basically automates what goes on in this [CIAO thread](https://cxc.cfa.harvard.edu/ciao/threads/fluxes_multiobi/). See the Quirks section for things that can go wrong with the automated fine correction.
* `reference_obsid`
  * The OBSID to correct other observations to when making fine astrometry corrections. Required to be set for the code to run - if you're running the code for the first time and don't know whether you'll need to make fine astrometric corrections, just set it to one of the observations you're merging and don't worry about it - you can always change it afterwards and rerun the code.
* `obsids_to_correct`
  * List of OBSIDs to make and apply fine astrometric corrections to. If left blank, the code will correct all observations other than `reference_obsid`.
* `force_continue`
  * Options of "True" or "False". Normally, if the code detects fine correction is needed, it will do the fine correction, create the files, and then stop so you can rerun the code with `fine_correct` set to 'use' instead of 'make'. However, this can sometimes lead to erroneous results (see the Quirks section), so this option just forces the code to continue the merging process using the uncorrected images. Should be used sparingly, as it can give bad results if used wrong.

The user will want to make sure that the `band_check` and `filter_check` parameters match which repro directory in the OBSIDs they want to look in. 

Additionally, the user should have copy-pasted the `params.config` file from their run of `chandra_xray_analysis` to the code directory, as the merging code makes use of several parameters from that config file. We list the parameters out here, though these are parameters that the user is likely to not have to change for most runs of either `chandra_xray_analysis` or `chandra_merge_obs` (and thus won't have to care about).

`coord_frame`, `wavdet_analysis_scales`, `fluximage_psfecf`, `analysis_psfecf`, `analysis_energy`, `bg_radius_factor`, `srcflux_PhoIndex`, `round_dec`, `flux_ref`.


A Sample Walkthrough of a Code Run
---

Let's say you have three observations (OBSIDs 1, 2, 3) that you want to analyze. OBSIDs 1 and 3 are also observations of the same target, so you'd like to combine them for analysis. 

You first use the `chandra_xray_analysis` code on all observations. After this is done, create a new, empty directory and put the chandra_merge_obs/ folder in it. Then copy-paste OBSIDs 1 and 3 (the OBSIDs you want to merge) to this directory, so that the only folders in this directory are the code folder and the to-be-merged OBSID folders. Also copy-paste the `params.config` file (that was used for the `chandra_xray_analysis` code run) to the code folder.

Next go into the code folder and adjust the parameters in `params_merged.config`. Make sure that the `distance` variable is set to the correct distance for your galaxy, and that the `filter_check` and `band_check` variables match what you used for the `chandra_xray_analysis` code run. You'll also need to set `reference_obsid` to either 1 or 3 (one of the OBSIDs you're merging) for the code to proceed.

Finally, open a terminal in the code directory and initialize CIAO, then run the code as 

        $./merge_obs.py

As the code runs, it will print out where it is in the analysis to the terminal.


Understanding the Output
---

The code will create an analysis/ folder in the same directory as the OBSID and chandra_merge_obs folders - this contains all the files and folders that the code creates and can, for the most part, be safely ignored. The analysis/ folder follows the same naming schema as in `chandra_xray_analysis`, and depends on the `filter_check` and `band_check` variables. E.g. if `filter_check` = '0.5-7' and `band_check` = '2-10', the resulting analysis folder will be named 'analysis_2-10-band__0.5-7-filter'.

The results (source fluxes, etc.) will be printed to a .txt file located in the code folder (following the same naming schema as the analysis/ folder). There will only be one results file per run of the code (in comparison to `chandra_xray_analysis`, which creates multiple output files).


Rerunning the Code
---

This follows the same methodology as in `chandra_xray_analysis`, in that time-intensive CIAO tasks (such as `merge_obs`, `wavdetect`, and `srcflux`) are skipped if the corresponding folder is found to exist. Folder names and the tasks they correspond to are (all located within the relevant analysis/ folder):
* merged/ - `merge_obs` (code run with `fine_correct` = 'make')
* merged_finecor/ - `merge_obs` for finely-corrected images (code run with `fine_correct` = 'use'
* wavdet/ - `wavdetect` (code run with `fine_correct` = 'make')
* wavdet_finecor/ - `savdetect` for finely-corrected images (code run with `fine_correct` = 'use'
* flux/   - `srcflux`

For more detail on this method, see the 'Rerunning the Code' section of the `chandra_xray_analysis` repo. In general, if you want to completely rerun the code (from scratch) using the same `filter_check` and `band_check` variables, you should delete the relevant analysis/ folder and then run `merge_obs.py`.

For this code, there is also a finecor/ folder where all the information related to the fine astrometric corrections is located.

Quirks
---

* This code is meant to be run after running `chandra_xray_analysis` on the individual OBSIDs (see the Setup and Sample Walkthrough sections for more detail).
* The automated fine astrometry correction assumes that there are actual sources present in all the to-be-merged observations. If you have several observations, and only one of them has a detected source in it, the code will diligently search out source centroids in all the observations and try to correct the observations to `reference_obsid`. This can throw problems, since if there's no detected source in an observation, then the source centroid doesn't really exist, and will change location depending on what precise area of the image you're looking for it in. This can present as running through the code with `fine_correct` = 'make', making the fine corrections, rerunning the code with `fine_correct` = 'use', and lo and behold you're still detecting offsets between sources in different observations. Even after the fine corrections have been applied. If you want to check whether this is happening, you can go and look at the actual (unmerged) images, where the code is saying that the source centroids are, and see if it actually looks like a source is there. If something like this is actually what is happening, you can run the code with `fine_correct` = 'make' and `force_continue` = "True", which will make the code merge the original, non-fine-corrected images and do analysis on that.


If something went wrong
---

* The code stopped saying it couldn't find the event files (`evt2.fits`)
  * This generally happens when there's a mismatch between the merging variables `filter_check` and `band_check`, and the repro/ directory inside the OBSID folder. You'll either need to update the merging-specific variables, or rerun `chandra_xray_analysis` on the individual OBSIDs using the filter band and flux band you want.
* The code's terminal output says that 'Galaxy regions had differing RAs/Decs/radii above user-specified limit set by max_std'
  * The code finds the final galaxy region by averaging the galaxy regions from all the different to-be-merged OBSIDs. The code should have printed out all the RAs/Decs/radii that exceeded the limit.
* The code says 'Sources ... had differing per-obsid radii above user-specified limit set by max_std'
  * The code finds final source regions by calculating source regions on a per-OBSID basis, then combining them to find the radii to use for the final source region (the position is gotten from the merged wavdetect source). These generally shouldn't differ too much, assuming the galaxy of interest was at the aimpoint of the S3 chip for all observations. By default, `max_std` is set fairly low, so if the source radii don't look that different (< ~0.1 pixels) you're probably fine. If the radii have differences of ~1+ pixels, that's probably an indication that the methods used in this merging code might not work best for your case.
* An error involving:   obsids_to_correct.remove(reference_obsid)   ValueError: list.remove(x): x not in list
  * `reference_obsid` needs to be set to an OBSID that is one of the OBSIDs you're trying to merge. This error occurs when the code can't find `reference_obsid`.
