Merging and Analysis of Chandra Observations
===
This code merges different observations of the same celestial object together, and detects and analyzes sources based on the merged data. This code is meant to be run after the observations have already been analyzed using the code in the `chandra_xray_analysis` repo.


Setup
---

You will need to download and install the current version of [CIAO](https://cxc.cfa.harvard.edu/ciao/). Additionally, this script was meant to be run after the user has already run `xray_flux.py` from the `chandra_xray_analysis` repo on all the individual observations. To this end the user should first run `chandra_xray_analysis` on all observations (without caring about merging), then copy-paste the relevant OBSID folders to the same directory as the chandra_merge_obs folder. The user will also need to copy-paste the `params.config` file they used in the run of `chandra_xray_analysis` into the chandra_merge_obs folder.



Brief Rundown of the Code
---


* We look in the appropriate repro directory (as dictated by the `filter_check` and `band_check` variables) in each OBSID for the cleaned, astrometrically corrected event file created by `xray_flux.py`. This will be the `new_evt2.fits` file. We then merge these event files using the CIAO function `merge_obs`. All output and interim files created by the code will be placed in the created analysis/ folder, which is labelled according to the `filter_check` and `band_check` variables in the same manner as for the `chandra_xray_analysis` code.
* We define our 'galaxy region' using the already-made galaxy region files located on the OBSID repro directories. For each galaxy region file, we retrieve the RA, Dec, and radius, then average them together to get the final coordinates and radius we'll use for our galaxy region. If we find that the individual regions' coordinates/radii differ above a user-specified limit (set by the `max_std` variable), we print this to the screen.
* We run the CIAO function `wavdetect` on the merged data products to obtain a list of detected sources. We then focus in on only the sources that lie within the galaxy region.
* We check the astrometry alignment between the different observations. Using the positions of the (merged) wavdetect sources as a starting point, we search for each source's centroid in each separate `new_evt2.fits` file and compare them. If the separations exceed a user-specified limit (set by the `max_sep` variable) we print the offsets of each source between each pair of OBSIDs, then quit the script, as merging observations with too large an astrometric offset can lead to incorrect results. See the relevant [CIAO thread](https://cxc.cfa.harvard.edu/ciao/threads/fluxes_multiobi/) for how to manually correct the astrometry of observations.
* We next create our source regions using a similar method as in `xray_flux.py`. Unfortunately this method will not work if used directly on the merged images, so instead we create one source region per source using each OBSIDs' `new_evt2.fits` file (at the coordinates of the merged wavdetect source), then combine them to find the final source region radii.
* Finally, we get fluxes using `srcflux` and print out the results to a text file specific to the filtering band and flux band (`filter_check` and `band_check` variables respectively). This text file is located in the code folder.



Variables the User Will Have to Care About
---

Merging-specific parameters, as detailed in params_merged.config:
 
* `distance`
        * The distance to the galaxy, in Mpc. You'll have to fill this in again, instead of just pulling the distance from the copy-pasted `params.config` file.
* `band_check`
        * The same as it is for the `chandra_xray_analysis` repo, this is the energy band (in keV) you want to find source fluxes in, e.g. 0.5-2, 2-10, etc.
* `filter_check`
        * Same as for `chandra_xray_analysis`, band to filter the image by. 

The user will want to make sure that the `band_check` and `filter_check` parameters match which repro directory in the OBSIDs they want to look in. 

Additionally, the user should have copy-pasted the `params.config` file from their run of `chandra_xray_analysis` to the code directory, as the merging code makes use of several parameters from that config file. We list the parameters out here, though these are parameters that the user is likely to not have to change for most runs of either `chandra_xray_analysis` or `chandra_merge_obs` (and thus won't have to care about).

`coord_frame`, `wavdet_analysis_scales`, `fluximage_psfecf`, `analysis_psfecf`, `analysis_energy`, `bg_radius_factor`, `srcflux_PhoIndex`, `round_dec`, `flux_ref`.



A Sample Walkthrough of a Code Run
---

Let's say you have three observations (OBSIDs 1, 2, 3) that you want to analyze. OBSIDs 1 and 3 are also observations of the same target, so you'd like to combine them for analysis. 

You first use the `chandra_xray_analysis` code on all observations. After this is done, create a new, empty directory and put the chandra_merge_obs folder in it. Then copy-paste OBSIDs 1 and 3 (the obsids you want to merge) to this directory, so that the only folders in this directory are the code folder and the to-be-merged OBSID folders. Also copy-paste the `params.config` file (that was used for the `chandra_xray_analysis` code run) to the code folder.

Next go into the code folder and adjust the parameters in `params_merged.config`. Make sure that the `distance` variable is set to the correct distance for your galaxy, and that the `filter_check` and `band_check` variables match what you used for the `chandra_xray_analysis` code run. 

Finally, open a terminal in the code directory and initialize CIAO, then run the code as 

        $./merge_obs.py

As the code runs, it will print out where it is in the analysis to the terminal.


Understanding the Output
---

The code will create an analysis/ folder in the same directory as the OBSID and chandra_merge_obs folders - this contains all the files and folders that the code creates and can, for the most part, be safely ignored. The analysis/ folder follows the same naming schema as in `chandra_xray_analysis`, and depends on the `filter_check` and `band_check` variables. E.g. if `filter_check` = '0.5-7' and `band_check` = '2-10', the resulting analysis folder will be named analysis_2-10-band__0.5-7-filter.

The results (source fluxes, etc.) will be printed to a .txt file located in the code folder (following the same naming schema as the analysis/ folder). There will only be one results file per run of the code (in comparison to `chandra_xray_analysis`, which creates multiple output files).


Rerunning the Code
---

This follows the same methodology as in `chandra_xray_analysis`, in that time-intensive CIAO tasks (such as `merge_obs`, `wavdetect`, and `srcflux`) are skipped if the corresponding folder is found to exist. Folder names and the tasks they correspond to are (all located within the relevant analysis.../ folder):
* merged/ - `merge_obs`
* wavdet/ - `wavdetect`
* flux/   - `srcflux`


Quirks
---


If something went wrong
---

* The code stopped saying it couldn't find the event files (`evt2.fits`)
* The code's terminal output says that 'Galaxy regions had differing RAs/Decs/radii above user-specified limit set by max_std`
* The code says 'Excessive offsets found!' after running `wavdetect` and stops there
* The code says 'Sources ... had differing per-obsid radii above user-specified limit set by max_std'


