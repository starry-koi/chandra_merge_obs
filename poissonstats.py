"""Collection of stats that I tend to use.  Including:
   Poisson confidence intervals from gehrels 1986 and Kraft et al. 1991.


   err_low, err_high = gehrelscl(n, cl) --- output Gehrels confidence limits (cl=0,1,2,3 for 2-sided 68, 90, 95, and 99% errors)
                               
    err_low, err_high = kraftcl(n, bg, cl) --- output Kraft condiecnce limits (cl=1,2,3 --> 90, 95, 99% condience [no 68%]).
    

"""

import numpy as np
from scipy.stats import stats
from astropy.stats import bootstrap
from astropy.io import ascii
import sys
import os

    
def gehrelscl(n, cl):
    """  tmperr_lower, tmperr_upper = gehrelscl(n,cl)

   Returns Confidence Limits with Poisson Statistics from Gehrels et al. 1986, ApJ, 303, 336.  Designed for Assigning confidence limits to X-ray detections, when ignoring the background (use mystats.kraftcl() if background is significant and n<10).  Can return 68%, 90%, or 95% confidence intervals.
   Use Analytic expression in Gehrels, Eq 9 for Upper-limit, and Eq 14 for lower limits (also see Table 3).  
   **NOTE: Gehrels gives one-sided upper and lower- confidence intervals, and this routine will return 2-sided limits.  So, asking for cl=68% will use Gehrels values for cl=0.84.  Asking for cl=90% here will use Gehrels values for cl=0.95.  Asking for cl=95% will use Gehrels values for cl=0.975.  Aksing for cl=99% will use Gehrels values for cl=0.995.**
   IN:'
     n= total number of counts (integer not required)
     cl=*2-SIDED* confidence interval.  0=68%, 1=90%, 2=95%, 3=99%
   OUT:'
     err=2-d array.  err[0]=minimum counts.  err[1]=maximum counts in confidence interval.
     *note* error on counts is the $(n-B)^{+(err[1]-(n-b))}_{-((n-b)-err[0])}
"""

    # Get parameters from Table 3 of Gehrels86 (needed for low-limit) for cl = 0,1,2,3
    # Set up a dictionary to return a tuple (s,beta,gamma)
    # S=number of sigma
    # beta, gamma = params from Geherels
    cl_dict = {0:(1.00,  0.0,   1.0),     # (s, beta, gamma) for cl=0; 68%-->cl=0.84
               1:(1.645, 0.031, -2.5),  # cl=1; 90%-->cl=0.95
               2:(1.960, 0.058, -2.22), # cl=2; 95%-->cl=0.975
               3:(2.576, 0.141, -2.00)} # cl=3, 99%-->cl=0.995
    s, beta, gamma = cl_dict[cl]

    # get upper and lower-limits
    tmperr_upper = (n+1)*(1. - (1./(9.*(n+1))) + (s/(3.*np.sqrt(n+1.))))**3.   # Equation (9)
    tmperr_lower = n*(1. - (1./(9.*n)) - (s/(3.*np.sqrt(n))) + (beta*n**gamma))**3. # Equation 14
    
    return tmperr_lower, tmperr_upper



def kraftcl(n, bg, cl):
    """ err_lo, err_high = mystats.kraftcl(n,b,cl)

    Returns Confidence limits from Bayesian Formulism of Kraft et al. 1991, ApJ, 374, 344.  Designed for Assigning confidence limits to X-ray detections.  If n total observed photons, and b background photons, returns minimum and maximum counts expected for 90%, 95%, or 99% confidence interval (cl).
    Reads in CL tables located in the /kraft91 subfolder
    IN:
        n  = total number of observed counts (including background).  Must be an integer <=10'
        b  = expected number of background counts.  Can be non-integer, but must be <=10'
        cl = confidence interval.  0=90%, 1=95%, 2=99%'
  OUT:'
    err_lo, err_hi = minimium, maximum counts in confidence interval.'
    *note*: total number of net counts is then (n-b)^  +(err_hi - (n-b))_ -(n-b-err_lo)'

    **Transcribed from an idl program I wrote (kraftcl.pro)
  """

    if cl not in [1,2,3]:
        print("Error: CL must = 1,2,3 for 90, 95, 99% confidence intervals, respectively")
        sys.exit(1)
    n = int(n)
    if n>10:
        print("Error: n must be integer <=10")
        sys.exit(1)
    if bg>10:
        print("Error: bg must be <=10")
        sys.exit(1)

    current_dir = os.getcwd()
    clrt = current_dir + '/kraft91/' #should work as long as you don't move the kraft91 folder
    # Define Confidence table to use
    infile_key = {1:clrt + 'kraft91_tab1_90cl.txt',
                  2:clrt + 'kraft91_tab2_95cl.txt',
                  3:clrt + 'kraft91_tab3_99cl.txt'}
    infile = infile_key[cl]

    #print("**Using n={0:d}, bg={1:.1f}, and confidence table {2:s}**".format(n,bg,infile.split('/')[-1]))

    # Read in cl table, and reformat to get lower and upper limits (each BKG level has 2 lines, top for lower limits, bottom for upper limits).
    cldata = ascii.read(infile, data_start=2)
    nlines = len(cldata) /2       # number of pairs of LL/UL confdence limits  
    wlower = np.arange(nlines)*2  # indices of lower CLs per BG
    wupper = wlower + 1           # indices of upper CLs per BG
    wlower = wlower.astype(int) #have to turn into integers for python3 to deal with it
    wupper = wupper.astype(int)

    # Extract the correct column for number of counts N=0-->10 (col2-->col12)   
    usecol = {0:'col2',
              1:'col3',
              2:'col4',
              3:'col5',
              4:'col6',
              5:'col7',
              6:'col8',
              7:'col9',
              8:'col10',
              9:'col11',
              10:'col12'}
    #Extract columns for BG counts and # of lower and upper counts
    #wlower = wlower.astype(int)
    #wupper = wupper.astype(int)
    nbg   = cldata['col1'].data[wlower]
    nmin  = cldata[usecol[n]].data[wlower]
    nmax  = cldata[usecol[n]].data[wupper]

    err_lo = np.interp(bg, nbg, nmin)
    err_hi = np.interp(bg, nbg, nmax)    

    return err_lo, err_hi

