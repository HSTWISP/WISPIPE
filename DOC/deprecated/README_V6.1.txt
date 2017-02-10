# WISPIPE README
 Based on Ivano Baronchelli's howtolaunch.txt and
 Sophia Dai's wispipe-dr-guideline.txt
 Marc Rafelski October 2016

New update: Now use the bash shell for everything. (Both Euereka and astroconda)

 -------------------------------
 SYSTEM REQUIREMENTS / INSTALL TIPS
 -------------------------------
IDL version: Tested with: 8.5.1 (Ivano) 8.2 (Sophia), 8.4 (Marc)
Ureka version: Tested with 1.5.2 (Ivano & Marc) 1.4.1.2 (Sophia)
Astroconda version: 4.1.11 (Ivano) 4.1.1 (Sophia) (latest Oct 25, Marc)

 The pipeline also depends on the latest version of the IDLastro
 library and pyregion (see BEFORE STARTING / SETUP below)

How to install Anaconda:

Download Anaconda with python 2.7 at https://www.continuum.io/downloads
 
Download command-line installer for Python 2.7 (64 bit), and then install it:
>> bash Anaconda2-4.1.1-MacOSX-x86_64.sh
 
install astroconda in python 
 
>> conda config --add channels http://ssb.stsci.edu/astroconda
 
>> conda create -n astroconda python=2.7 iraf pyraf stsci

Then whenever you want it, type: source activate astroconda (in BASH only!!)

How to install Eureka:

Download this link:
OS X < 10.12 (or linux)
http://ssb.stsci.edu/ureka/1.5.2/install_ureka_1.5.2

OS X >= 10.12
http://ssb.stsci.edu/ureka/1.5.2/install_ureka_1.5.2_sierra

Run the installer, e.g. sh install_ureka_1.5.2

To activate Ureka, type: ur_setup

**** Remember: Do not have both Ureka and astroconda active at the same
time. You can do ur_forget or source deactivate astroconda to remove
one or the other.

 -------------------------------
 BEFORE STARTING / SETUP
 -------------------------------

 1) For pyraf, make sure you have an iraf directory in ~/. In that
 directory, make sure you have a login.cl file. An example login.cl
 file is located in the DOC directory. Make sure to have it correctly
 point to your user information!


 2) Before beeing able to start, create a folder structure where you
 will reduce your data, e.g:
 $PATH/aXe
 $PATH/data
 The first folder (aXe) correspond to the final destination of the reduced
 data. The second one (data) correspond to the field where the
 original raw data downloaded from MAST are saved. Note that once the
 pipeline has started, the original files in the "data" folder will
 be modified.  So you should save a copy of the data you downloaded
 in case you want to run the pipeline again. One option is to create
 a second folder called:
 $PATH/download
 In the "data" folder, the raw data should be saved in folders named
 like (example for field 356):
 $PATH/data/Par356

 
 3) In your .bashrc and .cshrc files, you have to have specific paths set. 
 .cshrc:
 
 setenv WISPIPE /Users/mrafelski/WISPIPE/
 setenv WISPDATA /Users/mrafelski/data/wisp/
 setenv iref ~/data/iref/
 setenv crrefer ~/data/synphot/
 setenv mtab ~/data/synphot/mtab/
 setenv crotacomp ~/data/synphot/comp/ota/
 setenv crwfc3comp ~/data/synphot/comp/wfc3/
 
 .bashrc:
 
 export WISPIPE=/Users/mrafelski/WISPIPE/
 export WISPDATA=/Users/mrafelski/data/wisp/
 export iref=~/data/iref/
 export mtab=~/data/synphot/mtab/
 export crrefer=~/data/synphot/
 export crotacomp=~/data/synphot/comp/ota/
 export crwfc3comp=~/data/synphot/comp/wfc3/

Note: We have now switched completely to bash, but I have not tested
if it is not necessary to set the paths of the .cshrc or not. 

 4) Make sure you have synphot installed
 http://www.stsci.edu/institute/software_hardware/stsdas/synphot
 More info:
 http://www.stsci.edu/hst/observatory/crds/throughput.html
 Download files here (or get them from someone else):
 ftp://archive.stsci.edu/pub/hst/pysynphot/
 Check for latest Pysynphot files:
 http://www.stsci.edu/hst/observatory/crds/SIfileInfo/pysynphottables/current_tmc_html
 Make sure have the latest tab, mtab, files. Can download latest
 here: ftp.stsci.edu/cdbs
 (check mtab directory, comp/wfc3 directory)

 5) Make sure you have the files for the $iref directory. Best to get
 these from someone, but can also be downloaded here:
 ftp.stsci.edu/cdbs/iref
 Note: You want to get the latest darks from Marc Rafelski

 6) tweakreg needs to have interactive = False set as
 your default before running the pipeline. You can do this by
 entering pyraf, importing astrodrizzle, and then eparing tweakreg,
 and setting interative to False.
 MR NOTE: Do we want to change this in the pipeline to specify this so
 this doesn't depend on a system default???

 7) While we use Ureka, you need to install pyregion seperatly. e.g.:
 While Ureka is loaded, type:
 pip install pyregion

 8) Make sure you have the IDLastro library installed and up to date:
 https://github.com/wlandsman/IDLAstro
 Needs things like modfits, etc.

 -------------------------------
 DOWNLOAD DATA
 -------------------------------
 Download the raw data from HST MAST (http://archive.stsci.edu/hst/search.php)
1) Choose the proposal ID and Obset ID (orbits), refer to the google
 drive excel sheet for information on individual field:
 (https://docs.google.com/spreadsheets/d/1QWagshUlPFoGPTFxSIwl6dVKE-MWEck2jC2DF398iJ4/edit?usp=sharing)

2) Download the 'calibrated' data for IR exposures, and 'uncalibrated'+'reference files' for UVIS exposures.  
 (If space allows, I would just download everything to simplify the procedure.)

 -------------------------------
 PREPARATION
 -------------------------------
 1) Before running everything, check the kind of data you want to
 process. At this moment we are only able to process data with both
 the IR filters and grisms observed. UVIS data may or may not be present.

 2) Check input exposures using wispipe_initialcheck.sh
  (or use findf.pro in IDL directory). 
 The single exposures of all the filters will be shown in different
 ds9 windows for different filters/grism.
 Single exposures can be removed creating an
 appropriate folder and moving the contaminated (or showing other
 problems) ones in it. e.g.

cd $WISPIPE
source wispipe_initialcheck.sh Par302
cd $WISPDATA/Par302
mkdir badframes
mv icdxd3l2q* badframes

Notes on what programs to use when (but read below):

For fields with both IR filters and no UVIS, use wispipe.sh
For field with 2 UVIS filters, use wispipe_uvis.sh
For field with 1 UVIS filter, use wispipe_uvis_single.sh (may be obsolete as the only difference is in tweakreg_uvis_single.py and make_uvis_helpfile.py) 
For field with only G141 filter (no G102), use wispipe_g141.sh
For field with only G141 filter (no G102) and only 1 exposure (no tweakreg needed), use wispipe_g141_single.sh
For field with F140 + G102 + G141 filters, use wispipe_F140.sh


 -------------------------------
 CASE 1: UVIS DATA ARE NOT PRESENT 
 -------------------------------

 All in  a bash shell, Eureka environment
 PROGRAM TO USE: wispipe.sh

 If uvis data are not present, the pipeline can be run as follows:
 1) cd to the path where the WISPIPE folder is installed
    > cd $WISPIPE
 2) run the current version of the pipeline
    > source wispipe.sh Par > LOG/log_Par_6.1.log

 Example:
    > cd $WISPIPE
    > source wispipe.sh Par302 > LOG/log_Par302_6.1.log

 The data reduction can be launched without interruptions among one
 field and the next, using the program "multiple_par.sh", saved in
 this same folder. Read in this file for more informations.
 NOTE that the multiple par mode can't be used to preprocess uvis
 data. In this case the preprocess must be run using a bash terminal
 and an astroconda environment.

Note: Can use multiple_par.sh to run multiple at once; edit as appropriate from DOC

 -------------------------------
 CASE 2: UVIS DATA ARE PRESENT 
 -------------------------------

 Preprocessing in bash, astroconda  environment,
 Reduction in  a bash shell, Eureka environment
 PROGRAM TO USE:, uvis_preprocess.pro or multiple_uvis_preprocess.pro,
 wispipe_uvis.sh

 FIRST STEP: preprocess uvis files (bash, astroconda)
 To preprocess the uvis data, you can run it either individually on a
 specific field. Multiple fields can be pre-processed in a row if they
 have the same characteristics (ex: all observations before 2012,
 and both uvis filters present).

 Run IDL program to run mutiple fields, e.g.:
 multiple_uvis_preprocess, ['Par302', 'Par303', 'Par304']

 or run a single field:
 uvis_preprocess, 'Par302'

IMPORTANT NOTE: Data before around October 2012 need to be run with
the flag, /nopostflash in the idl program. If you do not do this,
things will crash. It depends if there are post-flash darks or not. If
you are processing data round this time change, check with Marc if
unsure.

Note: You need to have all the calibration files. You can run the code
with flag of /darksonly to have it tell you what darks you need. Talk
to Marc about this if you have questions.

Note: You can also run just calwf3 with this (e.g. if you ran CTE
corrections, but didn't end up having the darks to calibrate them. In
that case, use the /calwf3only flag. 

 SECOND STEP: reduction (tcsh, Eureka)
 The reduction can be run as in the case 1 (no uvis data present),
 but using the wispipe_uvis.sh program instad. Again, the
 reduction of multiple fields in a row can be obtained using
 "multiple_par.sh", but ONLY AFTER the uvis data are already
 preprocessed in a bash shell and astroconda environment.

    > source wispipe_uvis.sh Par302 > LOG/log_Par302_6.1.log

Note: May use  multiple_uvis_preprocess,multiple_par.sh
 (if desired, in DOC)

 -------------------------------
 CHECK RESULTS 
 -------------------------------

Check the log for any errors. 

check the output:

this program displays drizzled IR direct+grism, and UVIS data + pdf file for 1D spectra
IDL > qacheck,'258'

if w/UVIS observation, 
IDL > qacheck,'258',/uvis


==========  Tips ==========
to reduce several fields in a row, generate a shell script dr.sh with content:
./wispipe.sh ParXXX > & log-206-160310-1810.txt 
./wispipe.sh ParXXX > & log-216-160310-1900.txt 
...
...
...

then in shell source dr.sh:
% ./dr.sh

