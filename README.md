# actin-tracking
Extraction of accurate cytoskeletal actin velocity distributions from noisy measurements

### Included in this repository:
1. README
2. MainScript.m 	- the main script to fit the data in sample_data to the analytical jump model ( https://doi.org/10.1101/2020.08.13.247304 )
3. sample_data.zip 	- a folder of sample data, containing 3 cells, can be downloaded from the release page (https://github.com/caylamiller/actintracking/releases)
4. matlab_functions 	- a folder of supporting scripts needed to run MainScript.m
5. MainScript_Weibull.m	- an alternative to MainScript.m which fits the data to a Weibull distribution

### System requirements: 
Matlab with the following toolboxes:
1. optimization
2. statistics and machine learning
This code was run successfully on both Matlab 2015a and 2020b on Windows 10.  

### To run:
1. Download all files.
2. Add matlab_functions directory to the matlab path.
3. Unzip sample_data.zip. Move sample_data folder to the directory in which MainScript.m is stored.
4. Run MainScript.m from the directory in which MainScript is stored. (Otherwise input and output directories in the script, lines 4 - 17 would need to be modified appropriately.)

For the first run, you may wish to modify certain run parameters to run more quickly (at the cost of loss of information):
1. imgnums (line 20) = the image numbers (as names in sample_data) that will be analyzed. You may wish to run only one image (e.g. 18) during the first run to ensure all toolboxes are properly installed and all scripts run. 
2. datamax (line 42)  = the maximum number of points to fit at each timescale. It is currently set to 5000; running at a low number (e.g. 100) may be useful to test that code runs completely, though the output will not be accurate.
3. fitts (line 44) = the timescales on which to fit the data; set conservatively now to [10 20 30 40] (seconds). Running only 1 timescale (e.g. 20) will similarly speed analysis for quick testing. 

As written, this script ran in about 20 minutes on a desktop computer.

### Outputs:
The results in the paper rely on a much larger dataset (n = 32 cells across multiple experimental replicates). While the results from this small number of cells may roughly match those in the paper, it is expected that results will have greater errors from the smaller sample size. 

Output files are stored in matlab_output if run as provided. 

##### .mat files created:
1. [date]_[imgnum]_tracks_n2v.mat - one for each image analyzed, contains structure named "results," with each line corresponding to 1 tracked speckle from that image, as well as variables and parameters used for calculations
2. [date]_all_tracks_tracks_n2v.mat - combines all of the tracking results from the above files into one, contains structure named "allresults" which contains all tracks from all images analyzed
3. [date]_all_displacements.mat - saved workspace variables after calculating displacements over varying timescales, includes structure "stepresults" which contains these displacements
4. [date]_analytical_jump_model_fit_params.mat - includes variables [pop]paramsMLE (containing best fit parameters) and [pop]cisMLE (containing the 95% confidence intervals on those parameters)

##### .fig files created:
1. [date]_allpops_measured_displacements.fig - showing the measured displacement distributions, pooled across all populations and all cells
2. [date]_[pop]_fit_displacements.fig - one for each subcellular population (4 total), showing the measured displacement distribution, pooled across all cells, and the best fit to that population by the analytical jump model
3. [date]_[pop]_velocity_distributions.fig - one for each subcellular population (4 total), showing the best fit velocity distributions from the analytical jump model
4. [date]_analytical_jump_model_fit_params.fig - showing the best fit parameters for the four populations, as a function of timescale

### Expanding to fit other models:
This code fits the provided data to the analytical 1-D jump model (without reversals). Other models in the text were tested by modifying the "Fig1DAnalyticalJumpModel," "AnalyticalJumpModelWithNoise," and "AnalyticalJumpModel" codes. 

As an example, "FitWbl" and "WblWithNoise" are included. These are called in a script called "MainScript_Weibull.m". The primary modifications to fit this model are:
1. Creating the FitWbl function. This function is directly derived from the Fit1DAnalyticalJumpModel function, modifying the following:
	a. The functions "myfunLSQ" and "myfunMLE" are updated to use WblWithNoise, and fit parameters redefined accordingly
	b. The "bestfit" vector is redefined to use WblWithNoise
2. Creating the "WblWithNoise" function to generate a distribution of measured displacements as defined by Bayes theorum, from an underlying Weibull distribution of true displacements, with localiation error defined by the parameters in noiseparams. This function is directly derived from the "AnalyticalJumpModelWithNoise" function, modifying the following:
	a. pSj is redefined to call "wblpdf" instead of "AnalyticalJumpModel"
	b. fd is simplified such that it does not include the discrete probability of moving zero jumps
3. Creating "MainScript_Weibull" to call the above functions and properly save and display outputs. This script is directly derived from the "MainScript", modifying the following:
	a. paramseed is redefined to seed appropriate Weibull shape and scale parameters
	b. vrange begins at 0 rather than 0.001, because the Weibull is not discontinuous at 0 as the analytical jump model is
	c. fit calls FitWbl (line 120)
	d. figures and the fit parameter .mat file are renamed (lines 124, 150, 183, 185)
	e. velocity plots call the wblpdf function (line 136) and do not include a discontinuity at 0
	f. fit parameters are labeled appropriately for the function (lines 169, 178)

These functions can be modified similarly to fit any desired underlying distribution. 
