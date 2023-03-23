ANALYSES

cell_counts_cutoff: Counts the number of cells with various relaxation times
    Used to generate minor subplot of Figure 3

glm_model_stim_specific_v3.m: Contains the equations used to model the temporal
	response profile. Called by the ML_fit_params files.

loop_process_ML_results.m: Takes in all of the outputs from the maximum
	likelihood stuff, and generates most (Non-LDA/Time Decoding) of 
	the numbers used in the paper. 

ML_fit_init_v2.m: Initialization script for ML_fit_params_v2.m.
	handles things like version control. Used to determine which cells
	were Temporal Context Cells before adding in Even/Odd requirement.
	
ML_fit_init_v2_EO: Similar to ML_fit_init_v2.m, but set up to determine
	if cells are significantly fit for both even and odd trials.

ML_fit_init_v4.m: Similar to ML_fit_init_v2.m, but used to fit the ideal
	time cells.
	
ML_fit_init_v5.m: Similar to ML_fit_init_v2.m, but used to fit the ideal
	temporal context cells.
	
ML_fit_params_v2.m: Uses maximum likelihood approach to determine
	if a cell has a temporal response profile. Used before adding in Even/Odd.
	Runs in Matlab/2016a.
	
ML_fit_params_v2_EO: Similar to ML_fit_params_v2.m, but set up to determine
	if cells are significantly fit for both even and odd trials. Runs in Matlab/2016a

ML_fit_params_v4.m: Similar to ML_fit_params_v2.m, but used to fit the ideal
	time cells. Runs in Matlab/2016a
	
ML_fit_params_v5.m: Similar to ML_fit_params_v2.m, but used to fit the ideal
	temporal context cells. Runs in Matlab/2016a

process_ideal.m: Used to process the results of ML_fit_params_v4.m and
	ML_fit_params_v5.m for the ideal cells.
	
process_ML_results2.m: Generates rasters and PSTH's.

LDATemporalDecoder_not_looped_FINAL.m: Runs a LDA analysis.  Current options
    are configured to run on simulated data 
    (simulated_exponential_decay_cells_20Hz_max20Hz.mat and 
    ),
    but can also be set to run on cleaned_data_v1.mat
    Used to generate Figure S4.

LDATemporalDecoder_leave_bins_out_FINAL.m: Runs a LDA analysis multiple times,
    each time dropping more time bins.  Used to generate Figure S5.
    Runs on cleaned_data_v1.mat

LDATemporalDecoder_leave_cells_out_FINAL.m: Runs a LDA analysis multiple times,
    each time dropping cells with longer relaxation times.
    Used to generate Figure 3
    Runs on cleaned_data_v1.mat and uses finalparameters.csv to identify the
    relaxation times
    
 vec2struct.m: Converts a vector into a structure


DATA 

(Some of the output files from the maximum likelihood fits were too large to upload to git. 
	Those files are marked with a "*" in this document and can be downloaded at the following 
	link: https://drive.google.com/drive/folders/1WJqmPdezQp9el54Kj7qIyieLDELitjnb?usp=sharing)


cleaned_data_v1.mat: Reformatted version of "theData180329.mat" that is in
	the correct format to be used by the maximum likelihood code. Only
	incorporates the fields still needed from "theData180329.mat."

f_average_buffalo.mat: Trial averaged firing rates for cells for all trials,
	novel only, and repeat only.

finalparameters.csv: Contains the parameter values for the 109 temporal context cells.

lag_cos_df.csv: Used to calculate the error bars in the lag analysis.

lag_cos_means_CIs.csv: Holds the error bar info for the lag analysis.

MaxFiringRate.mat: Contains the max average firing rate for each cell. Used for
	normalization in "loop_process_ML_results.m"

ML_Results*: Contains the results of running ML_fit_params_v2 on each cell.

ML_Results_EO*: Contains the results of running ML_fit_params_v2_EO on only those
	with a significant temporal response (before even/odd requirement).

ML_Results_TC*: Contains the results of running ML_fit_params_v4 on the 
	ideal time cells.
	
ML_Results_TCC*: Contains the results of running ML_fit_params_v5 on the
	ideal temporal context cells.

ParameterCorrelations.jasp: Contains the results of the Bayesian Kendall's 
	Tau analysis on cell parameters.

simulated_exponential_decay_cells_20Hz_max20Hz: Contains simulated ideal temporal
    context cells data.  Used by LDATemporalDecoder_not_looped_FINAL.m to
    generate Figure S4.

simulated_exponential_decay_cells_400Hz.mat: Contains spike times for
	idealized temporal context cells. Used by ML_fit_params_v5.m.

simulated_ideal_time_cells_samepeaks_00Hz_max400Hz.mat: Contains spike times for
	idealized time cells. Used by ML_fit_params_v4.m.

simulated_ideal_time_cells_1Hz_max40Hz: Contain simulated ideal time cell
    data.  Used by LDATemporalDecoder_not_looped_FINAL.m to generate Figure
    S4.

theData180329.mat: Original form of the EC data. Contains data from both monkeys and session information
