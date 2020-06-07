%% initialize

%Adds folders that the code needs access to.
addpath(genpath('/project/ecog-eeg/stim_specific_time_cells/Buffalo/export_fig'));
addpath(genpath('/projectnb/ecog-eeg/stim_specific_time_cells/Buffalo'));
%Do I want to also plot things?
params.plot_on = 0;
%Sets where the results will be saved
params.folder_results = 'results_ML';
params.folder_mat = '/projectnb/ecog-eeg/stim_specific_time_cells/Buffalo';
%p-value for cells
params.max_p = 0.05/55;
%Version control 
params.plot_version = 102;

%This is the data to be fit
params.data = 'simulated_ideal_time_cells_samepeaks_00Hz_max400Hz.mat'; % the data structure with the following fields: 

%This is how many times I want to fail to improve
%the previous best fit before ending the code.
params.N=3;
%Set to 1 if running on cluster
params.submit =1;
if params.submit % submit as a job
    params.plot_on = 0;
else
    if params.plot_on
        
        %params.faverage = 'f_average_single_cell.mat';
        %params.faverage = 'f_average_ss_Hpc_CorrCond_LaterHalf.mat'; 
        % use preprocess_data_cromer_v1.m to get f_average_single_cell.mat
    end
end
