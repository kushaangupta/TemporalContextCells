%% initialize
addpath(genpath('ecog-eeg/stim_specific_time_cells/Buffalo/export_fig'));
addpath(genpath('ecog-eeg/stim_specific_time_cells/Buffalo'));
params.plot_on = 0; % set to zero when submitted as a job
params.folder_results = 'results_ML';
%params.folder_mat = '/projectnb/ecog-eeg/stim_specific_time_cells/Cromer/mat_ML';
%params.folder_mat = '/projectnb/ecog-eeg/stim_specific_time_cells/BrincatMiller/mat_ML';
params.folder_mat = '../ML_Results_TCC';
params.max_p = 0.05/70;
%v19 complete results all trials, including correct and incorrect 
params.plot_version = 201;

%params.data = 'dataConvHpc_correct_LaterHalf_ss.mat'; % the data structure with the following fields: 
params.data = 'Data/simulated_exponential_decay_cells_400Hz.mat'; % the data structure with the following fields: 

%This is how many times I want to fit the data.
params.N=5; %number of trials for PSO
params.submit =1;
if params.submit % submit as a job
    params.plot_on = 0;
else
    if params.plot_on

        params.faverage = 'Data/f_average_single_cell.mat';
        %params.faverage = 'f_average_ss_Hpc_CorrCond_LaterHalf.mat'; 
        % use preprocess_data_cromer_v1.m to get f_average_single_cell.mat
    end
end
