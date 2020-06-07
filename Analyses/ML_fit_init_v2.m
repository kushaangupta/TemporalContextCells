%% initialize
%Plot Version
%Is it set to submit
%Is Plot off
%What is N at

%Adds folders where stuff will be saved
addpath(genpath('/project/ecog-eeg/stim_specific_time_cells/Buffalo/export_fig'));
addpath(genpath('/projectnb/ecog-eeg/stim_specific_time_cells/Buffalo'));
%Do I want to generate a plot while running
params.plot_on = 1;
%Sets where the results will be saved
params.folder_results = 'results_ML';
params.folder_mat = '/projectnb/ecog-eeg/stim_specific_time_cells/Buffalo';
%Significance level (Bonferoni correction)
params.max_p = .05/349;
%Use this for version control 
params.plot_version = 50;
%This is the data that I'm going to fit
params.data = 'cleaned_data_v1.mat'; 

%This is how many times I want to fail to improve
%the previous best fit before ending the code.
params.N=20;
%params.submit must = 1 if submited as job
params.submit =0;
%Checks if job was submitted
if params.submit
    %Makes sure no plot is generated
    params.plot_on = 0;
else
    if params.plot_on
        %Gives filename of smoothed firing rate data
        params.faverage = 'f_average_buffalo.mat';
    end
end
