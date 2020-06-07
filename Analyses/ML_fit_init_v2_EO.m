%% initialize
%Plot Version
%Is it set to submit
%Is Plot off
%What is N at

%Adds some folders that the code needs to access.
addpath(genpath('/project/ecog-eeg/stim_specific_time_cells/Buffalo/export_fig'));
addpath(genpath('/projectnb/ecog-eeg/stim_specific_time_cells/Buffalo'));
%Do I want to generate a plot while running
params.plot_on = 1;
%Sets where the plots are saved
params.folder_results = 'results_ML';
params.folder_mat = '/projectnb/ecog-eeg/stim_specific_time_cells/Buffalo';
%Significance level (Bonferoni Corrected)
params.max_p = .05/349;
%Version control
params.plot_version = 220;
%The data that we want to analyze
params.data = 'cleaned_data_v1.mat'; 

%This is how many times I want to fail to improve
%the previous best fit before ending the code.
params.N=20;
%params.submit must = 1 if submited as job
params.submit =1;
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
