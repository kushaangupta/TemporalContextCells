%% This does preprocessing to complete rest of paper.
warning('off','all')
%These will hold the parameters for the convolution model
%Holds all of the unit sigmas
EMGsigs = [];
%Holds all of the unit mus
EMGmus = [];
%Holds all of the unit taus
EMGtaus = [];
%Holds all of the unit indices
EMGindex = [];
%Holds all of the unit peak values
EMGpeak = [];
%Holds all of the unit baseline firing values
EMGcons = [];
%Counts how many units were significantly fit by the
%convolution model
EMGcount = 1;

%These will hold the parameters for the Gaussian model
%Holds all of the unit sigmas
Gsigs = [];
%Holds all of the unit mus
Gmus = [];
%Holds all of the unit taus
%This will just be NaN's
Gtaus = [];
%Holds all of the unit indices
Gindex = [];
%Holds all of the unit peak values
Gpeak = [];
%Holds all of the unit baseline firing values
Gcons = [];
%Counts how many units were significantly fit by the
%Gaussian
Gcount = 1;



%Loads up the initialization script that generated the fits.
ML_fit_init_v2

%Finds all the parameters for the EC units
%Units 238-245 are hippocampal
for i = [1:237 246:357]
    try
        %Prints what unit is it on
        cell_no=i
        
        %Loads the fit for the given unit i
        load(sprintf('%s/glm_cell_%i_v_%i.mat',params.folder_mat,cell_no,params.plot_version))
        params.folder_mat = '../ML_Results';
        
        %Checks if a version of the convolution model fit the spikes better
        %than the constant
        if  (max(T) > 0.003) && ((pValue < .05/349) || (pValue_I < .05/349))
            %Checks if negative or positive EMG was the better model
            if (LL_I < LL_T)
                %Positive EMG is winner
                if xT(5) > 0.001
                    %Stores unit sigma
                    EMGsigs(EMGcount) = xT(2);
                    %Stores unit mu
                    %Need to move t=0 back to actual 0.
                    EMGmus(EMGcount) = xT(3) + 2.25;
                    %Stores unit tau
                    %Report as 1/tau in model
                    EMGtaus(EMGcount) = 1/xT(4);
                    %Stores baseline firing rate
                    EMGcons(EMGcount) = xT(1);
                    %Stores unit peak height
                    EMGpeak(EMGcount) = xT(5);
                    %Stores the unit index
                    EMGindex(EMGcount) = i;
                    %Increases iterator
                    EMGcount = EMGcount + 1;
                end
            else
                %Negative EMG is the winner
                if xI(5) < -0.001
                    %Stores unit sigma
                    EMGsigs(EMGcount) = xI(2);
                    %Stores unit mu
                    %Need to move t=0 back to actual 0.
                    EMGmus(EMGcount) = xI(3) + 2.25;
                    %Stores unit tau
                    %Report as 1/tau in model
                    EMGtaus(EMGcount) = 1/xI(4);
                    %Stores baseline firing rate
                    EMGcons(EMGcount) = xI(1);
                    %Stores the unit peak height
                    EMGpeak(EMGcount) = xI(5);
                    %Stores the unit index
                    EMGindex(EMGcount) = i;
                    %Increases iterator
                    EMGcount = EMGcount + 1;
                end
            end
        end
        
        %Checks if a version of the Gaussian model fit the spikes better
        %than the constant
        
        if max(G) > 0.003 && (pValue_G < .05/349 || pValue_IG < .05/349)
            %Checks if negative or positive Gauss was the better model
            if (LL_IG < LL_G)
                %Positive Gauss is winner
                if xG(5) > 0.001
                    %Stores unit sigma
                    Gsigs(Gcount) = xG(2);
                    %Stores unit mu
                    %Need to move t=0 back to actual 0.
                    Gmus(Gcount) = xG(3) + 2.25;
                    %Stores unit tau
                    %Unit tau does not exist for Gaussian
                    Gtaus(Gcount) = NaN;
                    %Stores baseline firing rate
                    Gcons(Gcount) = xG(1);
                    %Stores unit peak height
                    Gpeak(Gcount) = xG(5);
                    %Stores the unit index
                    Gindex(Gcount) = i;
                    %Increases iterator
                    Gcount = Gcount + 1;
                end
            else
                %Negative Gaussian is the winner
                if xIG(5) < -0.001
                    %Stores unit sigma
                    Gsigs(Gcount) = xIG(2);
                    %Stores unit mu
                    %Need to move t=0 back to actual 0.
                    Gmus(Gcount) = xIG(3) + 2.25;
                    %Stores unit tau
                    %Unit tau does not exist for Gaussian
                    Gtaus(Gcount) = 1/xIG(4);
                    %Stores baseline firing rate
                    Gcons(Gcount) = xIG(1);
                    %Stores the unit peak height
                    Gpeak(Gcount) = xIG(5);
                    %Stores the unit index
                    Gindex(Gcount) = i;
                    %Increases iterator
                    Gcount = Gcount + 1;
                end
            end
        end
    catch ME
        fprintf('File not found for unit %d\n', cell_no);
        continue;
    end
end


%Will hold if cells were better with the Gaussian or EMG
nopass = [];

G_betterthan_EMG = [];
pvalues_G_EMG = [];

iter=1;
iter1 = 1;
%Check for any cells where the EMG was not better than the Gaussian
for i = EMGindex
    
    %Check if the given unit was also significantly fit
    if(ismember(i, Gindex))
        
        %Get location of the cells to be able to compare
        EMG_ind = find(i == EMGindex);
        G_ind = find(i == Gindex);
        
        %Loads the fit for the given unit i
        load(sprintf('%s/glm_cell_%i_v_%i.mat',params.folder_mat,i,params.plot_version))
        params.folder_mat = '../ML_Results';
        
        %Determine if EMG is positive or negative
        if(EMGpeak(EMG_ind) > 0)
            %EMGpeak is positive
            if (pValue_GE >= .05/349)
                %Gaussian is winner
                G_betterthan_EMG(iter) = i;
                pvalues_G_EMG(iter)=pValue_GE;
                iter = iter+1;
                
            end
        else
            %EMGpeak is negative
            if (pValue_IG >= .05/349)
                %Gaussian is winner
                G_betterthan_EMG(iter) = i;
                pvalues_G_EMG(iter)=pValue_IG;
                iter = iter+1;
            end
        end
    end
    
    %Checks if the fit for guassian is better fit
    if(ismember(i, G_betterthan_EMG))
        %Gaussian is better
        nopass(iter1) = 1;
    else
        %EMG is better
        nopass(iter1) = 0;
    end
    
    %Increase
    iter1 = iter1+1;
    
end

emgcelltau =[];
gauscelltau =[]';

iter = 1;
iter1 = 1;

for i = 1:length(EMGpeak)
    
    if(nopass(i) == 0)
        emgcelltau(iter) =EMGtaus(i);
        iter = iter+1;
    else
        gauscelltau(iter1) =EMGtaus(i);
        iter1=iter1+1;
    end
    
    
end


%Cells Better fit by EMG
G_Better_EMG = length(emgcelltau);





%Holds all the parameters for the cells
parameters = [EMGindex; EMGmus; EMGsigs; EMGtaus; EMGpeak].';
parameters = array2table(parameters, 'VariableNames', {'Unit', 'Mu', 'Sigma', 'Tau', 'Peak'});
%Used by JASP
writetable(parameters,'Data/finalparameters.csv');

%Loads up the initialization script that generated the fits.
ML_fit_init_v2_EO

%Loads up the final parameters, even though they should already be loaded
%Under its new name
finalparams = readtable('Data/finalparameters.csv','HeaderLines',1);
finalparams = finalparams{:,:};
%Track which cells pass
pass = zeros(128,5);
%loop through the cells
for i = 1:128
    %Set the cell number
    cell_no = finalparams(i,1);
    %Stores the cell number
    pass(i,5) = cell_no;
    try
        %Loads the fit for the given unit i
        load(sprintf('%s/glm_cell_%i_v_%i.mat',params.folder_mat,cell_no,params.plot_version))
        params.folder_mat = '../ML_Results';
    catch
    end
    %goes through the steps in The Salz and Tiganj
    %papers for even odd fits.
    
    if pValue_O < .01
        pass(i,1) = 1;
    end
    
    if pValueE < .01
        pass(i,2) = 1;
    end
    
    Evens = reshape(TE,5500,[]);
    Evens = mean(Evens,2);
    
    Odds = reshape(TO,5500,[]);
    Odds = mean(Odds,2);
    
    
    [h p] = corr(Odds, Evens);
    
    %Store the rho
    pass(i,4) = h;
    %Checks
    if pass(i,1) == 1 & pass(i,2) == 1 & pass(i,4) >= 0.4
        pass(i,3) = 1;
    end
    
end

units = pass(:,5);

for i = 1:length(units)
    
    targetrowfinal = find(finalparams(:,1) == units(i));
    
    if pass(i,3) == 1
        
        finalparams(targetrowfinal,6) = pass(i,4);
    else
        
        finalparams(targetrowfinal,:) = [];
        
    end
    
end

parameters = array2table(finalparams, 'VariableNames', {'Unit', 'Mu', 'Sigma', 'Tau', 'Peak','rho'});
%Used by JASP
writetable(parameters,'Data/finalparameters.csv');

%% Paragraph on basic fit results
%This is the number of cells (out of 349) that were fit by the EMG model.
number_of_fit_cells = length(finalparams)
%Number of cells that increase firing;
increasing_cells = sum(finalparams(:,5) > 0)
%Number of cells that decrease firing rate
decreasing_cells = number_of_fit_cells - increasing_cells

%% Paragraph on parameter distributions
%Median value of mu
median_mu = median(finalparams(:,2))
%Gives the interquatile range of mu
low_quartile_mu = quantile(finalparams(:,2),.25)
high_quartile_mu = quantile(finalparams(:,2), .75)
%Gives last decile of mu
mu_90_percent = quantile(finalparams(:,2), .9)


EMGtaus_90 = finalparams(:,4);
%Median value of tau
median_tau = median(EMGtaus_90)
%Gives the interquatile range of mu
low_quartile_tau = quantile(EMGtaus_90,.25)
high_quartile_tau = quantile(EMGtaus_90, .75)
%Gives last decile of mu
tau_90_percent = quantile(EMGtaus_90, .9)

%Median value of sigma
median_sigma = median(finalparams(:,3))
%Gives the interquatile range of mu
low_quartile_sigma = quantile(finalparams(:,3),.25)
high_quartile_sigma = quantile(finalparams(:,3), .75)
%Gives last decile of mu
sigma_90_percent = quantile(finalparams(:,3), .9)


%% Paragraph on Unit Correlation to Stimulus

%Loads up dataset with stimulus information
load('Data/theData180329.mat')
%renames the data
data = olddata;
%Holds the unit correlations
all_correlations = [];
shuf_correlations = [];
%iterator for storing unit values
iter = 1;

%Runs through all the EC units
for cell = [1:237 246:357]
    
    
    %iterator for images
    iter1 = 1;
    
    %Holds spikes for first and second presentations
    novel_spikes = [];
    snovel_spikes= [];
    repeat_spikes = [];
    
    %Loops through all first presentation images
    for image = 1:length(data.novelTrlI{cell,1})
        
        %checks if there is a corresponding second presentaiton image
        if length(find(data.imgID{cell,1} == data.novelTrlI{cell,1}(image))) == 2
            
            %Stores trials where the image appears
            img_neuron_id = find(data.imgID{cell,1} == data.novelTrlI{cell,1}(image));
            %Stores number of spikes during first presentation
            novel_spikes(iter1) = nansum(data.spks{cell,1}(img_neuron_id(1), 501:5500));
            %Stores number of spikes during second presentation
            repeat_spikes(iter1) = nansum(data.spks{cell,1}(img_neuron_id(2), 501:5500));
            %add to iterator
            iter1 = iter1+1;
        end
        
    end
    
    forw = 1;
    
    if ~mod(length(novel_spikes), 2)
        maximum = length(novel_spikes);
    else
        maximum = length(novel_spikes)-1;
    end
    
    snovel_spikes = novel_spikes;
    
    for s = 1:maximum
        
        if forw
            snovel_spikes(s+1) = novel_spikes(s);
            forw = 0;
        else
            snovel_spikes(s-1) = novel_spikes(s);
            forw = 1;
        end
        
    end
    
    %Checks if there was a complete block
    if iter1 >= 30
        
        %Calculates the Correlation
        all_correlations(iter) = corr(novel_spikes.', repeat_spikes.', 'Type', 'Kendall');
        shuf_correlations(iter) =  corr(snovel_spikes.', repeat_spikes.', 'Type', 'Kendall');
        %increases the correlation
        iter = iter+1;
    end
    
end

%Perform t-test
[h_corr_all p_corr_all CI_corr_all t_corr_all] = ttest(all_correlations);
%Report results of t test
t_corr_all.tstat
t_corr_all.df
p_corr_all
corr_cohens_d_all = mean(all_correlations)/std(all_correlations)

%Perform t-test
[h_corr_shuf p_corr_shuf CI_corr_shuf t_corr_shuf] = ttest(shuf_correlations);
%Report results of t test
t_corr_shuf.tstat
t_corr_shuf.df
p_corr_shuf
corr_cohens_d_shuf = mean(shuf_correlations)/std(shuf_correlations)

%Perform t-test
[h_corr_alls p_corr_alls CI_corr_alls t_corr_alls] = ttest(all_correlations, shuf_correlations);
%Report results of t test
t_corr_alls.tstat
t_corr_alls.df
p_corr_alls
corr_cohens_d_alls = mean(all_correlations-shuf_correlations)/std(all_correlations-shuf_correlations)

%Perform Wilcoxon signed rank test
[p_WSRT_all h_WSRT_all stats_WSRT_all] = signrank(all_correlations);
%Report results
stats_WSRT_all.signedrank
p_WSRT_all

%Perform Wilcoxon signed rank test
[p_WSRT_shuf h_WSRT_shuf stats_WSRT_shuf] = signrank(shuf_correlations);
%Report results
stats_WSRT_shuf.signedrank
p_WSRT_shuf


%Holds the unit correlations for Temporal context cells
TCC_correlations = [];
TCC_shufcorrelations = [];
%iterator for storing unit values
iter = 1;

EMGindex = finalparams(:,1);
%Runs through Temporal Context Cells
for temp = 1:length(EMGindex)
    
    cell = EMGindex(temp);
    
    %iterator for images
    iter1 = 1;
    
    %Holds spikes for first and second presentations
    novel_spikes = [];
    snovel_spikes= [];
    repeat_spikes = [];
    
    %Loops through all first presentation images
    for image = 1:length(data.novelTrlI{cell,1})
        
        %checks if there is a corresponding second presentaiton image
        if length(find(data.imgID{cell,1} == data.novelTrlI{cell,1}(image))) == 2
            
            %Stores trials where the image appears
            img_neuron_id = find(data.imgID{cell,1} == data.novelTrlI{cell,1}(image));
            %Stores number of spikes during first presentation
            novel_spikes(iter1) = nansum(data.spks{cell,1}(img_neuron_id(1), 501:5500));
            %Stores number of spikes during second presentation
            repeat_spikes(iter1) = nansum(data.spks{cell,1}(img_neuron_id(2), 501:5500));
            %add to iterator
            iter1 = iter1+1;
        end
        
    end
    
    forw = 1;
    
    if ~mod(length(novel_spikes), 2)
        maximum = length(novel_spikes);
    else
        maximum = length(novel_spikes)-1;
    end
    
    snovel_spikes = novel_spikes;
    
    for s = 1:maximum
        
        if forw
            snovel_spikes(s+1) = novel_spikes(s);
            forw = 0;
        else
            snovel_spikes(s-1) = novel_spikes(s);
            forw = 1;
        end
        
    end
    
    %Checks if there was a complete block
    if iter1 >= 30
        
        %Calculates the Correlation
        TCC_correlations(iter) = corr(novel_spikes.', repeat_spikes.', 'Type', 'Kendall');
        TCC_shufcorrelations(iter) =  corr(snovel_spikes.', repeat_spikes.', 'Type', 'Kendall');
        %increases the correlation
        iter = iter+1;
    end
    
end

%Perform t-test
[h_corr_TCC p_corr_TCC CI_corr_TCC t_corr_TCC] = ttest(TCC_correlations);
%Report results of t test
t_corr_TCC.tstat
t_corr_TCC.df
p_corr_TCC
corr_cohens_d_TCC = mean(TCC_correlations)/std(TCC_correlations)

%Perform t-test
[h_corr_sTCC p_corr_sTCC CI_corr_sTCC t_corr_sTCC] = ttest(TCC_shufcorrelations);
%Report results of t test
t_corr_sTCC.tstat
t_corr_sTCC.df
p_corr_sTCC
corr_cohens_d_TCC = mean(TCC_shufcorrelations)/std(TCC_shufcorrelations)

%Perform t-test
[h_corr_TCCs p_corr_TCCs CI_corr_TCCs t_corr_TCCs] = ttest(TCC_correlations, TCC_shufcorrelations);
%Report results of t test
t_corr_TCCs.tstat
t_corr_TCCs.df
p_corr_TCCs
corr_cohens_d_TCCs = mean(TCC_correlations-TCC_shufcorrelations)/std(TCC_correlations-TCC_shufcorrelations)


%Perform Wilcoxon signed rank test
[p_WSRT_TCC h_WSRT_TCC stats_WSRT_TCC] = signrank(TCC_correlations);
%Report results
stats_WSRT_TCC.signedrank
p_WSRT_TCC

%Perform Wilcoxon signed rank test
[p_WSRT_sTCC h_WSRT_sTCC stats_WSRT_sTCC] = signrank(TCC_shufcorrelations);
%Report results
stats_WSRT_sTCC.signedrank
p_WSRT_sTCC


%Holds the unit correlations for non-Temporal context cells
NTCC_correlations = [];
NTCC_shufcorrelations = [];
NTCC_Num = [];
%iterator for storing unit values
iter = 1;


%Runs through non-Temporal Context Cells
for cell = [1:237 246:357]
    
    %Make sure cell is not Temporal Context
    if ~(ismember(cell, EMGindex))
        %iterator for images
        iter1 = 1;
        
        %Holds spikes for first and second presentations
        novel_spikes = [];
        snovel_spikes= [];
        repeat_spikes = [];
        
        %Loops through all first presentation images
        for image = 1:length(data.novelTrlI{cell,1})
            
            %checks if there is a corresponding second presentaiton image
            if length(find(data.imgID{cell,1} == data.novelTrlI{cell,1}(image))) == 2
                
                %Stores trials where the image appears
                img_neuron_id = find(data.imgID{cell,1} == data.novelTrlI{cell,1}(image));
                %Stores number of spikes during first presentation
                novel_spikes(iter1) = nansum(data.spks{cell,1}(img_neuron_id(1), 501:5500));
                %Stores number of spikes during second presentation
                repeat_spikes(iter1) = nansum(data.spks{cell,1}(img_neuron_id(2), 501:5500));
                %add to iterator
                iter1 = iter1+1;
            end
            
        end
        
        forw = 1;
        
        if ~mod(length(novel_spikes), 2)
            maximum = length(novel_spikes);
        else
            maximum = length(novel_spikes)-1;
        end
        
        snovel_spikes = novel_spikes;
        
        for s = 1:maximum
            
            if forw
                snovel_spikes(s+1) = novel_spikes(s);
                forw = 0;
            else
                snovel_spikes(s-1) = novel_spikes(s);
                forw = 1;
            end
            
        end
        
        %Checks if there was a complete block
        if iter1 >= 30
            
            %Calculates the Correlation
            NTCC_correlations(iter) = corr(novel_spikes.', repeat_spikes.', 'Type', 'Kendall');
            NTCC_shufcorrelations(iter) =  corr(snovel_spikes.', repeat_spikes.', 'Type', 'Kendall');
            %increases the correlation
            iter = iter+1;
        end
    end
end
    
    %Perform t-test
    [h_corr_NTCC p_corr_NTCC CI_corr_NTCC t_corr_NTCC] = ttest(NTCC_correlations);
    %Report results of t test
    t_corr_NTCC.tstat
    t_corr_NTCC.df
    p_corr_NTCC
    corr_cohens_d_NTCC = mean(NTCC_correlations)/std(NTCC_correlations)
    
    %Perform t-test
    [h_corr_sNTCC p_corr_sNTCC CI_corr_sNTCC t_corr_sNTCC] = ttest(NTCC_shufcorrelations);
    %Report results of t test
    t_corr_sNTCC.tstat
    t_corr_sNTCC.df
    p_corr_sNTCC
    corr_cohens_d_sNTCC = mean(NTCC_shufcorrelations)/std(NTCC_shufcorrelations)
    
    
     %Perform t-test
    [h_corr_NTCCs p_corr_NTCCs CI_corr_NTCCs t_corr_NTCCs] = ttest(NTCC_correlations,NTCC_shufcorrelations);
    %Report results of t test
    t_corr_NTCCs.tstat
    t_corr_NTCCs.df
    p_corr_NTCCs
    corr_cohens_d_NTCCs = mean(TCC_correlations-TCC_shufcorrelations)/std(TCC_correlations-TCC_shufcorrelations)

    %Perform Wilcoxon signed rank test
    [p_WSRT_NTCC h_WSRT_NTCC stats_WSRT_NTCC] = signrank(NTCC_correlations);
    %Report results
    stats_WSRT_NTCC.signedrank
    p_WSRT_NTCC
    
    %Compare TCC to NTCC Cells
    [h_corr_NTCC_TCC p_corr_NTCC_TCC CI_corr_NTCC_TCC t_corr_NTCC_TCC] ...
        = ttest2(NTCC_correlations,TCC_correlations);
    %Report results of t test
    t_corr_NTCC_TCC.tstat
    t_corr_NTCC_TCC.df
    p_corr_NTCC_TCC
    pooledsd = sqrt(((length(EMGindex)-1)*var(TCC_correlations)+(349-length(EMGindex)-1)*var(NTCC_correlations))/(349-2));
    corr_cohens_d_TCC_NTCC=(mean(TCC_correlations)-mean(NTCC_correlations))/pooledsd
    %Perform Wilcoxon rank sum test
    [p_WSRT_NTCC_TCC h_WSRT_NTCC_TCC stats_WSRT_NTCC_TCC] = ranksum(NTCC_correlations,TCC_correlations);
    %Report results
    stats_WSRT_NTCC_TCC.ranksum
    p_WSRT_NTCC_TCC
    
    
    
    
    
    
    
    %% Population Similarity Analysis
    
    %Load up maximum firing rate information to normalize cells
    load('Data/MaxFiringRate.mat');
    %Holds when blocks start for each session
    start_index=[];
    %iterator
    iter = 0;
    %We fix the block length to 30 to keep things similar across sessions
    used_block_len = 30;
    %Store the actual block lengths
    true_block_len(1:17) = 30;
    true_block_len(18:43) = 60;
    %This will hold the cosine similarities at each lag
    cosine_lag = zeros(69, 2*used_block_len-1);
    %This is used to determine how many times each dot_lag is calculated
    cosine_lag_counter = zeros(69,2*used_block_len-1);
    
    %This saves the image numbers that start each block
    for i = [1:15 17]
        
        start_index(i,:) = [1 61 31];
        
    end
    
    for i = 18:43
        
        start_index(i,:) = [1 61 NaN];
        
    end
    
    %Skips Non LEC Units in session 29
    for ses_no = [1:28 30:length(data.session.name)]
        
        %Checks for empty session 16
        if ~isempty(data.session.imgID{ses_no}) && length(data.session.imgID{ses_no})>60
            
            %find neurons from that session
            neuron_count = 0;
            clear pop_vec
            for neuron_no = 1:length(data.spks)
                if strcmp(data.name{neuron_no}(1:8), data.session.name{ses_no})
                    neuron_count = neuron_count + 1;
                    %grand_neuron_count = grand_neuron_count+1;
                    ses_neurons{ses_no}(neuron_count) = neuron_no;
                end
            end
            % find population vector corresponding to the first and the second presentation of each image
            for img_no = 1:max(data.session.imgID{ses_no}) % loop through all images
                for neuron_no = 1:neuron_count % loop through all neurons
                    % find neurons that were recorded during both presentation of
                    % the image
                    if length(find(data.imgID{ses_neurons{ses_no}(neuron_no)} == img_no)) == 2
                        img_neuron_id = find(data.imgID{ses_neurons{ses_no}(neuron_no)} == img_no);
                        
                        %Second Presentation
                        pop_vec(img_no,2,neuron_no) = (nanmean(data.spks{ses_neurons{ses_no}(neuron_no)}(img_neuron_id(2), 1:5500)))/maxfiring(ses_neurons{ses_no}(neuron_no));
                        %First Presentation
                        pop_vec(img_no,1,neuron_no) = (nanmean(data.spks{ses_neurons{ses_no}(neuron_no)}(img_neuron_id(1), 1:5500)))/maxfiring(ses_neurons{ses_no}(neuron_no));
                        
                    end
                end
            end
            
            %Gets number of blocks in the session
            number_of_blocks(ses_no) = floor(length(data.session.imgID{ses_no})/(2*true_block_len(ses_no)));
            %Loops through blocks
            for block_no = 1:number_of_blocks(ses_no)
                %Loops through blocks in dot_lag and
                %dot_lag_counter.
                iter = iter + 1;
                %Sees what image number to start with.
                start_ind = start_index(ses_no, block_no);
                for img_no = start_ind:start_ind+used_block_len-1
                    if length(find(data.session.imgID{ses_no} == img_no)) == 2 % if the image was presented twice
                        for img_no_diff = start_ind:start_ind+used_block_len-1
                            cosine_lag(iter, img_no_diff-img_no+used_block_len) =...
                                cosine_lag(iter, img_no_diff-img_no+used_block_len) ...
                                + (1-pdist([reshape(pop_vec(img_no_diff,1,:), 1, neuron_count);reshape(pop_vec(img_no,2,:), 1, neuron_count)], 'cosine'));
                            %+ dot(pop_vec(img_no_diff,1,:),pop_vec(img_no,2,:));
                            cosine_lag_counter(iter, img_no_diff-img_no+used_block_len) = ...
                                cosine_lag_counter(iter, img_no_diff-img_no+used_block_len) + 1;
                        end
                    end
                end
            end
        end
    end
    
    
    %Store normalized cosine similarities
    cosine_lag_all = cosine_lag./cosine_lag_counter;
    
    %This will hold all the values that will go into the regression
    regression_df = nan(207, 3);
    %Iterator
    iter = 1;
    %Loop through blocks
    for i = 1:69
        %Loop through lags -1, 0, and 1
        for j = 29:31
            %Store block, lag, and cosine similarity
            regression_df(iter,:) = [i j-30 cosine_lag_all(i,j)];
            iter = iter+1;
        end
    end
    
    %Exports to be used by R
    csvwrite('lag_cos_df.csv', regression_df);
    
    %Plus 1 t-test
    [h_corr_Plus1 p_corr_Plus1 CI_corr_Plus1 t_corr_Plus1] = ttest(cosine_lag_all(:,30),cosine_lag_all(:,31));
    %Report results of t test
    t_corr_Plus1.tstat
    t_corr_Plus1.df
    p_corr_Plus1
    corr_cohens_d_Plus1=(nanmean(cosine_lag_all(:,30))-nanmean(cosine_lag_all(:,31)))/(nanstd(cosine_lag_all(:,30)-cosine_lag_all(:,31)))
    
    %-1 t-test
    [h_corr_Neg1 p_corr_Neg1 CI_corr_Neg1 t_corr_Neg1] = ttest(cosine_lag_all(:,30),cosine_lag_all(:,29));
    %Report results of t test
    t_corr_Neg1.tstat
    t_corr_Neg1.df
    p_corr_Neg1
    corr_cohens_d_Neg1=(nanmean(cosine_lag_all(:,30))-nanmean(cosine_lag_all(:,29)))/(nanstd(cosine_lag_all(:,30)-cosine_lag_all(:,29)))
    
    
    
    
    
    %Permutation test
    %Subset data to just lags -1 to 1
    subset_data = cosine_lag_all(:,29:31);
    %Get actual mean difference of Lag 0 to Lag 1
    one_back_mean_diff = nanmean(subset_data(:,2) - subset_data(:,1));
    %Get actual mean difference of Lag 0 to Lag -1
    one_forward_mean_diff = nanmean(subset_data(:,2) - subset_data(:,3));
    
    %Time to permute
    %Number of permutations
    permutations = 100000;
    %holds the mean differences for one back and one forward
    perm_diffs = zeros(permutations, 2);
    
    for i = 1:permutations
        
        %Get size of data
        [M, N] = size(subset_data);
        %Save row indices
        rowIndex = repmat((1:M)',[1 N]);
        %Gets randomize column indices
        [~, randomizedColIndex] = sort(rand(M,N),2);
        %Helps to created B
        newLinearIndex = sub2ind([M,N],rowIndex, randomizedColIndex);
        B = subset_data(newLinearIndex);
        
        %Get mean difference of Lag 0 to Lag 1
        perm_diffs(i,1) = nanmean(B(:,2) - B(:,1));
        %Get mean difference of Lag 0 to Lag -1
        perm_diffs(i,2) = nanmean(B(:,2) - B(:,3));
        
    end
    
    %Permuted value for lag -1
    perm_lag_neg_1 = sum(perm_diffs(:,1) >=one_back_mean_diff)
    %Permuted value for lag +1
    perm_lag_pos_1 = sum(perm_diffs(:,2) >=one_forward_mean_diff)
    
    %% NTCC Analysis
    
    %Load up maximum firing rate information to normalize cells
    load('Data/MaxFiringRate.mat');
    %Holds when blocks start for each session
    start_index=[];
    %iterator
    iter = 0;
    %We fix the block length to 30 to keep things similar across sessions
    used_block_len = 30;
    %Store the actual block lengths
    true_block_len(1:17) = 30;
    true_block_len(18:43) = 60;
    %This will hold the cosine similarities at each lag
    cosine_lag = zeros(69, 2*used_block_len-1);
    %This is used to determine how many times each dot_lag is calculated
    cosine_lag_counter = zeros(69,2*used_block_len-1);
    
    %This saves the image numbers that start each block
    for i = [1:15 17]
        
        start_index(i,:) = [1 61 31];
        
    end
    
    for i = 18:43
        
        start_index(i,:) = [1 61 NaN];
        
    end
    
    %Skips Non LEC Units in session 29
    for ses_no = [1:28 30:length(data.session.name)]
        
        %Checks for empty session 16
        if ~isempty(data.session.imgID{ses_no}) && length(data.session.imgID{ses_no})>60
            
            %find neurons from that session
            neuron_count = 0;
            clear pop_vec
            for neuron_no = 1:length(data.spks)
                if strcmp(data.name{neuron_no}(1:8), data.session.name{ses_no}) & ~ismember(neuron_no, EMGindex)
                    neuron_count = neuron_count + 1;
                    %grand_neuron_count = grand_neuron_count+1;
                    ses_neurons{ses_no}(neuron_count) = neuron_no;
                end
            end
            % find population vector corresponding to the first and the second presentation of each image
            for img_no = 1:max(data.session.imgID{ses_no}) % loop through all images
                for neuron_no = 1:neuron_count % loop through all neurons
                    % find neurons that were recorded during both presentation of
                    % the image
                    if length(find(data.imgID{ses_neurons{ses_no}(neuron_no)} == img_no)) == 2
                        img_neuron_id = find(data.imgID{ses_neurons{ses_no}(neuron_no)} == img_no);
                        
                        %Second Presentation
                        pop_vec(img_no,2,neuron_no) = (nanmean(data.spks{ses_neurons{ses_no}(neuron_no)}(img_neuron_id(2), 1:5500)))/maxfiring(ses_neurons{ses_no}(neuron_no));
                        %First Presentation
                        pop_vec(img_no,1,neuron_no) = (nanmean(data.spks{ses_neurons{ses_no}(neuron_no)}(img_neuron_id(1), 1:5500)))/maxfiring(ses_neurons{ses_no}(neuron_no));
                        
                    end
                end
            end
            
            %Gets number of blocks in the session
            number_of_blocks(ses_no) = floor(length(data.session.imgID{ses_no})/(2*true_block_len(ses_no)));
            %Loops through blocks
            for block_no = 1:number_of_blocks(ses_no)
                %Loops through blocks in dot_lag and
                %dot_lag_counter.
                iter = iter + 1;
                %Sees what image number to start with.
                start_ind = start_index(ses_no, block_no);
                for img_no = start_ind:start_ind+used_block_len-1
                    if length(find(data.session.imgID{ses_no} == img_no)) == 2 % if the image was presented twice
                        for img_no_diff = start_ind:start_ind+used_block_len-1
                            cosine_lag(iter, img_no_diff-img_no+used_block_len) =...
                                cosine_lag(iter, img_no_diff-img_no+used_block_len) ...
                                + (1-pdist([reshape(pop_vec(img_no_diff,1,:), 1, neuron_count);reshape(pop_vec(img_no,2,:), 1, neuron_count)], 'cosine'));
                            %+ dot(pop_vec(img_no_diff,1,:),pop_vec(img_no,2,:));
                            cosine_lag_counter(iter, img_no_diff-img_no+used_block_len) = ...
                                cosine_lag_counter(iter, img_no_diff-img_no+used_block_len) + 1;
                        end
                    end
                end
            end
        end
    end
    
    
    %Store normalized cosine similarities
    cosine_lag_all = cosine_lag./cosine_lag_counter;
    
    %This will hold all the values that will go into the regression
    regression_df = nan(207, 3);
    %Iterator
    iter = 1;
    %Loop through blocks
    for i = 1:69
        %Loop through lags -1, 0, and 1
        for j = 29:31
            %Store block, lag, and cosine similarity
            regression_df(iter,:) = [i j-30 cosine_lag_all(i,j)];
            iter = iter+1;
        end
    end
    
    %Exports to be used by R
    %csvwrite('lag_cos_df.csv', regression_df);
    
    %Plus 1 t-test
    [h_corr_Plus1 p_corr_Plus1 CI_corr_Plus1 t_corr_Plus1] = ttest(cosine_lag_all(:,30),cosine_lag_all(:,31));
    %Report results of t test
    t_corr_Plus1.tstat
    t_corr_Plus1.df
    p_corr_Plus1
    corr_cohens_d_Plus1=(nanmean(cosine_lag_all(:,30))-nanmean(cosine_lag_all(:,31)))/(nanstd(cosine_lag_all(:,30)-cosine_lag_all(:,31)))
    
    %-1 t-test
    [h_corr_Neg1 p_corr_Neg1 CI_corr_Neg1 t_corr_Neg1] = ttest(cosine_lag_all(:,30),cosine_lag_all(:,29));
    %Report results of t test
    t_corr_Neg1.tstat
    t_corr_Neg1.df
    p_corr_Neg1
    corr_cohens_d_Neg1=(nanmean(cosine_lag_all(:,30))-nanmean(cosine_lag_all(:,29)))/(nanstd(cosine_lag_all(:,30)-cosine_lag_all(:,29)))
    
    
    
    
    
    %Permutation test
    %Subset data to just lags -1 to 1
    subset_data = cosine_lag_all(:,29:31);
    %Get actual mean difference of Lag 0 to Lag 1
    one_back_mean_diff = nanmean(subset_data(:,2) - subset_data(:,1));
    %Get actual mean difference of Lag 0 to Lag -1
    one_forward_mean_diff = nanmean(subset_data(:,2) - subset_data(:,3));
    
    %Time to permute
    %Number of permutations
    permutations = 100000;
    %holds the mean differences for one back and one forward
    perm_diffs = zeros(permutations, 2);
    
    for i = 1:permutations
        
        %Get size of data
        [M, N] = size(subset_data);
        %Save row indices
        rowIndex = repmat((1:M)',[1 N]);
        %Gets randomize column indices
        [~, randomizedColIndex] = sort(rand(M,N),2);
        %Helps to created B
        newLinearIndex = sub2ind([M,N],rowIndex, randomizedColIndex);
        B = subset_data(newLinearIndex);
        
        %Get mean difference of Lag 0 to Lag 1
        perm_diffs(i,1) = nanmean(B(:,2) - B(:,1));
        %Get mean difference of Lag 0 to Lag -1
        perm_diffs(i,2) = nanmean(B(:,2) - B(:,3));
        
    end
    
    %Permuted value for lag -1
    perm_lag_neg_1 = sum(perm_diffs(:,1) >=one_back_mean_diff)
    %Permuted value for lag +1
    perm_lag_pos_1 = sum(perm_diffs(:,2) >=one_forward_mean_diff)
    
    
    %%
    %TCC Cells
    
    %Load up maximum firing rate information to normalize cells
    load('Data/MaxFiringRate.mat');
    %Holds when blocks start for each session
    start_index=[];
    %iterator
    iter = 0;
    %We fix the block length to 30 to keep things similar across sessions
    used_block_len = 30;
    %Store the actual block lengths
    true_block_len(1:17) = 30;
    true_block_len(18:43) = 60;
    %This will hold the cosine similarities at each lag
    cosine_lag = zeros(69, 2*used_block_len-1);
    %This is used to determine how many times each dot_lag is calculated
    cosine_lag_counter = zeros(69,2*used_block_len-1);
    
    %This saves the image numbers that start each block
    for i = [1:15 17]
        
        start_index(i,:) = [1 61 31];
        
    end
    
    for i = 18:43
        
        start_index(i,:) = [1 61 NaN];
        
    end
    
    %Skips Non LEC Units in session 29
    for ses_no = [1:28 30:length(data.session.name)]
        
        %Checks for empty session 16
        if ~isempty(data.session.imgID{ses_no}) && length(data.session.imgID{ses_no})>60
            
            %find neurons from that session
            neuron_count = 0;
            clear pop_vec
            for neuron_no = 1:length(data.spks)
                if strcmp(data.name{neuron_no}(1:8), data.session.name{ses_no}) & ismember(neuron_no, EMGindex)
                    neuron_count = neuron_count + 1;
                    %grand_neuron_count = grand_neuron_count+1;
                    ses_neurons{ses_no}(neuron_count) = neuron_no;
                end
            end
            if neuron_count < 2
                continue
            end
            % find population vector corresponding to the first and the second presentation of each image
            for img_no = 1:max(data.session.imgID{ses_no}) % loop through all images
                for neuron_no = 1:neuron_count % loop through all neurons
                    % find neurons that were recorded during both presentation of
                    % the image
                    if length(find(data.imgID{ses_neurons{ses_no}(neuron_no)} == img_no)) == 2
                        img_neuron_id = find(data.imgID{ses_neurons{ses_no}(neuron_no)} == img_no);
                        
                        %Second Presentation
                        pop_vec(img_no,2,neuron_no) = (nanmean(data.spks{ses_neurons{ses_no}(neuron_no)}(img_neuron_id(2), 1:5500)))/maxfiring(ses_neurons{ses_no}(neuron_no));
                        %First Presentation
                        pop_vec(img_no,1,neuron_no) = (nanmean(data.spks{ses_neurons{ses_no}(neuron_no)}(img_neuron_id(1), 1:5500)))/maxfiring(ses_neurons{ses_no}(neuron_no));
                        
                    end
                end
            end
            
            %Gets number of blocks in the session
            number_of_blocks(ses_no) = floor(length(data.session.imgID{ses_no})/(2*true_block_len(ses_no)));
            %Loops through blocks
            for block_no = 1:number_of_blocks(ses_no)
                %Loops through blocks in dot_lag and
                %dot_lag_counter.
                iter = iter + 1;
                %Sees what image number to start with.
                start_ind = start_index(ses_no, block_no);
                for img_no = start_ind:start_ind+used_block_len-1
                    if length(find(data.session.imgID{ses_no} == img_no)) == 2 % if the image was presented twice
                        for img_no_diff = start_ind:start_ind+used_block_len-1
                            cosine_lag(iter, img_no_diff-img_no+used_block_len) =...
                                cosine_lag(iter, img_no_diff-img_no+used_block_len) ...
                                + (1-pdist([reshape(pop_vec(img_no_diff,1,:), 1, neuron_count);reshape(pop_vec(img_no,2,:), 1, neuron_count)], 'cosine'));
                            %+ dot(pop_vec(img_no_diff,1,:),pop_vec(img_no,2,:));
                            cosine_lag_counter(iter, img_no_diff-img_no+used_block_len) = ...
                                cosine_lag_counter(iter, img_no_diff-img_no+used_block_len) + 1;
                        end
                    end
                end
            end
        end
    end
    
    
    %Store normalized cosine similarities
    cosine_lag_all = cosine_lag./cosine_lag_counter;
    
    %This will hold all the values that will go into the regression
    regression_df = nan(207, 3);
    %Iterator
    iter = 1;
    %Loop through blocks
    for i = 1:69
        %Loop through lags -1, 0, and 1
        for j = 29:31
            %Store block, lag, and cosine similarity
            regression_df(iter,:) = [i j-30 cosine_lag_all(i,j)];
            iter = iter+1;
        end
    end
    
    %Exports to be used by R
    %csvwrite('lag_cos_df.csv', regression_df);
    
    %Plus 1 t-test
    [h_corr_Plus1 p_corr_Plus1 CI_corr_Plus1 t_corr_Plus1] = ttest(cosine_lag_all(:,30),cosine_lag_all(:,31));
    %Report results of t test
    t_corr_Plus1.tstat
    t_corr_Plus1.df
    p_corr_Plus1
    corr_cohens_d_Plus1=(nanmean(cosine_lag_all(:,30))-nanmean(cosine_lag_all(:,31)))/(nanstd(cosine_lag_all(:,30)-cosine_lag_all(:,31)))
    
    %-1 t-test
    [h_corr_Neg1 p_corr_Neg1 CI_corr_Neg1 t_corr_Neg1] = ttest(cosine_lag_all(:,30),cosine_lag_all(:,29));
    %Report results of t test
    t_corr_Neg1.tstat
    t_corr_Neg1.df
    p_corr_Neg1
    corr_cohens_d_Neg1=(nanmean(cosine_lag_all(:,30))-nanmean(cosine_lag_all(:,29)))/(nanstd(cosine_lag_all(:,30)-cosine_lag_all(:,29)))
    
    
    
    
    
    %Permutation test
    %Subset data to just lags -1 to 1
    subset_data = cosine_lag_all(:,29:31);
    %Get actual mean difference of Lag 0 to Lag 1
    one_back_mean_diff = nanmean(subset_data(:,2) - subset_data(:,1));
    %Get actual mean difference of Lag 0 to Lag -1
    one_forward_mean_diff = nanmean(subset_data(:,2) - subset_data(:,3));
    
    %Time to permute
    %Number of permutations
    permutations = 100000;
    %holds the mean differences for one back and one forward
    perm_diffs = zeros(permutations, 2);
    
    for i = 1:permutations
        
        %Get size of data
        [M, N] = size(subset_data);
        %Save row indices
        rowIndex = repmat((1:M)',[1 N]);
        %Gets randomize column indices
        [~, randomizedColIndex] = sort(rand(M,N),2);
        %Helps to created B
        newLinearIndex = sub2ind([M,N],rowIndex, randomizedColIndex);
        B = subset_data(newLinearIndex);
        
        %Get mean difference of Lag 0 to Lag 1
        perm_diffs(i,1) = nanmean(B(:,2) - B(:,1));
        %Get mean difference of Lag 0 to Lag -1
        perm_diffs(i,2) = nanmean(B(:,2) - B(:,3));
        
    end
    
    %Permuted value for lag -1
    perm_lag_neg_1 = sum(perm_diffs(:,1) >=one_back_mean_diff)
    %Permuted value for lag +1
    perm_lag_pos_1 = sum(perm_diffs(:,2) >=one_forward_mean_diff)
    
    
    
    
    
    %% Generates the regular plots
    %Make sure to set params.plot_on = 1;
    %in ML_fit_init_v2.m before running this
    
    for i = EMGindex
        
        process_ML_results2(i,finalparams);
        
    end
    
    
    Bindex1 = sortrows(parameters,[4 6]);
    load('Data/f_average_buffalo.mat')
    % Heatplot sorted by tau
    font_size=16;
    clear f_time_cells_all f_time_cells_all_resid
    condition = 3;
    time_cells_indexes = table2array(Bindex1(:,1));
    peaks = table2array(Bindex1(:,5));
    for i=1:length(time_cells_indexes)
        if (peaks(i)>0)
            %f_time_cells_all(i,:)=squeeze(sum(f_average_single_cell(time_cells_indexes(i),1:12,5000-1000:5000+delay_length+1000)));%delay_length %3018
            %f_time_cells_all(i,:)=squeeze(f_average_single_cell(time_cells_indexes(i),ddi,1:end));%delay_length %3018
            f_time_cells_all(i,:)=squeeze(f_average_single_cell(time_cells_indexes(i),condition,1:5500));%delay_length %3018
            %f_time_cells_all(i,:)=f_time_cells_all(i,:)-min(f_time_cells_all(i,:));
            f_time_cells_all(i,:)=(f_time_cells_all(i,:)/max(f_time_cells_all(i,:)));
            %f_time_cells_all(i,:)=f_time_cells_all(i,:)/max(f_time_cells_all(i,1+1000:end-1000));
            
        else
            %f_time_cells_all(i,:)=squeeze(sum(f_average_single_cell(time_cells_indexes(i),1:12,5000-1000:5000+delay_length+1000)));%delay_length %3018
            %f_time_cells_all(i,:)=squeeze(f_average_single_cell(time_cells_indexes(i),ddi,1:end));%delay_length %3018
            f_time_cells_all(i,:)=squeeze(f_average_single_cell(time_cells_indexes(i),condition,1:5500));%delay_length %3018
            f_time_cells_all(i,:)=f_time_cells_all(i,:)+(1-max(f_time_cells_all(i,:)));
            f_time_cells_all(i,:)=-(f_time_cells_all(i,:))+1;
            f_time_cells_all(i,:)=((f_time_cells_all(i,:)/max(f_time_cells_all(i,:))));
            %f_time_cells_all(i,:)=f_time_cells_all(i,:)/max(f_time_cells_all(i,1+1000:end-1000));
        end
    end
    figure;
    imagesc(f_time_cells_all)
    colormap('parula')
    set(gca,'xlim',[0 5.5*1000])
    set(gca,'XTick',500:1000:5500)
    set(gca,'XTickLabel',{0, 1, 2, 3, 4, 5})
    set(gca, 'YTick', [])
    h=xlabel('Time [s]');
    set(h,'fontsize',font_size)
    h=ylabel('Neuron ID');
    set(h,'fontsize',font_size)
    set(gca,'fontsize',font_size)
    title('Normalized firing rate')
    set(gcf,'color',[1 1 1])
    colorbar('Ticks', [0 .5 1])
    % added for the revision
    %hold on, plot([600 600], [0 240], 'k','linewidth',1)
    saveas(gcf,sprintf('export_fig/final_heatplot_time_cells_all_tau_v%i.fig',params.plot_version));
    saveas(gcf,sprintf('export_fig/final_heatplot_time_cells_all_tau_v%i.jpeg',params.plot_version));
    
    figure()
    scatter(finalparams(:,2), finalparams(:,4), 75, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0])%, 'MarkerFaceAlpha', 0.25)
    xlim([0 5])
    set(gca, 'YTick', [0 5 10 15 20])
    set(gca,'xlim',[0 5])
    set(gca,'XTick',0:1:5)
    xlabel('Response Latency [s]');
    ylabel('Relaxation Time [s]');
    set(gca,'fontsize',font_size)
    box off
    
    figure()
    histogram(finalparams(:,4), 50);
    xlim([-.5 20]);
    %set(gca, 'xscale', 'log')
    xlabel('Relaxation Time');
    ylabel('Frequency');
    set(gca,'fontsize',font_size)
    box off
    
    figure()
    histogram(finalparams(:,4), 50, 'FaceColor', [.5 .5 .5]);
    xlim([-.5 20]);
    %set(gca, 'xscale', 'log')
    xlabel('Relaxation Time');
    ylabel('Frequency');
    set(gca,'fontsize',font_size)
    box off
    
    figure()
    histogram(finalparams(:,4), 50);
    xlim([-.5 20]);
    set(gca, 'xscale', 'log')
    xlabel('Relaxation Time');
    ylabel('Frequency');
    set(gca,'fontsize',font_size)
    box off
    
    figure()
    histogram(finalparams(:,4), 50, 'FaceColor', [.5 .5 .5]);
    xlim([-.5 20]);
    set(gca, 'xscale', 'log')
    xlabel('Relaxation Time');
    ylabel('Frequency');
    set(gca,'fontsize',font_size)
    box off
    
    quants = 0:.1:1;
    breaks = zeros(length(quants),1);
    
    for i = 2:length(quants)
        
        breaks(i) = quantile(finalparams(:,4),quants(i));
        
        
    end
    
    figure()
    histogram(finalparams(:,4),breaks)
    xlim([-.5 20]);
    ylim([0 15]);
    %set(gca, 'xscale', 'log')
    xlabel('Relaxation Time');
    ylabel('Frequency');
    set(gca,'fontsize',font_size)
    box off
    
    figure()
    histogram(finalparams(:,4),breaks, 'FaceColor', [.5 .5 .5])
    xlim([-.5 20]);
    ylim([0 15]);
    %set(gca, 'xscale', 'log')
    xlabel('Relaxation Time');
    ylabel('Frequency');
    set(gca,'fontsize',font_size)
    box off
    
    figure()
    histogram(finalparams(:,4),breaks)
    xlim([-.5 20]);
    ylim([0 15]);
    set(gca, 'xscale', 'log')
    xlabel('Log Relaxation Time');
    ylabel('Frequency');
    set(gca,'fontsize',font_size)
    box off
    
    figure()
    histogram(finalparams(:,4),breaks, 'FaceColor', [.5 .5 .5])
    xlim([-.5 20]);
    ylim([0 15]);
    set(gca, 'xscale', 'log')
    xlabel('Log Relaxation Time');
    ylabel('Frequency');
    set(gca,'fontsize',font_size)
    box off
    
    
    
    figure()
    histogram(all_correlations)
    xlim([-0.5 0.5])
    xlabel('Unit Correlation Coefficient')
    ylabel('Frequency')
    set(gca,'fontsize',font_size)
    box off
    
    figure()
    histogram(all_correlations, 'FaceColor', [.5 .5 .5])
    xlim([-0.5 0.5])
    xlabel('Unit Correlation Coefficient')
    ylabel('Frequency')
    set(gca,'fontsize',font_size)
    box off
    
    
    %Load in csv of CI's
    figure()
    errorbar([-1 0 1], [nanmean(cosine_lag_all(:,29)) nanmean(cosine_lag_all(:,30)) nanmean(cosine_lag_all(:,31))], ...
        [lagcosmeansCIs(1,8) lagcosmeansCIs(2,8) lagcosmeansCIs(3,8)], 'ko', 'MarkerSize', 6, ...
        'LineWidth', 1.5, 'MarkerFaceColor', 'black')
    xlim([-1.33 1.33])
    ylim([0.83 0.87])
    set(gca, 'XTick', -1:1)
    set(gca, 'YTick', 0.83:0.01:0.87)
    xlabel('Lag')
    ylabel('Cosine Similarity')
    
    %suppcells = [138 92 141 101 292 35 17 209 104 289 330 123];
