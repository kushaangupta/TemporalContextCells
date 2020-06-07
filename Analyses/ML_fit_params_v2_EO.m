function ML_fit_params_v2(ii)

%Loads up all the initialization stuff
ML_fit_init_v2_EO

%Pulls cell number
cell_no=ii;

%Loads the spike data
load(params.data);
%Variables that need to be passable between functions
global f_spikes f_spikes_evenodd t c train_test_flag with_T lb ub x trial_length st ut y

%particle swarming options
hybridopts = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000);
options = optimoptions('particleswarm','SwarmSize',50,'HybridFcn',{@fmincon,hybridopts});

%Gets length of trial
trial_length = data.trial_length;
%Gets time before trial
time_before = data.time_before;
%Gets size of data
n_cells=size(data.spikes,1);
%Gets significance level
max_p=params.max_p;

%Holds lower bounds
lb = [];
%Holds upper bounds
ub = [];
%Prints cell number
cell_no = cell_no
%Holds all of the spikes
f_spikes = [];
f_spikes_odd = [];
f_spikes_even = [];

%loops through the trials to reformat for solver
for j=1:data.number_of_trials(cell_no)
    %Will hold all the spikes in the time window I care about
    %Each millisecond bin either has a spike (1) or doesn't (0)
    f_trial = zeros(trial_length,1);
    %Gets the spikes in the given trial
    data.spikes{cell_no,j};
    temp1 = ans;
    %Only loads up the spikes during the trial
    f_trial(temp1(time_before <  temp1 & temp1 < time_before + trial_length)) = 1;
    %Stores the spikes
    f_spikes=[f_spikes; f_trial];
    
    if mod(j,2)
        f_spikes_odd=[f_spikes_odd; f_trial];
    else
        f_spikes_even=[f_spikes_even; f_trial];
    end
    
end

%Reformats time to repeat the same period over and over again for solver
t=repmat([1:trial_length], 1, length(f_spikes)/trial_length)';
%Transforms time to milliseconds
t=t/1000;
%Sets t=0 to the middle of trial for solver
t=t-mean(t);

%% GLM

%Reads in the parameters from ML_fit_params_v2
finalparams = readtable('finalparameters.csv','HeaderLines',1);
%Turns the table into a matrix
finalparams = finalparams{:,:};
%Get the row with the current cell
row_num = find(finalparams(:,1) == cell_no);

%Checks if cell was significant before. If it wasn't the code exits,
%because we don't need to check if even and odd fits were significant.

if isempty(row_num)
    exit
end

%Loads what equations are going to be used.
fun=@glm_model_stim_specific_v3;

%Constant model odd (No Temporal Response)
f_spikes = f_spikes_odd;
%Reformats time to repeat the same period over and over again for solver
t=repmat([1:trial_length], 1, length(f_spikes)/trial_length)';
t=t/1000;
t=t-mean(t);

with_T=0;
%number of terms in equation (Probability must be between 0 and 1)
n=1;
%Constant value bounds
%Lower bound
lb.o=10^-10;
%Upper bound
ub.o=1/n;
%All the bounds for constant model
%Lower bounds
lbC = lb;
%Upper bounds
ubC = ub;
%Runs the solver and stores all the information from the fit
[xCO, xCO_tmp, LL_CO, LL_CO_tmp, LL_CO_trial_2blocks, aic_CO, bic_CO, CO] = run_solver(params, fun, options);

%Constant model even (No Temporal Response)
f_spikes = f_spikes_even;
%Reformats time to repeat the same period over and over again for solver
t=repmat([1:trial_length], 1, length(f_spikes)/trial_length)';
t=t/1000;
t=t-mean(t);

with_T=0;
%number of terms in equation (Probability must be between 0 and 1)
n=1;
%Constant value bounds
%Lower bound
lb.o=10^-10;
%Upper bound
ub.o=1/n;
%All the bounds for constant model
%Lower bounds
lbC = lb;
%Upper bounds
ubC = ub;
%Runs the solver and stores all the information from the fit
[xCE, xCE_tmp, LL_CE, LL_CE_tmp, LL_CE_trial_2blocks, aic_CE, bic_CE, CE] = run_solver(params, fun, options);

if finalparams(row_num,5) > 0
    
    %EMG with positive peak Even
    with_T = 1;
    %number of terms in the model (baseline and EMG)
    n = 2;
    %Sigma value bounds
    %Lower bound
    lb.sig = .001;
    %Upper bound
    ub.sig = 1;
    %Mu value bounds
    %format puts t=0 at very center of trial
    %Lower bound of mu (corresponds to trial onset)
    lb.mu = -2.25;
    %Upper bound of mu (corresponds to trial end)
    ub.mu = 2.75;
    %Tau value bounds
    %Lower bound of tau (model uses 1/tau, so 1/.05 = 20 secs)
    lb.tau = .05;
    %Upper bound of tau (model uses 1/tau, so 1/100 = 1 ms)
    ub.tau = 100;
    %Baseline firing rate bounds
    %Lower bound of baseline firing rate
    lb.o = 10^-7;
    %Upper bound of baseline firing rate
    ub.o=1/n;
    %Peak height bounds
    %Lower bound of peak height
    lb.peak = 10^-7;
    %Upper bound of peak height
    ub.peak = 1/n;
    %Stores all of the lower bounds
    lbT = lb;
    %Stores all of the upper bounds
    ubT = ub;
    %Runs the solver and stores all the information from the fit.
    [xTE, xTE_tmp, LL_TE, LL_TE_tmp, LL_TE_trial_2blocks, aic_TE, bic_TE, TE] = run_solver(params, fun, options);
    
    f_spikes = f_spikes_odd;
    %Reformats time to repeat the same period over and over again for solver
    t=repmat([1:trial_length], 1, length(f_spikes)/trial_length)';
    t=t/1000;
    t=t-mean(t);
    
    %EMG with positive peak odd trials
    with_T = 1;
    %number of terms in the model (baseline and EMG)
    n = 2;
    %Sigma value bounds
    %Lower bound
    lb.sig = .001;
    %Upper bound
    ub.sig = 1;
    %Mu value bounds
    %format puts t=0 at very center of trial
    %Lower bound of mu (corresponds to trial onset)
    lb.mu = -2.25;
    %Upper bound of mu (corresponds to trial end)
    ub.mu = 2.75;
    %Tau value bounds
    %Lower bound of tau (model uses 1/tau, so 1/.05 = 20 secs)
    lb.tau = .05;
    %Upper bound of tau (model uses 1/tau, so 1/100 = 1 ms)
    ub.tau = 100;
    %Baseline firing rate bounds
    %Lower bound of baseline firing rate
    lb.o = 10^-7;
    %Upper bound of baseline firing rate
    ub.o=1/n;
    %Peak height bounds
    %Lower bound of peak height
    lb.peak = 10^-7;
    %Upper bound of peak height
    ub.peak = 1/n;
    %Stores all of the lower bounds
    lbT = lb;
    %Stores all of the upper bounds
    ubT = ub;
    %Runs the solver and stores all the information from the fit.
    [xTO, xTO_tmp, LL_TO, LL_TO_tmp, LL_TO_trial_2blocks, aic_TO, bic_TO, TO] = run_solver(params, fun, options);
    
    
    
else
    %EMG with negative peak Even
    with_T = 1;
    %number of terms in the model (baseline and EMG)
    n = 2;
    %Sigma value bounds
    %Lower bound
    lb.sig = .001;
    %Upper bound
    ub.sig = 1;
    %Mu value bounds
    %format puts t=0 at very center of trial
    %Lower bound of mu (corresponds to trial onset)
    lb.mu = -2.25;
    %Upper bound of mu (corresponds to trial end)
    ub.mu = 2.75;
    %Tau value bounds
    %Lower bound of tau (model uses 1/tau, so 1/.05 = 20 secs)
    lb.tau = .05;
    %Upper bound of tau (model uses 1/tau, so 1/100 = 1 ms)
    ub.tau = 100;
    %Baseline firing rate bounds
    %Lower bound of baseline firing rate
    lb.o = 10^-7;
    %Upper bound of baseline firing rate
    ub.o=1/n;
    %Peak height bounds
    %Lowest value of peak height (Cannot go below zero)
    lb.peak = -xCE;
    %Highest value of peak height
    ub.peak = -(10^-7);
    %Stores all of the lower bounds
    lbT = lb;
    %Stores all of the upper bounds
    ubT = ub;
    %Runs the solver and stores all the information from the fit.
    [xTE, xTE_tmp, LL_TE, LL_TE_tmp, LL_TE_trial_2blocks, aic_TE, bic_TE, TE] = run_solver(params, fun, options);
    
    f_spikes = f_spikes_odd;
    %Reformats time to repeat the same period over and over again for solver
    t=repmat([1:trial_length], 1, length(f_spikes)/trial_length)';
    t=t/1000;
    t=t-mean(t);
    
    %EMG with positive peak odd trials
    with_T = 1;
    %number of terms in the model (baseline and EMG)
    n = 2;
    %Sigma value bounds
    %Lower bound
    lb.sig = .001;
    %Upper bound
    ub.sig = 1;
    %Mu value bounds
    %format puts t=0 at very center of trial
    %Lower bound of mu (corresponds to trial onset)
    lb.mu = -2.25;
    %Upper bound of mu (corresponds to trial end)
    ub.mu = 2.75;
    %Tau value bounds
    %Lower bound of tau (model uses 1/tau, so 1/.05 = 20 secs)
    lb.tau = .05;
    %Upper bound of tau (model uses 1/tau, so 1/100 = 1 ms)
    ub.tau = 100;
    %Baseline firing rate bounds
    %Lowest value of peak height (Cannot go below zero)
    lb.peak = -xCO;
    %Highest value of peak height
    ub.peak = -(10^-7);
    %Stores all of the lower bounds
    lbT = lb;
    %Stores all of the upper bounds
    ubT = ub;
    %Runs the solver and stores all the information from the fit.
    [xTO, xTO_tmp, LL_TO, LL_TO_tmp, LL_TO_trial_2blocks, aic_TO, bic_TO, TO] = run_solver(params, fun, options);
    
end


%%
%Clears out data to make the size of saved file smaller
clear data;

save(sprintf('%s/glm_cell_%i_v_%i.mat', params.folder_mat,cell_no,params.plot_version));

%Tests the fit of the EMG w/ positive peak compared to constant
%Has four more parameters
try
[hhE,pValueE,statE] = lratiotest(LL_TE,LL_CE,4);
catch
end
%Tests the fit of the EMG w/ negative peak compared to constant
%Has four more parameters
try
[hh_O,pValue_O,stat_O] = lratiotest(LL_TO,LL_CO,4);
catch
end


%Saves results of script
save(sprintf('%s/glm_cell_%i_v_%i.mat', params.folder_mat,cell_no,params.plot_version));

if params.plot_on
    %Can be used to generate the plot
end

%Checks if it was submitted on cluster
if params.submit
    %Closes program on cluster
    exit
end

end


function [xC, xC_tmp, LL_C, LL_C_tmp, LL_C_trial_2blocks, aic_C, bic_C, C] =...
    run_solver(params, fun, options)
global train_test_flag f_spikes lb ub
%global f_spikes t train_test_flag with_T lb ub trial_length
%create vectors from upper and lower bounds structures since
%particleswarm requires vector input
lb_cell=struct2cell(lb);
lb_vec=[lb_cell{:}];
ub_cell=struct2cell(ub);
ub_vec=[ub_cell{:}];
train_test_flag=1; %training - return the fit sum
tic
stop_loop=0;
c=1;
LL_C_tmp_max=-Inf;
i=0;
while stop_loop==0
    %Helps to handle occasions when probability is not between 0 and 1.
    try
        i=i+1;
        xC_tmp(i,:)=particleswarm(fun,length(lb_vec),lb_vec,ub_vec,options);
        LL_C_tmp{i}=-fun(xC_tmp(i,:));
        if LL_C_tmp_max>=LL_C_tmp{i}
            c=c+1
        else
            LL_C_tmp_max=LL_C_tmp{i};
            xC_tmp_max=xC_tmp(i,:);
            c=1
        end
        if c>=params.N
            stop_loop=1;
        end
    catch
    end
    
end
LL_C=LL_C_tmp_max;
xC=xC_tmp_max;
if isfield(lb,'t') || length(lb_vec) == 1
    train_test_flag=2; %compute LL per trial - divide trials into 2 blocks
    LL_C_trial_2blocks=-fun(xC);
else
    LL_C_trial_2blocks = [];
end
[aic_C,bic_C]=aicbic(LL_C,length(lb_vec),length(f_spikes));
toc
train_test_flag=0; %testing - rerun the fit
C=fun(xC);
clear f_spikes t train_test_flag with_T lb ub trial_length
end




%Make sure to compile with matlab/2016a
%compile localy
%mcc -mv -o stim_spec_time_cells ML_fit_params.m

%compile on the cluster
%mcc -mv -o stim_spec_time_cells_16a ML_fit_params.m -R -singleCompThread -R -nodisplay
