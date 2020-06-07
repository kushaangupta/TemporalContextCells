function ML_fit_params_v5(ii)

ML_fit_init_v5 % run initialization script - set all parameters there 

if params.submit 
    cell_no=str2double(ii);
else
    cell_no=ii;
    if params.plot_on
        load(params.faverage)
    end
end


load(params.data);

global f_spikes f_spikes_evenodd t c train_test_flag with_T lb ub x trial_length st ut

%particle swarming options
hybridopts = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000);
options = optimoptions('particleswarm','SwarmSize',50,'HybridFcn',{@fmincon,hybridopts}); 

%load('data_cromer.mat')

trial_length = data.trial_length;
time_before = 0; 
n_cells=size(data.spikes,1);
max_p=params.max_p;
    
lb = [];
ub = [];
cell_no = cell_no
f_spikes = [];
f_spikes_odd = [];
f_spikes_even = [];
for j=1:data.number_of_trials(cell_no) % loop through trials
    f_trial = zeros(trial_length,1);
    f_trial(ceil(data.spikes{cell_no,j}...
        (find(time_before <  data.spikes{cell_no,j} &...
        data.spikes{cell_no,j} < time_before + trial_length)))-time_before) = 1;
    f_spikes=[f_spikes; f_trial];
    if mod(j,2)
        f_spikes_odd=[f_spikes_odd; f_trial];
    else
        f_spikes_even=[f_spikes_even; f_trial];
    end 
end
f_spikes_evenodd = [f_spikes_even; f_spikes_odd];
t=repmat([1:trial_length], 1, length(f_spikes)/trial_length)';
t=t/1000;
t=t-mean(t);

%% GLM    

fun=@glm_model_stim_specific_v3;

%shuffle the trials
%rng(1)
%[f_spikes, f_spikes_i, f_spikes_j]=randblock(f_spikes,trial_length);

%with  no time

with_T=0;
n=1; %number of terms
%fit just the constant term
lb.o=10^-10; %offset coefficient - lower bound
ub.o=1/n; %offset coefficient - upper bount
lbC = lb;
ubC = ub;
[xC, xC_tmp, LL_C, LL_C_tmp, LL_C_trial_2blocks, aic_C, bic_C, C] = run_solver(params, fun, options);

with_T = 1;
n = 2;
lb.sig = .001;
ub.sig = 1;
lb.mu = -2.5;
ub.mu = 2.5;
lb.tau = .05;
ub.tau = 100;
lb.o = 10^-7;
ub.o=1/n;
lb.peak = 10^-7;
ub.peak = 1/n;
lbT = lb;
ubT = ub;

[xT, xT_tmp, LL_T, LL_T_tmp, LL_T_trial_2blocks, aic_T, bic_T, T] = run_solver(params, fun, options);




%%
%[hh_trial_2blocks,pValue_trial_2blocks,stat_trial_2blocks] = lratiotest(LL_T_trial_2blocks,LL_C_trial_2blocks,4);
%if sum(pValue_trial_2blocks<max_p)==2 %by requiring 2 out of 2 to be significant (0.01) we are 100 times less likely to get the significance by chance
 %   two_blocks_pass=1;
%else
%    two_blocks_pass=0;
%end

save(sprintf('%s/glm_examp_exp_cell_%i_v_%i.mat', params.folder_mat,cell_no,params.plot_version));

%Tests the fit of the EMG w/ positive peak compared to constant
%Has four more parameters
[hh,pValue,stat] = lratiotest(LL_T,LL_C,4);
%Tests the fit of the EMG w/ negative peak compared to constant
%Has four more parameters


%Clears out data to make the size of saved file smaller
clear data;
%Saves results of script
save(sprintf('%s/glm_examp_exp_cell_%i_v_%i.mat', params.folder_mat,cell_no,params.plot_version));



if params.plot_on
    figure
    fs = 1000;
    C_tmp = reshape(C,trial_length,length(C)/trial_length);
    C_vec = mean(C_tmp');
    hold on, plot(fs*C_vec)
    T_tmp = reshape(T,trial_length,length(T)/trial_length);
    T_vec = mean(T_tmp');
    plot(fs*T_vec)
    plot(fs*squeeze(f_average_single_cell(cell_no,1,time_before:end)))
    set(gca,'xlim',[0 trial_length])
    set(gcf,'color','w')
    title(sprintf('cell %i, time cell: %i',cell_no, time_cell))
    legend('Constat','Time','Firing rate')
    ylabel('Hz')
    xlabel('ms')
    export_fig(sprintf('%s/cell_%i_v_%i.pdf',params.folder_results,cell_no,params.plot_version));
    saveas(gcf,sprintf('%s/cell_%i_v_%i.fig',params.folder_results,cell_no,params.plot_version));
end
            
%save(sprintf('%s/glm_cell_%i_v_%i.mat', params.folder_mat,cell_no,params.plot_version), 't', 'dAIC', 'dBIC', 'xC_tmp', 'xT_tmp', 'LL_C_tmp', 'LL_T_tmp', 'aic_C', 'aic_T', 'bic_C', 'bic_T', 'hh','pValue','stat','hh_trial_2blocks','pValue_trial_2blocks','stat_trial_2blocks','hh_trial_3blocks','pValue_trial_3blocks','stat_trial_3blocks','three_blocks_pass','two_blocks_pass','one_block_pass','LL_T_trial_2blocks','LL_T_trial_3blocks','LL_C_trial_2blocks','LL_C_trial_3blocks','LL_C','LL_T','xT','xC','T','C','f_spikes','f_fr_vec_per_cell','time_cell', 'lb', 'ub');
clear data f_spikes_even f_spikes_evenodd f_spikes_odd C_tmp CS13 CS14 f_average_single_cell RCS12_ RCS13_ RCS14_ T_tmp T_vec 
%save(sprintf('%s/glm_cell_%i_v_%i.mat', params.folder_mat,cell_no,params.plot_version));
if params.submit
    exit %must have when submitted as a job 
end

end
%dAIC=aic_P(:,2)-aic_TP :if possitive time contributes more than the second
%position term with probability exp(d/2)https://en.wikipedia.org/wiki/Akaike_information_criterion
%The AIC penalizes the number of parameters less strongly than does the BIC
% so we should use AIC to be more conservative
%dBIC=bic_P(:,2)-bic_TP


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
   try
       i=i+1;
       xC_tmp(i,:)=particleswarm(fun,length(lb_vec),lb_vec,ub_vec,options);
       LL_C_tmp{i}=-fun(xC_tmp(i,:));
       if LL_C_tmp_max>=LL_C_tmp{i}
           c=c+1;
       else
           LL_C_tmp_max=LL_C_tmp{i};
           xC_tmp_max=xC_tmp(i,:);
           c=1;
       end
       if c>=params.N
           stop_loop=1;
       end
   catch
       warning('Had an issue, but its probably cool? Lets try again and pretend it didnt happen.')
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



        
       
%compile localy
%mcc -mv -o stim_spec_time_cells ML_fit_params.m

%compile on the cluster
%mcc -mv -o stim_spec_time_cells ML_fit_params.m -R -singleCompThread -R -nodisplay
%mcc -mv -o stim_spec_time_cells_16a ML_fit_params.m -R -singleCompThread -R -nodisplay
