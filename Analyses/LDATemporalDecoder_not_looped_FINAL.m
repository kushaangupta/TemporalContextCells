%{
This analysis was used to generate Figure S4 c.
In this analysis, spiking data is binned and sorted then used to train
and test a Linear Discriminate Classifier

%50, 100, 250, 500
bins_all=bins_max-glob_itr; %11, 22, 55 110

variables that can be set to adjust the performance of this
analysis include:
bin_size_constant
    Bin size might be set to 50, 100, 250, 500 to get 110, 55, 22, 11 bins 
    (respectively).  It is currently set to 250 to get reasonable 
    classification accuracy while still having a reasonable number of bins. 

permutations
    Controls how many permutations are used in the permutation testing of
    classifier performance.  1000 permutations provides a strong statistical strength,
    but is slower to run.  

trials
    number of trials used for training classifier and the number of trials
    used for testing the classifier.  A higher number does get slightly
    better performance, but also runs much slower.

Other variable configuration options not used in the final implementation,
these configuration options are not fully bug-tested and may be incomplete: 

c_bins, f_bins
    bins used for classification category (c_bins) and bins used as input 
    (f_bins) can be set seperately although there is little use case for 
    this.  Currently, they are set to be the same value and set
    using bin_size_constant

traingroup_size, testgroup_size
    trials can be clustered together to improve classifier performance,
    This option is currently not used (a cluster of 1)

left_ol, right_overlap
    Sets the training and testing bins to overlap on the left and right
    sides.  This introduces a mathematical problem with the independence of
    the bins and was thus not used (no overlap, i.e. overlap set to 0).

glob_itr
    Not used for this analysis.  Used for supplementary figure S5c.  Used to
    cut earlier time bins.
%}

%cross temporal classifier of time
%filename='cleaned_data_v1.mat'; 

filename='simulated_exponential_decay_cells_20Hz_max20Hz';
%filename='simulated_ideal_time_cells_1Hz_max40Hz';
bin_size_constant=250;
load(filename)
bins_max=(data.trial_length-500)/bin_size_constant;

glob_itr=0;

%clearvars
clearvars -except bins_max glob_itr bin_size_constant filename
%load('simulated_ideal_time_cells_k4_constant_50.mat')
%close all
save=0;
%number of classifier category bins
%50, 100, 250, 500
bins_all=bins_max-glob_itr; %11, 22, 55 110
c_bins=bins_all;
%bins of firing rate
f_bins=bins_all;

%%load the data and average firing rate in bins for training trials
%randomly subsample without replacement training trials and testing trials
load(filename)
load('params.mat')
[unit_count,trial_count]=size(data.spikes);

%comment out omission of certain cells when running on simuated data
%{
omit=find(data.number_of_trials<30);
data.spikes(omit,:)=[];
data.number_of_trials(omit)=[];
unit_count=unit_count-length(omit);
%}

it=500+glob_itr*bin_size_constant; %initial time cut out in ms
f_bin_size=bin_size_constant;
c_bin_size=bin_size_constant;

%500 trials used for training and 500 trials used for testing
trials=500;
traingroup_size=1;
testgroup_size=1;

fr_all_train=zeros(f_bins*trials,unit_count);
fr_all_test=zeros(f_bins*trials,unit_count);

left_ol=1;
right_ol=left_ol-1;

%% Bin the spiking rate data and sort into training and testing trials
for trial=1:trials
    unit_index=0;
for unit=1:unit_count
%for unit=units
    unit_index=unit_index+1;
    for bin=1:f_bins
        
        train_trial=randsample(1:2:data.number_of_trials(unit),traingroup_size);   %train odd   
        test_trial=randsample(2:2:data.number_of_trials(unit),testgroup_size);    % test even
        spikes_count_train=0;
        spikes_count_test=0;
        for tt=train_trial;
        spikes_count_train=spikes_count_train+sum(((bin-left_ol)*f_bin_size+it) < data.spikes{unit,tt} & data.spikes{unit,tt} < ((bin+right_ol)*f_bin_size)+it);
        end
        for tt=test_trial
        spikes_count_test=spikes_count_test+sum(((bin-left_ol)*f_bin_size+it) < data.spikes{unit,tt} & data.spikes{unit,tt} < ((bin+right_ol)*f_bin_size)+it);
        end
        fr_all_train((trial-1)*f_bins+bin,unit_index)=(spikes_count_train)/((min(bin,left_ol)+min(f_bins-bin,right_ol))*f_bin_size)/traingroup_size+rand*10^-13;
        fr_all_test((trial-1)*f_bins+bin,unit_index)=(spikes_count_test)/((min(bin,left_ol)+min(f_bins-bin,right_ol))*f_bin_size)/testgroup_size+rand*10^-13;
    end
end
end

%set the true category according to the bin count
%works best if f_bins is multiple of c_bins
cat_ind_train(f_bins*trials)=0;
for trial=1:trials
for bin=1:f_bins
    cat_ind_train((trial-1)*f_bins+bin)=ceil((bin/f_bins)*c_bins);
end
end

%% run the classifier.  See the matlab help on classify for options
%'linear' uses a linear discriminate analysis
[result,err,posterior,logp,coef] = classify(fr_all_test,fr_all_train,categorical(cat_ind_train)','linear');
%% Analyze and plot the classifier performance
posterior_averaged=zeros(f_bins,c_bins);
expected=ones(f_bins,c_bins)/f_bins;
DoF=f_bins*(c_bins-1);

%Intialize values
rmean(f_bins)=0;
abserrormean(f_bins)=0;
abserrorstd(f_bins)=0;
for bin=1:f_bins
    posterior_averaged(bin,:)=mean(posterior(bin:f_bins:(f_bins*(trials-1)+bin),:),1);
    rerror(bin)=std(grp2idx(result(bin:f_bins:(f_bins*(trials-1)+bin))));
    rmean(bin)=mean(grp2idx(result(bin:f_bins:(f_bins*(trials-1)+bin))));
    abserrormean(bin)=mean(abs(cat_ind_train(bin)-grp2idx(result(bin:f_bins:(f_bins*(trials-1)+bin)))));
    abserrorstd(bin)=std(abs(cat_ind_train(bin)-grp2idx(result(bin:f_bins:(f_bins*(trials-1)+bin)))));
end

figure; imagesc(posterior_averaged);
colormap(gray)
colorbar
xlabel('Decoded Time Bin')
ylabel('Actual Time Bin')
title('Posterior Distribution')
%saveas(gcf,sprintf('%s_PosteriorAveraged_bin%ims_binOL%i_grouping%i-%i.fig',savename,f_bin_size,left_ol+right_ol,traingroup_size,testgroup_size))
%saveas(gcf,sprintf('%s_PosteriorAveraged_bin%ims_binOL%i_grouping%i-%i.png',savename,f_bin_size,left_ol+right_ol,traingroup_size,testgroup_size))

figure; hold on;
imagesc(log(1/(100)+posterior_averaged),[log(.02) log(.5)]); colormap(gray);
yl1=ylim;
xl1=xlim;

for bin_i=1:c_bins
   max_f_bin=find(max(posterior_averaged(:,bin_i))==posterior_averaged(:,bin_i));
   plot(bin_i,max_f_bin,'r.','MarkerSize',16);
end
set(gca,'YDir','reverse')
%xlim(xl1);
%ylim(yl1);
xlim([.5 bins_max+.5])
ylim([.5 bins_max+.5])
%colorbar
ylabel('Decoded Time (s)')
xlabel('Actual Time (s)')
%title('Log of Posterior Distribution (.01 Threshold)')
set(gca,'xtick',[1:4:f_bins])
set(gca,'xticklabel',{.001*[0:(4*f_bin_size):data.trial_length]})
set(gca,'ytick',[1:4:c_bins])
set(gca,'yticklabel',{.001*[0:(4*c_bin_size):data.trial_length]})
pbaspect([1 1 1])
colorbar
if save
saveas(gcf,sprintf('%s_PosteriorLog_bin%ims_binOL%i_grouping%i-%i.fig',savename,f_bin_size,left_ol+right_ol,traingroup_size,testgroup_size))
saveas(gcf,sprintf('%s_PosteriorLog_bin%ims_binOL%i_grouping%i-%i.png',savename,f_bin_size,left_ol+right_ol,traingroup_size,testgroup_size))
end

%% final figure
mdl = fitlm(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),.001*c_bin_size*abserrormean,'y ~ x1')

%figure; plot(cat_ind_train(1:f_bins),rmean,'-')
figure; hold on;
%plot(mdl)
plot(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),mdl.Coefficients{1,1}+mdl.Coefficients{2,1}*.001*f_bin_size*cat_ind_train(1:f_bins),'k--')
%plot(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),1.25*ones(size(.001*f_bin_size*(cat_ind_train(1:f_bins)-1))),'k--')
plot(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),.001*c_bin_size*abserrormean,'Color','red','Marker','.','LineStyle','-','MarkerSize',26);
xlabel('Time (s)')
xlim([-.1 4.9])
yl=ylim;
%ylim([0 (yl(2)+.1)])
ylim([-.05 .8])
%ylim([0 2.5])
ylabel('Error In Decoding (s)')

if save
saveas(gcf,sprintf('%s_AbsError_bin%ims_binOL%i_grouping%i-%i.fig',savename,f_bin_size,left_ol+right_ol,traingroup_size,testgroup_size))
saveas(gcf,sprintf('%s_AbsError_bin%ims_binOL%i_grouping%i-%i.png',savename,f_bin_size,left_ol+right_ol,traingroup_size,testgroup_size))
end

figure; hold on;
plot(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),mdl.Coefficients{1,1}+mdl.Coefficients{2,1}*.001*f_bin_size*cat_ind_train(1:f_bins),'r-')
%errorbar(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),.001*c_bin_size*abserrormean,.001*c_bin_size*abserrorstd./sqrt(trials),'ro','MarkerSize',7);
plot(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),.001*c_bin_size*abserrormean,'r.-','MarkerSize',16);
xlabel('Time (s)')
%yl=ylim;
%ylim([-.05 .8]);
xl=xlim;
xlim([-.1 4.9])
ylabel('Error In Decoding (s)')
%title('Absolute Value of Decoding Error')

diagonal_bins(f_bins)=0;
for i=1:f_bins
    diagonal_bins(i)=posterior_averaged(i,i);
end

up_to=it-500
%The bin are not indepedent, bins that are closer together are likely to be
%correlated, so a chi squared test is not the actually the right
%statistical test
p_value_chi2=1-chi2cdf(sum(sum((trials*(posterior_averaged-expected).^2)./expected)),DoF)
[~, p_value_proptest]=ztest(diagonal_bins,expected(1,1),sqrt(expected(1,1)*(1-expected(1,1))/trials),'Tail','right')

overall_averaged_error=mean(.001*c_bin_size*abserrormean)
