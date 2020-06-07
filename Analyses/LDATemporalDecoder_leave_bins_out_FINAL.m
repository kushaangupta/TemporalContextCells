%{
This analysis was used to generate Figure S5.
In this analysis, spiking data is binned and sorted then used to train
and test a Linear Discriminate Classifier
This process is repeated as bins from earlier in the interval are dropped.
This demonstrates that temporal decoding can be performed across the
interval, even as earlier, easier to decode time bins are dropped


variables that can be set to adjust the performance of this
analysis include:
bin_size_constant
    Bin size might be set to 50, 100, 250, 500 to get 110, 55, 22, 11 bins 
    (respectively, for the first iteration with no bins dropped).  
    It is currently set to 250 to get reasonable 
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
%}
%% Intialization
filename='cleaned_data_v1.mat'; %skip 238-245 

bin_size_constant=250;
permutations=1000;
meanabserrormean_perm(permutations)=0;

load(filename)
bins_max=(data.trial_length-500)/bin_size_constant;
significant_permutations(bins_max-1)=0;
%glob_itr=0;
%for glob_itr=0:(bins_max-2)
%% Repeat for different numbers of bins cut out of the analysis
for glob_itr=0:(bins_max-2)
%clear the work space each loop to keep memory usage down and to avoid
%intializing variables with non-zero values
clearvars -except bins_max glob_itr bin_size_constant filename permutations meanabserrormean_perm sig_flag_1 sig_flag_2 sig_flag_3 sig_itrs tstats p_ significant_permutations

save=0;
%number of classifier category bins
%50, 100, 250, 500
bins_all=bins_max-glob_itr; %11, 22, 55 110
c_bins=bins_all;
%bins of firing rate
f_bins=bins_all;

load(filename)
load('params.mat')
[unit_count,trial_count]=size(data.spikes);
%omit cells with low number of trials and other problems
omit=find(data.number_of_trials<30);
omit=union(omit,238:245);
omit=union(omit,Bindex1(:,1));

data.spikes(omit,:)=[];
data.number_of_trials(omit)=[];

%randcells=datasample(1:(unit_count-length(omit)),129,'Replace',false);
%data.spikes=data.spikes(randcells,:);
%data.number_of_trials=data.number_of_trials(randcells);
unit_count=unit_count-length(omit);
%unit_count=129;
%}
data.spikes=data.spikes;

it=500+glob_itr*bin_size_constant; %initial time cut out in ms
f_bin_size=bin_size_constant;
c_bin_size=bin_size_constant;


trials=100;
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
%% Aanalyze and plot the classifier performance
posterior_averaged=zeros(f_bins,c_bins);

rerror(f_bins)=0;
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

figure(3+glob_itr); hold on;
bigger_matrix=zeros(bins_max,bins_max);
%upper_bound=min(.5,max(max(posterior_averaged)));
upper_bound=.5;
bigger_matrix(:,:)=log(upper_bound);
bigger_matrix((bins_max-c_bins+1):(bins_max),(bins_max-f_bins+1):(bins_max))=log(1/(100)+posterior_averaged);

%imagesc(bigger_matrix,[log(1/(c_bins)) log(upper_bound)]); colormap(gray);
imagesc(bigger_matrix,[log(.02) log(.5)]); colormap(gray);

yl1=ylim;
xl1=xlim;
line([glob_itr+.5, glob_itr+.5],[.5, bins_max],'Color','r','LineWidth',3.0,'LineStyle','--')
line([.5, bins_max],[glob_itr+.5, glob_itr+.5],'Color','r','LineWidth',3.0,'LineStyle','--')
for bin_i=1:c_bins
   max_f_bin=find(max(posterior_averaged(:,bin_i))==posterior_averaged(:,bin_i));
   plot(bin_i+glob_itr,glob_itr+max_f_bin,'r.','MarkerSize',16);
end

set(gca,'YDir','reverse')
xlim([.5 bins_max+.5])
ylim([.5 bins_max+.5])
colorbar
ylabel('Decoded Time (s)')
xlabel('Actual Time (s)')
%title('Log of Posterior Distribution (.01 Threshold)')
set(gca,'xtick',[1:4:f_bins])
set(gca,'xticklabel',{.001*[0:(4*bin_size_constant):data.trial_length]})
set(gca,'ytick',[1:4:c_bins])
set(gca,'yticklabel',{.001*[0:(4*bin_size_constant):data.trial_length]})

pbaspect([1 1 1])
saveas(gcf,sprintf('PNGs/Posteriors/heatmap_linear_cut_%.2i.png',glob_itr));

meanabserrormean=mean(abserrormean);
expected = zeros(size(posterior_averaged));
bias=mean(posterior_averaged,1);
for i=1:c_bins
    expected(i,:)=bias;
end


%diagonal_bins(f_bins)=0;
%for i=1:f_bins
%    diagonal_bins(i)=posterior_averaged(i,i);
%end

meanabserrormean_perm(permutations)=0;

%% Test using permutations of the labels.
for perm=1:permutations
fr_all_train=fr_all_train(randperm(f_bins*trials),:);
[result,err,posterior,logp,coef] = classify(fr_all_test,fr_all_train,categorical(cat_ind_train)','linear');
posterior_averaged=zeros(f_bins,c_bins);
rerror(f_bins)=0;
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
meanabserrormean_perm(perm)=mean(abserrormean);
end
significant_permutations(glob_itr+1)=sum(meanabserrormean_perm<meanabserrormean);
up_to=it-500;
Actual_mean=meanabserrormean*bin_size_constant;
Mean_permutted_mean=mean(meanabserrormean_perm)*bin_size_constant;
if glob_itr==0
figure(1); hold on
histogram(meanabserrormean_perm*bin_size_constant/1000,16)
yl2=ylim;
line([Actual_mean/1000, Actual_mean/1000],[0, yl2(2)],'Color','k','LineWidth',3.0)
xl=xlim;
xlim([xl(1)-.05 xl(2)])
xlabel('Mean Absolute Value of Decoding Error (s)')
ylabel('Frequency')
end
Mean_permutted_std=std(meanabserrormean_perm)*bin_size_constant;
%}
%close 3
end

figure(2)
plot(bin_size_constant/1000*(0:(bins_max-2)),significant_permutations,'Marker','.','Color','black','LineStyle','-','Markersize',18);