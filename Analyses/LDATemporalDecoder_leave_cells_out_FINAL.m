%{
This analysis was used to generate Figure 3.
In this analysis, spiking data is binned and sorted then used to train
and test a Linear Discriminate Classifier
This process is repeated as cell with faster relaxation times are removed.
This demonstrates that the population as a whole contains temporal
information.


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
    Not used for Figure 3.  Used for supplementary figure S5c.  Used to
    cut earlier time bins.
%}
%cross temporal classifier of time

%use cleaned_data_v1.mat to access the actual data set
%%Configuration and Setup
filename='cleaned_data_v1.mat'; %skip 238-245 

%the desired size of bins for analysis
bin_size_constant=250;
load(filename)
bins_max=(data.trial_length-500)/bin_size_constant;

%not used for this implementation
glob_itr=0;

%clear the workspaces of unneeded variables
clearvars -except bins_max glob_itr bin_size_constant filename
permutations=1000;

%save is set to 0, because the plots appearance needs to be manually
%adjusted before 
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
paramtable=readtable('finalparameters.csv');

figure(1); hold on;
line([0 5], [1.65 1.65],'Color','black')
xlabel('Actual Time')
ylabel('Decoding Error')
tau_thresh=0:.125:2.5;
overall_averaged_error(length(tau_thresh))=0;
overall_averaged_std(length(tau_thresh))=0;
meanabserrormean_perm(permutations,length(tau_thresh))=0;
%% Main loop, repeating the analysis on successively less cells
for itr=1:length(tau_thresh)
    unit_list_indexes=find(table2array(paramtable(:,4))>=tau_thresh(itr));
    unit_list=table2array(paramtable(unit_list_indexes,1));

%Remove units with low numbers of trials
omit=find(data.number_of_trials<30);
%omit=union(omit,238:245);
data.spikes(omit,:)=[];
data.number_of_trials(omit)=[];
unit_list=setdiff(unit_list,omit);

it=500+glob_itr*bin_size_constant; %initial time cut out in ms
f_bin_size=bin_size_constant;
c_bin_size=bin_size_constant;

%How many trials are cells with low number of trials bootstrapped too
%And how many trials are cells with high higher numbers of trial are
%subsampled from
trials=100;


traingroup_size=1;
testgroup_size=1;

%Initialize Training and testing data, trial number by number of cells in size
fr_all_train=zeros(f_bins*trials,length(unit_list));
fr_all_test=zeros(f_bins*trials,length(unit_list));

%Left_ol right
left_ol=1;
right_ol=left_ol-1;

%% Bin the spiking rate data and sort into training and testing trials
for trial=1:trials
    unit_index=0;
for unit_index=1:length(unit_list)
%for unit=units
    unit=unit_list(unit_index);
    for bin=1:f_bins
        
        train_trial=randsample(1:2:data.number_of_trials(unit),traingroup_size);   %train odd   
        test_trial=randsample(2:2:data.number_of_trials(unit),testgroup_size);    % test even
        
        %{
        odds=1:2:data.number_of_trials(unit);
        evens=2:2:data.number_of_trials(unit);
        train_trial=odds(trial);
        test_trial=evens(trial);
        %}
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
%intialize the posterior distribution
posterior_averaged=zeros(f_bins,c_bins);
expected=ones(f_bins,c_bins)/f_bins;
DoF=f_bins*(c_bins-1);

%more intializiation
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
unpermuttedmean=abserrormean;
unpermuttedstd=abserrorstd;
%figure; imagesc(posterior_averaged);
%colormap(gray)
%colorbar
%xlabel('Decoded Time Bin')
%ylabel('Actual Time Bin')
%title('Posterior Distribution')
%saveas(gcf,sprintf('%s_PosteriorAveraged_bin%ims_binOL%i_grouping%i-%i.fig',savename,f_bin_size,left_ol+right_ol,traingroup_size,testgroup_size))
%saveas(gcf,sprintf('%s_PosteriorAveraged_bin%ims_binOL%i_grouping%i-%i.png',savename,f_bin_size,left_ol+right_ol,traingroup_size,testgroup_size))

%create the Figure for Posterior Distribution (Figure 3a)
figure(3); hold on;
imagesc(log(1/(100)+posterior_averaged),[log(.02) log(.5)]); colormap(gray);

yl1=ylim;
xl1=xlim;
for bin_i=1:c_bins
   max_f_bin=find(max(posterior_averaged(:,bin_i))==posterior_averaged(:,bin_i));
   plot(bin_i,max_f_bin,'r.','MarkerSize',16);
end
line([itr/2, itr/2],[.5, bins_max],'Color','r','LineWidth',3.0,'LineStyle','--')
line([.5, bins_max],[itr/2, itr/2],'Color','r','LineWidth',3.0,'LineStyle','--')

set(gca,'YDir','reverse')
%xlim(xl1);
%ylim(yl1);
xlim([.5 bins_max+.5])
ylim([.5 bins_max+.5])

ylabel('Decoded Time (s)')
xlabel('Actual Time (s)')
title(sprintf('Log of Posterior Distribution (%0.3f)',tau_thresh(itr)))
set(gca,'xtick',[1:4:f_bins])
set(gca,'xticklabel',{.001*[0:(4*f_bin_size):data.trial_length]})
set(gca,'ytick',[1:4:c_bins])
set(gca,'yticklabel',{.001*[0:(4*c_bin_size):data.trial_length]})

pbaspect([1 1 1])
saveas(gcf,sprintf('Video_decodingheatmap3/heatmap_thresh_%.2i_%0.3f.png',itr,tau_thresh(itr)));
clf(3)

%% final figure
%mdl = fitlm(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),.001*c_bin_size*abserrormean,'y ~ x1');

%figure; plot(cat_ind_train(1:f_bins),rmean,'-')
%plot(mdl)
%plot(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),mdl.Coefficients{1,1}+mdl.Coefficients{2,1}*.001*f_bin_size*cat_ind_train(1:f_bins),'k-')
%plot(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),1.25*ones(size(.001*f_bin_size*(cat_ind_train(1:f_bins)-1))),'k--')



overall_averaged_std(itr)=std(abserrormean)*.001*c_bin_size;
overall_averaged_error(itr)=mean(.001*c_bin_size*abserrormean);


if save
saveas(gcf,sprintf('%s_AbsError_bin%ims_binOL%i_grouping%i-%i.fig',savename,f_bin_size,left_ol+right_ol,traingroup_size,testgroup_size))
saveas(gcf,sprintf('%s_AbsError_bin%ims_binOL%i_grouping%i-%i.png',savename,f_bin_size,left_ol+right_ol,traingroup_size,testgroup_size))
end

%brown time cells: 179=#B3, 141=#8D, 107=#6B [179/255 141/255 107/255]
%purple exponentially decaying cells: 191=#BF, 129=#81, 255=#FF
%figure; hold on;
%plot(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),mdl.Coefficients{1,1}+mdl.Coefficients{2,1}*.001*f_bin_size*cat_ind_train(1:f_bins),'k-')
%errorbar(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),.001*c_bin_size*abserrormean,.001*c_bin_size*abserrorstd./sqrt(trials),'ko','MarkerSize',7);
%xlabel('Time (s)')
%yl=ylim;
%ylim([-.05 1.65]);
%xl=xlim;
%xlim([-.1 4.9])
%ylabel('Error In Decoding (s)')
%title('Absolute Value of Decoding Error')

%plot of values along diagonal of posterior matrix.  Not used.
diagonal_bins(f_bins)=0;
for i=1:f_bins
    diagonal_bins(i)=posterior_averaged(i,i);
end
Cells_cut_threshold=tau_thresh(itr)
cells_used=length(unit_list)

%% Permutation Test of Decoder Performance
for perm=1:permutations
    %shuffle bins
    fr_all_train=fr_all_train(randperm(f_bins*trials),:);
    %run the classifier
    [result,err,posterior,logp,coef] = classify(fr_all_test,fr_all_train,categorical(cat_ind_train)','linear');
    %intialize values
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
        %abserrors((1+(bin-1)*trials):(1+bin*trials))=abs(cat_ind_train(bin)-grp2idx(result(bin:f_bins:(f_bins*(trials-1)+bin))));
    end
    meanabserrormean_perm(perm,itr)=mean(abserrormean);
end
close 1
%{
figure(1); hold on
plot(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),.001*c_bin_size*unpermuttedmean,'Marker','.','LineStyle','-','MarkerSize',10);
errorbar(.001*f_bin_size*(cat_ind_train(1:f_bins)-1),abserrormean*bin_size_constant/1000,1.9623*abserrorstd/sqrt(1000)*bin_size_constant/1000);
xlabel('Time (s)')
xlim([-.1 4.9])
ylim([-.1 2])

ylabel('Error In Decoding (s)')
saveas(gcf,sprintf('PNGs/Error_auto/decoding_error_%.2i_%0.3f.png',itr,tau_thresh(itr)));
%}
figure (1);

end
%legend('Chance Performance','All TC cells (128 Cells)','Tau > .5 (42 Cells)','Tau > 1 (29 Cells)','Tau > 1.5 (17 Cells)','Tau > 2 (15 Cells)', 'Tau > 2.5 (14 Cells)')
figure(2); hold on;
%errorbar(tau_thresh,overall_averaged_error,overall_averaged_std,'bo-');
plot(tau_thresh,overall_averaged_error,'Marker','.','Color','black','LineStyle','-','Markersize',18);
%1.9623 is 95% confidence interval with degrees of freedom equal to 999
%errorbar(tau_thresh,mean(meanabserrormean_perm(:,:),1)*bin_size_constant/1000,1.9623*std(meanabserrormean_perm(:,:)*bin_size_constant/1000));
errorbar(tau_thresh,median(meanabserrormean_perm(:,:),1)*bin_size_constant/1000,(median(meanabserrormean_perm(:,:),1)-quantile(meanabserrormean_perm(:,:),.025))*bin_size_constant/1000,(-median(meanabserrormean_perm(:,:),1)+quantile(meanabserrormean_perm(:,:),.975))*bin_size_constant/1000,'Color','black');
%line([0 tau_thresh(end)], [1.65 1.65],'Color','Black')
xlabel('Tau > Threshold')
ylabel('Decoding Error')
legend('Average Decoding Error','Permutted Decoding Error')
ylim([-.1 2])
xlim([-.1 2.6])

%make video
workingDir='Video_decodingheatmap2';
outputVideo = VideoWriter(fullfile(workingDir,'Heatmap_video_out.avi'));
outputVideo.FrameRate = 10;
open(outputVideo);

imageNames = dir(fullfile(workingDir,'*.png'));
imageNames = {imageNames.name}';

for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,imageNames{ii}));
   writeVideo(outputVideo,img)
end

close(outputVideo)