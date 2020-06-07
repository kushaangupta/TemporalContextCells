%% This does preprocessing to complete rest of paper.
warning('off','all')
%These will hold the parameters for the convolution model
%Holds all of the unit sigmas
Expsigs = [];
%Holds all of the unit mus
Expmus = [];
%Holds all of the unit taus
Exptaus = [];
%Holds all of the unit indices
Expindex = [];
%Holds all of the unit peak values
Exppeak = [];
%Holds all of the unit baseline firing values
Expcons = [];
%Counts how many units were significantly fit by the
%convolution model
Expcount = 1;

%Loads up the initialization script that generated the fits.
ML_fit_init_v5

%Finds all the parameters for exp cells
for i = 1:70
    
    %Prints what unit is it on
    cell_no=i
    try
    %Loads the fit for the given unit i
    load(sprintf('%s/glm_examp_exp_cell_%i_v_%i.mat',params.folder_mat,cell_no,params.plot_version))
    
    %Checks if a version of the convolution model fit the spikes better
    %than the constant
    if  (pValue < .05/349)

            %Stores unit sigma
            Expsigs(Expcount) = xT(2);
            %Stores unit mu
            %Need to move t=0 back to actual 0.
            Expmus(Expcount) = xT(3) + 2.5;
            %Stores unit tau
            %Report as 1/tau in model
            Exptaus(Expcount) = 1/xT(4);
            %Stores baseline firing rate
            Expcons(Expcount) = xT(1);
            %Stores unit peak height
            Exppeak(Expcount) = xT(5);
            %Stores the unit index
            Expindex(Expcount) = i;
            %Increases iterator
            Expcount = Expcount + 1;
            
    end
    catch
    end
    
    
    
end

[R P] = corr(Expmus', Exptaus')

%Loads up the initialization script that generated the fits.
ML_fit_init_v5

smoothFRS = zeros(50,5000);

for i = 2:51
    
    cell_no=i
    
    %Loads the fit for the given unit i
    load(sprintf('%s/glm_examp_exp_cell_%i_v_%i.mat',params.folder_mat,cell_no,params.plot_version))
    
    tempp = reshape(f_spikes, [5000,1000]);
    smoothFRS((i-1),:) = smooth(mean(tempp,2));
    
    
end

figure()
imagesc(smoothFRS)
set(gca, 'XTick', [])
set(gca, 'YTick', [])

ML_fit_init_v4

smoothFRS = zeros(50,5000);

for i = 2:51
    
    cell_no=i
    
    %Loads the fit for the given unit i
    load(sprintf('%s/glm_examp_time_cell_%i_v_%i.mat',params.folder_mat,cell_no,params.plot_version))
    
    tempp = reshape(f_spikes, [5000,1000]);
    smoothFRS((i-1),:) = smooth(mean(tempp,2));
    
    
end

figure()
imagesc(smoothFRS)
set(gca, 'XTick', [])
set(gca, 'YTick', [])

figure()
hold on
iter = 0;
for i = [63 53 43 33 23 13]

%Prints what unit is it on
    cell_no=i
    
    %Loads the fit for the given unit i
    load(sprintf('%s/glm_examp_exp_cell_%i_v_%i.mat',params.folder_mat,cell_no,params.plot_version))
    
    tempp = reshape(f_spikes, [5000,1000]);
    temppp = smooth(mean(tempp,2));
    plot(1:5000, temppp, 'Color', [.6-iter .6-iter .6-iter], 'LineWidth', 1)
    
    
    iter = iter + .1;
end

ylim([0 .42])
set(gca, 'XTick', 0:1000:5000)
set(gca, 'XTickLabel',{0, 1, 2, 3, 4, 5})
set(gca, 'YTick', 0)
set(gca, 'YTickLabel','')
ylabel('Activation')
xlabel('Time [s]')


%These will hold the parameters for the convolution model
%Holds all of the unit sigmas
TCsigs = [];
%Holds all of the unit mus
TCmus = [];
%Holds all of the unit taus
TCtaus = [];
%Holds all of the unit indices
TCindex = [];
%Holds all of the unit peak values
TCpeak = [];
%Holds all of the unit baseline firing values
TCcons = [];
%Counts how many units were significantly fit by the
%convolution model
TCcount = 1;

%Loads up the initialization script that generated the fits.
ML_fit_init_v4


%Finds all the parameters for exp cells
for i = 1:45
    try
    %Prints what unit is it on
    cell_no=i
    
    %Loads the fit for the given unit i
    load(sprintf('%s/glm_examp_time_cell_%i_v_%i.mat',params.folder_mat,cell_no,params.plot_version))
    
    %Checks if a version of the convolution model fit the spikes better
    %than the constant
    if  (pValue < .05/349)

            %Stores unit sigma
            TCsigs(TCcount) = xT(2);
            %Stores unit mu
            %Need to move t=0 back to actual 0.
            TCmus(TCcount) = xT(3) + 7.5;
            %Stores unit tau
            %Report as 1/tau in model
            TCtaus(TCcount) = 1/xT(4);
            %Stores baseline firing rate
            TCcons(TCcount) = xT(1);
            %Stores unit peak height
            TCpeak(TCcount) = xT(5);
            %Stores the unit index
            TCindex(TCcount) = i;
            %Increases iterator
            TCcount = TCcount + 1;
            
    end
    catch
        
        
    end
    
    
    
end

mdl = fitlm(TCmus',TCtaus');



%Loads up the initialization script that generated the fits.
ML_fit_init_v4

figure()
hold on
iter = 0;
for i = [45 36 27 18 9 1]

%Prints what unit is it on
    cell_no=i
    
    %Loads the fit for the given unit i
    load(sprintf('%s/glm_examp_time_cell_%i_v_%i.mat',params.folder_mat,cell_no,params.plot_version))
    
    tempp = reshape(f_spikes, [15000,1000]);
    temppp = smooth(mean(tempp,2));
    plot(1:5000, temppp(1:5000), 'Color', hsv2rgb([0 .4+iter, 1]), 'LineWidth', 1.75)
    set(gca, 'XTick', 0:1000:5000)
    set(gca, 'XTickLabel',{0, 1, 2, 3, 4, 5})
    set(gca, 'YTick', 0)
    set(gca, 'YTickLabel','')
    xlabel('Time [s]')
    ylabel('Activation')
    ylim([0 .42])
    
    
    iter= iter+.12;
    
end




figure()
plot(.001:.001:5, .7991*(.001:.001:5)-.4011, 'color', hsv2rgb([0 .4, 1]))
hold on
iter = 0;
for i = [63 53 43 33 23 13]

%Prints what unit is it on
    cell_no=i
    
    
    scatter(Expmus(i),Exptaus(i), 80, [.6-iter .6-iter .6-iter],'filled')
    xlim([0 5])
    ylim([0 5])
    
    iter = iter + .1;
end

iter = 0;
for i = [45 36 27 18 9 1]

%Prints what unit is it on
    cell_no=i
    
    
    scatter(TCmus(i),TCtaus(i), 80, hsv2rgb([0 .4+iter 1]),'filled')
    xlim([0 5])
    ylim([0 5])
    
    iter = iter + .12;
end

xlabel('\mu')
ylabel('\tau')

set(gca, 'YTick', 0:1:5)
set(gca, 'YTickLabel',{0, 1, 2, 3, 4, 5})
set(gca, 'box', 'off')

