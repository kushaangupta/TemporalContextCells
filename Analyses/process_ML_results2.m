% Process ML results and display rasters for cells with index ind. 
% If ind is not given then process all the cells.

function process_ML_results2(ind, finalparameters)

ML_fit_init_v2  % run initialization script - set all parameters there 

%Loads the trial averaged firing rates
load(params.faverage)
%Loads the data that ML_fit_params_v2 used.
load(params.data)
%Loads the data that contains session info
load('theData180329.mat')

fs = 1000;
%Need this for code to work properly
global f_spikes t c train_test_flag with_T lb ub trial_length ut st pValue T xT


if ~exist('ind','var')
    ind = 1:size(data.spikes,1);
end


for cell_index=1:length(ind)
    cell_no = ind(cell_index);
    cell_no=cell_no
    if exist(sprintf('%s/glm_cell_%i_v_%i.mat',params.folder_mat,cell_no,params.plot_version), 'file')==2 
        load(sprintf('%s/glm_cell_%i_v_%i.mat',params.folder_mat,cell_no,params.plot_version))
        %Set up to plot any cell regardless of its significance
        if (pValue > -1)
        % make raster
        font_size=30;
        temp = figure;
        %Makes it so that the figure doesn't appear
        set(temp, 'Visible', 'off');
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        subplot(2,1,[1])
        for trial=1:data.number_of_trials(cell_no)
            spikes_trial_ind = data.spikes{cell_no,trial};
            trial_spks(trial) = length(spikes_trial_ind);
            %spikes_trial_ind = spikes_trial_ind(find(spikes_trial_ind>data.time_before...
                %& spikes_trial_ind<data.time_before+data.trial_length))-data.time_before; 
            hold on, plot(spikes_trial_ind,zeros(size(spikes_trial_ind))+trial,'k.','markersize',4);
            %plot(olddata.msViewing{cell_no,1}(trial)+500,data.number_of_trials(cell_no)-trial+1, 'r.','markersize',6);
        end
        set(gca,'xlim',[0 5500])
        set(gca,'ylim',[0 data.number_of_trials(cell_no)+1])
        % add the 0.6s line
        hold on, %plot([600 600], [0 data.number_of_trials(cell_no)+1], 'c:','linewidth',1)
        h=ylabel(sprintf('%i Trials',data.number_of_trials(cell_no)));
        set(h,'fontsize',font_size)
        k=xlabel('');
        set(k, 'fontsize',font_size)
        set(h,'fontsize',font_size)
        set(gca,'XTick',[])
        set(gca,'XTickLabel',[])
        set(gca,'YTick',[])
        set(gca,'YTickLabel',[])
        set(gca, 'linewidth', 4)
        %set(gca,'Ytick',(data.number_of_trials(cell_no)-(floor( data.number_of_trials(cell_no) / 60) * 60)):60:data.number_of_trials(cell_no)-30)
        %set(gca,'Yticklabel',(floor(data.number_of_trials(cell_no) / 60) * 60):-60:60)
        
        %pos=get(h,'Pos');
        %set(h,'Pos',[xpos pos(2) pos(3)])
        box on
        
        % plot average fit and average firing rate
        subplot(2,1,[2])
        C_tmp = reshape(C,data.trial_length,length(C)/data.trial_length);
        C_vec = mean(C_tmp');
        %hold on, %plot(fs*C_vec,'b:','linewidth',4)
        
        if (LL_I < LL_T)
            peak = xT(5);
            baseline = zeros(1,5500) +xT(1);
            hold on, plot(fs*baseline, 'black:', 'linewidth', 4)
            T_tmp = reshape(T,data.trial_length,length(T)/data.trial_length);
            T_vec = mean(T_tmp');
            plot(fs*squeeze(f_average_single_cell(cell_no,3,1:end)),'black','linewidth',4)
            plot(fs*T_vec,'magenta--','linewidth',4)
            set(gca,'xlim',[0 data.trial_length])
            set(gcf,'color','w')
            %legend('Constant term only','Time and constant term','Firing rate')
            h=ylabel('Firing rate [Hz]');
            set(h,'fontsize',font_size)
            h=xlabel('Time [s]');
            set(h,'fontsize',font_size)
            set(gca,'XTick', [])
            %set(gca,'XTickLabel',{0:1:(data.time_before+data.trial_length+data.time_after)/1000})
            set(gca,'fontsize',font_size) 
            set(gcf,'color',[1 1 1]) 
            set(gca, 'linewidth', 4)
            %max_y = round(max(max(fs*squeeze(f_average_single_cell(cell_no,:,1:end))')))+1;
            %set(gca,'ylim',[0 max_y])
        
        % add the 0.6sset(gca,'fontsize',font_size)  line
        %hold on, %plot([600 600], [0 max_y], 'c:','linewidth',1)
        
        %subplot(2,1,1), title(sprintf('cell %i Increasing Firing',cell_no))
        
        else
            peak = xT(5);
            baseline = zeros(1,5500) +xI(1);
            hold on, plot(fs*baseline, 'black:', 'linewidth', 4)
            T_tmp = reshape(I,data.trial_length,length(I)/data.trial_length);
            T_vec = mean(T_tmp');
            plot(fs*squeeze(f_average_single_cell(cell_no,3,1:end)),'black','linewidth',4)
            plot(fs*T_vec,'magenta--','linewidth',4)
            set(gca,'xlim',[0 data.trial_length])
            set(gcf,'color','w')
            %legend('Constant term only','Time and constant term','Firing rate')
            h=ylabel('Firing rate [Hz]');
            set(h,'fontsize',font_size)
            h=xlabel('Time [s]');
            set(h,'fontsize',font_size)
            set(gca,'XTick',[])
            %set(gca,'XTickLabel',{0:1:(data.time_before+data.trial_length+data.time_after)/1000})
            set(gca,'fontsize',font_size) 
            set(gcf,'color',[1 1 1])
            max_y = round(max(max(fs*squeeze(f_average_single_cell(cell_no,:,1:end))')))+1;
            set(gca,'ylim',[0 max_y])
            set(gca, 'linewidth', 4)
            % add the 0.6s line
            hold on, %plot([600 600], [0 max_y], 'c:','linewidth',1)
            
            %subplot(2,1,1), title(sprintf('cell %i Increasing Firing',cell_no))
            %subplot(2,1,1), title(sprintf('Cell 268, Session mp140117'))
            
            
        end
        %Saves the plot
        
        subplot(2,1,1);
        %n = title(sprintf('Relaxation time = %.2f s',finalparameters(find(finalparameters(:,1)==cell_no),4)));
        %set(n,'fontsize',font_size)
        
        %if (abs(peak) > .001)
            %saveas(gcf,sprintf('/projectnb/ecog-eeg/stim_specific_time_cells/Buffalo/finalcell_%i_v_%i.jpeg',cell_no,params.plot_version));
            saveas(gcf,sprintf('/projectnb/ecog-eeg/stim_specific_time_cells/Buffalo/corrcells_%i_v_%i.jpeg',cell_no,params.plot_version));
        
            %end
        end
        
    else
        warning('MISSING CELL %i',cell_no)
    end
end

