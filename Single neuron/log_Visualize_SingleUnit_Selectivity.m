%% Log_Visualize_SingleUnit_Selectivity
%  This script produces two plots for a given neuron to illustrate its
%  response profile. First, it computes firing rate of individual neuron 
%  under different behavioral conditions. 
% 1st plot: a proportional histogram of firing rates for a given behavior, compared to rest. 
% 2nd plot: smooth spike trace with labeled time points according to a set of behaviors.
% C. Testard, July 2022

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,15:16,18];
session_range_with_partner=[1:6,11:13,15:16,18];

%Set parameters
with_partner =0;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 10*temp_resolution;%set the smoothing window size (sigma)
threat_precedence =1;
exclude_sq=1;
plot_toggle=0;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end


s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/Example_units'];

    %% Load data

    %% Get data with specified temporal resolution and channels
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

        disp('Data Loaded')

        Spike_count_raster = zscore(Spike_rasters');

        %Extract behavior labels
        behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure

        %Set colormap
        %uisetcolor([0.6 0.8 1])
        Cmap = [[1 0 0];...%Aggression; red
            [1 0.4 0.1];...%Approach; dark orange
            [0 0 0];...%But sniff; NA
            [0.3 0.7 1];...%Drinking; light blue
            [0 0.7 0];...%Foraging; dark green
            [1 0 1];...%Groom sollicitation; magenta
            [0 1 1];...%Groom partner; cyan
            [0 0 1];...%Getting groomed; dark blue
            [0.8 0 0];...%Threat to partner; dark red
            [1 0 0];...%Threat to subject; red
            [0.9 0.9 0];...%leave; dark yellow
            [0 0 0];...%Lipsmack
            [0.2 0.9 0.76];...%Masturbating; turquoise
            [0.7 0 1];...%Mounting; light purple
            [0.9 0.5 0];...%Other monkeys vocalize; orange
            [1 0.8 0.1];...%Travel; yellow orange
            [0 0 0];...%Proximity; NA
            [0 0 0];...%Rowdy room; NA
            [0.6314 0.5059 0.0118];...%Scratch; maroon
            [0.5 0.2 0.5];...%Self-groom; dark purple
            [ 1 0.07 0.65];...%Submission; dark pink
            [0 0.4 0.5];...%Vocalzation; blue green
            [0 0 0];...%Yawning; NA
            [0.8 0.8 0.8]];%Rest; grey

        %Exclude axctivity at the very beginning and very end of the
        %session
        session_with_buffer = [250*temp_resolution:size(Spike_count_raster,1)-250*temp_resolution];
        activity_all = Spike_count_raster(session_with_buffer,:)';
        activity_vlpfc = Spike_count_raster(session_with_buffer,strcmp(brain_label,'TEO'))';
        activity_teo = Spike_count_raster(session_with_buffer,strcmp(brain_label,'vlPFC'))';

        %% Plot actvity for one session in relation to behavior

        % Sort activity using rastermap.
        [isort1_all, isort2_all, amap_all] = mapTmap(activity_all); %all units considered
        [isort1_teo, isort2_teo, amap_teo] = mapTmap(activity_teo); %just teo units
        [isort1_vlpfc, isort2_vlpfc, amap_vlpfc] = mapTmap(activity_vlpfc); %just vlpfc units
        sorted_activity_all = activity_all(isort1_all,:);
        sorted_activity_byArea = [activity_vlpfc(isort1_vlpfc,:); activity_teo(isort1_teo,:)];

        if plot_toggle %plot rastermap
            figure; imagesc(amap_all); xline(block_times.start_time(2)*temp_resolution,'LineWidth',2); xline(block_times.start_time(3)*temp_resolution,'LineWidth',2); colorbar; caxis([-2 2])
            figure; imagesc(sorted_activity_all); xline(block_times.start_time(2)*temp_resolution,'LineWidth',2); xline(block_times.start_time(3)*temp_resolution,'LineWidth',2); colorbar; caxis([-2 2])

            figure; imagesc([amap_vlpfc;amap_teo]); xline(block_times.start_time(2)*temp_resolution,'LineWidth',2); xline(block_times.start_time(3)*temp_resolution,'LineWidth',2); colorbar; caxis([-2 2]); yline(size(amap_vlpfc,1),'LineWidth',2,'LineStyle','-')
            figure; imagesc(sorted_activity_byArea); xline(block_times.start_time(2)*temp_resolution,'LineWidth',2); xline(block_times.start_time(3)*temp_resolution,'LineWidth',2); colorbar; caxis([-2 2]); yline(size(amap_vlpfc,1),'LineWidth',2)

            figure; imagesc(Spike_count_raster(session_with_buffer,:)'); xline(block_times.start_time(2)*temp_resolution,'LineWidth',2); xline(block_times.start_time(3)*temp_resolution,'LineWidth',2); colorbar; caxis([-2 2])
            figure; imagesc(behavior_labels(session_with_buffer,:)'); xline(block_times.start_time(2)*temp_resolution,'LineWidth',4); xline(block_times.start_time(3)*temp_resolution,'LineWidth',4);colormap(Cmap); colorbar

            y=mean(Spike_count_raster(session_with_buffer,:)');
            y_std=std(Spike_count_raster(session_with_buffer,:)');
            figure; hold on
            upper_lim=y+y_std;
            lower_lim=y-y_std;
            p = fill([1:length(y) length(y):-1:1],[upper_lim flip(lower_lim)],'red');
            p.FaceColor = [0.9 0.7 0.12];
            p.FaceAlpha = 0.3;
            p.EdgeColor = 'none';
            plot(y,'Color',[0.9 0.7 0.12],'LineWidth',2)
            ylim([-2 4])

        end

        %% Pre-process behavior

        %Simplify behavioral categories
        behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

        %Lump all aggressive interactions together
        behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
        behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");

        %Lump all travel together
        behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
        behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

        
    %% Plot firing rate profile across macaques behavioral repertoire

    %Set parameters
        unqLabels = [7,8,5,1];%1:max(behavior_labels)-1; %Get unique behavior labels (exclude rest)
        n_neurons(s) = size(Spike_rasters,1); %Get number of neurons
        n_behav = length(unqLabels); %Get number of unique behavior labels
   
    n= 138;%randsample(unit_count(3),1);%42;138;218;222
 
    %Plot distribution of firing rate across behaviors for example neuron n
    figure; hold on; set(gcf,'Position',[150 250 1500 300]); i=1;
    for b = unqLabels
        subplot(1, length(unqLabels), i); hold on
        idx = find(behavior_labels == b); %get idx where behavior b occurred
        idx_rest = find(behavior_labels == length(behav_categ));

        %plot firing rate distribution for example unit
        histogram(Spike_rasters(n, idx),'BinWidth',1,'Normalization','pdf', 'FaceColor',Cmap(b,:))
        histogram(Spike_rasters(n, idx_rest),'BinWidth',1,'Normalization','pdf', 'FaceColor',[0.5 0.5 0.5])
        legend({behav_categ(b),'Rest'})
        title([behav_categ(b)])

        set(gca,'FontSize',15);
        xlabel('Firing rate (Hz)'); ylabel('Normalized Frequency');
        ylim([0 0.2])
        i=i+1;
    end
    %saveas(gcf, [savePath '/ExampleUnitSelectivity.eps'])


    %Plot smoothed firing rate over the length of the session and highlight
    %behavior tipes

    %Get behavior indices
    idx_groom = find(behavior_labels == 7); 
    idx_getgroom = find(behavior_labels == 8); 
    idx_selfgroom = find(behavior_labels == 24);
    idx_forage = find(behavior_labels == 5);
    idx_agg = find(behavior_labels == 1);
    idx_travel = find(behavior_labels == 2 | behavior_labels == 12 |behavior_labels == 18);

    unit=n;
    figure; hold on; set(gcf,'Position',[150 250 1500 500]);
    xline([block_times.end_time_round(1), block_times.end_time_round(2)], "-",["Block 1 end", "Block 2 end"], "LineWidth",2);
    session_with_buffer = [120:session_length-120];
    
    %Full spike rate trace
    plot(1:session_length, Spike_rasters(unit,1:session_length));

    %Rest
    to_plot=nan(1,session_length); to_plot(idx_rest)=Spike_rasters(unit,idx_rest);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0.5 0.5 0.5])

    %Getting groomed
    to_plot=nan(1,session_length); to_plot(idx_getgroom)=Spike_rasters(unit,idx_getgroom);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",'b')

    %Grooming
    to_plot=nan(1,session_length); to_plot(idx_groom)=Spike_rasters(unit,idx_groom);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",'c')

    %Self-groom
    to_plot=nan(1,session_length); to_plot(idx_selfgroom)=Spike_rasters(unit,idx_selfgroom);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0.8 0 0.8])

    %Foraging
    to_plot=nan(1,session_length); to_plot(idx_forage)=Spike_rasters(unit,idx_forage);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0 0.7 0])

    %Aggression
    to_plot=nan(1,session_length); to_plot(idx_agg)=Spike_rasters(unit,idx_agg);
    plot(1:session_length, to_plot, "LineWidth",2, "Color","r")

    %Travel
    to_plot=nan(1,session_length); to_plot(idx_travel)=Spike_rasters(unit,idx_travel);
    plot(1:session_length, to_plot, "LineWidth",2, "Color","y")

    legend({'','','Smooth spike trace', 'Rest','Getting groomed','Grooming','Self-groom','Foraging','Aggression','Travel','Scratch'},...
        'Location','bestoutside')

    set(gca,'FontSize',15);
    ylabel('Firing rate (Hz)'); xlabel('Time (s)');
    title(['Example Unit from ' brain_label(n)])
    %saveas(gcf, [savePath '/ExampleUnitLabeledSmoothTrace.eps'])
    
    close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Old code
%     %Plot example spike trains
%     %figure; hold on
%     spiketrains = 50;
%     length_bout = 10;
%     %xlim([0 length_bout]); ylim([0 spiketrains])
%     spike_times_behav = []; spike_times_rest = []; trials_behav = []; trials_rest = [];
%     for st=1:spiketrains
%         Hz_behav =Spike_rasters(n,randsample(idx,length_bout));
%         Hz_rest =Spike_rasters(n,randsample(idx_rest,length_bout));
% 
%         total_spikes_behav = 0; total_spikes_rest = 0;
%         for lb = 1:length_bout
%             ticks_behav = rand(Hz_behav(lb),1);
%             spike_times_behav = [spike_times_behav; lb-1+ticks_behav];
%             total_spikes_behav = total_spikes_behav + Hz_behav(lb);
% 
%             ticks_rest = rand(Hz_rest(lb),1);
%             spike_times_rest=[spike_times_rest; lb-1+ticks_rest];
%             total_spikes_rest = total_spikes_rest+Hz_rest(lb);
%         end
% 
%         trials_rest = [trials_rest, ones(1,total_spikes_rest)*st];
%         trials_behav = [trials_behav, ones(1,total_spikes_behav)*st];
%     end
%     figure; subplot(2,1,1); spikeRasterPlot(seconds(spike_times_behav), trials_behav,'ColorOrder',[0 0 1])
%     subplot(2,1,2); spikeRasterPlot(seconds(spike_times_rest), trials_rest,'ColorOrder',[1 0 0])