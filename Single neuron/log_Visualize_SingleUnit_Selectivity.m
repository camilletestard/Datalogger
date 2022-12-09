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
agg_precedence =1;

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
        %Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma, agg_precedence );
        end

        disp('Data Loaded')

        Spike_count_raster = zscore(Spike_rasters');

        %Extract behavior labels
        behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure

        %Set colormap
        Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0.3 0.7 1];[0 0.7 0];[1 0 1];[0 1 1];...
            [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0.7 0 1];...
            [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0.8 0 0];[1 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
            [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];


        session_with_buffer = [250*temp_resolution:size(Spike_count_raster,1)-250*temp_resolution];

        figure; imagesc(Spike_count_raster(session_with_buffer,:)'); xline(block_times.start_time(2),'LineWidth',2); xline(block_times.start_time(3),'LineWidth',2); colorbar; caxis([-2 2])
        figure; imagesc(behavior_labels(session_with_buffer,:)'); xline(block_times.start_time(2),'LineWidth',4); xline(block_times.start_time(3),'LineWidth',4);colormap(Cmap); colorbar

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

%         %Hierarchical clustering for visualization
%         Y=squareform(pdist(Spike_count_raster(session_with_buffer,:)','cityblock'));
%         Z=linkage(Y,"average");
%         dendrogram(Z)
%         c = cophenet(Z,Y)
%         T = cluster(Z,"cutoff",1.152); figure; hist(T)
%         [~, idx_order]=sort(T);
%         idx_order = [find(strcmp(brain_label,'TEO')), find(strcmp(brain_label,'vlPFC'))];
%         figure; imagesc(Spike_count_raster(session_with_buffer,idx_order)'); colorbar; caxis([-3 3])

        %% Pre-process behavior
  
        %Simplify behavioral categories
        behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

        %Lump all aggressive interactions together
    behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
    behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");

    %Lump all travel together
    behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
    behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");


    %% Set parameters
    unqLabels = [1,5,7,8,24,18];%1:max(behavior_labels)-1; %Get unique behavior labels (exclude rest)
    n_neurons(s) = size(Spike_rasters,1); %Get number of neurons
    n_behav = length(unqLabels); %Get number of unique behavior labels


   

    %% Plot firing rate profile across macaques behavioral repertoire

    %Set color scheme
   color_scheme = {'r','','','',[0 0.7 0],'','c','b','','','','','','','','','','y','','','','',[0.8 0.6 0],[0.8 0 0.8]};
%     Cmap = [[0 0 0];[1 0.4 0.1];[0 0 0];[0 0.6 0.8];[0 0.7 0];[1 0 0.65];[0 1 1];...
%             [0 0 1];[0.5 0 0];[1 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];...
%             [0 0 0];[0 0 0];[0.9 0.7 0.12];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0.8 0.6 0];[0.8 0 0.8];...
%             [0 0 0];[0 0 0];[0.2 0.2 1]];

    n= 138;%randsample(unit_count(3),1);%42;138;218;222
    cohend_per_neuron = cohend(n,:);
    [~, orderIdx] = sort(cohend_per_neuron, 'ascend');

    %Plot distribution of firing rate across behaviors for example neuron n
    figure; hold on; set(gcf,'Position',[150 250 1500 300]); i=1;
    for b = unqLabels(orderIdx)
        subplot(1, length(orderIdx), i); hold on
        idx = find(behavior_labels == b); %get idx where behavior b occurred

        %plot firing rate distribution for example unit
        histogram(Spike_rasters(n, idx),'BinWidth',1,'Normalization','pdf', 'FaceColor',color_scheme{b})
        histogram(Spike_rasters(n, idx_rest),'BinWidth',1,'Normalization','pdf', 'FaceColor',[0.5 0.5 0.5])
        legend({behav_categ(b),'Rest'})
        title([behav_categ(b)])

        set(gca,'FontSize',15);
        xlabel('Firing rate (Hz)'); ylabel('Normalized Frequency');
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