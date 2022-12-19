%% Log_Visualize_SingleUnit_Selectivity
%  This script computes firing rate of individual neuron under different
%  behavioral conditions. Then, it computes a cohen's d (or effect size)
%  difference between the distribution of firing rates during behavior X
%  with a baseline firing rate.

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range=[1,2,11,12];

%Set parameters
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/'];

    %clearvars -except savePath filePath is_mac s home

    %% Load data

    %Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);

    session_length = size(Spike_rasters,2); % get session length

    %Extract behavior labels
    behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

    %Set parameters
    unqLabels = 1:max(behavior_labels)-1; %Get unique behavior labels (exclude rest)
    n_neurons = size(Spike_rasters,1); %Get number of neurons
    n_behav = length(unqLabels); %Get number of unique behavior labels

    %Estimate "baseline" neural firing distribution.
    idx_rest=find(behavior_labels==length(behav_categ));%Get idx of "rest" epochs.
    mean_baseline = mean(Spike_rasters(:,idx_rest),2);
    std_baseline = std(Spike_rasters(:,idx_rest),0,2);

    %Check visually that baseline is taken from epochs throughout the session
    y=zeros(1, session_length); y(idx_rest)=1;
    figure; plot(1:session_length, y); ylim([-0.5, 1.5])
    yticks([0 1]); yticklabels(["Behavior", "Rest/baseline"])
    xlabel('Time in s'); title('Baseline epochs')
    set(gca,'FontSize',15);
%     saveas(gcf, [savePath '/Baseline_epochs.png']); close all

    %Standerdize Unit rasters
    Spike_raster_relative2baseline = Spike_rasters./mean_baseline; %Divide by baseline firing
    Spike_raster_smooth = smoothdata(Spike_raster_relative2baseline,'gaussian',20);

    idx_groom = find(behavior_labels == 7); idx_getgroom = find(behavior_labels == 8);
    idx_agg = find(behavior_labels == 9 | behavior_labels == 10 |behavior_labels == 1);
    unit=209;
    figure; hold on; set(gcf,'Position',[150 250 1500 500]);
    xline([block_times.end_time_round(1), block_times.end_time_round(2)], "-",["Block 1 end", "Block 2 end"], "LineWidth",2);
    plot(1:session_length, Spike_rasters(unit,1:session_length)); 
    plot(idx_groom, Spike_rasters(unit,idx_groom), "LineWidth",2)
    plot(idx_getgroom, Spike_rasters(unit,idx_getgroom), "LineWidth",2)
    plot(idx_agg, Spike_rasters(unit,idx_agg), "LineWidth",2)
    plot(idx_rest, Spike_rasters(unit,idx_rest), "LineWidth",2)

    yline(1)


end
