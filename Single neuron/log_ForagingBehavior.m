%% Log_SingleNeuronTuning_zscore_batch
%  This script computes the mean z-scored firing rate of individual neuron under different
%  behavioral conditions. 
%  C. Testard July 2022

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
plot_toggle = 0;
select_behav=0;
with_partner = 0;
temp_resolution = 1; %Temporal resolution of firing rate. 1: 1sec; 10:100msec; 0.1: 10sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
min_occurrence =30*temp_resolution;
cohend_cutoff=0.3; p_cutoff=0.01;%Set thresholds
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size in sec (sigma)
null=0;%Set whether we want the null
threat_precedence =0;
exclude_sq=1;

%Initialize session batch variables:
n_behav = 26;
mean_cohend_per_behav = nan(length(sessions), n_behav);
median_cohend_per_behav = nan(length(sessions), n_behav);
std_cohend_per_behav = nan(length(sessions), n_behav);
se_cohend_per_behav = nan(length(sessions), n_behav);
prop_selective_per_behav = nan(length(sessions), n_behav);
num_selective_behav_per_neuron=cell(1,length(sessions));
n_per_behav = nan(length(sessions),n_behav);

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/Subject_behav'];
    

    %% Load data

    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);


    session_length = size(Spike_rasters,2); % get session length

    %Extract behavior labels
    behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure
    
    %Pool these behaviors with rest (effectively excluding them)
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").
    behavior_labels(behavior_labels==find(behav_categ=="Rowdy Room"))=length(behav_categ)+1; %exclude rowdy room for now
    behavior_labels(behavior_labels==find(behav_categ=="Other monkeys vocalize"))=length(behav_categ)+1; %exclude rowdy room for now
    
    %Pool travel, approach and leave
    behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel"); %Consider 'approach' to be 'Travel'.
    behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel"); %Consider 'leave' to be 'Travel'.

    if null
        %Simulate fake labels
        [sim_behav] = GenSimBehavior(behavior_labels,behav_categ, temp_resolution, plot_toggle);
        behavior_labels = sim_behav;
    end

    %% Set parameters
    unqLabels = find(behav_categ=="Foraging");%1:length(behav_categ)-1; %Get unique behavior labels (exclude rest)
    n_neurons(s) = size(Spike_rasters,1); %Get number of neurons
    n_behav = length(unqLabels); %Get number of unique behavior labels

    %Estimate "baseline" neural firing distribution.
    Spike_rasters_zscore = zscore(Spike_rasters,0,2);


    %% Compute cohen's d

    cohend = nan(n_neurons(s),n_behav);
    cohend_shuffle = nan(n_neurons(s),n_behav);
    mean_beh = nan(n_neurons(s), n_behav);
    mean_beh_shuffle = nan(n_neurons(s), n_behav);
    std_beh = nan(n_neurons(s), n_behav);
    std_beh_shuffle = nan(n_neurons(s), n_behav);
    p = nan(n_neurons(s), n_behav);
    p_rand = nan(n_neurons(s), n_behav);

    for n = 1:n_neurons(s)

        for b = 5
            idx = find(behavior_labels == b); %get idx where behavior b occurred
            n_per_behav(s)=length(idx);

            if n_per_behav(s)>min_occurrence

                mean_beh(n)=mean(Spike_rasters_zscore(n, idx),2);
                std_beh(n)=std(Spike_rasters_zscore(n, idx),0,2);
            
            end

        end
    end

    save_meanBehav{s}=mean_beh;
    save_meanBehav_teo{s}=mean_beh(strcmp(brain_label, "TEO"),:);
    save_meanBehav_vlpfc{s}=mean_beh(strcmp(brain_label, "vlPFC"),:);
    
    figure; hold on
    histogram(save_meanBehav_teo{s}, 20)
    histogram(save_meanBehav_vlpfc{s}, 20)
    legend({'TEO','vlPFC'})
    xlabel('Z-scored firing rate')
    ylabel('Neuron count')
    title(sessions(s).name)

    median_diff(s) = median(save_meanBehav_teo{s}) - median(save_meanBehav_vlpfc{s})


end

close all


