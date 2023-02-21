%% Log_PopulationSparseness_batch
%  This script computes the mean z-scored firing rate of individual neuron under different
%  behavioral conditions. 
%  C. Testard Jan 2023

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
agg_precedence =1;


%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/Subject_behav'];
    

    %% Load data
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
            log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma, agg_precedence);
    end

    session_length = size(Spike_rasters,2); % get session length

    %Extract behavior labels
    behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").
 

    n_neurons = size(Spike_rasters,1); %Get number of neurons

    %Estimate "baseline" neural firing distribution.
    idx_rest = find(behavior_labels==length(behav_categ));
    baseline = Spike_rasters(:,idx_rest);
    baseline_std = std(baseline');
    baseline_mean = mean(baseline');

    %Get population sparsness
    behav_firing = Spike_rasters(:,behavior_labels~=length(behav_categ));
    for ti = 1:size(behav_firing,2)
        active(:,ti)=(behav_firing(:,ti)>(baseline_mean+baseline_std)' | behav_firing(:,ti)<(baseline_mean-baseline_std)');
    end
    
    num_units_active{s} = sum(active);
    prop_units_active{s} = sum(active)./n_neurons;

end

