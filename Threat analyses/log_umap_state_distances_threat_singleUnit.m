%% log_umap_state_distances_threat
% This script applies UMAP to the data and computes the distance between
% threat data points to a baseline state (centre of mass for rest).

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
%session_range_no_partner=[1:3,11:13];
%session_range_no_partner=[4:6,15:16,18];
session_range_with_partner=[1:6,11:13,15:16,18];

%Set parameters
with_partner =0;
temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=1; %lump similar behavioral categories together
threat_precedence=0; % 1: aggression takes precedence; 0: Threat to partner and subject states take precedence
time_before_threat = 10*temp_resolution;
time_after_threat = 10*temp_resolution;
exclude_sq=0;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1; chan=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];


    %% Get data with specified temporal resolution and channels
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label, behavior_log]= ...
            log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);
    end

    disp('Data Loaded')

    %Raw data
    Spike_count_raster = Spike_rasters';%zscore(Spike_rasters');

    %% Select behaviors to visualize

    %Extract behavior labels and frequency
    behavior_labels = cell2mat({labels{:,3}}');
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

    %Lump all travel together
    behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
    behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

    %Extract block labels
    block_labels = cell2mat({labels{:,13}}');


    %Only consider indices with behavior of interest
    Spike_count_raster_final = Spike_count_raster;%Only keep timepoints where the behaviors of interest occur in spiking data
    behavior_labels_final = behavior_labels;%Same as above but in behavior labels
    block_labels_final =  block_labels;


    %% Get center of mass of rest epochs.

    idx_rest=find(behavior_labels==length(behav_categ)); block_rest=block_labels(behavior_labels==length(behav_categ));
    idx_rest_paired=idx_rest(block_rest==1); idx_rest_alone=idx_rest(block_rest==0);

    %get equal representation of rest during paired and alone blocks.
    idx_equalBlocks=[randsample(idx_rest_paired,min(length(idx_rest_paired),length(idx_rest_alone))); randsample(idx_rest_alone,min(length(idx_rest_paired),length(idx_rest_alone)))];

    rest_com = mean(Spike_count_raster(idx_rest,:));

 
    %% Find threat to subject onset times

    threat_to_subject_onset{s} = behavior_log.start_time_round(find(strcmp(behavior_log.Behavior, "HIS")));

    for event = 1:length(threat_to_subject_onset{s})

        threat_idx = threat_to_subject_onset{s}(event)+1 : threat_to_subject_onset{s}(event)+time_after_threat;
        baseline_idx = threat_to_subject_onset{s}(event)-time_before_threat+1 : threat_to_subject_onset{s}(event);

        block_threat =  block_labels_final(threat_idx);

        for unit=1:size(Spike_count_raster_final,2)
            threat_response(unit, event)=mean(Spike_count_raster_final(threat_idx,unit));
            baseline_activity(unit) = rest_com(unit);%mean(Spike_count_raster_final(baseline_idx,unit));
            responseIdx{s}(unit, event) = abs((threat_response(unit, event)-baseline_activity(unit))/baseline_activity(unit));
        end

        block{s}(event)=unique(block_threat);

    end

    meanResponse{s}=mean(responseIdx{s});
    meanResponse_paired{s}=mean(responseIdx{s}(:,block{s}==1),2);
    meanResponse_alone{s}=mean(responseIdx{s}(:,block{s}==0),2);

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(s)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end %end of session for loop

cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/UMAP_results']);
% save('NeuralDistancesThreat.mat','distance_to_baseline_umap_paired','distance_to_baseline_umap_alone',...
%     'distance_to_baseline_pca_paired','distance_to_baseline_pca_alone','time_before_threat','time_after_threat')
% load('NeuralDistancesThreat.mat')

limits=[0 5];
paired=cell2mat(cat(1,meanResponse_paired(a_sessions))');
alone=cell2mat(cat(1,meanResponse_alone(a_sessions))');
figure; hold on
scatter(paired,alone,'filled','r')  
xlabel('Absolute response index when paired')
ylabel('Absolute response index when alone')
xlim(limits); ylim(limits)

paired=cell2mat(cat(1,meanResponse_paired(h_sessions))');
alone=cell2mat(cat(1,meanResponse_alone(h_sessions))');
scatter(paired,alone,'filled','b')  
xlabel('Absolute response index when paired')
ylabel('Absolute response index when alone')
xlim(limits); ylim(limits)
plot(limits,limits)
legend({'Amos','Hooke','Diagonal'})