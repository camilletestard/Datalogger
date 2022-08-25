%% Log_transition_prob
%This script finds behavioral transtions, computes a transition matrix
%and plots a transition probability graph.
%Testard C. Feb 2022

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,15:16];
session_range_with_partner=[1:3,11:13];

%Set parameters
with_partner =0;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)


%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

behav_rest=nan(1,length(sessions)); % rest
behav_single=nan(1,length(sessions)); % single behavior recorded
behav_coWithProx=nan(1,length(sessions)); % behvior co-occurs with proximity or RR
behav_cooccur=nan(1,length(sessions)); % two key behaviors co-occurring
all_seconds=nan(1,length(sessions));

s=1;
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
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    end

    disp('Data Loaded')

    %Get probability of co-occrrence
    co_occurrence = cell2mat({labels{:,5}}');
    %Important notes:
    %0: no behavior recorded (rest)
    %1: only one behavior at a time
    %2: Behavior co-occurs with proximity or RR
    %3: Proximity and RR co-occur
    %4: Co-occur with groom present
    %5: Co-occur with other monkeys vocalize
    %6: Two 'key' behaviors co-occur

    behav_rest(s) = length(find(co_occurrence==0)); % rest
    behav_single(s) = length(find(co_occurrence==1)); % single behavior recorded
    behav_coWithProx(s) = length(find(co_occurrence==1|co_occurrence==2)); % behvior co-occurs with proximity or RR
    behav_cooccur(s) = length(find(co_occurrence>=4)); % two key behaviors co-occurring
    all_seconds(s) = size(labels,1);

end %end of session loop

nansum(behav_rest)/nansum(all_seconds)
nansum(behav_cooccur)/nansum(all_seconds)

