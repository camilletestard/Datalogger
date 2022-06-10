%% Ron_DoSmth
% Have fun!

%Set session list
home = '~'; % set home directory
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);

%Set parameters
is_mac = 0; %For loading the data
with_partner =1;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: Noise cluster (NC) is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well-isolated units

for s = 1 %For now focus on session Amos_2021-07-29

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/partner_vs_subject'];

    %% Load data

    %Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
    % Output:
    %       1. Spike_rasters: n_neuron X time in seconds (if temp_resolution = 1)
    %       2. labels: time X 11 - check the function description for
    %       details. Essentially subject behavior is column 3.
    %       3. labels_partner: time X 11. Partner behavior is column 3.
    %       4. behav_categ: string vector; behavioral category do the numbers in
    %       labels correspond to.
    %       5. block_times: specifies the block times and ID
    %       6. monkey: subjecet monkey ID
    %       7. reciprocal_set: set of behaviors that are deeemed "reciprocal"
    %       8. social_set: set of behaviors that are deeemed "social"
    %       9. ME_final: Motion energy values
    %       10. unit_count: number of units in the session
    %       11. groom_labels_all: labels for grooming. Check function
    %       description for details

    session_length = size(Spike_rasters,2); % get session length
    Spike_count_raster = Spike_rasters';

    %Extract behavior labels for subject and partner
    behavior_labels_subject_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels_partner_init = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
    
    % Set proximity as rest
    behavior_labels_subject_init(behavior_labels_subject_init==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").
    behavior_labels_partner_init(behavior_labels_partner_init==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

end