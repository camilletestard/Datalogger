%% Loading data from example session
% C. Testard August 2022

%Set path
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:); s=1;

%Set parameters
with_partner =1; %Include partner behavior data 1:yes; 0:no
temp = 1; temp_resolution = 1;
channel_flag = "all"; %"all", "vlPFC" or "TEO"
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units


%Set path
filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];


%% Get data with specified temporal resolution and channels
%Note: Thescription of inpout and output arguments detailed in te
%function's documentation.
[Spike_rasters, SpikeData, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= log_GenerateDataToRes_Example_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
disp('Data Loaded')

