%% Running Notes
%Script for testing out function to generate simulated behavior for any
%particular session based on the statistics of that session.

%Currently thinking to "fit" p value for bernoulli for each behavior and
%some term for how long each behavior tends to last.  All "coins"
%for each behavior are flipped and winning behavior begins (set some random
%way to pick between ties).  Then draw a sample for the duration of that
%bout of the behavior.  Once bout is over, flip coins again to select new
%behavior and repeat until the session is over.

%% Load in data

%Parameters for setting Path
is_mac = 0; %For loading the data
is_ron = 1; %For setting path



%Parameters for setting sessions to use and brain regions to look at
S_list = [1]; %List of session to pull; %For now focus on first Amost sessions
BRs = ["TEO", "vlPFC", "all"]; %Channels considered (sets brain region TEO vlPFC or all)

%Parameters for neural data to be used across sessions
with_NC =1; %0: Noise cluster (NC) is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well-isolated units
temp_resolution = 1; %Temporal resolution of firing rate. 1sec %%Note this will change if we introduce smoothing
smooth_fr = 0; %Toggle to smooth the spike rasters (currently not setting up this code as may require working with ms neural data and then downsampling to get back)

s = 1;
br = 3;

if is_ron

    sessions = dir('C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_data'); sessions = sessions(3:end,:);
    filePath = ['C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_data\' sessions(s).name];
    savePath = ['C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_results'];

else


     %Set session list
home = '~'; % set home directory
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/partner_vs_subject'];
end


channel_flag = BRs(br;


[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_Ron(filePath, temp_resolution, channel_flag, is_mac,is_ron, with_NC, isolatedOnly);

%% Figure out any manipulation that need to be done before/in the function given the data that needs to be passed.

