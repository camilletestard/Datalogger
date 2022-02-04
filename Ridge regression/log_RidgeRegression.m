%% Log_RidgeRegression
% Run a ridge regression (code from Musall et. al)


%% Load data

%Set path
is_mac = 0;
if is_mac
    cd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
else
    cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
end
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)

if is_mac
    cd('~/Dropbox (Penn)/Datalogger/Results/')
else
    cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Results/')
end
savePath = uigetdir('', 'Please select the result directory');

clearvars -except savePath filePath is_mac

%Set temporal resolution, channel
temp = 1; temp_resolution = 1; chan = 1; channel_flag = "all";
%for temp_resolution = [1/30, 1/20,  1/10, 1/5, 1/2, 1, 2, 5, 10]
%1 for second resolution, 10 for 100msec resolution, 100 for 10msec resolution, 1000 for msec resolution. etc.
%0.1 for 10sec resolution, 1/5 for 5sec resolution

%Set channels: 'TEO', 'vlPFC' or 'all'
%for channel_flag = ["vlPFC", "TEO", "all"]

%Get data with specified temporal resolution and channels
[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac);
%filePath is the experimental data path
%Temp_resolution is the temporal resolution at which we would like to
%analyze the data
%Channel_flag specifies with channels to include: only TEO array, only
%vlPFC array or all channels
%is_mac is whether a mac or a pc is being used
disp('Data Loaded')
Spike_count_raster = Spike_rasters';

%% Preprocessing, extract and format behavioral events

%Extract behavior info:
behavior_labels_subject = cell2mat({labels{:,3}}'); %Extract unique behavior info
behavior_labels_partner = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
block_labels = cell2mat({labels{:,10}}'); %Extract block info

%Compute the %of timepoints with unidentified behavior
%Important note to keep in mind: I still think that the timepoints where
%there were no behaviors assigned should be excluded...
no_behavs_timepoints = length(find(behavior_labels_subject == 28))/length(behavior_labels_subject);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
fprintf('Percent timepoints without a behavior assigned: %s \n', num2str(no_behavs_timepoints));
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%Create subject behavior regressors:
subject_behav_reg = zeros(size(Spike_count_raster,1), length(behav_categ)-1); %initialize
for b = 1:length(behav_categ)-1
    idx = find(behavior_labels_subject == b);
    subject_behav_reg(idx,b) = 1;
end

%Create partner behavior regressors:
partner_behav_reg = zeros(size(Spike_count_raster,1), length(behav_categ)-1); %initialize
for b = 1:length(behav_categ)-1
    idx = find(behavior_labels_partner == b);
    partner_behav_reg(idx,b) = 1;
end

%Create block regressor:
block_reg = zeros(size(Spike_count_raster,1), length(unique(block_labels))); %initialize
for b = 1:length(unique(block_labels))
    idx = find(block_labels == b);
    block_reg(idx,b) = 1;
end

% %Create behavior change regressor (exclude for now because we do not seem
% to be able to detect behavior shifts irrespective of the shift)
% subject_behav_change_reg = zeros(size(Spike_count_raster,1), 1); %initialize
% subject_behav_change_reg(find(diff(behavior_labels_subject)~=0)+1) = 1;

%Combine behavioral events
Behavior_events=[subject_behav_reg, partner_behav_reg, block_reg];

disp('Preprocessing - adding events done')

%% Set up event types

%Get variable names
BehaviorEventNames = [behav_categ(1:end-1), append('partner.',behav_categ(1:end-1)),string(block_times{:,"Behavior"})'];

%Set up event type here.
allEventTypes=ones(size(BehaviorEventNames))*3;
% For now I will assign all behavioral events as event type 3
%Event Type 1 = Whole trial
%Event Type 2 = from stimulus onset to the rest of the trial
%Event Type 3 = from before movement onset to after movement onset (pre & post time)
allEventsInfo = [BehaviorEventNames; allEventTypes]';

disp('Preprocessing - Assigning event types done')

%% Preprocessing - Selecting Analog Tracking
%For now we don't need to do anything here...


%% Setup Design Matrix - Regressor labels

regLabels = BehaviorEventNames';
regLabels(length(allEventsInfo)+1) = {'Motion Energy L'};
regLabels(length(allEventsInfo)+2) = {'Motion Energy R'};


disp('Setup Design Matrix - Regressor labels done')

%% Preprocessing - Grouping regressors for later cross-validation
%For now only include behavior and motion energy

regGroups = {'BehavioralEvents' 'Movements'};
regGroups{2,1} = BehaviorEventNames';
regGroups{2,2} = {'Motion Energy L' 'Motion Energy R'} ;

disp('Preprocessing - Grouping regressors for later cross-validation done')

%% Setup Design Matrix - Behavioral Events


  %Creates task regressors with the time varying kernels as described in Churchland
    [behavR, behavIdx] = makeDesignMatrix_ST_CT_RWD(temptaskEvents, taskEventType, opts);
    
%end
