%% Log GenerateLabelsPerSecond
% Format the raw data to have two elements:
% 1. Neural data matrix size [Time (in sec) x #neurons]
% 2. Label vector which describes the behavior at time t [Time (in sec) x 1]
% Camille Testard - Sept. 2021

%% Load data
is_mac = 1;

if is_mac
    cd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
else
    cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
end
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)

session = filePath(end-9:end);
monkey = filePath(end-14:end-10);

behavior_log = readtable(['EVENTLOG_restructured_',monkey,session,'.csv']);% Load behavioral data
load(['Neural_data_' session '.mat']) % Load neural data; array1 is in TEO and array2 is in vlPFC
length_recording = size(Unit_rasters,2); %Unit rasters in second resolution

%% Preprocessing: behavioral log and neural data at specified resolution

%Set temporal resolution
temp_resolution = 10; %1 for second resolution, 10 for 100msec resolution, 100 for 10msec resolution, 1000 for msec resolution. etc.
                      %0.1 for 10sec resolution, 1/5 for 5sec resolution

%Behavioral log
behavior_log{:,'start_time_round'}=round(behavior_log{:,'start_time'}*temp_resolution);
behavior_log{:,'end_time_round'}=round(behavior_log{:,'end_time'}*temp_resolution);
behavior_log{:,'duration_round'}=behavior_log{:,'end_time_round'}*temp_resolution-behavior_log{:,'start_time_round'}*temp_resolution;

%Neural data
Chan_name = fieldnames(SpikeData); %Identify channel names
C = regexp(Chan_name,'\d*','Match');
C_char = cellfun(@char, C{:}, 'UniformOutput', false);
Chan_num = str2num(C_char{1, 1});

% % % %Separate channels by array
% % % array1_chan = [1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,...
% % %     42,45,46,49,50,53,54,57,58,61,62,65,66,69,70,73,74,77,78,81,82,85,86,...
% % %     89,90,93,94,97,98,101,102,105,106,109,110,113,114,117,118,121,122,125,126];
% % % 
% % % array2_chan = [3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39,40,...
% % %     43,44,47,48,51,52,55,56,59,60,63,64,67,68,71,72,75,76,79,80,83,84,...
% % %     87,88,91,92,95,96,99,100,103,104,107,108,111,112,115,116,119,120,123,124,127,128];
% % % 
% % % chan_idx_array1 = find(ismember(Chan_num,array1_chan))';
% % % chan_idx_array2 = find(ismember(Chan_num,array2_chan))';

unit=1;
for i = 1:length(fields(SpikeData)) %For all channels
    
    if ~isempty(SpikeData.(Chan_name{i})) %If there are sorted units on this channel
        for j = 1:length(SpikeData.(Chan_name{i})) %For all units
            
            Spike_rasters(unit,:) = zeros(1,round(length_recording*temp_resolution)+1); %Fill the line with zeros to initiate raster for that trial
            ticks = round(SpikeData.(Chan_name{i}){j}*temp_resolution);
            Spike_counts = hist(ticks, round(length_recording*temp_resolution)+1);
            Spike_rasters(unit, :) = Spike_counts; %Fill in spikes in the raster
            clear ticks Spike_counts
            
            unit = unit+1;
        end
    end
    
end

length_recording = size(Spike_rasters,2);

%% Get behavior label vector for each second

%Create event intervals:
start_times = behavior_log{:,'start_time_round'};
end_times = behavior_log{:,'end_time_round'};
Intervals = [start_times end_times];

%Create behavior key
behav_categ = unique(behavior_log{:,'Behavior'}); %Get all the unique behaviors
behav_categ{20}='Rest'; %Add rest as a behavior (no defined behavior ongoing)
double_behav_set = [find(matches(behav_categ,'Proximity')), find(matches(behav_categ,"RR"))];%, find(matches(behav_categ,"HIS")), find(matches(behav_categ,"HIP"))]; %For behaviors that often co-occur with other behaviors
omv = find(matches(behav_categ,'Other monkeys vocalize'));

%Create behavior label vector (label every second of the session)
% This cell matrix will have three columns. The first column is the full
% name of the behavior label
labels = cell(length_recording,3); %initialize dataframe
for s = 1:length_recording %for all secs in a session
    % this finds the index of he rows(2) that have x in between
    idx = find(s > Intervals(:,1) & s < Intervals(:,2)); %find if this second belong to any interval
    %IMPORTANT note: interval exclude boundaries as is.
    if ~isempty(idx) %if yes
        labels{s,1} = behavior_log{idx,'Behavior'}; %add behavior id (in plain english)
        labels{s,2} = find(matches(behav_categ,labels{s,1})); %add behavior number 
        if length(labels{s,2})>1 %If one behavior co-occurs with proximity or RR
            labels{s,4} = 'co-occur';
            labels{s,3} = setdiff(labels{s,2}, double_behav_set); % only consider the other behavior (it that takes precedence over proximity and RR)
        else %If only one behavior happens in that sec
            labels{s,3} = labels{s,2};
            labels{s,4} = 'single';
        end
        if length(labels{s,3})~=1 %If two behaviors are co-occurring which do not include proximity or RR
            if any(labels{s,3}==omv) % if one of the behavior includes other monkey vocalize
                labels{s,3}=omv; %Keep OMV
            else %Otherwise choose the second behavior
%                 error('More than one behavior simultansouly')
%                 return
                labels{s,3}= labels{s,3}(2); %2nd behavior (HIP/HIS take precedence over aggression)
            end
        end
    else %if not
        labels{s,1} = NaN; labels{s,2} = 20; labels{s,3} = 20; %Set behavior category to "NaN" and 
    end
end

%length(find(strfind(labels{:,4},'co-occur')))

    % %% Save variables
    save([filePath '/Data_' num2str(1000/temp_resolution) 'msec_res.mat'],'labels','behav_categ','Spike_rasters');