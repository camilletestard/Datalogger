%% Log GenerateLabelsPerSecond
% Format the raw data to have two elements:
% 1. Neural data matrix size [Time (in sec) x #neurons]
% 2. Label vector which describes the behavior at time t [Time (in sec) x 1]
% Camille Testard - Sept. 2021

%% Load data
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)

session = filePath(end-9:end);
monkey = filePath(end-14:end-10);

behavior_log = readtable(['EVENTLOG_restructured_',monkey,session,'.csv']);% Behavioral data
load('Neural_data.mat') % Neural data; array1 is in TEO and array2 is in vlPFC
session_length = size(Unit_rasters,2);
Spike_count_raster = Unit_rasters';

%Preprocessing: round times in behavioral log
behavior_log{:,'start_time_round'}=round(behavior_log{:,'start_time'});
behavior_log{:,'end_time_round'}=round(behavior_log{:,'end_time'});
behavior_log{:,'duration_round'}=behavior_log{:,'end_time_round'}-behavior_log{:,'start_time_round'};

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
labels = cell(session_length,3); %initialize dataframe
for s = 1:session_length %for all secs in a session
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
            else %Otherwise choose the first behavior
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
    save([filePath '/Labels_per_sec.mat'],'labels','behav_categ');