function [Spike_rasters, labels, behav_categ, block_times] = log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, monkey)
        %Log GenerateDataToRes_function
        % This function formats the raw data to have two elements:
        % 1. Neural data matrix, size [Time (to chosen resolution) x #neurons]
        % 2. Label vector which describes the behavior at time t [Time (to chosen resolution) x 4]
            %1st column includes all behaviors in "plain english"
            %2nd column behavior number code
            %3rd column unique behavior code (when two occur
            %simultaneously, we chose one, see below for details)
            %4th column whether behavior happens in isolation or co-occurs
            %with another.
            %5th column indicates the block in which we are
            %(Paired,monkey1; Paired monkey2 or Alone)
        %filePath is the experimental data path
        %Temp_resolution is the temporal resolution at which we would like to
        %analyze the dat
        %Channel_flag specifies with channels to include: only TEO array, only
        %vlPFC array or all channels

% Camille Testard - Nov. 2021

%% Load data
cd(filePath)

session = filePath(end-9:end);
monkey = filePath(end-14:end-10);

behavior_log = readtable(['EVENTLOG_restructured_',monkey,session,'.csv']);% Load behavioral data
load(['Neural_data_' session '.mat']) % Load neural data; array1 is in TEO and array2 is in vlPFC
length_recording = size(Unit_rasters,2); %Unit rasters in second resolution

%% Preprocessing: behavioral log and neural data at specified resolution

%Behavioral log
behavior_log{:,'start_time_round'}=round(behavior_log{:,'start_time'}*temp_resolution);
behavior_log{:,'end_time_round'}=round(behavior_log{:,'end_time'}*temp_resolution);
behavior_log{:,'duration_round'}=behavior_log{:,'end_time_round'}-behavior_log{:,'start_time_round'};

%Eliminate behaviors that do not meet the minimum length
%Note that this will be an issue only for time resolution >1sec
min_length = 1/temp_resolution;
idx = find(behavior_log{:,'duration_s'}<min_length);
behavior_log{idx,'start_time_round'} = 0; behavior_log{idx,'end_time_round'} = 0;

%Get block times (at the end of the EventLog_Restructured)
block_times = behavior_log(end-2:end,:);
behavior_log(end-2:end,:) = [];

%Neural data
Chan_name = fieldnames(SpikeData); %Identify channel names
C = regexp(Chan_name,'\d*','Match');
C_char = cellfun(@char, C{:}, 'UniformOutput', false);
Chan_num = str2num(C_char{1, 1});

%Separate channels by array
if monkey == "Hooke"
array1_chan = [1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,...
    42,45,46,49,50,53,54,57,58,61,62,65,66,69,70,73,74,77,78,81,82,85,86,...
    89,90,93,94,97,98,101,102,105,106,109,110,113,114,117,118,121,122,125,126];

array2_chan = [3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39,40,...
    43,44,47,48,51,52,55,56,59,60,63,64,67,68,71,72,75,76,79,80,83,84,...
    87,88,91,92,95,96,99,100,103,104,107,108,111,112,115,116,119,120,123,124,127,128];
else

end

chan_idx_array1 = find(ismember(Chan_num,array1_chan))';
chan_idx_array2 = find(ismember(Chan_num,array2_chan))';

%Select channels
if strcmp(channel_flag,'TEO')
    channels = chan_idx_array1;
elseif strcmp(channel_flag,'vlPFC')
    channels = chan_idx_array2;
elseif strcmp(channel_flag,'all')
    channels = 1:length(fields(SpikeData)); %all channels
end

%Create spike matrix structure
unit=1;
for i = channels %For all channels
    
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

%% Get behavior label vector for each time bin at specified resolution

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
    if s<=block_times{1,'end_time_round'}
        labels{s,5} = string(block_times{1,'Behavior'});
        labels{s,6} = 1;
    elseif s>block_times{1,'end_time_round'} && s<=block_times{2,'end_time_round'}
        labels{s,5} = string(block_times{2,'Behavior'});
        labels{s,6} = 2;
    elseif s>block_times{2,'end_time_round'}
        labels{s,5} = string(block_times{3,'Behavior'});
        labels{s,6} = 3;
    end
end

behav_categ{7}='Threat to partner';
behav_categ{8}='Threat to subject';
behav_categ{11}='Travel';
behav_categ{13}='Rowdy Room';
behav_categ{14}='Squeeze partner';
behav_categ{15}='Squeeze Subject';

end