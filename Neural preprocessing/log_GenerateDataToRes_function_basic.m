function [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final, unit_count, groom_labels_all] = log_GenerateDataToRes_function_basic(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly)
%Log GenerateDataToRes_function
% Input data: 
%   1. Behavior of the subject: "EVENTLOG_restructured.csv"
%   2. Behavior of the partner: "EVENTLOG_restructured_partner.csv"
%   3. Spike sorted neural data: "Neural_data_*sessionname*.mat" OR
%   "Neural_data_*sessionname*_isolatedunits.mat"
%   4. Motion energy measure for each top view video: "ME.csv"

% This function formats the raw input data to have the following elements:
% 1. Neural data matrix, size [Time (to chosen resolution) x #neurons]
% 2. Label vector which describes the behavior at time t [Time (to chosen resolution) x 4]
%       1st column includes all behaviors in "plain english"
%       2nd column behavior number code
%       3rd column unique behavior code (when two occur
%       simultaneously, we chose one, see below for details)
%       4th column whether behavior happens in isolation or co-occurs
%       with another.
%       5th column indicates the block in which we are
%       (Paired,monkey1; Paired monkey2 or Alone)
%       6th column gives a corresponding numerical value to the block

% Arguments: 
%   filePath is the experimental data path
%   Temp_resolution is the temporal resolution at which we would like to
%   analyze the dat
%   Channel_flag specifies with channels to include: only TEO array, only
%   vlPFC array or all channels
%   is_mac: specifies whether the code is run on a mac or pc
%   with_NC: specifies whether the "noise cluster" (or the first cell of
%   every channel) is included (1) or not (0). If with_NC=2 then we only
%   inlcude the noise cluster (not the other neurons).
%   isolatedOnly: specifies if only the well isolated units are considered.

% Camille Testard - Nov. 2021

%% Load data
cd(filePath)

if is_mac
    split = '/';
else
    split = '\';
end
split_file_name = strsplit(filePath,split); full_session_name = split_file_name{end}; session_name_split = strsplit(full_session_name,'_');
session = session_name_split{2};
monkey = session_name_split{1};

%Load behavioral data
behavior_log = readtable('EVENTLOG_restructured.csv');% for subject
labels_partner = [];
groom_labels_all = [];% Load grooming context labeling

%Load neural data
if isolatedOnly==1
    load(['Neural_data_' session '_IsolatedUnits.mat'])%only including well isolated units
else
    load(['Neural_data_' session '.mat']) % Load neural data; array1 is in TEO and array2 is in vlPFC
end
length_recording = size(Unit_rasters,2); %Unit rasters in second resolution

% Load motion energy
ME_final = [];%readtable('ME.csv');

cd('~/Dropbox (Penn)/Datalogger/Results/All_sessions')
num_unit_allsessions = readtable('Session_log_num_units.csv');% Load number of unit data
session_idx = find(~cellfun(@isempty,(strfind(num_unit_allsessions.session_name,session))));
unit_count = [num_unit_allsessions.num_units_vlPFC(session_idx), num_unit_allsessions.num_units_TEO(session_idx), num_unit_allsessions.num_units(session_idx)];


%% Preprocessing: behavioral log and neural data at specified resolution

%Round times and get duration for behavioral logs (subject & partner)
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


%% Neural data

Chan_name = fieldnames(SpikeData); %Identify channel names
C = regexp(Chan_name,'\d*','Match');
C_char = cellfun(@char, C{:}, 'UniformOutput', false);
Chan_num = str2num(C_char{1, 1});

%Separate channels by array
%IMPORTANT NOTE: the channel mapping is reversed for each monkey
if strcmp(monkey,'Hooke')
    
    TEO_chan = [1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,...
        42,45,46,49,50,53,54,57,58,61,62,65,66,69,70,73,74,77,78,81,82,85,86,...
        89,90,93,94,97,98,101,102,105,106,109,110,113,114,117,118,121,122,125,126];
    
    vlPFC_chan = [3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39,40,...
        43,44,47,48,51,52,55,56,59,60,63,64,67,68,71,72,75,76,79,80,83,84,...
        87,88,91,92,95,96,99,100,103,104,107,108,111,112,115,116,119,120,123,124,127,128];
    
elseif strcmp(monkey,'Amos')
    
    vlPFC_chan = [1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,...
        42,45,46,49,50,53,54,57,58,61,62,65,66,69,70,73,74,77,78,81,82,85,86,...
        89,90,93,94,97,98,101,102,105,106,109,110,113,114,117,118,121,122,125,126];
    
    TEO_chan = [3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39,40,...
        43,44,47,48,51,52,55,56,59,60,63,64,67,68,71,72,75,76,79,80,83,84,...
        87,88,91,92,95,96,99,100,103,104,107,108,111,112,115,116,119,120,123,124,127,128];
    
end

chan_idx_TEO = find(ismember(Chan_num,TEO_chan))';
chan_idx_vlPFC = find(ismember(Chan_num,vlPFC_chan))';

%Select channels
if strcmp(channel_flag,'TEO')
    channels = chan_idx_TEO;
elseif strcmp(channel_flag,'vlPFC')
    channels = chan_idx_vlPFC;
elseif strcmp(channel_flag,'all')
    channels = 1:length(fields(SpikeData)); %all channels
end

%Create spike matrix structure

if with_NC==0 && isolatedOnly==0 %If don't include noise cluster
    unit=1;
    for i = channels %For all channels
        if length(SpikeData.(Chan_name{i}))>1 %If there are sorted units on this channel
            for j = 2:length(SpikeData.(Chan_name{i})) %For all units except Noise cluster

                Spike_rasters(unit,:) = zeros(1,round(length_recording*temp_resolution)); %Fill the line with zeros to initiate raster for that trial (IMPORTANT NOTE: removed +1)
                ticks = round(SpikeData.(Chan_name{i}){j}*temp_resolution);
                Spike_counts = hist(ticks, round(length_recording*temp_resolution));
                Spike_rasters(unit, :) = Spike_counts; %Fill in spikes in the raster
                clear ticks Spike_counts

                unit = unit+1;
            end
        end
    end
elseif with_NC==2 && isolatedOnly==0 %if ONLY include the noise cluster
    unit=1;
    for i = channels %For all channels
        if ~isempty(SpikeData.(Chan_name{i})) %If there are sorted units on this channel
            for j = 1 %For the first channel only

                Spike_rasters(unit,:) = zeros(1,round(length_recording*temp_resolution)); %Fill the line with zeros to initiate raster for that trial (IMPORTANT NOTE: removed +1)
                ticks = round(SpikeData.(Chan_name{i}){j}*temp_resolution);
                Spike_counts = hist(ticks, round(length_recording*temp_resolution));
                Spike_rasters(unit, :) = Spike_counts; %Fill in spikes in the raster
                clear ticks Spike_counts

                unit = unit+1;
            end
        end
    end
else %include noise cluster (or first channel)
    unit=1;
    for i = channels %For all channels
        if ~isempty(SpikeData.(Chan_name{i})) %If there are sorted units on this channel
            for j = 1:length(SpikeData.(Chan_name{i})) %For all units

                Spike_rasters(unit,:) = zeros(1,round(length_recording*temp_resolution)); %Fill the line with zeros to initiate raster for that trial (IMPORTANT NOTE: removed +1)
                ticks = round(SpikeData.(Chan_name{i}){j}*temp_resolution);
                Spike_counts = hist(ticks, round(length_recording*temp_resolution));
                Spike_rasters(unit, :) = Spike_counts; %Fill in spikes in the raster
                clear ticks Spike_counts

                unit = unit+1;
            end
        end
    end
end


length_recording = size(Spike_rasters,2);

%% Get behavior label vector for each time bin at specified resolution

%Create behavior key
%behav_categ = unique(behavior_log{:,'Behavior'}); %Get all the unique behaviors
%Set it constant across sessions. This makes coding consistent across
%sessions.
behav_categ = ["Aggression","Proximity","Groom Give", "HIP","Foraging", "Vocalization","SS", "Masturbating",...
    "Submission", "Approach","Yawning","Self-groom","HIS","Other monkeys vocalize", "Lip smack",...
    "Groom Receive","Leave","Drinking","SP","Pacing/Travel","Scratch","RR", "Butt sniff","Grm prsnt",...
    "Head Bobbing", "Swinging", "Object Manipulation"];
behav_categ = sort(behav_categ);

behav_categ{length(behav_categ)+1}='Rest'; %Add rest as a behavior (no defined behavior ongoing)

%For behaviors that often co-occur with other behaviors, determine priority
double_behav_set = [find(matches(behav_categ,'Proximity')), find(matches(behav_categ,"RR"))];% When co-occurring, other behaviors will take precedence over these ones %, find(matches(behav_categ,"HIS")), find(matches(behav_categ,"HIP"))];
omv = find(matches(behav_categ,'Other monkeys vocalize')); %When co-occurring with toher behaviors, vocalizations take precedence
reciprocal_set = [find(matches(behav_categ,'Proximity')), find(matches(behav_categ,"Groom Give")), find(matches(behav_categ,"Groom Receive")),...
    find(matches(behav_categ,"SS")), find(matches(behav_categ,"HIP")), find(matches(behav_categ,"HIS")), find(matches(behav_categ,"SP"))];
social_set = [find(matches(behav_categ,'Proximity')), find(matches(behav_categ,"Groom Give")), find(matches(behav_categ,"Groom Receive")),...
    find(matches(behav_categ,"Submission")), find(matches(behav_categ,"Approach")), find(matches(behav_categ,"Leave")), find(matches(behav_categ,"Butt sniff")),...
    find(matches(behav_categ,"Grm prsnt")), find(matches(behav_categ,"Aggression"))];


%% Create behavior label vector (label every window of the session)
% This cell matrix will have six columns. The first column is the full
% name of the behavior label

%Create event intervals:
%For subject monkey
start_times = behavior_log{:,'start_time_round'};
end_times = behavior_log{:,'end_time_round'};
Intervals = [start_times end_times];

%%%%%% Create labels vector for SUBJECT monkey %%%%%%
labels = cell(length_recording,11); %initialize dataframe
for s = 1:length_recording %for all secs in a session
    % this finds the index of the rows(2) that have x in between
    idx = find(s >= Intervals(:,1) & s < Intervals(:,2)); %find if this second belong to any interval
    %IMPORTANT note: interval includes lower bound but excludes upper boundary as is.
    if ~isempty(idx) %if yes
        labels{s,1} = behavior_log{idx,'Behavior'}; %add behavior label in [plain english]
        labels{s,2} = find(matches(behav_categ,labels{s,1})); %add behavior label in [number]
        if length(labels{s,2})>1 %If one behavior co-occurs with other behavior(s)
            if ~isempty(setdiff(labels{s,2}, double_behav_set)) %If behavior co-occurs with proximity or RR
                labels{s,3} = setdiff(labels{s,2}, double_behav_set); % only consider the other behavior (it takes precedence over proximity and RR)
                labels{s,4} = 'co-occur with prox or RR';
                labels{s,5} = 3;
            else %If proximity & RR co-occur
                labels{s,3}=find(matches(behav_categ,'Proximity')); % prioritize proximity
                labels{s,4} = 'co-occur';
                labels{s,5} = 2;
            end
        else %If only one behavior happens in that sec
            labels{s,3} = labels{s,2};
            labels{s,4} = 'single';
            labels{s,5} = 1;
        end
        if length(labels{s,3})~=1 %If two behaviors are co-occurring which do not include proximity or RR
            if any(labels{s,3}==omv) % if one of the behavior includes other monkey vocalize
                labels{s,3}=omv; %Keep OMV
                labels{s,4} = 'co-occur';
                labels{s,5} = 2;
            else %Otherwise just choose the second behavior
                %                 error('More than one behavior simultansouly')
                %                 return
                labels{s,3}= labels{s,3}(2); %2nd behavior (HIP/HIS take precedence over aggression)
                labels{s,4} = 'co-occur';
                labels{s,5} = 2;
            end
        end
    else %if not
        labels{s,1} = NaN; labels{s,2} = length(behav_categ); labels{s,3} = length(behav_categ); labels{s,4} = 'NA'; labels{s,5} = 0;%Set behavior category to "NaN" and label to rest
    end
    %Add behavior information: reciprocal vs non-reciprocal. Reciprocal
    %behavior is the exact reverse of the subject behavior (i.e. we can
    %100% predict the behavior label of the partner based on the behavior label
    %of the subject)
    if any(reciprocal_set == labels{s,3}) %if behavior is reciprocal
        labels{s,6} = "reciprocal";
        labels{s,7} = 1;
    else
        labels{s,6} = "non-reciprocal";
        labels{s,7} = 0;
    end
    %Add behavior information: social vs. non-social
    if any(social_set == labels{s,3})
        labels{s,8} = "social"; %if behavior is social
        labels{s,9} = 1;
    else
        labels{s,8} = "non-social";
        labels{s,9} = 0;
    end
    %Add block information
    if s<=block_times{1,'end_time_round'}
        labels{s,10} = string(block_times{1,'Behavior'});
        labels{s,11} = 1;
    elseif s>block_times{1,'end_time_round'} && s<=block_times{2,'end_time_round'}
        labels{s,10} = string(block_times{2,'Behavior'});
        labels{s,11} = 2;
    elseif s>block_times{2,'end_time_round'}
        labels{s,10} = string(block_times{3,'Behavior'});
        labels{s,11} = 3;
    end

end


%Rename behavior category to not have acronyms
behav_categ{find(matches(behav_categ,'HIP'))}='Threat to partner';
behav_categ{find(matches(behav_categ,'HIS'))}='Threat to subject';
behav_categ{find(matches(behav_categ,'Pacing/Travel'))}='Travel';
behav_categ{find(matches(behav_categ,'RR'))}='Rowdy Room';
behav_categ{find(matches(behav_categ,'SP'))}='Squeeze partner';
behav_categ{find(matches(behav_categ,'SS'))}='Squeeze Subject';


end