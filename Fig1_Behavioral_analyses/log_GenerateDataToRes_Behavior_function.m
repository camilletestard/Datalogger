function [labels, behav_categ, block_times, monkey,behavior_log, behav_categ_original] = ...
    log_GenerateDataToRes_Behavior_function(filePath, temp_resolution, is_mac,threat_precedence, exclude_sq)

%Log GenerateDataToRes_function
% Input data:
%   1. Behavior of the subject: "EVENTLOG_restructured.csv"
%   2. Behavior of the partner: "EVENTLOG_restructured_partner.csv"
%   3. Spike sorted neural data: "Neural_data_*sessionname*.mat" OR
%   "Neural_data_*sessionname*_isolatedunits.mat"
%   4. Block info: "session_block_schedule.csv"

% Arguments:
%   filePath: is the experimental data path
%   Temp_resolution: is the temporal resolution at which we would like to
%   analyze the dat
%   Channel_flag: specifies with channels to include: only TEO array, only
%   vlPFC array or all channels
%   is_mac: specifies whether the code is run on a mac or pc
%   with_NC: specifies whether the "noise cluster" (or the first cell of
%   every channel) is included (1) or not (0). If with_NC=2 then we only
%   include the noise cluster (not the other neurons).
%   isolatedOnly: specifies if only the well isolated units are considered.
%   smooth: is the data smoothed using function conv.
%   sigma: size of smoothing
%   threat_precedence: does aggression take precedence over threat events.

% This function formats the raw input data to have the following elements (output data):
% 1. Spike_rasters: Neural data matrix, size [Time (to chosen resolution) x #neurons]
% 2. labels: Label vector which describes the behavior at time t [Time (to chosen resolution) x 4]
%       1st column includes all behaviors in "plain english"
%       2nd column behavior number code
%       3rd column unique behavior code (when two occur
%       simultaneously, we chose one, see below for details)
%       4th column whether behavior happens in isolation or co-occurs
%       with another. String description.
%       5th column numerical code for the type of co-occurrence
%       6th column is the behavior "reciprocal" (i.e. partner behavior can
%       be 100% predictted by subject behavior and vice-versa)
%       7th column binary code reciprocal (1) vs. not (0)
%       8th column is the behavior "social" or not (i.e. done with a
%       conspecific)
%       9th column binary code for social (1) or not (0).
%       10th column indicates the block ID in which we are
%       (Paired, "female" neighbor; Paired, "male" neighbor or "alone")
%       11th column gives a corresponding numerical value to the block order in time (1st, 2nd 3rd).
%       12th column is a numerical version of block ID (social context).
%       13th column whether individual is paired (1) or alone (0)
% 3. labels_partner: same as above but for the partner
% 4. behav_categ: behavioral categories
% 5. block_times: Order and timing of blocks during the session
% 6. monkey: ID of the subject monkey
% 7. reciprocal_set: Set of reciprocal behaviors
% 8. social_set: Set of social behaviors
% 9.unit_count: Number of units per brain area
% 10.groom_labels_all: Label of grooming category. Columns correspond to:
%       1. Behavior label. If not grooming (7 groom give, 8 groom receive),
%       then all categories below will be 0. If grooming, the bout can be
%       qualified as:
%       2. Is it the start (1), end (2) or middle (3) of grooming bout
%       3. Grooming after a threat (2) or not (1)
%       4. Grooming reciprocated (2) or not (1)
%       5. Grooming initiated by subject (2) or not (1)
% 11.brain_label: brain area label for each unit
% 12.behavior_log: raw behavioral log
% 13.behav_categ_original: original category labeling

% Camille Testard - Created Nov. 2021
% Last update - March 2023

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
block_log = readtable('session_block_schedule.csv');% for block info

% Assign threat to partner and threat to sibject labels to squeeze events.
if exclude_sq
    index=find(strcmp('SS',behavior_log{:,'Behavior'}));
    behavior_log{index,'Behavior'}={'HIS'};
    index=find(strcmp('SP',behavior_log{:,'Behavior'}));
    behavior_log{index,'Behavior'}={'HIP'};
end


%% Preprocessing: behavioral log and neural data at specified resolution

%Round times and get duration for behavioral logs (subject & partner)
behavior_log{:,'start_time_round'}=round(behavior_log{:,'start_time'}*temp_resolution);
behavior_log{:,'end_time_round'}=round(behavior_log{:,'end_time'}*temp_resolution);
behavior_log{:,'duration_round'}=behavior_log{:,'end_time_round'}-behavior_log{:,'start_time_round'};

%Eliminate behaviors that do not meet the minimum length
%Note that this will be an issue only for time resolution >1sec
min_length = 1/temp_resolution;
idx = find(behavior_log{:,'duration_round'}<min_length);
behavior_log{idx,'start_time_round'} = 0; behavior_log{idx,'end_time_round'} = 0;

%Get block times (at the end of the EventLog_Restructured)
block_times = behavior_log(end-2:end,:);
behavior_log(end-2:end,:) = [];



%% Get behavior label vector for each time bin at specified resolution

%Create behavior key
%behav_categ = unique(behavior_log{:,'Behavior'}); %Get all the unique behaviors
%Set it constant across sessions. This makes coding consistent across
%sessions.
behav_categ = ["Aggression","Proximity","Groom Give", "HIP","Foraging", "Vocalization", "Masturbating",...
    "Submission", "Approach","Yawning","Self-groom","HIS","Other monkeys vocalize", "Lip smack",...
    "Groom Receive","Leave","Drinking","Pacing/Travel","Scratch","RR", "Butt sniff","Grm prsnt","Mounting","SS","SP"];
behav_categ = sort(behav_categ);

behav_categ{length(behav_categ)+1}='Rest'; %Add rest as a behavior (no defined behavior ongoing)

%For behaviors that often co-occur with other behaviors, determine priority
double_behav_set = [find(matches(behav_categ,'Proximity')), find(matches(behav_categ,"RR")), find(matches(behav_categ,'Other monkeys vocalize'))];% When co-occurring, other behaviors will take precedence over these ones

%annotation: omv = other monkeys vocalize
%            grmpr = groom presentation
grmpr = find(matches(behav_categ,'Grm prsnt'));
omv = find(matches(behav_categ,'Other monkeys vocalize'));

%Create set of reciprocal behaviors (i.e. behavior of partner can be 100%
%predicted by behavior of subject)
reciprocal_set = [find(matches(behav_categ,'RR')), find(matches(behav_categ,'Other monkeys vocalize')), ...
    find(matches(behav_categ,'Proximity')), find(matches(behav_categ,"Groom Give")), find(matches(behav_categ,"Groom Receive")),...
    find(matches(behav_categ,"HIP")), find(matches(behav_categ,"HIS"))];

%Create set of social behaviors (which happen with conspecifics)
social_set = [find(matches(behav_categ,'Proximity')), find(matches(behav_categ,"Groom Give")), find(matches(behav_categ,"Groom Receive")),...
    find(matches(behav_categ,"Submission")), find(matches(behav_categ,"Approach")), find(matches(behav_categ,"Leave")), find(matches(behav_categ,"Butt sniff")),...
    find(matches(behav_categ,"Grm prsnt")), find(matches(behav_categ,"Aggression"))];


%% Create behavior label vector (label every window of the session)
% This cell matrix will have six columns. The first column is the full
% name of the behavior label

%Create event time intervals:
%For subject monkey
start_times = behavior_log{:,'start_time_round'};
end_times = behavior_log{:,'end_time_round'};
Intervals = [start_times end_times];
Intervals(strcmp(behavior_log{:,'Behavior'},"Camera Sync"),2) =Intervals(strcmp(behavior_log{:,'Behavior'},"Camera Sync"),1);


%%%%%% Create labels vector for SUBJECT monkey %%%%%%
length_recording = max(behavior_log{:,'end_time_round'}+10);
labels = cell(length_recording,12); %initialize dataframe
for s = 1:length_recording %for all secs in a session
    % this finds the index of the rows(2) that have x in between
    idx = find(s >= Intervals(:,1) & s < Intervals(:,2)); %find if this second belong to any interval
    %IMPORTANT note: interval includes lower bound but excludes upper boundary as is.
    if ~isempty(idx)%if it belongs to an interval
        labels{s,1} = behavior_log{idx,'Behavior'}; %add behavior label in [plain english]
        labels{s,2} = find(matches(behav_categ,labels{s,1})); %add behavior label in [number]

        if all(~strcmp(labels{s,1},"Camera Sync")) %if not camera sync

            if length(labels{s,2})>1 %If one behavior co-occurs with other behavior(s)
                if ~isempty(setdiff(labels{s,2}, double_behav_set)) %If behavior co-occurs with proximity, other monkeys vocalize or RR
                    labels{s,3} = setdiff(labels{s,2}, double_behav_set); % only consider the other behavior (it takes precedence over proximity and RR)
                    labels{s,4} = 'co-occur with prox, omv or RR';
                    labels{s,5} = 2;
                elseif isempty(setdiff(labels{s,2}, double_behav_set(1:2))) && any(labels{s,3}==omv) %if RR or proximity co-occur with omv
                    labels{s,3}=omv; % prioritize Other monkeyz vocalize
                    labels{s,4} = 'omv co-occurs with prox & RR';
                    labels{s,5} = 3;
                else
                    isempty(setdiff(labels{s,2}, double_behav_set(1:2)));%If proximity and RR co-occur
                    labels{s,3}=find(matches(behav_categ,'RR')); % prioritize Rowdy Room
                    labels{s,4} = 'prox & RR co-occur';
                    labels{s,5} = 4;
                end
            else %If only one behavior happens in that sec OR if other types of co-occurrence happen
                labels{s,3} = labels{s,2};
                labels{s,4} = 'single';
                labels{s,5} = 1;
            end


            if length(labels{s,3})~=1 %If two behaviors are co-occurring which do not include omv, proximity or RR

                if any(labels{s,3}==grmpr) % if one of the behavior includes groom present
                    labels{s,3}=grmpr; %Keep groom present
                    labels{s,4} = 'grmpr co-occur';
                    labels{s,5} = 5;
                elseif any(labels{s,3}==find(matches(behav_categ,'Scratch'))) %if on of the behavior includes scratch
                    labels{s,3}=setdiff(labels{s,3}, find(matches(behav_categ,'Scratch'))); %Consider the other behavior
                    labels{s,4} = 'scratch co-occur';
                    labels{s,5} = 6;
                    if length(labels{s,3})~=1 %there are more than 2 behaviors that co-occur with scratch,
                        labels{s,3}= labels{s,3}(2); %choose 2nd behavior
                    end
                elseif any(labels{s,3}==find(matches(behav_categ,'Aggression'))) %If one of the behavior includes aggression
                    labels{s,3}=find(matches(behav_categ,'Aggression')); %keep aggression
                    labels{s,4} = 'aggression co-occur';
                    labels{s,5} = 7;
                elseif any(labels{s,3}==find(matches(behav_categ,'Submission'))) %If one of the behavior includes submission
                    labels{s,3}=find(matches(behav_categ,'Submission')); %keep submission
                    labels{s,4} = 'submission co-occur';
                    labels{s,5} = 8;
                elseif any(labels{s,3}==find(matches(behav_categ,'Masturbating'))) %If one of the behavior includes masturbation
                    labels{s,3}=find(matches(behav_categ,'Masturbating')); %keep masturbation
                    labels{s,4} = 'masturbating co-occur';
                    labels{s,5} = 9;
                elseif any(labels{s,3}==find(matches(behav_categ,'Mounting'))) %If one of the behavior includes Mounting
                    labels{s,3}=find(matches(behav_categ,'Mounting')); %keep mounting
                    labels{s,4} = 'mounting co-occur';
                    labels{s,5} = 10;
                elseif any(labels{s,3}==find(matches(behav_categ,'Approach'))) %If one of the behavior includes Approach
                    labels{s,3}=find(matches(behav_categ,'Approach')); %keep approach
                    labels{s,4} = 'approach co-occur';
                    labels{s,5} = 11;
                elseif any(labels{s,3}==find(matches(behav_categ,'Leave'))) %If one of the behavior includes Leave
                    labels{s,3}=find(matches(behav_categ,'Leave')); %keep leave
                    labels{s,4} = 'leave co-occur';
                    labels{s,5} = 12;
                else %Otherwise just choose the second behavior for now...
                    %                 error('More than one behavior simultansouly')
                    %                 return
                    labels{s,3}= labels{s,3}(2);
                    labels{s,4} = 'Other key behav co-occur';
                    labels{s,5} = 13;
                end

            end

            %%%%%%%%%%%%  THREAT PRECEDENCE CLAUSE %%%%%%%%%%%%
            if threat_precedence == 1
                if any(labels{s,2}==find(matches(behav_categ,'HIP'))) %If one of the behavior includes threat to partner
                    labels{s,3}=find(matches(behav_categ,'HIP')); %Keep HIP
                    labels{s,4} = 'HIP co-occur';
                    labels{s,5} = 14;
                elseif any(labels{s,2}==find(matches(behav_categ,'HIS'))) %If one of the behavior includes threat to subject
                    labels{s,3}=find(matches(behav_categ,'HIS')); %Keep HIS
                    labels{s,4} = 'HIS co-occur';
                    labels{s,5} = 15;
                end
            end

        else
            labels{s,1} = NaN; labels{s,2} = length(behav_categ); labels{s,3} = length(behav_categ); labels{s,4} = 'NA'; labels{s,5} = 0;%Set behavior category to "NaN" and label to rest
        end

    else %if second belongs to no interval, set it to "rest"
        labels{s,1} = NaN; labels{s,2} = length(behav_categ); labels{s,3} = length(behav_categ); labels{s,4} = 'NA'; labels{s,5} = 0;%Set behavior category to "NaN" and label to rest

    end

    %%%%%%%%%%%%%%%%%%%%%%%%
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

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %Add block information
    if s<=block_times{1,'end_time_round'}
        labels{s,10} = string(block_log{strcmp(full_session_name, block_log{:,'session_name'}),2}); %identity of block (female neighbor, male neighbor or alone block)
        labels{s,11} = 1; %block order

        if labels{s,10}=="female"
            labels{s,12} = 1;%numerical form of block identity
            labels{s,13} = 1; %paired block
        elseif labels{s,10}=="male"
            labels{s,12} = 2;
            labels{s,13} = 1;%paired block
        else
            labels{s,12} = 3;
            labels{s,13} = 0;%alone block
        end

    elseif s>block_times{1,'end_time_round'} && s<=block_times{2,'end_time_round'}
        labels{s,10} = string(block_log{strcmp(full_session_name, block_log{:,'session_name'}),3});
        labels{s,11} = 2;

        if labels{s,10}=="female"
            labels{s,12} = 1;%numerical form
            labels{s,13} = 1; %paired block
        elseif labels{s,10}=="male"
            labels{s,12} = 2;
            labels{s,13} = 1; %paired block
        else
            labels{s,12} = 3;
            labels{s,13} = 0;%alone block
        end

    elseif s>block_times{2,'end_time_round'}
        labels{s,10} = string(block_log{strcmp(full_session_name, block_log{:,'session_name'}),4});
        labels{s,11} = 3;

        if labels{s,10}=="female"
            labels{s,12} = 1;%numerical form
            labels{s,13} = 1; %paired block
        elseif labels{s,10}=="male"
            labels{s,12} = 2;
            labels{s,13} = 1; %paired block
        else
            labels{s,12} = 3;
            labels{s,13} = 0;%alone block
        end

    end

end

%Rename behavior category to not have acronyms
behav_categ_original = behav_categ;
behav_categ{find(matches(behav_categ,'HIP'))}='Threat to partner';
behav_categ{find(matches(behav_categ,'HIS'))}='Threat to subject';
behav_categ{find(matches(behav_categ,'Pacing/Travel'))}='Travel';
behav_categ{find(matches(behav_categ,'RR'))}='Rowdy Room';
behav_categ{find(matches(behav_categ,'Grm prsnt'))}='Groom sollicitation';
behav_categ{find(matches(behav_categ,'Groom Give'))}='Groom partner';
behav_categ{find(matches(behav_categ,'Groom Receive'))}='Getting groomed';


end