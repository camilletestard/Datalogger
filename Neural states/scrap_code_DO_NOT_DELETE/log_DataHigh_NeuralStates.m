
%% Load data

filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)

session = filePath(end-9:end);
monkey = filePath(end-14:end-10);

%Load behavioral data
behavior_log = readtable(['EVENTLOG_restructured_',monkey,session,'.csv']);% Behavioral data
behavior_log{:,'start_time_round'}=round(behavior_log{:,'start_time'});
behavior_log{:,'end_time_round'}=round(behavior_log{:,'end_time'});
behavior_log{:,'duration_round'}=behavior_log{:,'end_time_round'}-behavior_log{:,'start_time_round'};

%Load neural data
load(['Neural_data_' session '.mat']) % Neural data

%% Neural states during "awake" state. 1- Compare active behaviors states

% Select behaviors of a particular duration:
min_duration = 1;%in sec
idx = find(behavior_log{:,'duration_round'}>=min_duration);
new_log=behavior_log(idx,:);

% Create structure with chunked neural data per event types.
time_around_epoch = 0;%in sec
event_window=min_duration;%in sec
B = struct();
pseudo_trial = 1;
for behav = 1:size(new_log,1)
    
    total_time = new_log{behav,'start_time_round'}-time_around_epoch ...
        :new_log{behav,'end_time_round'}+time_around_epoch;
    
    divide=floor(length(total_time)./event_window);
    mode=mod(length(total_time),event_window);
    new_chunks = reshape(total_time(1:end-mode),divide,[]);
    for ii = 1:divide
        
        B(pseudo_trial).data = Unit_rasters(:,new_chunks(ii,:));
        
        B(pseudo_trial).condition = char(new_log{behav,'Behavior'});
        %     if strcmp(new_log{behav,'Behavior'},'Resting')
        %         B(behav).condition = 'Non-Resting';
        %     else
        %         B(behav).condition = 'Resting';
        %     end
        
        B(pseudo_trial).epochStarts=1;
        if strcmp(new_log{behav,'Behavior'},'HIP')
            B(pseudo_trial).epochColors=[1,0,0];
        elseif strcmp(new_log{behav,'Behavior'},'Aggression')
            B(pseudo_trial).epochColors=[0.6,0,1];
        elseif strcmp(new_log{behav,'Behavior'},'Approach')
            B(pseudo_trial).epochColors=[0.6,1,1];
        elseif strcmp(new_log{behav,'Behavior'},'Leave')
            B(pseudo_trial).epochColors=[0.2,0.2,1];
        elseif strcmp(new_log{behav,'Behavior'},'Groom Give')
            B(pseudo_trial).epochColors=[0,0,1];
        elseif strcmp(new_log{behav,'Behavior'},'Groom Receive')
            B(pseudo_trial).epochColors=[0,0.75,1];
        elseif strcmp(new_log{behav,'Behavior'},'Proximity')
            B(pseudo_trial).epochColors=[0,0.75,1];
        elseif strcmp(new_log{behav,'Behavior'},'Self-groom')
            B(pseudo_trial).epochColors=[0,1,0];
        elseif strcmp(new_log{behav,'Behavior'},'Scratch')
            B(pseudo_trial).epochColors=[1,0.2,0.3];
        elseif strcmp(new_log{behav,'Behavior'},'HIS')
            B(pseudo_trial).epochColors=[1,0.75,0];
        elseif strcmp(new_log{behav,'Behavior'},'Pacing/Travel')
            B(pseudo_trial).epochColors=[1,0,1];
        elseif strcmp(new_log{behav,'Behavior'},'RR')
            B(pseudo_trial).epochColors=[0,0.5,1];
        elseif strcmp(new_log{behav,'Behavior'},'SS')
            B(pseudo_trial).epochColors=[0,0.8,0.2];
        elseif strcmp(new_log{behav,'Behavior'},'SP')
            B(pseudo_trial).epochColors=[0.3,0.9,0.2];
            elseif strcmp(new_log{behav,'Behavior'},'Submission')
            B(pseudo_trial).epochColors=[1,0.2,1];
        elseif strcmp(new_log{behav,'Behavior'},'Foraging')
            B(pseudo_trial).epochColors=[0,0.5,0];
        elseif strcmp(new_log{behav,'Behavior'},'Drinking')
            B(pseudo_trial).epochColors=[0.5,0,1];
        end
        B(pseudo_trial).type='state';
        
        pseudo_trial=pseudo_trial+1;
    end
end

% Select a random subsample
behaviors = unique({B.condition});
B_subsample = B([],1);
subsample_size = 50
for b = 1:length(behaviors)
    idx=find(ismember({B.condition},behaviors(b)));
    if length(idx)>subsample_size
    B_subsample=[B_subsample;B(idx(randi(length(idx),50,1)))'];
    else
        B_subsample=[B_subsample;B(idx)'];
    end 
end

%Shuffling test
B_rand = B_subsample;
randomization = randperm(length(B));
for i = 1:length(randomization)
B_rand(i).data = B(randomization(i)).data;
end

% % Average "chunked" trials (for better visualization):
% table_B = struct2table(B); avg_B=struct();
% tabulate(table_B{:,'condition'});
% behav_list = unique(table_B{:,'condition'});
% num_events=3; avg_state=1;
% for b = 1:length(behav_list)
%     
%     idx=find(ismember(table_B{:,'condition'},behav_list{b}));
%     divide=floor(length(idx)./num_events);
%     mode=mod(length(idx),num_events);
%     average_over = reshape(idx(1:end-mode),divide,[]);
%     
%     for avg = 1:size(average_over,1)
%         avg_B(avg_state).data = mean(cat(3,B(average_over(avg,1)).data,B(average_over(avg,2)).data,...
%             B(average_over(avg,3)).data,B(average_over(avg,4)).data),3);
%         avg_B(avg_state).condition = behav_list{b};
%         avg_B(avg_state).epochStarts = 1;
%         avg_B(avg_state).epochColors = B(average_over(avg,1)).epochColors;
%         avg_B(avg_state).type='state';
%         
%         avg_state=avg_state+1;
%     end
% end

% Call DataHigh
DataHigh(B_subsample,'DimReduce')