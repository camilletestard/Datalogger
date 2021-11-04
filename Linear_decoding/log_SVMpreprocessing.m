%% Log SVM pre-processing
% Format the raw data to have two elements:
% 1. Neural data matrix size [Time (in sec) x #neurons]
% 2. Label vector which describes the behavior at time t [Time (in sec) x 1]
% Camille Testard - Sept. 2021

%% Load data
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)

%behav_log = readtable('behavioral_log_session1.csv');
behavior_log = readtable('EVENTLOG_restructured.csv');% Behavioral data
load('Neural_data.mat') % Neural data; array1 is in TEO and array2 is in vlPFC
session_length = size(Unit_rasters,2);
Spike_count_raster = Unit_rasters';

%Preprocessing: round times in behavioral log
behavior_log{:,'start_time_round'}=round(behavior_log{:,'start_time'});
behavior_log{:,'end_time_round'}=round(behavior_log{:,'end_time'});
behavior_log{:,'duration_round'}=behavior_log{:,'end_time_round'}-behavior_log{:,'start_time_round'};

%% Get behavior label vector for each second

%Create intervals:
start_times = behavior_log{:,'start_time_round'};
end_times = behavior_log{:,'end_time_round'};
Intervals = [start_times end_times];

%Create behavior key
behav_categ = unique(behavior_log{:,'Behavior'});
double_behav_set = [11,12];

%Create label vector
labels = cell(session_length,3);
for s = 1:session_length %for all secs in a session
    % this finds the index of he rows(2) that have x in between
    idx = find(s > Intervals(:,1) & s < Intervals(:,2));
    if ~isempty(idx)
        labels{s,1} = behavior_log{idx,'Behavior'};
        labels{s,2} = find(matches(behav_categ,labels{s,1}));
        if length(labels{s,2})>1
            labels{s,3} = setdiff(labels{s,2}, double_behav_set);
        else
            labels{s,3} = labels{s,2};
        end
        if length(labels{s,3})~=1
            labels{s,3}= labels{s,3}(1);
        end
    else
        labels{s,1} = NaN; labels{s,2} = 0; labels{s,3} = 0;
    end
end

    % %% Save variables
    save([filePath '\Labels.mat'],'labels','behav_categ');


behavior_labels = cell2mat({labels{:,3}}');
behav_freq_table = tabulate(behavior_labels);
behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=0,:); % Discard 0 (non-defined behaviors)

% Select behaviors with a minimum # of occurrences
min_occurrences = 30;
%behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%[3,4,5,6,7,8,13,14,15,16];
behav = [3,4,5,6,7,8,13,14,15,16]; %excluding proximity which is a source of confusion.

idx = find(ismember(behavior_labels,behav));
Spike_count_raster_final = Spike_count_raster(idx,:);
behavior_labels_final = behavior_labels(idx,:);
tabulate(behavior_labels);

%Time shift behaviors
shift_length = 5;%in sec
behavior_labels_shifted = behavior_labels_final(shift_length:end);
Spike_count_raster_shifted = Spike_count_raster_final(1:end-shift_length+1,:);

%% Run SVM over multiple iterations
num_iter = 100;

for iter = 1:num_iter
    
    clearvars -except behav_categ behavior_labels_final behavior_labels_shifted Spike_count_raster_final Spike_count_raster_shifted num_iter iter hitrate hitrate_shuffled C C_shuffled
    
    %Balance number of trials per class
    Labels = behavior_labels_final;
    Input_matrix = Spike_count_raster_final;
    
    uniqueLabels = unique(Labels); %IDentify unique labels (useful when not numbers)
    NumOfClasses = length(uniqueLabels); % Total number of classes
    numericLabels = 1:NumOfClasses; %Numeric name of labels
    
    labels_temp = Labels;
    for i=1:NumOfClasses,
        idx = Labels == uniqueLabels(i);
        labels_temp(idx) = numericLabels(i);
        labels_id{i,1} = uniqueLabels(i); labels_id{i,2}=behav_categ{uniqueLabels(i)} ;
    end
    Labels = labels_temp;
    
    num_trials = hist(Labels,numericLabels); %number of trials in each class
    minNumTrials = min(num_trials); %find the minimum one
    chosen_trials = [];
    for i = 1:NumOfClasses %for each class
        idx = find(Labels == numericLabels(i)); %find indexes of trials belonging to this class
        rand_i = randsample(length(idx), minNumTrials); %Select a random n number of them
        chosen_trials = [chosen_trials; idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
    end
    Input_matrix = Input_matrix(chosen_trials, :);
    Labels = Labels(chosen_trials, :);
    Labels_shuffled = Labels(randperm(length(Labels)));
    tabulate(Labels);
    
    % % %% Save variables
    % save([filePath '\SVM_input.mat'],'Input_matrix', 'Labels','Labels_shuffled');
    
    % Run svm
    [hitrate(iter), C{iter}] = log_SVM_basic_function(Input_matrix, Labels, 5, 0, 0);
    [hitrate_shuffled(iter), C_shuffled{iter}] = log_SVM_basic_function(Input_matrix, Labels_shuffled, 5, 0, 0);
end

mean_hitrate = mean(hitrate)
sd_hitrate = std(hitrate)
mean_hitrate_shuffled = mean(hitrate_shuffled)
sd_hitrate_shuffled = std(hitrate_shuffled)

C_concat=cat(3,C{:});
confusion_mat_avg=round(mean(C_concat,3)*100)
rowNames = {labels_id{:,2}}; colNames = {labels_id{:,2}};
C_table = array2table(confusion_mat_avg,'RowNames',rowNames,'VariableNames',colNames)

figure(1);
confusionchart(confusion_mat_avg)

figure(2);
fig = uifigure;
uit = uitable(fig,'Data',C_table);


