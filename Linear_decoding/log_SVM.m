%% Log_SVM
%% Run a linear decoder on a the neural activity over all or a subset of behaviors

%% Load data
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)

load('Labels_per_sec.mat')
load('Neural_data.mat') % Neural data; array1 is in TEO and array2 is in vlPFC
session_length = size(Unit_rasters,2);
Spike_count_raster = Unit_rasters';

%Compute freq of behavior for the session
behavior_labels = cell2mat({labels{:,3}}');
behav_freq_table = tabulate(behavior_labels);
behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=0,:); % Discard 0 (non-defined behaviors)

% Select behaviors with a minimum # of occurrences
% % min_occurrences = 20;
% % behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%[3,4,5,6,7,8,13,14,15,16];
% % behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
% % behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
behav = [5,6] ;%[1:6,9,17];%[1:6,9:11,16,17]; %manually select behaviors of interest 
behav_categ(behav)

idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
behavior_labels_final = behavior_labels(idx,:);%Same as above but in behavior labels
tabulate(behavior_labels_final);

% % % %Time shift behaviors
% % % shift_length = 5;%in sec
% % % behavior_labels_shifted = behavior_labels_final(shift_length:end);
% % % Spike_count_raster_shifted = Spike_count_raster_final(1:end-shift_length+1,:);

% Run SVM over multiple iterations
num_iter = 1;

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
    %tabulate(Labels);
    
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

% % % figure(1);
% % % confusionchart(confusion_mat_avg)
% % % 
% % % figure(2);
% % % fig = uifigure;
% % % uit = uitable(fig,'Data',C_table);
% % % 


