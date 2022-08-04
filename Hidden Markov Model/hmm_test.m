

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,15:16];
session_range_with_partner=[1:3,11:13];
session_range_with_neighbor=11;

%Set parameters
with_partner =0;% 0: no partner; 1: with partner; 2: with neighbor
temp = 1; temp_resolution = 500;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 50;%Number of SVM iterations

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
elseif with_partner ==2
    session_range = session_range_with_neighbor;
    a_sessions = nan; h_sessions = 11;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=1;
%Set path
filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];

chan = 1;
[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function_hmm(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);

disp('Data Loaded')

Spike_count_raster = zscore(Spike_rasters');
behavior_labels = cell2mat({labels{:,3}}');

%remove rest epochs
idx_keep = find(behavior_labels~=length(behav_categ)); %1:length(behavior_labels);%
Spike_count_raster_final = Spike_count_raster(idx_keep,:);
behavior_labels_final = behavior_labels(idx_keep);

%reset labels
unique_lbl = unique(behavior_labels_final);
labels_final = behavior_labels_final;
for lb = 1:length(unique_lbl)
    labels_final(behavior_labels_final == unique_lbl(lb)) = lb;
end
[1:length(unique_lbl)'; unique_lbl(1:length(unique_lbl))']

%% HMM estimation

%%%%%%%%%%%%%
%Chen Matlab toolbox
num_states = 4;
[model, llh] = hmmEm(labels_final',num_states);
figure; heatmap(model.A) %transition probabilities
figure; heatmap(model.E) %observation probability in each state
figure; plot(llh)

z = hmmViterbi(labels_final',model);
figure; plot(z)
figure; hist(z)

%%%%%%%%%%%%%
%Matlab toolbox
seq = labels_final';
[TRANS,EMIS] = hmmtrain(seq, model.A, model.E,'Tolerance',1e-10, 'Maxiterations',1e5);
figure; heatmap(TRANS)
STATES = hmmviterbi(seq,TRANS,EMIS)
figure; hold on; plot(STATES); plot(labels_final)
count=hist(STATES);

%%%%%%%%%%%%%%
%Decode
num_iter = 100;
min_occurrences = min(count);
state_labels_final = STATES;
for iter = 1:num_iter

    %             clearvars -except savePath behav_categ behavior_labels_final behavior_labels_shifted Spike_count_raster_final...
    %                 Spike_count_raster_shifted num_iter iter hitrate hitrate_shuffled C C_shuffled temp_resolution...
    %                 channel_flag filePath chan temp mean_hitrate sd_hitrate mean_hitrate_shuffled C_table behavs_eval behav is_mac min_occurrences

    %Balance number of trials per class
    Labels = state_labels_final';
    Input_matrix = Spike_count_raster_final;

    uniqueLabels = unique(Labels); %IDentify unique labels (useful when not numbers)
    NumOfClasses = length(uniqueLabels); % Total number of classes
    numericLabels = 1:NumOfClasses; %Numeric name of labels

    labels_temp = Labels;
    for i=1:NumOfClasses
        idx = Labels == uniqueLabels(i);
        labels_temp(idx) = numericLabels(i);
        labels_id{i,1} = uniqueLabels(i); labels_id{i,2}=behav_categ{uniqueLabels(i)} ;
    end
    Labels = labels_temp;

    num_trials = hist(Labels,numericLabels); %number of trials in each class
    minNumTrials = min_occurrences;%min(num_trials); %find the minimum one %CT change to have 30 of each class
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

    disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
end
mean(hitrate)
mean(hitrate_shuffled)