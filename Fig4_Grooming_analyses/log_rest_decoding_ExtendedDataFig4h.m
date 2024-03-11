%% Log_rest_decoding_ExtendedDatFig4h.m
% This script runs a linear decoder on resting bouts. The goal is to
% compare neural tracking of time between grooming and resting.
% C. Testard July 2023

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range=[1:6,11:13,15:16,18];
a_sessions = 1:6; h_sessions = [11:13,15:16,18];

%Set parameters
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_MU =1;%0: MU cluster is excluded; 1:MU cluster is included; 2:ONLY multi-unit cluster
randomsample=0;
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 50;%Number of SVM iterations
smooth= 0; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=0;%lump similar behavioral categories together to increase sample size.
threat_precedence =1;
exclude_sq=1;

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];


    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_MU, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')


    Spike_count_raster{s} = Spike_rasters';
    session_length(s) = size(Spike_count_raster{s},1);
    behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ);
    block_labels = cell2mat({labels{:,12}}');


    %% Extract rest bouts
    RestBouts = zeros(size(behavior_labels));
    RestBouts(behavior_labels== length(behav_categ))=1;

    rest_bout_start = [1; find(diff(RestBouts)==1)+1];
    rest_bout_end = find(diff(RestBouts)==-1);

    if length(rest_bout_end)<length(rest_bout_start) %can happen if grooming went until very end of session
        rest_bout_end(length(rest_bout_start))=length(RestBouts);
    end
    rest_duration = rest_bout_end-rest_bout_start;
    rest_bout=[rest_bout_start, rest_bout_end, rest_duration];
    rest_bout_final = rest_bout;%(rest_bout(:,3)>9,:); %Only consider rest bouts which last at least 9 sec


    %% Extract labels
    neural_data = zscore(Spike_count_raster{s});
    time = 1:length(behavior_labels);

    bouts_to_consider = 1:size(rest_bout_final,1);
    rest_idx=[]; boutid_all=[];
    for b=1:length(bouts_to_consider)
        idx = rest_bout_final(bouts_to_consider(b),1):rest_bout_final(bouts_to_consider(b),2);
        bout_id = ones(size(idx))*b;

        rest_idx = [rest_idx, idx];
        boutid_all = [boutid_all, bout_id];

    end

    neural_data_final=neural_data(rest_idx,:);
    behavior_labels_final = round(rescale(boutid_all)*10);

    tabulate(behavior_labels_final)


    %% Decoding

        disp('Start running SVM...')
        for iter = 1:num_iter

            %subsample to match number of neurons across brain areas
            Labels = behavior_labels_final;
            Input_matrix = neural_data_final;

            %Balance number of trials per class
            uniqueLabels = unique(Labels); %IDentify unique labels (useful when not numbers)
            NumOfClasses = length(uniqueLabels); % Total number of classes
            numericLabels = 1:NumOfClasses; %Numeric name of labels

            labels_temp = Labels;
            for i=1:NumOfClasses
                idx = Labels == uniqueLabels(i);
                labels_temp(idx) = numericLabels(i);
            end
            Labels = labels_temp;

            num_trials = hist(Labels,numericLabels); %number of trials in each class
            minNumTrials = min(num_trials); %use the minimum # of instances

            chosen_trials = [];
            for i = 1:NumOfClasses %for each class
                idx = find(Labels == numericLabels(i)); %find indexes of trials belonging to this class
                rand_i = randsample(length(idx), minNumTrials); %Select a random n number of them
                chosen_trials = [chosen_trials, idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
            end
            Input_matrix = Input_matrix(chosen_trials, :);
            Labels = Labels(chosen_trials);
            Labels_shuffled = Labels(randperm(length(Labels)));

            % Run svm
            [hitrate(iter), C{iter}] = log_SVM_basic_function(Input_matrix, Labels', 5, 0, 0);
            [hitrate_shuffled(iter), C_shuffled{iter}] = log_SVM_basic_function(Input_matrix, Labels_shuffled', 5, 0, 0);

            if mod(iter,10)==1
                disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
            end
        end %end of SVM iterations

        mean_hitrate(s) = mean(hitrate);
        mean_hitrate_shuffled(s) = mean(hitrate_shuffled);

end

figure; hold on
data = mean_hitrate(session_range)';
data_shuffle = mean_hitrate_shuffled(session_range)';
bp = bar([mean(data(:,:)); mean(data_shuffle(:,:))]','FaceAlpha',0.2);
sp1 = scatter(ones(size(data,1))*1,data(:,1), 'filled','b');
sp1 = scatter(ones(size(data_shuffle,1))*2,data_shuffle(:,1), 'filled','r');
