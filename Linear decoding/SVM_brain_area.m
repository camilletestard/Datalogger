%% Log_SVM_subject_batch
% Run a linear decoder on the neural activity for the subject's behavior
% (only including behaviors with a minimum # occurrence)
% Batch version of script
% March 2022 - Camille Testard

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


%Set parameters
with_partner =0;
temp = 1; temp_resolution = 1;
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
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];


    %% Get data with specified temporal resolution and channels
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
    end

    disp('Data Loaded')

    Spike_count_raster_final = Spike_rasters;
    behavior_labels_final = brain_label;

    tabulate(removecats(categorical(behavior_labels_final)))


    %% Run SVM over multiple iterations

    disp('Start running SVM...')
    for iter = 1:num_iter


        Labels = behavior_labels_final;
        Input_matrix = Spike_count_raster_final;

        %Balance number of trials per class
        uniqueLabels = unique(Labels); %IDentify unique labels (useful when not numbers)
        NumOfClasses = length(uniqueLabels); % Total number of classes
        numericLabels = 1:NumOfClasses; %Numeric name of labels

        labels_temp = Labels;
        for i=1:NumOfClasses
            idx = strcmp(Labels, uniqueLabels(i));
            labels_temp(idx) = numericLabels(i);
            labels_id{i,1} = uniqueLabels(i);
        end
        Labels = str2double(labels_temp);

        num_trials = hist(Labels,numericLabels); %number of trials in each class
        minNumTrials = min(num_trials); %30; %find the minimum one %CT change to have 30 of each class
        chosen_trials = [];
        for i = 1:NumOfClasses %for each class
            idx = find(Labels == numericLabels(i)); %find indexes of trials belonging to this class
            rand_i = randsample(length(idx), minNumTrials); %Select a random n number of them
            chosen_trials = [chosen_trials; idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
        end
        if size(chosen_trials) ~= [minNumTrials*2,1]
            chosen_trials = reshape(chosen_trials,[minNumTrials*2,1]);
        end
        Input_matrix = Input_matrix(chosen_trials, :);
        Labels = Labels(chosen_trials)';
        Labels_shuffled = Labels(randperm(length(Labels)));

        % Run svm
        [hitrate(iter), C{iter}] = log_SVM_basic_function(Input_matrix, Labels, 5, 0, 0);
        [hitrate_shuffled(iter), C_shuffled{iter}] = log_SVM_basic_function(Input_matrix, Labels_shuffled, 5, 0, 0);

        if mod(iter,10)==1
            disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
        end
    end %End of SVM for loop

    mean_hitrate(s) = mean(hitrate)
    sd_hitrate(s) = std(hitrate);
    mean_hitrate_shuffled(s) = mean(hitrate_shuffled)
    sd_hitrate_shuffled = std(hitrate_shuffled);

    %         C_concat=cat(3,C{:}); %Get confusion matrix
    %         confusion_mat_avg=round(mean(C_concat,3)*100); %Average over SVM iterations
    %         rowNames{s} = {labels_id{:,1}}; colNames{s} = {labels_id{:,1}}; %Get behavior names
    %         C_table{s} = array2table(confusion_mat_avg,'RowNames',rowNames{s},'VariableNames',colNames{s});

    clear labels_id

end %End of session for loop

%% Plot all sessions results

%Change savePath for all session results folder:
cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);

%Plot decoding accuracy for all sessions, separated by monkey
figure;  hold on;
for s = a_sessions
    y = mean_hitrate(s);
    scatter(s, y, 'r', 'filled','MarkerFaceAlpha',0.7)
end
for s = h_sessions
    y = mean_hitrate(s);
    scatter(s, y, 'b', 'filled','MarkerFaceAlpha',0.7)
end
chance_level = 1/2;
yline(chance_level,'--','Chance level', 'FontSize',16)
xlim([0.8 16.2]); ylim([0.4 1])
ax = gca;
ax.FontSize = 14;
ylabel('Deconding accuracy','FontSize', 18); xlabel('Session#','FontSize', 18)
title('Decoding accuracy for brain area','FontSize', 14)

saveas(gcf,['Decoding_BrainArea_allSessions.png'])
close all

