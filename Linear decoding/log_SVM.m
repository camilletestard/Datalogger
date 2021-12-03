%% Log_SVM
%% Run a linear decoder on a the neural activity over all or a subset of behaviors

%% Load data

%Set path
is_mac = 1;
if is_mac
    cd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
else
    cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
end
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)

if is_mac
    cd('~/Dropbox (Penn)/Datalogger/Results/')
else
    cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Results/')
end
savePath = uigetdir('', 'Please select the result directory');

%Set temporal resolution
temp = 1; temp_resolution = 1;
for temp_resolution = [1, 5, 10, 100] %5sec, 1sec, 100msec
    %temp_resolution = [1/5, 1/2, 1, 5, 10, 100] %5sec, 1sec, 100msec
    %1 for second resolution, 10 for 100msec resolution, 100 for 10msec resolution, 1000 for msec resolution. etc.
    %0.1 for 10sec resolution, 1/5 for 5sec resolution

    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all";
    for channel_flag = ["vlPFC", "TEO", "all"]

        %Get data with specified temporal resolution and channels
        [Spike_rasters, labels, behav_categ]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag);
        %filePath is the experimental data path
        %Temp_resolution is the temporal resolution at which we would like to
        %analyze the dat
        %Channel_flag specifies with channels to include: only TEO array, only
        %vlPFC array or all channels
        disp('Data Loaded')

        %Previously loaded files individually. Keep for now.
        % % % load(['Data_' num2str(1000/temp_resolution) 'msec_res.mat'])%Load data
        % % % load('Labels.mat')%Behavioral labels
        % % % load('Neural_data.mat') % Neural data; array1 is in TEO and array2 is in vlPFC

        Spike_count_raster = Spike_rasters';

        %% Select behaviors to decode
        %Compute freq of behavior for the session
        behavior_labels = cell2mat({labels{:,3}}');
        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=0,:); % Discard 0 (non-defined behaviors)

        % Select behaviors with a minimum # of occurrences
        % min_occurrences = 50;
        % behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%[3,4,5,6,7,8,13,14,15,16];
        % behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
        % behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
        behav = [3,4,5,6,7,8,17]; %[1:6,9:11,16,17]; %manually select behaviors of interest
        behavs_eval = behav_categ(behav);

        idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx,:);%Same as above but in behavior labels
        tabulate(behavior_labels_final);

        %% Time shift behaviors
        % % % shift_length = 5;%in sec
        % % % behavior_labels_shifted = behavior_labels_final(shift_length:end);
        % % % Spike_count_raster_shifted = Spike_count_raster_final(1:end-shift_length+1,:);


        %% Run SVM over multiple iterations
        num_iter = 1000;
        
        disp('Start running SVM...')
        for iter = 1:num_iter

            clearvars -except savePath behav_categ behavior_labels_final behavior_labels_shifted Spike_count_raster_final...
                Spike_count_raster_shifted num_iter iter hitrate hitrate_shuffled C C_shuffled temp_resolution...
                channel_flag filePath chan temp mean_hitrate sd_hitrate mean_hitrate_shuffled C_table behavs_eval

            %Balance number of trials per class
            Labels = behavior_labels_final;
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
            minNumTrials = 50; %min(num_trials); %find the minimum one %CT change to have 200 of each class
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
%             [hitrate_shuffled(iter), C_shuffled{iter}] = log_SVM_basic_function(Input_matrix, Labels_shuffled, 5, 0, 0);
            
            disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
        end

        disp('****************************************************************************')
        disp([num2str(1000/temp_resolution) 'msec resolution, channels: ' channel_flag '. DONE'])
        disp('****************************************************************************')

        mean_hitrate(temp, chan) = mean(hitrate)
        sd_hitrate(temp, chan) = std(hitrate);
%         mean_hitrate_shuffled(temp, chan) = mean(hitrate_shuffled)
%         sd_hitrate_shuffled = std(hitrate_shuffled);

        C_concat=cat(3,C{:});
        confusion_mat_avg=round(mean(C_concat,3)*100);
        rowNames = {labels_id{:,2}}; colNames = {labels_id{:,2}};
        C_table{temp, chan} = array2table(confusion_mat_avg,'RowNames',rowNames,'VariableNames',colNames);

        chan = chan +1;

        clearvars -except savePath mean_hitrate sd_hitrate mean_hitrate_shuffled C_table temp_resolution channel_flag filePath chan temp behavs_eval
    end
    temp = temp+1;
    
end

%rowNames = ["5sec", "2sec", "1sec", "500msec", "100msec", "10msec"]; colNames = ["vlPFC","TEO","all"];
%rowNames = ["1sec", "500msec", "100msec", "10msec"]; colNames = ["vlPFC","TEO","all"];
rowNames = ["1sec"]; colNames = ["vlPFC","TEO","all"];
result_hitrate = array2table(mean_hitrate,'RowNames',rowNames,'VariableNames',colNames)
result_sdhitrate = array2table(sd_hitrate,'RowNames',rowNames,'VariableNames',colNames)

% save([savePath '\SVM_results_' num2str(length(behavs_eval)) 'behav.mat'], 'mean_hitrate', 'sd_hitrate', 'C_table', 'result_hitrate', 'result_sdhitrate', 'behavs_eval')
% writetable(result_hitrate,[savePath '\SVM_results_' num2str(length(behavs_eval)) 'behav.csv'],'WriteRowNames',true,'WriteVariableNames',true); 

figure; hold on
for b = 1:size(mean_hitrate,1)
    y = mean_hitrate(b,:);
    std_dev = sd_hitrate(b,:);
    errorbar(y,std_dev,'s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
    %plot(x,y,'Color','k')
end
chance_level = 1/length(behavs_eval);
yline(chance_level,'--','Chance level')
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 1])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14; 
ylabel('Deconding accuracy','FontSize', 18); xlabel('Brain area','FontSize', 18)
title('Decoding accuracy by brain area (window size = 1sec)','FontSize', 18)



% x = 1:3; x2 = [x, fliplr(x)];
% cmap = cool(size(mean_hitrate,1));
% cmap_sd = cool(size(mean_hitrate,1));
% for b = 1:size(mean_hitrate,1)
%     y = mean_hitrate(b,:);
%     std_dev = sd_hitrate(b,:);
%     curve1 = y + std_dev;
%     curve2 = y - std_dev;
%     inBetween = [curve1, fliplr(curve2)];
%     h = fill(x2, inBetween,cmap(b,:));
%     h.FaceAlpha = 0.5; 
%     plot(x,y,'Color','k')
% end
% chance_level = 1/length(behavs_eval);
% yline(chance_level,'--','Chance level')
% xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 1])
% xticklabels({'','vlPFC','TEO','all',''})
% ylabel('Deconding accuracy'); xlabel('Brain area')
% title('Decoding accuracy by window size and brain area')
%leg = legend("1sec",'',"500msec",'',"100msec",'',"10msec");
%title(leg,'Window size')

cd(savePath)
saveas(gcf,['SVM_results_' num2str(length(behavs_eval)) 'behav.png'])
close all