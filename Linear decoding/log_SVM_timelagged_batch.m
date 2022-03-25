%% Log_SVM_timelagged
%% Run a linear decoder on a the neural activity over all or a subset of behaviors


%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,15:18];
session_range_with_partner=[1:3,11:13];


%Set parameters
with_partner =0;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 100;%Number of SVM iterations
min_occurrences = 30; %Min number of occurrence for a particular behavior

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:18];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];
    
    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1;
    for channel_flag = ["vlPFC", "TEO", "all"]

        pred =1; predict =1;
        for predict =[1,2] %Forward or backward lag. I.e. prediction of future behavior (predict =1) or impact of past behavior (predict=2)?

            lag=1; lag_length = 0;
            for lag_length = [0,2,5,10,20] %in window size
                %IMPORTANT NOTE: at 1 sec, thre are very few non-overlapping
                %windows which makes the comparison difficult...

                % Get data with specified temporal resolution and channels
                if with_partner ==1
                    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
                else
                    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
                end

                disp('Data Loaded')

                Spike_count_raster_init = Spike_rasters';
                behavior_labels_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject


                %% Time shift behaviors

                %lag_length = 0; %Set length of time lag (in window size). Will include in a for loop eventually

                if predict == 1 %If predicting future behavior
                    behavior_labels = behavior_labels_init(1+lag_length:end); %Time shift behavior regressor
                    Spike_count_raster = Spike_count_raster_init(1:end-lag_length,:); %Truncate neural data

                    %check the amount of labels that differ in the lag vs non-lagged versions
                    behavior_labels_nonlagged = behavior_labels_init(1:end-lag_length);
                    diff_behav_vectors = behavior_labels - behavior_labels_nonlagged;
                    perc_diff_before = round(length(find(diff_behav_vectors~=0))/length(diff_behav_vectors)*100);
                    fprintf('Percent difference between-time lagged and regular labels BEFORE selecting behaviors: %s \n', num2str(perc_diff_before));

                elseif predict == 2  %If predicting past behavior
                    behavior_labels = behavior_labels_init(1:end-lag_length); %Time shift behavior regressor
                    Spike_count_raster = Spike_count_raster_init(1+lag_length:end,:);  %Truncate neural data

                    %check the amount of labels that differ in the lag vs non-lagged versions
                    behavior_labels_nonlagged = behavior_labels_init(1+lag_length:end,:);
                    diff_behav_vectors = behavior_labels - behavior_labels_nonlagged;
                    perc_diff_before = round(length(find(diff_behav_vectors~=0))/length(diff_behav_vectors)*100);
                    fprintf('Percent difference between-time lagged and regular labels BEFORE selecting behaviors: %s \n', num2str(perc_diff_before));
                end


                %% Select behaviors to decode

                %Compute freq of behavior for the session
                behav_freq_table = tabulate(behavior_labels);
                behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

                % Select behaviors with a minimum # of occurrences
                behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
                behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
                behav = behav(behav~=find(matches(behav_categ,'Rowdy Room')));%excluding Rowdy room which is a source of confusion.
                behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
                behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.

                % OR select behaviors manually
                %behav = [4,6,7,8,17,23]; %manually select behaviors of interest
                %Behaviors selecte on 1 window size difference, to make sure we
                %have equivalent # categories across sessions

                %Only keep the behaviors of interest
                idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
                Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
                behavior_labels_final = behavior_labels(idx,:);%Same as above but in behavior labels
                tabulate(behavior_labels_final);

                %check the amount of labels that differ in the lag vs non-lagged
                %versions after selecting the behaviors of interest.
                if lag_length~=0
                    behavior_labels_final_nonlagged = behavior_labels_nonlagged(idx,:);
                    idx_diff_behav = find((behavior_labels_final - behavior_labels_final_nonlagged)~=0);
                    perc_diff_after = round(length(idx_diff_behav)/length(diff_behav_vectors)*100);
                    fprintf('Percent difference between-time lagged and regular labels AFTER selecting behaviors: %s \n', num2str(perc_diff_after));

                    %Only consider windows where the behaviors of the lagged and
                    %non-lagged versions do not overlap
                    Spike_count_raster_final = Spike_count_raster_final(idx_diff_behav,:);%Only keep timepoints where the behaviors of interest occur in spiking data
                    behavior_labels_final = behavior_labels_final(idx_diff_behav,:);%Same as above but in behavior labels
                    behav_freq_table = tabulate(behavior_labels_final);

                    %Again, only consider behaviors we have more than X
                    %occurrences of
                    behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);
                    idx = find(ismember(behavior_labels_final,behav)); %find the indices of the behaviors considered
                    Spike_count_raster_final = Spike_count_raster_final(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
                    behavior_labels_final = behavior_labels_final(idx,:);
                    tabulate(behavior_labels_final);

                    %Print behaviors selected
                    behavs_eval = behav_categ(behav);
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                    fprintf('Behaviors evaluated are: %s \n', behavs_eval);
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


                end


                %% Run SVM over multiple iterations
               
                disp('Start running SVM...')
                for iter = 1:num_iter

                    %                 clearvars -except savePath behav_categ behavior_labels_final behavior_labels_shifted Spike_count_raster_final...
                    %                     Spike_count_raster_shifted num_iter iter hitrate hitrate_shuffled C C_shuffled temp_resolution...
                    %                     channel_flag filePath chan temp mean_hitrate sd_hitrate mean_hitrate_shuffled C_table behavs_eval behav is_mac min_occurrences

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
                    minNumTrials = min_occurrences; %30; %find the minimum one %CT change to have 30 of each class
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

                channel = char(channel_flag);
                disp('****************************************************************************')
                disp([num2str(lag_length) 's window lag, channels: ' channel '. Prediction:' num2str(predict) '. DONE'])
                disp(['Mean hitrate: ' num2str(round(mean(hitrate)*100,2)) '%, suffled:' num2str(round(mean(hitrate_shuffled)*100,2)) '%'])
                disp('****************************************************************************')
                pause(1)

                mean_hitrate{s}(chan, lag, pred) = mean(hitrate);
                sd_hitrate{s}(chan, lag, pred) = std(hitrate);
                mean_hitrate_shuffled{s}(chan, lag, pred) = mean(hitrate_shuffled);
                sd_hitrate_shuffled = std(hitrate_shuffled);

                C_concat=cat(3,C{:});
                confusion_mat_avg=round(mean(C_concat,3)*100);
                rowNames{s} = {labels_id{:,2}}; colNames = {labels_id{:,2}};
                C_table{s,chan,lag,pred} = array2table(confusion_mat_avg,'RowNames',rowNames{s},'VariableNames',colNames);

                lag = lag+1;

                clear labels_id
            end %end of lag length loop

            pred = pred + 1;
        end %End of past or future events to be decoded

        chan = chan +1;
    end %End of channel loop
    %     temp = temp+1;

% mean_hitrate = permute( mean_hitrate , [1,3,2] ); sd_hitrate = permute( sd_hitrate , [1,3,2] );
% mean_hitrate_shuffled = permute( mean_hitrate_shuffled , [1,3,2] );

% % colNames = ["0sec","2sec","5sec","10sec","20sec"]; rowNames = ["vlPFC","TEO","all"];
% % result_hitrate_future = array2table(mean_hitrate(:,:,1),'RowNames',rowNames,'VariableNames',colNames)
% % result_hitrate_future_shuffled = array2table(mean_hitrate_shuffled(:,:,1),'RowNames',rowNames,'VariableNames',colNames)
% % result_sdhitrate_future = array2table(sd_hitrate(:,:,1),'RowNames',rowNames,'VariableNames',colNames);
% % result_hitrate_past = array2table(mean_hitrate(:,:,2),'RowNames',rowNames,'VariableNames',colNames)
% % result_sdhitrate_past = array2table(sd_hitrate(:,:,2),'RowNames',rowNames,'VariableNames',colNames);
% % 
% % % save([savePath '\SVM_results_social_context.mat'], 'mean_hitrate', 'sd_hitrate', 'C_table', 'result_hitrate', 'result_sdhitrate', 'behavs_eval')
% % % writetable(result_hitrate,[savePath '\SVM_results_social_context.csv'],'WriteRowNames',true,'WriteVariableNames',true);
% % save([savePath '\SVM_timelagged_results_' num2str(length(behav)) 'behav.mat'], 'mean_hitrate', 'mean_hitrate_shuffled','sd_hitrate', 'C_table', 'behavs_eval')
% % writetable(result_hitrate_future,[savePath '\SVM_results_' num2str(length(behav)) 'future_behav.csv'],'WriteRowNames',true,'WriteVariableNames',true);
% % writetable(result_hitrate_past,[savePath '\SVM_results_' num2str(length(behav)) 'past_behav.csv'],'WriteRowNames',true,'WriteVariableNames',true);

%Plotting results for decoding future behaviors
figure; hold on; set(gcf,'Position',[150 250 1000 500])
cmap = hsv(size(mean_hitrate,2));
for b = 1:size(mean_hitrate{s},1) %for all channels
    y = mean_hitrate{s}(b,:,1)./mean_hitrate_shuffled{s}(b,:,1);
    std_dev = sd_hitrate{s}(b,:,1);
    errorbar(y,std_dev,'s','MarkerSize',10,...
        'MarkerEdgeColor',cmap(b,:),'MarkerFaceColor',cmap(b,:))
    %plot(x,y,'Color','k')
end
leg = legend(['vlPFC',"TEO","All"],'Location','best');
title(leg,'Brain area')
chance_level = 1;
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0.8 1 2 3 4 5 5.2]); xlim([0.8 5.2]); ylim([0 10])
xticklabels({'','0sec','2sec','5sec','10sec','20sec',''})
ax = gca;
ax.FontSize = 14;
ylabel('Deconding accuracy relative to chance','FontSize', 18); xlabel('Time lag length (#windows)','FontSize', 18)
title('Decoding accuracy for future behaviors','FontSize', 20)
cd(savePath)
saveas(gcf,['SVM_results_FutureBehav.png'])
close all

%Plotting results for decoding past behaviors
figure; hold on; set(gcf,'Position',[150 250 1000 500])
cmap = hsv(size(mean_hitrate{s},2));
for b = 1:size(mean_hitrate{s},1)
    y = mean_hitrate{s}(b,:,2)./mean_hitrate_shuffled{s}(b,:,2);
    std_dev = sd_hitrate{s}(b,:,2);
    errorbar(y,std_dev,'s','MarkerSize',10,...
        'MarkerEdgeColor',cmap(b,:),'MarkerFaceColor',cmap(b,:))
    %plot(x,y,'Color','k')
end
leg = legend(['vlPFC',"TEO","All"],'Location','best');
title(leg,'Brain area')
chance_level = 1;
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0.8 1 2 3 4 5 5.2]); xlim([0.8 5.2]); ylim([0 10])
xticklabels({'','0sec','2sec','5sec','10sec','20sec',''})
ax = gca;
ax.FontSize = 14;
ylabel('Deconding accuracy','FontSize', 18); xlabel('Time lag length (#windows)','FontSize', 18)
title('Decoding accuracy for past behaviors','FontSize', 20)
cd(savePath)
saveas(gcf,['SVM_results_PastBehav.png'])
close all

end %End of session loop