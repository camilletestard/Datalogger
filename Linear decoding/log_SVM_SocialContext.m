%% Log_SVM
%% Run a linear decoder on a the neural activity over all or a subset of behaviors
% This script allows to deocde a combination of factors:
% 1. Subject behavior
% 2. Partner's behavior
% 3. Predict future behavior (with specified lag time in terms of windows)
% 4. Detect past behaviors (with specified lag time in terms of windows)
% 5. Behavior shifts (irrespective of what the shift is and the TYPE of shift)
% 6. Window size used fo decoding behaviors can be manipulated
% 7. Brain area

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

clearvars -except savePath filePath is_mac

%Set parameters
temp = 1; temp_resolution = 1;
chan = 1; channel_flag = "all";
randomsample=0;

for channel_flag = ["vlPFC", "TEO", "all"]

    %Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final, unit_count]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac);
    %filePath is the experimental data path
    %Temp_resolution is the temporal resolution at which we would like to
    %analyze the data
    %Channel_flag specifies with channels to include: only TEO array, only
    %vlPFC array or all channels
    %is_mac is whether a mac or a pc is being used
    disp('Data Loaded')

    Spike_count_raster = Spike_rasters';
    behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    block_labels = cell2mat({labels{:,11}}');

    %% Check which behaviors occur in different blocks

%     unq_behav = unique(behavior_labels);
%     behav_in_block = zeros(length(unq_behav), 3);
%     for b = 1:length(unq_behav)
%         for bl = 1:3
%             behav_in_block(b, bl) = length(intersect(find(behavior_labels == unq_behav(b)), find(block_labels==bl)));
%         end
%     end
% 
%     figure; hold on; set(gcf,'Position',[150 250 1500 500])
%     bar(behav_in_block, 'stacked')
%     xticks(1:length(unq_behav))
%     xticklabels(behav_categ(unq_behav));
%     xtickangle(45)
%     leg = legend(block_times.Behavior{:});
%     title(leg,'Block')

    %% Select behaviors to decode

    % Select behaviors which occur in multiple blocks and the
    %behav =[find(matches(behav_categ,'Groom Give')), find(matches(behav_categ,'Groom Receive'))]; %manually select behaviors of interest
    %behav =[find(matches(behav_categ,'Threat to subject')), find(matches(behav_categ,'Threat to partner'))];
    behav = 28;

    %Only keep the behaviors of interest
    idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
    Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
    behavior_labels_final = block_labels(idx,:);%Same as above but in behavior labels
    tabulate(behavior_labels_final);

    %% Run SVM over multiple iterations
    num_iter = 1000;

    disp('Start running SVM...')
    for iter = 1:num_iter

        %         clearvars -except savePath behav_categ behavior_labels_final behavior_labels_shifted Spike_count_raster_final...
        %             Spike_count_raster_shifted num_iter iter hitrate hitrate_shuffled C C_shuffled temp_resolution...
        %             channel_flag filePath chan temp mean_hitrate sd_hitrate mean_hitrate_shuffled C_table behavs_eval behav is_mac min_occurrences

        %subsample to match number of neurons across brain areas
        Labels = behavior_labels_final;
        if randomsample==1 && channel_flag~="all"
            Input_matrix = Spike_count_raster_final(:,randsample(unit_count(chan), min(unit_count)));
        else
            Input_matrix = Spike_count_raster_final;
        end


        %Balance number of trials per class
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
        minNumTrials = 30;%min(num_trials); %30; %find the minimum one %CT change to have 30 of each class
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
    disp([num2str(1000/temp_resolution) 'msec resolution, channels: ' channel '. DONE'])
    disp('****************************************************************************')

    mean_hitrate(temp, chan) = mean(hitrate)
    sd_hitrate(temp, chan) = std(hitrate);
    mean_hitrate_shuffled(temp, chan) = mean(hitrate_shuffled)
    sd_hitrate_shuffled = std(hitrate_shuffled);

    C_concat=cat(3,C{:});
    confusion_mat_avg=round(mean(C_concat,3)*100);
    rowNames = {labels_id{:,2}}; colNames = {labels_id{:,2}};
    C_table{temp, chan} = array2table(confusion_mat_avg,'RowNames',rowNames,'VariableNames',colNames);

    chan = chan +1;

    clearvars -except is_mac savePath mean_hitrate sd_hitrate mean_hitrate_shuffled C_table temp_resolution channel_flag filePath chan temp behavs_eval behav randomsample onlysingle
end
% temp = temp+1;

cd(savePath)

%Plot confusion matrices
vlpfc = table2array(C_table{:,1}); TEO = table2array(C_table{1,2}); all_units = table2array(C_table{1,3});
D = vlpfc-TEO;
CustomAxisLabels = string(C_table{1,1}.Properties.VariableNames)
figure; pfc = heatmap(vlpfc,'Colormap', jet); pfc.XDisplayLabels = CustomAxisLabels; pfc.YDisplayLabels = CustomAxisLabels; title(pfc,'vlpfc confusion matrix'); caxis([0, 100]);
figure; teo = heatmap(TEO,'Colormap', jet); teo.XDisplayLabels = CustomAxisLabels; teo.YDisplayLabels = CustomAxisLabels; title(teo,'TEO confusion matrix'); caxis([0, 100]);
figure; all = heatmap(all_units,'Colormap', jet); all.XDisplayLabels = CustomAxisLabels; all.YDisplayLabels = CustomAxisLabels; title(all,'All units confusion matrix'); caxis([0, 100]);
figure; h = heatmap(D,'Colormap', jet); h.XDisplayLabels = CustomAxisLabels; h.YDisplayLabels = CustomAxisLabels; title(h,'vlpfc-TEO confusion matrix');caxis([-50, 50]);

saveas(pfc,['SVM_results_' num2str(length(behavs_eval)) 'behav_NOsubsample_unique_confmatPFC.png'])
saveas(teo,['SVM_results_' num2str(length(behavs_eval)) 'behav_NOsubsample_unique_confmatTEO.png'])
saveas(all,['SVM_results_' num2str(length(behavs_eval)) 'behav_NOsubsample_unique_confmatALL.png'])
saveas(h,['SVM_results_' num2str(length(behavs_eval)) 'behav_NOsubsample_unique_confmatPFC-TEO.png'])

close all

rowNames = ["1sec"]; colNames = ["vlPFC","TEO","all"];
result_hitrate = array2table(mean_hitrate,'RowNames',rowNames,'VariableNames',colNames)
result_hitrate_shuffled = array2table(mean_hitrate_shuffled,'RowNames',rowNames,'VariableNames',colNames)
result_sdhitrate = array2table(sd_hitrate,'RowNames',rowNames,'VariableNames',colNames)

% save([savePath '\SVM_results_social_context.mat'], 'mean_hitrate', 'sd_hitrate', 'C_table', 'result_hitrate', 'result_sdhitrate', 'behavs_eval')
% writetable(result_hitrate,[savePath '\SVM_results_social_context.csv'],'WriteRowNames',true,'WriteVariableNames',true);
save([savePath '\SVM_results_' num2str(length(behav)) 'behav_subsample_unique.mat'], 'mean_hitrate', 'sd_hitrate', 'C_table', 'result_hitrate', 'result_sdhitrate', 'behavs_eval')
writetable(result_hitrate,[savePath '\SVM_results_' num2str(length(behav)) 'behav_subsample_unique.csv'],'WriteRowNames',true,'WriteVariableNames',true);

%Plotting results decoding accuracy for all behaviors at 1sec and lower resolution
figure; hold on; set(gcf,'Position',[150 250 700 500])
cmap = hsv(size(mean_hitrate,1));
for b = 1:size(mean_hitrate,1)
    y = mean_hitrate(b,:);
    std_dev = sd_hitrate(b,:);
    errorbar(y,std_dev,'s','MarkerSize',10,...
        'MarkerEdgeColor',cmap(b,:),'MarkerFaceColor',cmap(b,:))
    %plot(x,y,'Color','k')
end
%leg = legend([rowNames, 'Chance']);
%title(leg,'Window size')
chance_level = 1/length(behav);
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 1])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14;
ylabel('Deconding accuracy','FontSize', 18); xlabel('Brain area','FontSize', 18)
title('Decoding accuracy for subject current behavioral states (unique behaviors, subsample neurons)','FontSize', 14)

cd(savePath)
saveas(gcf,['SVM_results_' num2str(length(behavs_eval)) 'behav_subsample_unique.png'])
close all

