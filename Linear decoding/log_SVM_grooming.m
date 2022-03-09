%% Log_SVM_grooming
%% Run a linear decoder on a the neural activity for different grooming contexts
% This script allows to deocde grooming:
% 1. Start vs. end
% 2. Post-threat or not
% 3. Reciprocated or not
% 4. Initiated or not
%Camille Testard, March 2022


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
with_NC=1;
isolatedOnly=0;

%for channel_flag = ["vlPFC", "TEO", "all"]

%Get data with specified temporal resolution and channels
%[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final, unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final, unit_count, groom_labels_all]= log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);

%filePath is the experimental data path
%Temp_resolution is the temporal resolution at which we would like to
%analyze the data
%Channel_flag specifies with channels to include: only TEO array, only
%vlPFC array or all channels
%is_mac is whether a mac or a pc is being used
disp('Data Loaded')

Spike_count_raster = Spike_rasters';
behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject

%% Add labelling for grooming
%Beginning vs. end of grooming bout
%Grooming post-threat vs not.
%Grooming after reciprocation or not
%Groom received after grm Prsnt or not

groom_categ_label = {'Star.vs.end', 'Post-threat','Reciprocated','Initiated'};

b = 1;
for behav = [7,8]
    groom_behav={'Give','Receive'};

    for groom_categ = 1:4
        %Note: after threat, almost always groom RECEIVE, not given by
        %subject. Also, grooming after groom present is only for groom
        %RECEIVE.

        groom_labels = groom_labels_all(:,groom_categ+1);


        %% Select behaviors to decode

        % Select behaviors with a minimum # of occurrences
        %behav =[find(matches(behav_categ,'Groom Give'))]; %find(matches(behav_categ,'Groom Give')),


        %Only keep the behaviors of interest
        idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = groom_labels(idx,:);%Same as above but in behavior labels

        %If groom label is start vs. end
        if groom_categ==1
            idx_epoch = find(~ismember(behavior_labels_final,3)); %remove middle chunk
            Spike_count_raster_final = Spike_count_raster_final(idx_epoch,:);%Only keep timepoints where the behaviors of interest occur in spiking data
            behavior_labels_final = behavior_labels_final(idx_epoch,:);%Same as above but in behavior labels
        end

        behav_size=tabulate(behavior_labels_final);
        disp('########################')
        tabulate(behavior_labels_final)
        disp('########################')

        channel = char(channel_flag);
        disp('****************************************************************************')
        disp(['Groom categ: ' groom_categ_label{groom_categ} ', Channels: ' channel ', Behavior: Groom ' groom_behav{b}])
        disp('****************************************************************************')

        %pause(5)
        if all(behav_size(:,2)>=30) && length(behav_size(:,2))>1
            %% Run SVM over multiple iterations
            num_iter = 100;

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
                if min(num_trials)<50 %If there are less than 50 instances for a given behavior
                    minNumTrials = min(num_trials); %use the minimum # of instances
                else
                    minNumTrials = 50; %impose 50 occurrences per category
                end
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

                %disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
            end

            mean_hitrate(chan,b,groom_categ) = mean(hitrate);
            sd_hitrate(chan,b,groom_categ) = std(hitrate);
            mean_hitrate_shuffled(chan,b,groom_categ) = mean(hitrate_shuffled);
            sd_hitrate_shuffled(chan,b,groom_categ) = std(hitrate_shuffled);

            %             C_concat=cat(3,C{:});
            %             confusion_mat_avg=round(mean(C_concat,3)*100);
            %             rowNames = {labels_id{:,2}}; colNames = {labels_id{:,2}};
            %             C_table{chan,b,groom_categ} = array2table(confusion_mat_avg,'RowNames',rowNames,'VariableNames',colNames);
            %
        else
            mean_hitrate(chan,b,groom_categ) = nan;
            sd_hitrate(chan,b,groom_categ) = nan;
        end
    end
    b = b+1;
end
chan = chan +1;
clearvars -except is_mac savePath mean_hitrate sd_hitrate mean_hitrate_shuffled sd_hitrate_shuffled C_table temp_resolution channel_flag filePath chan temp behavs_eval behav randomsample onlysingle b groom_labels_all isolatedOnly with_NC groom_categ_label
%end

squeeze(mean_hitrate)
squeeze(sd_hitrate)
(squeeze(mean_hitrate) - squeeze(mean_hitrate_shuffled))./squeeze(mean_hitrate_shuffled)

cd(savePath)

rowNames = ["1sec"]; colNames = ["vlPFC","TEO","all"];
result_hitrate = array2table(mean_hitrate,'RowNames',rowNames,'VariableNames',colNames)
result_hitrate_shuffled = array2table(mean_hitrate_shuffled,'RowNames',rowNames,'VariableNames',colNames)
result_sdhitrate = array2table(sd_hitrate,'RowNames',rowNames,'VariableNames',colNames)

% save([savePath '\SVM_results_social_context.mat'], 'mean_hitrate', 'sd_hitrate', 'C_table', 'result_hitrate', 'result_sdhitrate', 'behavs_eval')
% writetable(result_hitrate,[savePath '\SVM_results_social_context.csv'],'WriteRowNames',true,'WriteVariableNames',true);
save([savePath '\SVM_results_' num2str(length(behav)) 'behav_NOsubsample_unique.mat'], 'mean_hitrate', 'sd_hitrate', 'C_table', 'result_hitrate', 'result_sdhitrate', 'behavs_eval')
writetable(result_hitrate,[savePath '\SVM_results_' num2str(length(behav)) 'behav_NOsubsample_unique.csv'],'WriteRowNames',true,'WriteVariableNames',true);

groom_give_hitrate = reshape(mean_hitrate(:,1,:), [3,4]);
groom_give_hitrate_shuffled = reshape(mean_hitrate_shuffled(:,1,:), [3,4]);
groom_receive_hitrate = reshape(mean_hitrate(:,2,:), [3,4]);
groom_receive_hitrate_shuffled = reshape(mean_hitrate_shuffled(:,2,:), [3,4]);

%Plotting results decoding accuracy for grooming give
figure; hold on; set(gcf,'Position',[150 250 700 500])
cmap = hsv(size(groom_give_hitrate,1));
for b = 1:size(groom_give_hitrate,1)
    y = groom_give_hitrate(b,[1,3]);
    %     std_dev = sd_hitrate(b,:);
    %     errorbar(y,std_dev,'s','MarkerSize',10,...
    %         'MarkerEdgeColor',cmap(b,:),'MarkerFaceColor',cmap(b,:))
    scatter(1:2,y,40,'filled','Color',cmap(b,:))
end
leg = legend(["vlPFC", "TEO","All"]);
title(leg,'Brain Area')
chance_level = 1/2;
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0.8 1 2 2.2]); xlim([0.8 2.2]); ylim([0.4 1])
xticklabels({'','Start vs. end','Reciprocal or not',''})
ax = gca;
ax.FontSize = 14;
ylabel('Deconding accuracy','FontSize', 18); xlabel('Grooming context','FontSize', 18)
title('Decoding accuracy for the context of grooming give','FontSize', 14)

cd(savePath)
saveas(gcf,['Decoding grooming given context.png'])

%Plotting results decoding accuracy for grooming receive
figure; hold on; set(gcf,'Position',[150 250 700 500])
cmap = hsv(size(groom_receive_hitrate,1));
for b = 1:size(groom_receive_hitrate,1)
    y = groom_receive_hitrate(b,:);
    y_shuffle = groom_receive_hitrate_shuffled(b,:);
    %     std_dev = sd_hitrate(b,:);
    %     errorbar(y,std_dev,'s','MarkerSize',10,...
    %         'MarkerEdgeColor',cmap(b,:),'MarkerFaceColor',cmap(b,:))
    scatter(1:4,y,40,'filled','Color',cmap(b,:))
    %     scatter(1:4,y_shuffle,40,'filled','Color','k')
end
leg = legend(["vlPFC", "TEO","All"]);
title(leg,'Brain Area')
chance_level = 1/2;
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0.8 1 2 3 4 4.2]); xlim([0.8 4.2]); ylim([0.4 1])
xticklabels({'','Start vs. end','Post-threat or not','Reciprocal or not', 'Subject-initiated vs. not',''})
ax = gca;
ax.FontSize = 14;
ylabel('Deconding accuracy','FontSize', 18); xlabel('Grooming context','FontSize', 18)
title('Decoding accuracy for the context of grooming received','FontSize', 14)

cd(savePath)
saveas(gcf,['Decoding grooming received context.png'])

close all

