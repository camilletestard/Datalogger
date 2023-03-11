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
session_range_no_partner=[1:6,11:13,15:16,18];
session_range_with_partner=[1:6,11:13,15:16,18];


%Set parameters
with_partner =0;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 10;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null 
simplify=0;%lump similar behavioral categories together to increase sample size.
threat_precedence =0;
exclude_sq=1;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];

    chan = 1;
    for channel_flag = ["vlPFC", "TEO", "all"]


        %% Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);
        end

        disp('Data Loaded')

        %Raw data
        Spike_count_raster = Spike_rasters';

        %Low-pass filter data
        %Spike_count_raster = lowpass(Spike_rasters',0.05,1);

        behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

        co_occurrence = cell2mat({labels{:,5}}');

        if null
            %Simulate fake labels
            [sim_behav] = GenSimBehavior(behavior_labels,behav_categ, temp_resolution);
            behavior_labels = sim_behav;
        end

        if unq_behav==1%Select epochs where only one behavior happens at any given time (i.e. no co-occurrence).
            idx_single = find(co_occurrence<4); %(no co-occurrence, or with RR or with proximity)
            Spike_count_raster = Spike_count_raster(idx_single,:);
            behavior_labels = behavior_labels(idx_single,:);
        end

        if simplify
            %Simplify behavioral catagories
            %Lump all aggressive interactions together
            behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
            behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");

            %Lump all travel together
            behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
            behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

            %Lump foraging and drinking
            behavior_labels(behavior_labels==find(behav_categ=="Drinking"))=find(behav_categ=="Foraging");

            %Lump grooming together
            behavior_labels(behavior_labels==find(behav_categ=="Getting groomed"))=find(behav_categ=="Groom partner");
        end

        %% Select behaviors to decode

        %Compute freq of behavior for the session
        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

        % Select behaviors with a minimum # of occurrences
        min_occurrences = 30;
        behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences

        % Remove behaviors that are ill-defined
        behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Other monkeys vocalize')));%excluding scratch which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Vocalization')));%excluding scratch which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rowdy Room')));%excluding scratch which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.


        % OR select behaviors manually
        %behav = [5,18];%[4,5,17];% [7,8]%[5,7:10,21];%[4,5,7:10];%[4:8,17]; %[1:6,9:11,16,17]; %manually select behaviors of interest

        %Print behaviors selected
        behavs_eval = behav_categ(behav);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        fprintf('Behaviors evaluated are: %s \n', behavs_eval);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        %Only keep the behaviors of interest
        idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx,:);%Same as above but in behavior labels
        tabulate(removecats(categorical(behavior_labels_final)));


        %% Run SVM over multiple iterations

        disp('Start running SVM...')
        for iter = 1:num_iter


            %subsample to match number of neurons across brain areas
            Labels = behavior_labels_final;
            if randomsample==1
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

            % Run svm
            [hitrate(iter), C{iter}] = log_SVM_basic_function(Input_matrix, Labels, 5, 0, 0);
            [hitrate_shuffled(iter), C_shuffled{iter}] = log_SVM_basic_function(Input_matrix, Labels_shuffled, 5, 0, 0);

            if mod(iter,10)==1
                disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
            end
        end %End of SVM for loop

        channel = char(channel_flag);
        disp('****************************************************************************')
        disp([num2str(1000/temp_resolution) 'msec resolution, channels: ' channel '. DONE'])
        disp('****************************************************************************')

        mean_hitrate{s}(chan) = mean(hitrate)
        sd_hitrate{s}(chan) = std(hitrate);
        mean_hitrate_shuffled{s}(chan) = mean(hitrate_shuffled)
        sd_hitrate_shuffled = std(hitrate_shuffled);

        C_concat=cat(3,C{:}); %Get confusion matrix
        confusion_mat_avg{s, chan}=round(mean(C_concat,3)*100); %Average over SVM iterations
        rowNames{s} = {labels_id{:,2}}; colNames{s} = {labels_id{:,2}}; %Get behavior names
        C_table{s, chan} = array2table(confusion_mat_avg{s, chan},'RowNames',rowNames{s},'VariableNames',colNames{s});

        chan = chan +1;
        clear labels_id

    end %end of channel for loop

    cd(savePath)

    %% Plot confusion matrices

% %     figure; set(gcf,'Position',[150 250 1500 500]);
% %     subplot(1,2,1)
% %     hp=heatmap(confusion_mat_avg{s, 1},'Colormap', jet); caxis([0 100]);
% %     hp.XDisplayLabels = rowNames{s}; hp.YDisplayLabels = colNames{s}; title(['vlPFC'])
% %     ylabel('Real'); xlabel('Predicted')
% %     ax = gca;
% %     ax.FontSize = 14;
% %     
% %     subplot(1,2,2)
% %     hp=heatmap(confusion_mat_avg{s, 2},'Colormap', jet); caxis([0 100]);
% %     hp.XDisplayLabels = rowNames{s}; hp.YDisplayLabels = colNames{s}; title(['TEO'])
% %     ylabel('Real'); xlabel('Predicted')
% %     ax = gca;
% %     ax.FontSize = 14;
% %     saveas(gcf,['ConfusionMatrix.pdf'])


%     vlpfc = table2array(C_table{s,1,1}); TEO = table2array(C_table{s,1,2}); all_units = table2array(C_table{s,1,3});
%     D{s} = vlpfc-TEO;
%     CustomAxisLabels = string(C_table{s,1,1}.Properties.VariableNames); figure; set(gcf,'Position',[150 250 1500 800]);
%     subplot(2,2,1); pfc = heatmap(vlpfc,'Colormap', hot); pfc.XDisplayLabels = CustomAxisLabels; pfc.YDisplayLabels = CustomAxisLabels; title(pfc,'vlpfc confusion matrix'); caxis([0, 100]);
%     subplot(2,2,2); teo = heatmap(TEO,'Colormap', hot); teo.XDisplayLabels = CustomAxisLabels; teo.YDisplayLabels = CustomAxisLabels; title(teo,'TEO confusion matrix'); caxis([0, 100]);
%     subplot(2,2,3); all = heatmap(all_units,'Colormap', hot); all.XDisplayLabels = CustomAxisLabels; all.YDisplayLabels = CustomAxisLabels; title(all,'All units confusion matrix'); caxis([0, 100]);
%     subplot(2,2,4); h = heatmap(D{s},'Colormap', cool); h.XDisplayLabels = CustomAxisLabels; h.YDisplayLabels = CustomAxisLabels; title(h,'vlpfc-TEO confusion matrix');caxis([-50, 50]);
    %close all

    %% Plot session result

%     rowNames = ["1sec"]; colNames = ["vlPFC","TEO","all"];
%     result_hitrate = array2table(mean_hitrate{s},'RowNames',rowNames,'VariableNames',colNames);
%     result_hitrate_shuffled = array2table(mean_hitrate_shuffled{s},'RowNames',rowNames,'VariableNames',colNames);
%     result_sdhitrate = array2table(sd_hitrate{s},'RowNames',rowNames,'VariableNames',colNames);

    %     save([savePath '\SVM_results_' num2str(length(behav)) 'behav_NOsubsample_unique.mat'], 'mean_hitrate', 'sd_hitrate', 'C_table', 'result_hitrate', 'result_sdhitrate', 'behavs_eval')
    %     writetable(result_hitrate,[savePath '\SVM_results_' num2str(length(behav)) 'behav_NOsubsample_unique.csv'],'WriteRowNames',true,'WriteVariableNames',true);

% % %     %Plotting results decoding accuracy for all behaviors at 1sec and lower resolution
% % %     figure; hold on; set(gcf,'Position',[150 250 700 500])
% % %     cmap = hsv(size(mean_hitrate,1));
% % %     for b = 1:size(mean_hitrate,1)
% % %         y = mean_hitrate{s}(b,:);
% % %         std_dev = sd_hitrate{s}(b,:);
% % %         errorbar(y,std_dev,'s','MarkerSize',10,...
% % %             'MarkerEdgeColor',cmap(b,:),'MarkerFaceColor',cmap(b,:))
% % %         %plot(x,y,'Color','k')
% % %     end
% % %     chance_level = 1/length(behav);
% % %     yline(chance_level,'--','Chance level', 'FontSize',16)
% % %     xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 1])
% % %     xticklabels({'','vlPFC','TEO','all',''})
% % %     ax = gca;
% % %     ax.FontSize = 14;
% % %     ylabel('Deconding accuracy','FontSize', 18); xlabel('Brain area','FontSize', 18)
% % %     title('Decoding accuracy for subject current behavioral states','FontSize', 14)


end %End of session for loop

%% Plot all sessions results

%Change savePath for all session results folder:
cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);
save('SVM_results_subjectBehav.mat', "mean_hitrate","sd_hitrate","mean_hitrate_shuffled","behav","a_sessions","h_sessions","behav_categ")
load('SVM_results_subjectBehav.mat')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bar plot decoding accuracy

figure; hold on
data = cell2mat(mean_hitrate');
data_shuffle = cell2mat(mean_hitrate_shuffled');
bp = bar([mean(data(:,:)); mean(data_shuffle(:,:))],'FaceAlpha',0.2);

sp1 = scatter(ones(size(data,1))*0.77,data(:,1), 'filled','b');
sp1 = scatter(ones(size(data,1)),data(:,2), 'filled','r');
sp1 = scatter(ones(size(data,1))*1.22,data(:,3), 'filled','y');

sp1 = scatter(ones(size(data_shuffle,1))*1.77,data_shuffle(:,1), 'filled','b');
sp1 = scatter(ones(size(data_shuffle,1))*2,data_shuffle(:,2), 'filled','r');
sp1 = scatter(ones(size(data_shuffle,1))*2.22,data_shuffle(:,3), 'filled','y');

legend(bp,{'vlPFC','TEO','all'},'Location','best')

ylabel('Decoding Accuracy'); ylim([0 1])
xticks([1 2]); xticklabels({'Real', 'Shuffled'}); xlim([0.25 2.75])
ax = gca;
ax.FontSize = 16;

saveas(gcf,['SVM_results_subjectBehav.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %Confusion matrices
% % % % % vlpfc = table2array(C_table{s,1,1}); TEO = table2array(C_table{s,1,2}); all_units = table2array(C_table{s,1,3});
% % % % % D{s} = vlpfc-TEO;
% % % % % CustomAxisLabels = string(C_table{s,1,1}.Properties.VariableNames); figure; set(gcf,'Position',[150 250 1500 800]);
% % % % % subplot(2,2,1); pfc = heatmap(vlpfc,'Colormap', hot); pfc.XDisplayLabels = CustomAxisLabels; pfc.YDisplayLabels = CustomAxisLabels; title(pfc,'vlpfc confusion matrix'); caxis([0, 100]);
% % % % % subplot(2,2,2); teo = heatmap(TEO,'Colormap', hot); teo.XDisplayLabels = CustomAxisLabels; teo.YDisplayLabels = CustomAxisLabels; title(teo,'TEO confusion matrix'); caxis([0, 100]);
% % % % % subplot(2,2,3); all = heatmap(all_units,'Colormap', hot); all.XDisplayLabels = CustomAxisLabels; all.YDisplayLabels = CustomAxisLabels; title(all,'All units confusion matrix'); caxis([0, 100]);
% % % % % subplot(2,2,4); h = heatmap(D,'Colormap', cool); h.XDisplayLabels = CustomAxisLabels; h.YDisplayLabels = CustomAxisLabels; title(h,'vlpfc-TEO confusion matrix');caxis([-50, 50]);
