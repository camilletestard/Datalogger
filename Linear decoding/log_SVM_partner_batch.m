%% Log_SVM_partner
%% Run a linear decoder on a the neural activity for the partner behaviors

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range=[1,2,11,12];

%Set parameters
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
onlysingle=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 100;%Number of SVM iterations

s=1;

for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];

    chan = 1;

    for channel_flag = ["vlPFC", "TEO", "all"]

        %Get data with specified temporal resolution and channels
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac);
        %filePath is the experimental data path
        %Temp_resolution is the temporal resolution at which we would like to
        %analyze the data
        %Channel_flag specifies with channels to include: only TEO array, only
        %vlPFC array or all channels
        %is_mac is whether a mac or a pc is being used
        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';
        behavior_labels_subject_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        behavior_labels_partner_init = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
        block_labels = cell2mat({labels{:,10}}'); %Extract block info

        %SANITY CHECK: Compute overlap between partner behavior and subject behavior
        perc_overlap_allbehav = length(find(behavior_labels_subject_init == behavior_labels_partner_init))/length(behavior_labels_subject_init);
        overlap = tabulate(behavior_labels_subject_init(find(behavior_labels_subject_init == behavior_labels_partner_init)));
        subject_behav = tabulate(behavior_labels_subject_init); partner_behav = tabulate(behavior_labels_partner_init);
        perc_overlap_perbehav = [behav_categ', overlap(:,2), overlap(:,2)./subject_behav(:,2)*100];
        %Important notes:
        %1. There are some discrepancies when comparing the partner and
        %subject labels. Most notably in proximity, but also RR, HIP, HIS,
        %SP, SS which should be the same.
        %2. I just realized that depending on what other behavior is
        %co-occurring, the label can change. This may explain discrepancies
        %in RR and proximity
        %3. it is often the case that self-groom events co-occur for Amos and Lovelace. I expect
        %this to be the case for foraging in Hooke-pair.


        behavior_labels = behavior_labels_partner_init;

        %% Select behaviors to decode

        %Compute freq of behavior for the session
        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

        % Select behaviors with a minimum # of occurrences
        min_occurrences = 50;
        behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
        behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rowdy Room')));%excluding rest which is a source of confusion.

        % OR select non-reciprocal behaviors (particularly important
        % for decoding partner behavior)
        behav = setdiff(behav, reciprocal_set);

        % OR select behaviors manually
        behav = [4,5,17,23,25];%manually select behaviors of interest
        %Select behaviors manually to ensure that the same
        %behaviors are considered for the partner and subject comparisons.
        %This list could change from session to session.. I'll have to
        %think of a way to automatize this.

        %Print behaviors selected
        behavs_eval = behav_categ(behav);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        fprintf('Behaviors evaluated are: %s \n', behavs_eval);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        %Only keep the behaviors of interest
        idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx,:);%Same as above but in behavior labels
        tabulate(behavior_labels_final);

        %check the amount of labels that differ in the partner vs. subject
        %labels after selecting the behaviors of interest.
        if partner==1
            subject_behav_after_selection = behavior_labels_subject_init(idx);
            partner_behav_after_selection = behavior_labels_partner_init(idx);
            block_after_selection = block_labels(idx);

            overlap_partner_subject_after_selection = length(find((partner_behav_after_selection - subject_behav_after_selection)==0))/length(idx);
            fprintf('Percent overlap between subject and partner labels AFTER selecting behaviors: %s \n', num2str(overlap_partner_subject_after_selection));
            alone_block_obs = length(find(block_after_selection==3))/length(idx);
            fprintf('Percent observations in "alone" block AFTER selecting behaviors: %s \n', num2str(alone_block_obs));

            %Only consider windows where the behaviors of subject and
            %partner do not overlap
            diff_idx = find((partner_behav_after_selection - subject_behav_after_selection)~=0); %find the indices where subject and partner behavior do not overlap
            Spike_count_raster_final = Spike_count_raster_final(diff_idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
            behavior_labels_final = behavior_labels_final(diff_idx,:);%Same as above but in behavior labels
            tabulate(behavior_labels_final);
            %IMPORTANT NOTE: There will be much less instances

            %             %Only consider behaviors during the alone block
            %             block_after_selection_overlap_out = block_after_selection(diff_idx);
            %             alone_idx = find(block_after_selection_overlap_out==3);
            %             Spike_count_raster_final = Spike_count_raster_final(alone_idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
            %             behavior_labels_final = behavior_labels_final(alone_idx,:);%Same as above but in behavior labels
            %             tabulate(behavior_labels_final);
        end


        %% Run SVM over multiple iterations
        num_iter = 500;

        disp('Start running SVM...')
        for iter = 1:num_iter

            %             clearvars -except savePath behav_categ behavior_labels_final behavior_labels_shifted Spike_count_raster_final...
            %                 Spike_count_raster_shifted num_iter iter hitrate hitrate_shuffled C C_shuffled temp_resolution...
            %                 channel_flag filePath chan temp mean_hitrate sd_hitrate mean_hitrate_shuffled C_table behavs_eval behav is_mac min_occurrences

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
            minNumTrials = 30;%min(num_trials); %find the minimum one %CT change to have 30 of each class
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
        disp(['ID type ' num2str(partner) ', channels: ' channel '. DONE'])
        disp('****************************************************************************')

        mean_hitrate(chan, id) = mean(hitrate)
        sd_hitrate(chan, id) = std(hitrate);
        mean_hitrate_shuffled(chan, id) = mean(hitrate_shuffled)
        sd_hitrate_shuffled = std(hitrate_shuffled);

        C_concat=cat(3,C{:});
        confusion_mat_avg=round(mean(C_concat,3)*100);
        rowNames = {labels_id{:,2}}; colNames = {labels_id{:,2}};
        C_table{chan, id} = array2table(confusion_mat_avg,'RowNames',rowNames,'VariableNames',colNames);

        id = id+1;

    
    end% End of channel for loop
    chan = chan +1;
end% End of session for loop


colNames = ["subject", "partner"]; rowNames = ["vlPFC","TEO","all"];
result_hitrate = array2table(mean_hitrate,'RowNames',rowNames,'VariableNames',colNames)
result_sdhitrate = array2table(sd_hitrate,'RowNames',rowNames,'VariableNames',colNames);

% save([savePath '\SVM_results_social_context.mat'], 'mean_hitrate', 'sd_hitrate', 'C_table', 'result_hitrate', 'result_sdhitrate', 'behavs_eval')
% writetable(result_hitrate,[savePath '\SVM_results_social_context.csv'],'WriteRowNames',true,'WriteVariableNames',true);
save([savePath '\SVM_results_' num2str(length(behav)) 'behav.mat'], 'mean_hitrate', 'sd_hitrate', 'C_table', 'result_hitrate', 'result_sdhitrate', 'behavs_eval')
writetable(result_hitrate,[savePath '\SVM_results_' num2str(length(behav)) 'behav.csv'],'WriteRowNames',true,'WriteVariableNames',true);

%Plotting results for multiple behaviors at 1sec and lower resolution
figure; hold on; set(gcf,'Position',[150 250 1000 500])
cmap = hsv(size(mean_hitrate,1));
for b = 1:size(mean_hitrate,1)
    y = mean_hitrate(b,:);
    std_dev = sd_hitrate(b,:);
    errorbar(y,std_dev,'s','MarkerSize',10,...
        'MarkerEdgeColor',cmap(b,:),'MarkerFaceColor',cmap(b,:))
    %plot(x,y,'Color','k')
end
leg = legend([rowNames, 'Chance']);
title(leg,'Brain area')
chance_level = 1/length(behav);
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0.8 1 2 2.2]); xlim([0.8 2.2]); ylim([0 1])
xticklabels({'','subject','partner',''})
ax = gca;
ax.FontSize = 14;
ylabel('Deconding accuracy','FontSize', 18); xlabel('ID','FontSize', 18)
title('Decoding accuracy for subject vs. partner behavior','FontSize', 20)
cd(savePath)
saveas(gcf,['SVM_results_SubjectVsPartner.png'])
close all
