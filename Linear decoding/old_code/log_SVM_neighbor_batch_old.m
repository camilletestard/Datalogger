%% Log_SVM_neighbor_batch
% Run a linear decoder on a the neural activity for the partner behaviors
% Batch version to run across sessions
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
session_range_with_neighbor=11;

%Set parameters
with_partner =2;% 0: no partner; 1: with partner; 2: with neighbor
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 50;%Number of SVM iterations
min_occurrences = 30;%Minimum number of occurrence per behavior
alone_block=0; %1: during alone block; 0:during paired blocks; anything else: all blocks.
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)
null=0;%Set whether we want the null 
simplify=0;%lump similar behavioral categories together to increase sample size.

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

s=11;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];

    chan = 1;

    for channel_flag = ["vlPFC", "TEO", "all"]

        
        %% Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        elseif with_partner ==2
            [Spike_rasters, labels, labels_partner, labels_neighbor, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function_neighbor(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        end

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';
        behavior_labels_subject_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        behavior_labels_partner_init = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
        behavior_labels_neighbor_init = cell2mat({labels_neighbor{:,3}}');
        block_labels = cell2mat({labels{:,11}}'); %Extract block info
        alone_block_id = find(strcmp(block_times{:,1},"Alone.block"));

        %SANITY CHECK: Compute overlap between partner behavior and subject behavior
        perc_overlap_subject_partner = length(find(behavior_labels_subject_init == behavior_labels_partner_init))/length(behavior_labels_subject_init);
        perc_overlap_subject_neighbor = length(find(behavior_labels_subject_init == behavior_labels_neighbor_init))/length(behavior_labels_subject_init);
        perc_overlap_partner_neighbor = length(find(behavior_labels_partner_init == behavior_labels_neighbor_init))/length(behavior_labels_subject_init);   
    
        %Only consider pair or alone blocks
        if alone_block==1
            behav_freq_table = tabulate(behavior_labels_subject_init(block_labels== alone_block_id));
            [~, max_beh]= max(behav_freq_table(:,2));

            idx = find(block_labels== alone_block_id & behavior_labels_subject_init==max_beh);
            behavior_labels_subject_init = behavior_labels_subject_init(idx);
            behavior_labels_partner_init = behavior_labels_partner_init(idx);
            behavior_labels_neighbor_init = behavior_labels_neighbor_init(idx);
            Spike_count_raster = Spike_count_raster(idx,:);
            block_labels=block_labels(idx);

        elseif alone_block==0
            behavior_labels_subject_init = behavior_labels_subject_init(block_labels~= alone_block_id);
            behavior_labels_partner_init = behavior_labels_partner_init(block_labels~= alone_block_id);
            behavior_labels_neighbor_init = behavior_labels_neighbor_init(block_labels~= alone_block_id);
            Spike_count_raster = Spike_count_raster(block_labels~= alone_block_id,:);
            block_labels=block_labels(block_labels~= alone_block_id,:);

        elseif alone_block==3
            idx = find(behavior_labels_subject_init==7);
            behavior_labels_subject_init = behavior_labels_subject_init(idx);
            behavior_labels_partner_init = behavior_labels_partner_init(idx);
            behavior_labels_neighbor_init = behavior_labels_neighbor_init(idx);
            Spike_count_raster = Spike_count_raster(idx,:);
            block_labels=block_labels(idx);

        end


        behavior_labels = behavior_labels_neighbor_init;

        %% Select behaviors to decode

        %Compute freq of behavior for the session
        behav_freq_table = tabulate(behavior_labels);
        %behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

        % Select behaviors with a minimum # of occurrences
        behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
%         behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
%         behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.

        % Then select non-reciprocal behaviors
        behav = setdiff(behav, reciprocal_set);

        % OR select behaviors manually
        %behav = [4,5,17,23,25];%manually select behaviors of interest
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
        tabulate(removecats(categorical(behavior_labels_final)));

    
        %check the amount of labels that differ in the partner vs. subject
        %labels after selecting the behaviors of interest.
        subject_behav_after_selection = behavior_labels_subject_init(idx);
        partner_behav_after_selection = behavior_labels_partner_init(idx);
        neighbor_behav_after_selection = behavior_labels_neighbor_init(idx);
        block_after_selection = block_labels(idx);


        overlap_neighbor_subject_after_selection = length(find(neighbor_behav_after_selection == subject_behav_after_selection))/length(idx);
        fprintf('Percent overlap between subject and neighbor labels AFTER selecting behaviors: %s \n', num2str(overlap_neighbor_subject_after_selection))
        
        overlap_neighbor_partner_after_selection = length(find(neighbor_behav_after_selection == partner_behav_after_selection))/length(idx);
        fprintf('Percent overlap between partner and neighbor labels AFTER selecting behaviors: %s \n', num2str(overlap_neighbor_partner_after_selection))

        alone_block_obs = length(find(block_after_selection==3))/length(idx);
        fprintf('Percent observations in "alone" block AFTER selecting behaviors: %s \n', num2str(alone_block_obs))

        %Only consider windows where the behaviors of subject and
        %partner do not overlap
        diff_idx = find(neighbor_behav_after_selection ~= subject_behav_after_selection); %find the indices where subject and partner behavior do not overlap
        Spike_count_raster_final = Spike_count_raster_final(diff_idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels_final(diff_idx);%Same as above but in behavior labels
        subject_behav_after_selection = subject_behav_after_selection(diff_idx);
        block_after_selection_final = block_after_selection(diff_idx);
        behav = unique(behavior_labels_final);

        tabulate(removecats(categorical(behavior_labels_final)));
        tabulate(block_after_selection_final)
        crosstab(removecats(categorical(behavior_labels_final)), block_after_selection_final)
        %Note: not all behaviors are equally happening across blocks. It's
        %possible that what we decode is actually block ID and not the
        %behavior itself...

        %Display which behaviors will be decoded
        behavs_eval = behav_categ(behav);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        fprintf('Behaviors evaluated are: %s \n', behavs_eval);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        %             %Only consider behaviors during the alone block
        %             block_after_selection_overlap_out = block_after_selection(diff_idx);
        %             alone_idx = find(block_after_selection_overlap_out==3);
        %             Spike_count_raster_final = Spike_count_raster_final(alone_idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        %             behavior_labels_final = behavior_labels_final(alone_idx,:);%Same as above but in behavior labels
        %             tabulate(behavior_labels_final);


        %% Run SVM over multiple iterations

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

%         channel = char(channel_flag);
%         disp('****************************************************************************')
%         disp(['ID type ' num2str(partner) ', channels: ' channel '. DONE'])
%         disp('****************************************************************************')
% 
        mean_hitrate{s}(chan) = mean(hitrate)
        sd_hitrate{s}(chan) = std(hitrate);
        mean_hitrate_shuffled{s}(chan) = mean(hitrate_shuffled)
        sd_hitrate_shuffled = std(hitrate_shuffled);

        C_concat=cat(3,C{:});
        confusion_mat_avg=round(mean(C_concat,3)*100);
        rowNames = {labels_id{:,2}}; colNames = {labels_id{:,2}};
        C_table{s,chan} = array2table(confusion_mat_avg,'RowNames',rowNames,'VariableNames',colNames);

    chan = chan +1;
    end% End of channel for loop
    
    clear labels_id

end% End of session for loop


%Change savePath for all session results folder:
cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);
save('SVM_results_partnerBehav.mat', "mean_hitrate","mean_hitrate_shuffled","behav","a_sessions","h_sessions","behav_categ")


%Plot decoding accuracy for all sessions, separated by monkey
figure;  set(gcf,'Position',[150 250 700 700]);
subplot(2,1,1);hold on;
cmap = hsv(size(mean_hitrate,2));
for s = a_sessions
    y = mean_hitrate{s};
    std_dev = sd_hitrate{s};
    errorbar(y,std_dev,'s','MarkerSize',5,'MarkerFaceColor',cmap(s,:))
end
chance_level = 1/length(behav);
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 1])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14;
ylabel('Deconding accuracy','FontSize', 18); xlabel('Brain area','FontSize', 18)
title('Decoding accuracy for partner current behavioral states, Monkey A','FontSize', 14)

subplot(2,1,2);hold on;
cmap = hsv(size(mean_hitrate,2));
for s = h_sessions
    y = mean_hitrate{s};
    std_dev = sd_hitrate{s};
    errorbar(y,std_dev,'s','MarkerSize',5,'MarkerFaceColor',cmap(s,:))
end
chance_level = 1/length(behav);
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 1])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14;
ylabel('Decoding accuracy','FontSize', 18); xlabel('Brain area','FontSize', 18)
title('Decoding accuracy for partner current behavioral states, Monkey H','FontSize', 14)

if randomsample ==0 && onlysingle==1
    saveas(gcf,['SVM_results_accuracy_behav_NOsubsample_unique_allSessions_partner.png'])
elseif randomsample ==1 && onlysingle==1
    saveas(gcf,['SVM_results_accuracy_behav_subsample_unique_allSessions_partner.png'])
elseif randomsample ==0 && onlysingle==0
    saveas(gcf,['SVM_results_accuracy_behav_NOsubsample_NOTunique_allSessions_partner.png'])
elseif randomsample ==1 && onlysingle==0
    saveas(gcf,['SVM_results_accuracy_behav_subsample_NOTunique_allSessions_partner.png'])
end
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot decoding accuracy relative to chance for all sessions, separated by monkey
figure;  set(gcf,'Position',[150 250 700 700]);
subplot(2,1,1);hold on;
cmap = hsv(size(mean_hitrate,2));
for s = a_sessions
    y = mean_hitrate{s}./mean_hitrate_shuffled{s};
    scatter(1:3, y, 60,'filled','MarkerFaceAlpha',0.6)
end
legend({sessions(a_sessions).name},'Location','eastoutside')
chance_level = 1;
yline(chance_level,'--','At Chance', 'FontSize',16)
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 4])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14;
ylabel({'Multiple of chance level','hitrate/shuffled'},'FontSize', 18); xlabel('Brain area','FontSize', 18)
title('Decoding accuracy for partner current behavioral states, Monkey A','FontSize', 14)

subplot(2,1,2);hold on;
cmap = hsv(size(mean_hitrate,2));
for s = h_sessions
    y = mean_hitrate{s}./mean_hitrate_shuffled{s};
     scatter(1:3, y, 60,'filled','MarkerFaceAlpha',0.6)
end
legend({sessions(h_sessions).name},'Location','eastoutside')
chance_level = 1;
yline(chance_level,'--','At Chance', 'FontSize',16)
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 4])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14;
ylabel({'Multiple of chance level','hitrate/shuffled'},'FontSize', 18); xlabel('Brain area','FontSize', 18)
title('Decoding accuracy for partner current behavioral states, Monkey H','FontSize', 14)

if randomsample ==0 && unq_behav==1
    saveas(gcf,['SVM_results_behav_NOsubsample_unique_allSessions_partner.png'])
elseif randomsample ==1 && unq_behav==1
    saveas(gcf,['SVM_results_behav_subsample_unique_allSessions_partner.png'])
elseif randomsample ==0 && unq_behav==0
    saveas(gcf,['SVM_results_behav_NOsubsample_NOTunique_allSessions_partner.png'])
elseif randomsample ==1 && unq_behav==0
    saveas(gcf,['SVM_results_behav_subsample_NOTunique_allSessions_partner.png'])
end
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bar plot decoding accuracy

figure; hold on
data = cell2mat(mean_hitrate');
data_shuffle = cell2mat(mean_hitrate_shuffled');
bp = bar([mean(data(:,:)); mean(data_shuffle(:,:))],'FaceAlpha',0.2);

sp1 = scatter(ones(size(data,1))*0.77,data(:,1), 'filled','b');
sp1 = scatter(ones(size(data,1)),data(:,2), 'filled','r');
sp1 = scatter(ones(size(data,1))*1.22,data(:,3), 'filled','y');

sp1 = scatter(ones(size(data,1))*1.77,data_shuffle(:,1), 'filled','b');
sp1 = scatter(ones(size(data,1))*2,data_shuffle(:,2), 'filled','r');
sp1 = scatter(ones(size(data,1))*2.22,data_shuffle(:,3), 'filled','y');

legend(bp,{'vlPFC','TEO','all'},'Location','best')

ylabel('Decoding Accuracy'); ylim([0.1 0.6])
xticks([1 2]); xticklabels({'Real', 'Shuffled'}); xlim([0.25 2.75])
ax = gca;
ax.FontSize = 16;
saveas(gcf,['SVM_results_allSessions_allUnits_NEIGHBOR.png'])

