%% Log_SVM_transitions_batch
% Run a linear decoder to decode behavioral transitions from the neural activity
% Batch version
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
transition=1;%1: decode any behavior shift irrespective of what the shift is
             %2: decode SPECIFIC behavioral shifts 
time_around_shift = 4; %in sec

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

    chan = 1;
    for channel_flag = ["vlPFC", "TEO", "all"]


        %% Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        end

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';
        behavior_labels_subject_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        
        %% Extract transitions

        if transition==1
            %Get behavior shift times irrespective of what the shift is
            subject_behav_change_reg = ones(size(labels,1), 1); %initialize
            shift_times = find(diff(behavior_labels_subject_init)~=0);
            shift_idx=[];
            for st =1:length(shift_times)
                shifts = shift_times(st)-time_around_shift : shift_times(st)+time_around_shift;
                shift_idx = [shift_idx,shifts];
            end
            min_occurrences = 50;
            subject_behav_change_reg(shift_idx) = 2;
            behavior_labels = subject_behav_change_reg;
            behav = [1,2];
            %IMPORTANT NOTE: we are NOT able to decode behavior shifts vs.
            %non shifts
            
        elseif transition ==2
            %Get behavior shift with shift 'id'
            x=behavior_labels_subject_init(1:end-1); y=behavior_labels_subject_init(2:end);
            behavior_labels= [sscanf(sprintf('%d%d,',[x.';y.']),'%d,'); 0];
            shifts = x-y;
            
            shift_categ_table= tabulate(behavior_labels(shifts~=0));
            shift_categ_table = shift_categ_table(shift_categ_table(:,2)~=0,:);
            min_occurrences = 10;
            behav = shift_categ_table(shift_categ_table(:,2)>=min_occurrences,1);
            
            %IMPORTANT NOTE: We are able to decode above change what the
            %behavioral shift is when there is one! (that just happened or
            %is just about to happen)
        end
        
        
        %% Select behaviors to decode

        %Only keep the behaviors of interest
        idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx,:);%Same as above but in behavior labels
        tabulate(behavior_labels_final);
        
      
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
                labels_id{i,1} = uniqueLabels(i); %labels_id{i,2}=behav_categ{uniqueLabels(i)} ;
            end
            Labels = labels_temp;
            
            num_trials = hist(Labels,numericLabels); %number of trials in each class
            minNumTrials = min_occurrences;%min(num_trials); %30; %find the minimum one %CT change to have 30 of each class
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
        disp(['Transition type' num2str(transition) ', channels: ' channel '. DONE'])
        disp('****************************************************************************')
        
        mean_hitrate{s}(chan) = mean(hitrate)
        sd_hitrate{s}(chan) = std(hitrate);
        mean_hitrate_shuffled{s}(chan) = mean(hitrate_shuffled)
        sd_hitrate_shuffled = std(hitrate_shuffled);
        
        C_concat=cat(3,C{:});
        confusion_mat_avg=round(mean(C_concat,3)*100);
        rowNames{s}=cell2mat(labels_id);
        %rowNames = {labels_id{:,2}}; colNames = {labels_id{:,2}};
         C_table{s,chan} = array2table(confusion_mat_avg);
        %C_table{chan, trans} = array2table(confusion_mat_avg,'RowNames',rowNames,'VariableNames',colNames);

        chan = chan +1;
    end%end of channel loop

    clear labels_id
end %end of session loop

%% Plot decoding accuracy relative to chance for all sessions, separated by monkey

%Change savePath for all session results folder:
cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);


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
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 4.5])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14;
ylabel({'Multiple of chance level','hitrate/shuffled'},'FontSize', 18); xlabel('Brain area','FontSize', 18)
title('Decoding accuracy for behavioral transitions, Monkey A','FontSize', 14)

subplot(2,1,2);hold on;
cmap = hsv(size(mean_hitrate,2));
for s = h_sessions
    y = mean_hitrate{s}./mean_hitrate_shuffled{s};
     scatter(1:3, y, 60,'filled','MarkerFaceAlpha',0.6)
end
legend({sessions(h_sessions).name},'Location','eastoutside')
chance_level = 1;
yline(chance_level,'--','At Chance', 'FontSize',16)
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 4.5])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14;
ylabel({'Multiple of chance level','hitrate/shuffled'},'FontSize', 18); xlabel('Brain area','FontSize', 18)
title('Decoding accuracy for behavioral transitions, Monkey H','FontSize', 14)
 saveas(gcf,['SVM_results_state_transitions_allSessions.png'])