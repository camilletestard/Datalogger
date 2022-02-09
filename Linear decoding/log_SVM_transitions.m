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
is_mac = 0;
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

%Set temporal resolution
temp = 1; temp_resolution = 1;
%for temp_resolution = [1/30, 1/20,  1/10, 1/5, 1/2, 1, 2, 5, 10]
%1 for second resolution, 10 for 100msec resolution, 100 for 10msec resolution, 1000 for msec resolution. etc.
%0.1 for 10sec resolution, 1/5 for 5sec resolution

%Set channels: 'TEO', 'vlPFC' or 'all'
chan = 1; channel_flag = "all";
for channel_flag = ["vlPFC", "TEO", "all"]
    
    trans = 1; transition =1;
    for transition = [1,2]
        
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
        
        
        %% Extract transitions
        
        
        if transition==1
            %Get behavior shift times irrespective of what the shift is
            subject_behav_change_reg = ones(size(labels,1), 1); %initialize
            subject_behav_change_reg(find(diff(behavior_labels_subject_init)~=0)) = 2;
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
            min_occurrences = 15;
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
        num_iter = 1000;
        
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
            minNumTrials = 15;%min(num_trials); %30; %find the minimum one %CT change to have 30 of each class
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
        
        mean_hitrate(chan, trans) = mean(hitrate)
        sd_hitrate(chan, trans) = std(hitrate);
        mean_hitrate_shuffled(chan, trans) = mean(hitrate_shuffled)
        sd_hitrate_shuffled = std(hitrate_shuffled);
        
        C_concat=cat(3,C{:});
        confusion_mat_avg=round(mean(C_concat,3)*100);
        %rowNames = {labels_id{:,2}}; colNames = {labels_id{:,2}};
         C_table{chan, trans} = array2table(confusion_mat_avg);
        %C_table{chan, trans} = array2table(confusion_mat_avg,'RowNames',rowNames,'VariableNames',colNames);
        
        trans =trans+1;
        
        clearvars -except is_mac savePath mean_hitrate sd_hitrate mean_hitrate_shuffled C_table temp_resolution channel_flag filePath chan temp behavs_eval behav trans
    end
    chan = chan +1;
end
% temp = temp+1;
% end

colNames = ["Non-specific", "Specific"]; rowNames = ["vlPFC","TEO","all"];
result_hitrate = array2table(mean_hitrate,'RowNames',rowNames,'VariableNames',colNames)
result_sdhitrate = array2table(sd_hitrate,'RowNames',rowNames,'VariableNames',colNames);

save([savePath '\SVM_results_' num2str(length(behav)) 'behav.mat'], 'mean_hitrate', 'sd_hitrate', 'C_table', 'result_hitrate', 'result_sdhitrate', 'behavs_eval')
writetable(result_hitrate,[savePath '\SVM_results_' num2str(length(behav)) 'behav.csv'],'WriteRowNames',true,'WriteVariableNames',true); 

%Plotting results for non-specific state transitions 
figure; hold on; set(gcf,'Position',[150 250 1000 500])
    y = mean_hitrate(:,1);
    std_dev = sd_hitrate(:,1);
    errorbar(y,std_dev,'s','MarkerSize',10,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% chance_level = 1/length(behav);
% yline(chance_level,'--','Chance level', 'FontSize',16, 'Color', cmap(2))
chance_level = 1/2;
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 1])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14; 
ylabel('Deconding accuracy','FontSize', 18); xlabel('Brain Area','FontSize', 18)
title('Decoding non-specific state transitions','FontSize', 20)
cd(savePath)
saveas(gcf,['Decoding non-specific state transitions.png'])
close all


%Plotting results for specific state transitions 
figure; hold on; set(gcf,'Position',[150 250 1000 500])
y = mean_hitrate(:,2);
std_dev = sd_hitrate(:,2);
errorbar(y,std_dev,'s','MarkerSize',10,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
chance_level = 1/length(behav);
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 1])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14;
ylabel('Deconding accuracy','FontSize', 18); xlabel('Brain Area','FontSize', 18)
title('Decoding specific state transitions','FontSize', 20)
cd(savePath)
saveas(gcf,['Decoding specific state transitions.png'])
close all
