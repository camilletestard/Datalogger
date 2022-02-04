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

%Set parameters:
decode_categ =2; %1: decode partner or subject behavior; 2. Decode subject behavior change; 3.Decode block
partner =0; %Predicting partner behavior (partner =1) or subject behavior(partner=0)
predict=1; %Forward or backward lag. I.e. prediction of future behavior (predict =1) or impact of past behavior (predict=0)?
        
%Set temporal resolution
temp = 1; temp_resolution = 1;
%for temp_resolution = [1/30, 1/20,  1/10, 1/5, 1/2, 1, 2, 5, 10]
    %temp_resolution = [1/100,1/50,1/30, 1/20,  1/10, 1/5, 1/2, 1, 2, 5, 10]
    %temp_resolution = [1, 2, 5, 10] %5sec, 2sec, 1sec,500msec, 100msec
    %1 for second resolution, 10 for 100msec resolution, 100 for 10msec resolution, 1000 for msec resolution. etc.
    %0.1 for 10sec resolution, 1/5 for 5sec resolution

    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all"; 
    %for channel_flag = ["vlPFC", "TEO", "all"]

        %Get data with specified temporal resolution and channels
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac);
        %filePath is the experimental data path
        %Temp_resolution is the temporal resolution at which we would like to
        %analyze the data
        %Channel_flag specifies with channels to include: only TEO array, only
        %vlPFC array or all channels
        %is_mac is whether a mac or a pc is being used
        disp('Data Loaded') 

        Spike_count_raster_init = Spike_rasters';
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
        %2 it is often the case that self-groom events co-occur. I expect
        %this to be the case for foraging in Hooke-pair as well.
        
        %% Select subject or partner
        
        if decode_categ==1
            
            if partner ==1
                behavior_labels_init = behavior_labels_partner_init;
            else
                behavior_labels_init = behavior_labels_subject_init;
            end
            
        elseif decode_categ==2
            
            subject_behav_change_reg = ones(size(labels,1), 1); %initialize
            subject_behav_change_reg(find(diff(behavior_labels_subject_init)~=0)) = 2;
            behavior_labels_init = subject_behav_change_reg;
            %IMPORTANT NOTE: we are NOT able to decode behavior shifts vs.
            %non shifts
            
            x=behavior_labels_subject_init(1:end-1); y=behavior_labels_subject_init(2:end);
            behavior_labels_init= [sscanf(sprintf('%d%d,',[x.';y.']),'%d,'); 0];
            shifts = x-y;
            
            shift_categ_table= tabulate(behavior_labels_init(shifts~=0));
            shift_categ_table = shift_categ_table(shift_categ_table(:,2)~=0,:);
            min_occurrences = 15;
            behav = shift_categ_table(shift_categ_table(:,2)>=min_occurrences,1);
            
            %IMPORTANT NOTE: We are able to decode above change what the
            %behavioral shift is when there is one! (that just happened or
            %is just about to happen)
            
        elseif decode_categ==3
            
            behavior_labels_init = block_labels;
        end
        
        %% Time shift behaviors

        lag_length = 0; %Set length of time lag (in window size). Will include in a for loop eventually
        
        if predict %If trying to predict behavior
            behavior_labels = behavior_labels_init(1+lag_length:end);
            Spike_count_raster = Spike_count_raster_init(1:end-lag_length,:);
            
            %check the amount of labels that differ in the lag vs non-lagged versions
            behavior_labels_nonlagged = behavior_labels_init(1:end-lag_length);
            diff_behav_vectors = behavior_labels - behavior_labels_nonlagged; 
            perc_diff_before = round(length(find(diff_behav_vectors~=0))/length(diff_behav_vectors)*100);
            fprintf('Percent difference between-time lagged and regular labels BEFORE selecting behaviors: %s \n', num2str(perc_diff_before));
            
        else  %If testing the influence of past behavior
            behavior_labels = behavior_labels_init(1:end-lag_length);
            Spike_count_raster = Spike_count_raster_init(1+lag_length:end,:);
            
            %check the amount of labels that differ in the lag vs non-lagged versions
            behavior_labels_nonlagged = behavior_labels_init(1+lag_length:end,:);
            diff_behav_vectors = behavior_labels - behavior_labels_nonlagged; 
            perc_diff_before = round(length(find(diff_behav_vectors~=0))/length(diff_behav_vectors)*100);
            fprintf('Percent difference between-time lagged and regular labels BEFORE selecting behaviors: %s \n', num2str(perc_diff_before));
        end
        

        %% Select behaviors to decode
        
        %Compute freq of behavior for the session
%         behav_freq_table = tabulate(behavior_labels);
%         behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

        % Select behaviors with a minimum # of occurrences
%         min_occurrences = 50;
%         behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
%         behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
%         behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
%         behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.
        
        % OR select non-reciprocal behaviors (particularly important
        % for decoding partner behavior)
        %behav = setdiff(behav, reciprocal_set);
        
        % OR select behaviors manually
        %behav = [1 2];%[4:10,16,23];%[4,5,17];% [7,8]%[5,7:10,21];%[4,5,7:10];%[4:8,17]; %[1:6,9:11,16,17]; %manually select behaviors of interest
        
%         %Print behaviors selected
%         behavs_eval = behav_categ(behav);
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%         fprintf('Behaviors evaluated are: %s \n', behavs_eval);
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

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
            tabulate(behavior_labels_final);
        end
        
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
            
            %Only consider behaviors during the alone block
            block_after_selection_overlap_out = block_after_selection(diff_idx);
            alone_idx = find(block_after_selection_overlap_out==3);
            Spike_count_raster_final = Spike_count_raster_final(alone_idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
            behavior_labels_final = behavior_labels_final(alone_idx,:);%Same as above but in behavior labels
            tabulate(behavior_labels_final);
        end
        

        %% Run SVM over multiple iterations
        num_iter = 100;
        
        disp('Start running SVM...')
        for iter = 1:num_iter

            clearvars -except savePath behav_categ behavior_labels_final behavior_labels_shifted Spike_count_raster_final...
                Spike_count_raster_shifted num_iter iter hitrate hitrate_shuffled C C_shuffled temp_resolution...
                channel_flag filePath chan temp mean_hitrate sd_hitrate mean_hitrate_shuffled C_table behavs_eval behav is_mac min_occurrences

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

        clearvars -except is_mac savePath mean_hitrate sd_hitrate mean_hitrate_shuffled C_table temp_resolution channel_flag filePath chan temp behavs_eval behav
%     end
    temp = temp+1;
    
% end

%rowNames = ["100sec","50sec","30sec","20sec","10sec","5sec","2sec","1sec","500msec","200msec","100ms"]; colNames = ["vlPFC","TEO","all"];
rowNames = ["30sec","20sec","10sec","5sec","2sec","1sec","500msec","200msec","100ms"]; colNames = ["vlPFC","TEO","all"];
%rowNames = ["1sec", "500msec","200msec", "100msec"]; colNames = ["vlPFC","TEO","all"];
%rowNames = ["1sec"]; colNames = ["vlPFC","TEO","all"];
result_hitrate = array2table(mean_hitrate,'RowNames',rowNames,'VariableNames',colNames)
result_sdhitrate = array2table(sd_hitrate,'RowNames',rowNames,'VariableNames',colNames)

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
title(leg,'Window size')
chance_level = 1/length(behav);
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 1])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14; 
ylabel('Deconding accuracy','FontSize', 18); xlabel('Brain area','FontSize', 18)
%title('Decoding accuracy for behavioral states','FontSize', 20)
title('Decoding accuracy for social context','FontSize', 20)

cd(savePath)
%saveas(gcf,['SVM_results_SocialContext.png'])
%saveas(gcf,['SVM_results_grooming.png'])
saveas(gcf,['SVM_results_' num2str(length(behavs_eval)) 'behav.png'])
close all

%Plotting for all channels across many time windows - including long one.
%For grooming in vs. out.
figure; hold on; set(gcf,'Position',[150 250 1500 700])
cmap = hsv(size(mean_hitrate,1));
y = mean_hitrate(:,3); y = y(end:-1:1);
x = 1:size(mean_hitrate,1);
std_dev = sd_hitrate(:,3);
errorbar(x,y,std_dev,'-s','MarkerSize',10)
chance_level = 1/length(behav);
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([0 1:length(x) length(x)+1]);ylim([0 1]); xlim([0 10])
xticklabels({'','100msec','200msec','500msec','1sec','2sec','5sec','10sec','20sec','30sec',''})
ax = gca;
ax.FontSize = 14; 
ylabel('Deconding accuracy','FontSize', 18); xlabel('Window size','FontSize', 18)
title('Decoding accuracy for grooming by window size','FontSize', 20)

cd(savePath)
saveas(gcf,['SVM_results_grooming_allChannels.png'])
