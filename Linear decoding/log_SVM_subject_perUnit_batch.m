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
session_range_no_partner=[1:6,11:13,15:16];
session_range_with_partner=[1:6,11:13,15:16];


%Set parameters
with_partner =0;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 100;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null 
simplify=0;%lump similar behavioral categories together to increase sample size.
threat_precedence=1; % 1: aggression takes precedence; 0: Threat to partner and subject states take precedence
exclude_sq = 1;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=1;

for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];

    chan = 1;
    %for channel_flag = ["vlPFC", "TEO", "all"]


    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')

    Spike_count_raster = Spike_rasters';
        behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        co_occurrence = cell2mat({labels{:,5}}');

        if unq_behav==1%Select epochs where only one behavior happens at any given time (i.e. no co-occurrence).
            idx_single = find(co_occurrence==1 | co_occurrence==2 | co_occurrence==3); %(no co-occurrence, or with RR or with proximity)
            Spike_count_raster = Spike_count_raster(idx_single,:);
            behavior_labels = behavior_labels(idx_single,:);
        end

        %% Select behaviors to decode

        %Compute freq of behavior for the session
        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

        % Select behaviors with a minimum # of occurrences
        min_occurrences = 30;
        behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
        behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Other monkeys vocalize')));
        behav = behav(behav~=find(matches(behav_categ,'Rowdy Room')));%excluding scratch which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.


        % OR select behaviors manually
        % behav = [4:10,16,23];%[4,5,17];% [7,8]%[5,7:10,21];%[4,5,7:10];%[4:8,17]; %[1:6,9:11,16,17]; %manually select behaviors of interest

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
        uniqueLabels_session{s} = unique(behavior_labels_final);


        %% Run SVM over increasing numbers of units, over multiple iterations
        u = 1; 
        for unit = 1:size(Spike_count_raster_final,2)

            disp('Start running SVM...')
            for iter=1:num_iter

                %Select unit to run SVM
                Labels = behavior_labels_final;
                Input_matrix = Spike_count_raster_final(:,unit);


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

                % Run svm
                [hitrate(iter), C{iter}] = log_SVM_basic_function(Input_matrix, Labels, 5, 0, 0);
%                 [hitrate_shuffled(iter), C_shuffled{iter}] = log_SVM_basic_function(Input_matrix, Labels_shuffled, 5, 0, 0);
%                 hitrate_ratio(iter) = hitrate(iter)/hitrate_shuffled(iter);

                if mod(iter,50)==1
                    disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
                end
            end %End of SVM for loop

            disp(['SVM run for unit ' num2str(unit)])


            %         channel = char(channel_flag);
            %         disp('****************************************************************************')
            %         disp([num2str(1000/temp_resolution) 'msec resolution, channels: ' channel '. DONE'])
            %         disp('****************************************************************************')

            mean_C{s}(unit,1:25) = nan;
            mean_C{s}(unit,uniqueLabels) = diag(mean(cat(3,C{:}),3));
            mean_hitrate{s}(unit) = mean(hitrate);
            sd_hitrate{s}(unit) = std(hitrate);
            mean_C_aboveChance{s}(unit,1:25) = mean_C{s}(unit,1:25)>(1/NumOfClasses);
%             sd_hitrate_shuffled = std(hitrate_shuffled);
%             mean_hitrate_ratio{s,chan}(u) = mean(hitrate_ratio);
%             sd_hitrate_ratio{s,chan}(u) = std(hitrate_ratio);

            %chan = chan +1;
            clear labels_id

            u=u+1;
        end %End of unit number loop
        chan = chan+1;

    %end %end of channel for loop

    cd(savePath)

% % %     %Plotting results decoding accuracy for all behaviors at 1sec and lower resolution
% % %     figure; hold on; set(gcf,'Position',[150 250 700 500])
% % %     y1 = mean_hitrate{s,1};
% % %     %y2 = mean_hitrate{s,2};
% % %     std_dev1 = sd_hitrate{s,1};
% % %     %std_dev2 = sd_hitrate{s,2};
% % %     errorbar(y1,std_dev1,'s','MarkerSize',10)
% % %     %errorbar(y2,std_dev2,'s','MarkerSize',10)
% % %     chance_level = 1/behav;
% % %     yline(chance_level,'--','Chance level', 'FontSize',16)
% % %     xticks([1:length(unit_num_range)]); xlim([0.8 length(unit_num_range)+0.2]); ylim([0 1])
% % %     xticklabels(unit_num_range)
% % %     ax = gca;
% % %     ax.FontSize = 14;
% % %     ylabel('Decoding accuracy','FontSize', 18); xlabel('#units included','FontSize', 18)
% % %     title('Decoding accuracy for subject current behavioral states with increasing #units','FontSize', 14)
% % %     saveas(gcf, 'Decoding_subject_increasing_units.png')

    close all

end %End of session for loop

%% Plot all sessions results

%Change savePath for all session results folder:
cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);
save('SVM_per_Unit.mat')
%load('IncreasingUnits.mat')

hitrate_all = cell2mat(mean_hitrate); figure; histogram(hitrate_all)
C_all = cell2mat(mean_C');
C_final = C_all(:,~all(isnan(C_all))); 
behav_final = behav_categ(~all(isnan(C_all)));

figure; hold on; 
mean_dim = nanmean(C_final); [~, orderIdx] = sort(mean_dim);
boxchart(C_final(:,orderIdx))
xticklabels(behav_final(orderIdx))
ylabel('Decoding accuracy')

[orderOutput, orderIdx] = sort(sum(C_final>0.1)/size(C_final,1));
figure; hold on
scatter(1:length(orderOutput),orderOutput,40,'filled')
xticks(1:length(orderOutput))
xticklabels(behav_final(orderIdx))
ylabel('Prop decoding above chance')

test=C_final>0.1
figure; histogram(sum(test,2))
