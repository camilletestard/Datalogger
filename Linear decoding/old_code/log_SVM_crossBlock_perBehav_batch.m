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
num_iter = 5000;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=1;%lump similar behavioral categories together to increase sample size.
agg_precedence=0;

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
    %for channel_flag = ["vlPFC", "TEO", "all"]


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
                is_mac, with_NC, isolatedOnly, smooth, sigma, agg_precedence);
        end

        disp('Data Loaded')

        %Raw data
        Spike_count_raster = Spike_rasters';

        behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        %         behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Threat to partner");
        %         behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Threat to subject");
        %         behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

        co_occurrence = cell2mat({labels{:,5}}');

        blocks =  cell2mat({labels{:,12}}');

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
%             behavior_labels(behavior_labels==find(behav_categ=="Getting groomed"))=find(behav_categ=="Groom partner");
        end

        %% Select behaviors to decode

% % %         %Compute freq of behavior for the session
% % %         behav_freq_table = tabulate(behavior_labels);
% % %         behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)
% % % 
% % %         % Select behaviors with a minimum # of occurrences
% % %         min_occurrences = 50;
% % %         behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
% % % 
% % %         % Remove behaviors that are ill-defined
% % %         behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
% % %         behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
% % %         behav = behav(behav~=find(matches(behav_categ,'Other monkeys vocalize')));%excluding scratch which is a source of confusion.
% % %         behav = behav(behav~=find(matches(behav_categ,'Vocalization')));%excluding scratch which is a source of confusion.
% % %         behav = behav(behav~=find(matches(behav_categ,'Rowdy Room')));%excluding scratch which is a source of confusion.
% % %         behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.


        % OR select behaviors manually
        behav = [1,5,7,8,18,29];%[4,5,17];% [7,8]%[5,7:10,21];%[4,5,7:10];%[4:8,17]; %[1:6,9:11,16,17]; %manually select behaviors of interest

        %Print behaviors selected
        behavs_eval = behav_categ(behav);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        fprintf('Behaviors evaluated are: %s \n', behavs_eval);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        %Only keep the behaviors of interest
        idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx);%Same as above but in behavior labels
        block_final = blocks(idx);

        traindata = Spike_count_raster_final(block_final==1,:); trainlbls = behavior_labels_final(block_final==1);
        testdata = Spike_count_raster_final(block_final==2,:); testlbls = behavior_labels_final(block_final==2);

        block1_behav = tabulate(trainlbls);
        block2_behav = tabulate(testlbls);

        if length(unique(testlbls))==length(behav) && length(unique(trainlbls))==length(behav)

            if all(block1_behav(behav,2)>15) && all(block2_behav(behav,2)>15)


                %% Run SVM over multiple iterations

                disp('Start running SVM...')
                for iter = 1:num_iter


                    %Balance number of trials for train data
                    uniqueLabels = unique(trainlbls); %IDentify unique labels (useful when not numbers)
                    NumOfClasses = length(uniqueLabels); % Total number of classes
                    numericLabels = 1:NumOfClasses; %Numeric name of labels

                    labels_temp = trainlbls;
                    for i=1:NumOfClasses
                        idx = trainlbls == uniqueLabels(i);
                        labels_temp(idx) = numericLabels(i);
                        labels_id{i,1} = uniqueLabels(i); labels_id{i,2}=behav_categ{uniqueLabels(i)} ;
                    end
                    trainlbls = labels_temp;

                    num_trials = hist(trainlbls,numericLabels); %number of trials in each class
                    minNumTrials = min(num_trials); %find the minimum one %CT change to have 30 of each class
                    chosen_trials = [];
                    for i = 1:NumOfClasses %for each class
                        idx = find(trainlbls == numericLabels(i)); %find indexes of trials belonging to this class
                        rand_i = randsample(length(idx), minNumTrials); %Select a random n number of them
                        chosen_trials = [chosen_trials; idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
                    end
                    traindata = traindata(chosen_trials, :);
                    trainlbls = trainlbls(chosen_trials, :);
                    trainlbls_shuffled = trainlbls(randperm(length(trainlbls)));

                    %Balance number of trials for test data
                    uniqueLabels = unique(testlbls); %IDentify unique labels (useful when not numbers)
                    NumOfClasses = length(uniqueLabels); % Total number of classes
                    numericLabels = 1:NumOfClasses; %Numeric name of labels

                    labels_temp = testlbls;
                    for i=1:NumOfClasses
                        idx = testlbls == uniqueLabels(i);
                        labels_temp(idx) = numericLabels(i);
                        labels_id{i,1} = uniqueLabels(i); labels_id{i,2}=behav_categ{uniqueLabels(i)} ;
                    end
                    testlbls = labels_temp;

                    num_trials = hist(testlbls,numericLabels); %number of trials in each class
                    minNumTrials = min(num_trials); %find the minimum one %CT change to have 30 of each class
                    chosen_trials = [];
                    for i = 1:NumOfClasses %for each class
                        idx = find(testlbls == numericLabels(i)); %find indexes of trials belonging to this class
                        rand_i = randsample(length(idx), minNumTrials); %Select a random n number of them
                        chosen_trials = [chosen_trials; idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
                    end
                    testdata = testdata(chosen_trials, :);
                    testlbls = testlbls(chosen_trials, :);
                    testlbls_shuffled = testlbls(randperm(length(testlbls)));

                    % Run svm
                    [hitrate_baseline(iter), C{iter}] = log_SVM_basic_function(traindata,trainlbls, 5, 0, 0);
                    [hitrate(iter), C{iter}] = log_SVM_ccgp_function(trainlbls,traindata,...
                        testlbls, testdata, 0, 0);
                    [hitrate_shuffled(iter), C_shuffled{iter}] = log_SVM_ccgp_function(trainlbls,traindata,...
                        testlbls_shuffled, testdata, 0, 0);


                    if mod(iter,10)==1
                        disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
                    end
                end %End of SVM for loop

                channel = char(channel_flag);
                disp('****************************************************************************')
                disp([num2str(1000/temp_resolution) 'msec resolution, channels: ' channel '. DONE'])
                disp('****************************************************************************')

                mean_hitrate_baseline{s}(chan) = mean(hitrate_baseline)
                mean_hitrate{s}(chan) = mean(hitrate)
                sd_hitrate{s}(chan) = std(hitrate);
                mean_hitrate_shuffled{s}(chan) = mean(hitrate_shuffled)
                sd_hitrate_shuffled = std(hitrate_shuffled);

                C_concat=cat(3,C{:}); %Get confusion matrix
                confusion_mat_avg{s,chan}=round(mean(C_concat,3)*100); %Average over SVM iterations
                ccgp_per_behav{s,chan} = diag(confusion_mat_avg{s,chan});
%                 rowNames{s} = {labels_id{:,2}}; colNames{s} = {labels_id{:,2}}; %Get behavior names
%                 C_table{s, temp, chan} = array2table(confusion_mat_avg,'RowNames',rowNames{s},'VariableNames',colNames{s});

                chan = chan +1;
                clear labels_id

            end
        end
    %end %end of channel for loop

    cd(savePath)

end %End of session for loop

%% Plot all sessions results

%Change savePath for all session results folder:
cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);
% save('SVM_results_crossBlock.mat',"confusion_mat_avg","ccgp_per_behav", "mean_hitrate_baseline","mean_hitrate","sd_hitrate","mean_hitrate_shuffled","behav","a_sessions","h_sessions","behav_categ")
load('SVM_results_crossBlock.mat')
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Bar plot decoding accuracy
% % % 
% % % figure; hold on
% % % data_within = cell2mat(mean_hitrate_baseline');
% % % data_cross = cell2mat(mean_hitrate');
% % % data_shuffle = cell2mat(mean_hitrate_shuffled');
% % % bp = bar([mean(data_within(:,:)); mean(data_cross(:,:)); mean(data_shuffle(:,:))],'FaceAlpha',0.2);
% % % 
% % % sp1 = scatter(ones(size(data_within,1))*0.78,data_within(:,1), 'filled','b');
% % % sp1 = scatter(ones(size(data_within,1)),data_within(:,2), 'filled','r');
% % % sp1 = scatter(ones(size(data_within,1))*1.22,data_within(:,3), 'filled','y');
% % % 
% % % sp1 = scatter(ones(size(data_cross,1))*1.78,data_cross(:,1), 'filled','b');
% % % sp1 = scatter(ones(size(data_cross,1))*2,data_cross(:,2), 'filled','r');
% % % sp1 = scatter(ones(size(data_cross,1))*2.22,data_cross(:,3), 'filled','y');
% % % 
% % % sp1 = scatter(ones(size(data_shuffle,1))*2.78,data_shuffle(:,1), 'filled','b');
% % % sp1 = scatter(ones(size(data_shuffle,1))*3,data_shuffle(:,2), 'filled','r');
% % % sp1 = scatter(ones(size(data_shuffle,1))*3.22,data_shuffle(:,3), 'filled','y');
% % % 
% % % legend(bp,{'vlPFC','TEO','all'},'Location','best')
% % % 
% % % ylabel('Decoding Accuracy'); ylim([0.1 1])
% % % xticks([1 2 3]); xticklabels({'Within-context','Cross-context', 'Shuffled'}); xlim([0.25 3.75])
% % % ax = gca;
% % % ax.FontSize = 16;
% % % %saveas(gcf,['SVM_results_crossBlocks.pdf'])
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Plot cross-context generalizability performance per behavior
% % % % ccgp_vlpfc=cat(3,confusion_mat_avg{:,1});
% % % % ccgp_teo=cat(3,confusion_mat_avg{:,2});
% % % % ccgp_all=cat(3,confusion_mat_avg{:,3});
ccgp_all=cat(3,confusion_mat_avg{:,1});

data = ccgp_all;
behav_select = [1,5,7,8,18,29];
colormap={'r','g','c','b','y','k'};
figure; hold on
[sorted_data, idx]=sort(diag(median(data,3))); behav_sorted = behav_select(idx);
scatter(1:length(behav_select),sorted_data,400,'_')
scatter(ones(size(data,3))*1,reshape(data(idx(1),idx(1),:),1,size(data,3)),80, 'filled', colormap{idx(1)});
scatter(ones(size(data,3))*2,reshape(data(idx(2),idx(2),:),1,size(data,3)),80, 'filled', colormap{idx(2)});
scatter(ones(size(data,3))*3,reshape(data(idx(3),idx(3),:),1,size(data,3)),80, 'filled', colormap{idx(3)});
scatter(ones(size(data,3))*4,reshape(data(idx(4),idx(4),:),1,size(data,3)),80, 'filled', colormap{idx(4)});
scatter(ones(size(data,3))*5,reshape(data(idx(5),idx(5),:),1,size(data,3)),80, 'filled', colormap{idx(5)});
scatter(ones(size(data,3))*6,reshape(data(idx(6),idx(6),:),1,size(data,3)),80, 'filled', colormap{idx(6)});
ylabel('Decoding Accuracy'); ylim([0 100])
xticks(1:length(behav_select)); xticklabels(behav_categ(behav_sorted)); xlim([0.25 length(behav_select)+0.75])
ax = gca;
ax.FontSize = 16;
% % % %(gcf,['CCGP_crossBlocks.pdf'])

% % % data = cell2mat(ccgp_per_behav')';
% % % mean_data = nanmedian(data);
% % % figure; hold on
% % % [sorted_data, idx]=sort(mean_data); behav_sorted = behav_select(idx);
% % % scatter(1:length(behav_select),sorted_data,150,'_')
% % % scatter(ones(size(data,3))*1,data(:,idx(1)),80, 'filled', colormap{idx(1)});
% % % scatter(ones(size(data,3))*2,data(:,idx(2)),80, 'filled', colormap{idx(2)});
% % % scatter(ones(size(data,3))*3,data(:,idx(3)),80, 'filled', colormap{idx(3)});
% % % scatter(ones(size(data,3))*4,data(:,idx(4)),80, 'filled', colormap{idx(4)});
% % % scatter(ones(size(data,3))*5,data(:,idx(5)),80, 'filled', colormap{idx(5)});
% % % scatter(ones(size(data,3))*6,data(:,idx(6)),80, 'filled', colormap{idx(6)});
% % % % scatter(ones(size(data,3))*7,data(:,idx(7)),80, 'filled','b');
% % % % scatter(ones(size(data,3))*8,data(:,idx(8)),80, 'filled','b');
% % % ylabel('Decoding Accuracy'); %ylim([0 100])
% % % xticks([1:length(behav_select)]); xticklabels(behav_categ(behav_sorted)); xlim([0.25 length(behav_select)+0.75])
% % % ax = gca;
% % % ax.FontSize = 16;