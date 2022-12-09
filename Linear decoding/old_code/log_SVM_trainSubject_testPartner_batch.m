%% Log_SVM_trainSubject_testPartner_batch
% Based on the single unit responses, it seems that individual neurons
% respond in a similar direction when the subject is behaving compared to
% the partner. Is this also the case when we consider the neural population
% as a whole? Is neural population coding similar for partner vs. self?
% Train on partner behaving and test on subject behaving.
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
session_range_no_partner=[1:6,11:13,15:16,18];
session_range_with_partner=[1:6,11:13,15:16,18];

%Set parameters
with_partner =1;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 100;%Number of SVM iterations
min_occurrences = 30;%Minimum number of occurrence per behavior
alone_block=0; %1: during alone block; 0:during paired blocks; anything else: all blocks.
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)
null=0;%Set whether we want the null

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


        %% Load data

        %Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        end

        session_length = size(Spike_rasters,2); % get session length

        %Extract behavior labels for subject and partner
        behavior_labels_subject_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        behavior_labels_partner_init = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner

        %Consider sitting in proximity to be rest
        behavior_labels_subject_init(behavior_labels_subject_init==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").
        behavior_labels_partner_init(behavior_labels_partner_init==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

        %Lump approach and leave into travel
        behavior_labels_subject_init(behavior_labels_subject_init==find(behav_categ=="Approach"))=find(behav_categ=="Travel"); %Consider 'approach' to be 'Travel'.
        behavior_labels_partner_init(behavior_labels_partner_init==find(behav_categ=="Approach"))=find(behav_categ=="Travel"); %Consider 'approach' to be 'Travel'.
        behavior_labels_subject_init(behavior_labels_subject_init==find(behav_categ=="leave"))=find(behav_categ=="Travel"); %Consider 'leave' to be 'Travel'.
        behavior_labels_partner_init(behavior_labels_partner_init==find(behav_categ=="Leave"))=find(behav_categ=="Travel"); %Consider 'leave' to be 'Travel'.

        %Extract block info
        block_labels = cell2mat({labels{:,11}}');
        alone_block_id = find(strcmp(block_times{:,"Behavior"},"Alone.block"));

        %Only consider paired blocks
        if alone_block==1
            behavior_labels_subject_init = behavior_labels_subject_init(block_labels== alone_block_id);
            behavior_labels_partner_init = behavior_labels_partner_init(block_labels== alone_block_id);
            Spike_rasters = Spike_rasters(:,block_labels== alone_block_id);
            block_labels=block_labels(block_labels== alone_block_id,:);
        elseif alone_block==0
            behavior_labels_subject_init = behavior_labels_subject_init(block_labels~= alone_block_id);
            behavior_labels_partner_init = behavior_labels_partner_init(block_labels~= alone_block_id);
            Spike_rasters = Spike_rasters(:,block_labels~= alone_block_id);
            block_labels=block_labels(block_labels~= alone_block_id,:);
        end


        %% Get indices for which subject and partner do not behav similarly & don't engage in reciprocal behaviors
        % Select behaviors manually
        behav = [5,18]; %focus on a set of behaviors which happen enough time in both the subject and partner and are not reciprocal

        % Get indices where either:
        % 1. The subject is behaving and the partner is resting
        idx_sub = find(ismember(behavior_labels_subject_init,behav));
        idx_part = find(ismember(behavior_labels_partner_init,behav) &...
            ~ismember(behavior_labels_subject_init,behav));
        %ismember(behavior_labels_subject_init,length(behav_categ)));% &...
        %block_labels~=alone_block_id); %find the indices of the behaviors considered

        Spike_rasters_partner = Spike_rasters(:,idx_part)';%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_partner = behavior_labels_partner_init(idx_part);%Same as above but in behavior labels
        partner_behav = tabulate(behavior_labels_partner);

        Spike_rasters_subject = Spike_rasters(:,idx_sub)';
        behavior_labels_subject = behavior_labels_subject_init(idx_sub);block_labels = cell2mat({labels{:,11}}'); %Extract block info
        subject_behav = tabulate(behavior_labels_subject);

        if length(unique(behavior_labels_partner))==length(behav) && length(unique(behavior_labels_subject))==length(behav)
            if all(partner_behav(behav,2)>30) && all(subject_behav(behav,2)>30)


                %% Run SVM over multiple iterations

                disp('Start running SVM...')
                for iter = 1:num_iter

                    traindata = Spike_rasters_subject; trainlbls = behavior_labels_subject;
                    testdata = Spike_rasters_partner; testlbls = behavior_labels_partner;

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
                    [hitrate_subject(iter)] = log_SVM_basic_function(traindata, trainlbls, 5, 0, 0);
                    [hitrate_partner(iter)] = log_SVM_basic_function(testdata, testlbls, 5, 0, 0);
                    [hitrate_cross(iter), C{iter}] = log_SVM_ccgp_function(trainlbls,traindata,...
                        testlbls, testdata, 0, 0);
                    [hitrate_shuffled(iter), C_shuffled{iter}] = log_SVM_ccgp_function(trainlbls,traindata,...
                        testlbls_shuffled, testdata, 0, 0);

                    disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
                end

                %         channel = char(channel_flag);
                %         disp('****************************************************************************')
                %         disp(['ID type ' num2str(partner) ', channels: ' channel '. DONE'])
                %         disp('****************************************************************************')
                %

                mean_hitrate_partner{s}(chan) = mean(hitrate_partner)
                mean_hitrate_subject{s}(chan) = mean(hitrate_subject)
                mean_hitrate_cross{s}(chan) = mean(hitrate_cross)
                sd_hitrate{s}(chan) = std(hitrate_cross);
                mean_hitrate_shuffled{s}(chan) = mean(hitrate_shuffled)
                sd_hitrate_shuffled = std(hitrate_shuffled);

                C_concat=cat(3,C{:});
                confusion_mat_avg=round(mean(C_concat,3)*100);
                %             rowNames = {labels_id{:,2}}; colNames = {labels_id{:,2}};
                %             C_table{s,chan} = array2table(confusion_mat_avg,'RowNames',rowNames,'VariableNames',colNames);

                chan = chan +1;

            end %end of #behavior classes clause.
        end
    end% End of channel for loop

end% End of session for loop


%Change savePath for all session results folder:
cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);
save('SVM_results_TrainSubjectTestPatrner.mat', "mean_hitrate_subject","mean_hitrate_partner","mean_hitrate_cross","mean_hitrate_shuffled","behav","a_sessions","h_sessions","behav_categ")

load('SVM_results_TrainSubjectTestPatrner.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bar plot decoding accuracy

figure; hold on
data_subject = cell2mat(mean_hitrate_subject');
data_partner = cell2mat(mean_hitrate_partner');
data_cross = cell2mat(mean_hitrate_cross');
data_shuffle = cell2mat(mean_hitrate_shuffled');
bp = bar([mean(data_subject(:,:)); mean(data_partner(:,:)); ...
    mean(data_cross(:,:)); mean(data_shuffle(:,:))],'FaceAlpha',0.2);

sp1 = scatter(ones(size(data_subject,1))*0.78,data_subject(:,1), 'filled','b');
sp1 = scatter(ones(size(data_subject,1)),data_subject(:,2), 'filled','r');
sp1 = scatter(ones(size(data_subject,1))*1.22,data_subject(:,3), 'filled','y');

sp1 = scatter(ones(size(data_partner,1))*1.78,data_partner(:,1), 'filled','b');
sp1 = scatter(ones(size(data_partner,1))*2,data_partner(:,2), 'filled','r');
sp1 = scatter(ones(size(data_partner,1))*2.22,data_partner(:,3), 'filled','y');

sp1 = scatter(ones(size(data_cross,1))*2.78,data_cross(:,1), 'filled','b');
sp1 = scatter(ones(size(data_cross,1))*3,data_cross(:,2), 'filled','r');
sp1 = scatter(ones(size(data_cross,1))*3.22,data_cross(:,3), 'filled','y');

sp1 = scatter(ones(size(data_shuffle,1))*3.78,data_shuffle(:,1), 'filled','b');
sp1 = scatter(ones(size(data_shuffle,1))*4,data_shuffle(:,2), 'filled','r');
sp1 = scatter(ones(size(data_shuffle,1))*4.22,data_shuffle(:,3), 'filled','y');

legend(bp,{'vlPFC','TEO','all'},'Location','best')

ylabel('Decoding Accuracy'); ylim([0.4 1])
xticks([1 2 3 4]); xticklabels({'Subject','Partner','Cross', 'Shuffled'}); xlim([0.25 4.75])
ax = gca;
ax.FontSize = 16;
saveas(gcf,['SVM_results_trainSubjectTestPartner.pdf'])

