%% Log_SVM_subject_batch
% Run a linear decoder cross neuron populations (i.e. train based on a set
% of neurons and test based on another set of neurons)
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
channel_flag = "vlPFC";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 100;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)
null=0;%Set whether we want the null 
simplify=0;%lump similar behavioral categories together to increase sample size.

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

mean_hitrate_baseline = nan(1,max(session_range));
mean_hitrate = nan(1,max(session_range));
mean_hitrate_shuffled = nan(1,max(session_range));

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];


        %% Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        end

        disp('Data Loaded')

        %Raw data
        Spike_count_raster = Spike_rasters';

        %Low-pass filter data
        %Spike_count_raster = lowpass(Spike_rasters',0.05,1);

        behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
%         behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Threat to partner");
%         behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Threat to subject");
%         behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

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
            behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Aggression");
            behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Aggression");

            %Lump all travel together
            behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
            behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");
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
         
            idx_set1 = randsample(size(Spike_count_raster_final,2),round(size(Spike_count_raster_final,2)/2));
            idx_set2 = setdiff(1:size(Spike_count_raster_final,2),idx_set1)';

            traindata = Spike_count_raster_final(:,idx_set1); trainlbls = behavior_labels_final;
            testdata = Spike_count_raster_final(:,idx_set2); testlbls = trainlbls;
           
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
            [hitrate_baseline(iter)] = log_SVM_basic_function(traindata, trainlbls, 5, 0, 0);
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

        mean_hitrate_baseline(s) = mean(hitrate_baseline)
        mean_hitrate(s) = mean(hitrate)
        sd_hitrate(s) = std(hitrate);
        mean_hitrate_shuffled(s) = mean(hitrate_shuffled)
        sd_hitrate_shuffled = std(hitrate_shuffled);

        C_concat=cat(3,C{:}); %Get confusion matrix
        confusion_mat_avg=round(mean(C_concat,3)*100); %Average over SVM iterations
        rowNames{s} = {labels_id{:,2}}; colNames{s} = {labels_id{:,2}}; %Get behavior names
        C_table{s} = array2table(confusion_mat_avg,'RowNames',rowNames{s},'VariableNames',colNames{s});

        clear labels_id


    cd(savePath)

end %End of session for loop

%% Plot all sessions results

%Change savePath for all session results folder:
cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);
save('SVM_results_crossAreas.mat', "mean_hitrate","sd_hitrate","mean_hitrate_shuffled","behav","a_sessions","h_sessions","behav_categ")
%load('SVM_results_subjectBehav.mat')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bar plot decoding accuracy

figure; hold on
data_baseline = mean_hitrate_baseline';
data = mean_hitrate';
data_shuffle = mean_hitrate_shuffled';
bp = bar([nanmean(data_baseline); nanmean(data); nanmean(data_shuffle)],'FaceAlpha',0.2);

sp1 = scatter(ones(1,size(data,1)),data_baseline, 'filled','b');
sp1 = scatter(ones(1,size(data,1))*2,data, 'filled','b');
sp1 = scatter(ones(1,size(data,1))*3,data_shuffle, 'filled','b');

legend(bp,{'vlPFC','TEO','all'},'Location','best')  

ylabel('Decoding Accuracy'); ylim([0 1])
xticks([1 2 3]); xticklabels({'Baseline', 'Cross-area' , 'Shuffled'}); xlim([0.25 3.75])
ax = gca;
ax.FontSize = 16;

saveas(gcf,['SVM_results_crossAreas.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %Confusion matrices
% % % % % vlpfc = table2array(C_table{s,1,1}); TEO = table2array(C_table{s,1,2}); all_units = table2array(C_table{s,1,3});
% % % % % D{s} = vlpfc-TEO;
% % % % % CustomAxisLabels = string(C_table{s,1,1}.Properties.VariableNames); figure; set(gcf,'Position',[150 250 1500 800]);
% % % % % subplot(2,2,1); pfc = heatmap(vlpfc,'Colormap', hot); pfc.XDisplayLabels = CustomAxisLabels; pfc.YDisplayLabels = CustomAxisLabels; title(pfc,'vlpfc confusion matrix'); caxis([0, 100]);
% % % % % subplot(2,2,2); teo = heatmap(TEO,'Colormap', hot); teo.XDisplayLabels = CustomAxisLabels; teo.YDisplayLabels = CustomAxisLabels; title(teo,'TEO confusion matrix'); caxis([0, 100]);
% % % % % subplot(2,2,3); all = heatmap(all_units,'Colormap', hot); all.XDisplayLabels = CustomAxisLabels; all.YDisplayLabels = CustomAxisLabels; title(all,'All units confusion matrix'); caxis([0, 100]);
% % % % % subplot(2,2,4); h = heatmap(D,'Colormap', cool); h.XDisplayLabels = CustomAxisLabels; h.YDisplayLabels = CustomAxisLabels; title(h,'vlpfc-TEO confusion matrix');caxis([-50, 50]);
