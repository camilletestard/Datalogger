%% Log_SVM_SocialContext_batch
% Run a linear decoder on a the neural activity to decode whether the animal is paired or alone.
% Batch version
% March 2022 C. Testard


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
temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 100;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=1; %lump similar behavioral categories together
min_occurrence = 15;
threat_precedence = 1;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
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
    for channel_flag = ["vlPFC", "TEO", "all"]

        %Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence);
        end

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';
        behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        block_labels = cell2mat({labels{:,13}}');

 
        %% Select behaviors to decode

        if simplify == 1
            %Simplify behavioral catagories
            %Lump all aggressive interactions together
            behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
            behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");

            %Lump foraging and drinking
            behavior_labels(behavior_labels==find(behav_categ=="Drinking"))=find(behav_categ=="Foraging");

            %Lump all travel together
            behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
            behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

            behav = [1,5,18,29];
        else
            behav = [4,5,9,10,24,29];
        end

        % For all behaviors
        for b=1:length(behav)

            %Only keep the behaviors of interest
            idx = find(ismember(behavior_labels,behav(b))); %find the indices of the behaviors considered


            Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
            behavior_labels_final = block_labels(idx,:);%Same as above but in behavior labels
            behavfreq= tabulate(behavior_labels_final);

            if ~isempty(idx)

                if size(behavfreq,1)>1 && all(behavfreq(:,2)>=min_occurrence)

                    %% Run SVM over multiple iterations

                    disp('Start running SVM...')
                    for iter = 1:num_iter

                        Labels = behavior_labels_final;

                        %subsample to match number of neurons across brain area:
                        if randomsample==1 && channel_flag~="all"
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
                        end
                        Labels = labels_temp;

                        num_trials = hist(Labels,numericLabels); %number of trials in each class

                        minNumTrials = min_occurrence;%min(num_trials); %30; %find the minimum one %CT change to have 30 of each class
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

                        %disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
                    end

                    mean_hitrate{b}(s, chan) = mean(hitrate)
                    sd_hitrate{b}(s, chan) = std(hitrate);
                    mean_hitrate_shuffled{b}(s, chan) = mean(hitrate_shuffled)
                    sd_hitrate_shuffled{b}(s, chan) = std(hitrate_shuffled);


                    disp(['Behavior: ' behav_categ(behav(b)) ' done.'])

                end %end of min # occurrence clause
            end %end of if not empty loop

        end%end of behavior loop

        chan = chan +1;
        channel = char(channel_flag);
        disp('****************************************************************************')
        disp([num2str(1000/temp_resolution) 'msec resolution, channels: ' channel '. DONE'])
        disp('****************************************************************************')
    end %end of channel loop
    % temp = temp+1;


end %end of session loop

%% Plot all sessions results

%Change savePath for all session results folder:
cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);
save('SVM_allBehavs_AloneVSPaired.mat', "mean_hitrate","mean_hitrate_shuffled","behav","a_sessions","h_sessions","behav_categ")
load('SVM_allBehavs_AloneVSPaired.mat')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bar plot decoding accuracy


figure; hold on; set(gcf,'Position',[150 250 1500 500]); uplimit=20;

for b=1:length(behav)

    subplot(2,4,b); hold on

data = mean_hitrate{b}; data(data==0)=nan;
data_shuffled = mean_hitrate_shuffled{b}; data_shuffled(data_shuffled==0)=nan;
bp = bar([nanmean(data(:,:)); nanmean(data_shuffled(:,:))],'FaceAlpha',0.2);

sp1 = scatter(ones(size(data,1))*0.77,data(:,1), 'filled','b');
sp1 = scatter(ones(size(data,1)),data(:,2), 'filled','r');
sp1 = scatter(ones(size(data,1))*1.22,data(:,3), 'filled','y');

sp1 = scatter(ones(size(data,1))*1.77,data_shuffled(:,1), 'filled','b');
sp1 = scatter(ones(size(data,1))*2,data_shuffled(:,2), 'filled','r');
sp1 = scatter(ones(size(data,1))*2.22,data_shuffled(:,3), 'filled','y');

%legend(bp,{'vlPFC','TEO','all'},'Location','best')
title(behav_categ(behav(b)))

ylabel('Decoding Accuracy'); ylim([0.4 1])
xticks([1 2]); xticklabels({'Real', 'Shuffled'}); xlim([0.25 2.75])
ax = gca;
ax.FontSize = 16;

end
saveas(gcf,['SVM_results_allSessions_allBehavs_SocialContext.pdf'])
