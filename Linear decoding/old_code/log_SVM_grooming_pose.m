%% Log_SVM_grooming_batch
% Run a linear decoder on a the neural activity for different grooming contexts
% This script allows to decode grooming:
% 1. Start vs. end
% 2. Post-threat or not
% 3. Reciprocated or not
% 4. Initiated or not
% Batch version to run across sessions.
%Camille Testard, March 2022

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


%Set parameters
with_partner =0;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
randomsample=0;
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 100;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)
null=0;%Set whether we want the null 
simplify=0;%lump similar behavioral categories together to increase sample size.

%Initialize
% mean_hitrate = cell(length(sessions),3);
% sd_hitrate = cell(length(sessions),3);
% mean_hitrate_shuffled = cell(length(sessions),3);

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

    chan=1;
    for channel_flag = ["vlPFC", "TEO", "all"]


        %% Get data with specified temporal resolution and channels
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
            log_GenerateDataToRes_groom_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';
        behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        groom_labels = cell2mat({labels{:,14}}'); 

        %% Decode grooming pose for groom give and groom receive separately
        %Beginning vs. end of grooming bout
        %Grooming post-threat vs not.
        %Grooming after reciprocation or not
        %Groom received after grm Prsnt or not

        
        %% Select behaviors to decode

        % Select behaviors with a minimum # of occurrences
        behav = [8];

        %Only keep the behaviors of interest
        idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        groom_labels_final = groom_labels(idx,:);%Same as above but in behavior labels
        groom_labels_final=groom_labels_final(groom_labels_final~=0);

        pose = unique(groom_labels_final);


        behav_size=tabulate(groom_labels_final);
        disp('########################')
        tabulate(groom_labels_final)
        disp('########################')

        channel = char(channel_flag);
        %                 disp('****************************************************************************')
        %                 disp(['Groom categ: ' groom_categ_label{groom_categ} ', Channels: ' channel ', Behavior: Groom ' groom_behav{beh}])
        %                 disp('****************************************************************************')

        %pause(5)
        %Run SVM over multiple iterations

        disp('Start running SVM...')
        for iter = 1:num_iter

            Labels = groom_labels_final;
            Input_matrix = Spike_count_raster_final;

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
            if min(num_trials)<50 %If there are less than 50 instances for a given behavior
                minNumTrials = min(num_trials); %use the minimum # of instances
            else
                minNumTrials = 50; %impose 50 occurrences per category
            end
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
            [hitrate_shuffled(iter), C_shuffled{iter}] = log_SVM_basic_function(Input_matrix, Labels_shuffled, 5, 0, 0);

            if mod(iter,10)==1
                disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
            end
        end %end of SVM iterations

        mean_hitrate(chan,s) = mean(hitrate)
        sd_hitrate(chan,s) = std(hitrate);
        mean_hitrate_shuffled(chan,s) = mean(hitrate_shuffled)
        sd_hitrate_shuffled = std(hitrate_shuffled);

        clear labels_id

        chan = chan +1;
    end %End of channel for loop

    %% Plotting results decoding accuracy for grooming context

% % % %     figure;  set(gcf,'Position',[150 250 1300 400])
% % % %     subplot(1,2,1);hold on; %Groom Give
% % % %     for c = 1:3
% % % %         y = mean_hitrate{s,c}(1,:);
% % % %         std_dev = sd_hitrate{s,c}(1,:);
% % % %         scatter(1:4,y,60,'filled', 'MarkerFaceAlpha',0.7)
% % % %         %         errorbar(y,std_dev,'s')
% % % %     end
% % % %     leg = legend(["vlPFC","TEO","All"]);
% % % %     title(leg,'Brain Area')
% % % %     chance_level = 1/2;
% % % %     yline(chance_level,'--','Chance level', 'FontSize',16)
% % % %     xticks([1:4]); xlim([0.8 4.2]); ylim([0.4 0.85])
% % % %     xticklabels(groom_categ_label)
% % % %     ax = gca;
% % % %     ax.FontSize = 14;
% % % %     ylabel('Mean decoding accuracy','FontSize', 18); xlabel('Grooming Give Context','FontSize', 18)
% % % %     title('Decoding accuracy for the context of groom give','FontSize', 14)
% % % % 
% % % %     subplot(1,2,2);hold on %Groom Receive
% % % %     for c = 1:3
% % % %         y = mean_hitrate{s,c}(2,:);
% % % %         std_dev = sd_hitrate{s,c}(2,:);
% % % %         scatter(1:4,y,60,'filled', 'MarkerFaceAlpha',0.7)
% % % %         %         errorbar(y,std_dev,'s')
% % % %     end
% % % %     leg = legend(["vlPFC","TEO","All"]);
% % % %     title(leg,'Brain Area')
% % % %     chance_level = 1/2;
% % % %     yline(chance_level,'--','Chance level', 'FontSize',16)
% % % %     xticks([1:4]); xlim([0.8 4.2]); ylim([0.4 0.85])
% % % %     xticklabels(groom_categ_label)
% % % %     ax = gca;
% % % %     ax.FontSize = 14;
% % % %     ylabel('Mean decoding accuracy','FontSize', 18); xlabel('Grooming Receive Context','FontSize', 18)
% % % %     title('Decoding accuracy for the context of groom receive','FontSize', 14)
% % % % 
% % % %     cd(savePath)
% % % %     saveas(gcf,['Decoding grooming given context.png'])
% % % %     close all

end %End of session for loop

%Change savePath for all session results folder:
 cd('~/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/')
save('SVM_results_groomingCateg.mat', "mean_hitrate","mean_hitrate_shuffled","behav","a_sessions","h_sessions","behav_categ","home")
load('SVM_results_groomingCateg.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bar plot decoding accuracy

cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);
mean_hitrate(mean_hitrate==0)=nan;
mean_hitrate_shuffled(mean_hitrate_shuffled==0)=nan;

for groomcat = 2:4
    figure; hold on
    data = squeeze(mean_hitrate(1,groomcat,:,:))';
    data_shuffled = squeeze(mean_hitrate_shuffled(1,groomcat,:,:))';
    bp = bar([nanmean(data); nanmean(data_shuffled)],'FaceAlpha',0.2);

    sp1 = scatter(ones(size(data,1))*0.78,data(:,1), 'filled','b');
    sp1 = scatter(ones(size(data,1)),data(:,2), 'filled','r');
    sp1 = scatter(ones(size(data,1))*1.22,data(:,3), 'filled','y');

    sp1 = scatter(ones(size(data,1))*1.78,data_shuffle(:,1), 'filled','b');
    sp1 = scatter(ones(size(data,1))*2,data_shuffle(:,2), 'filled','r');
    sp1 = scatter(ones(size(data,1))*2.22,data_shuffle(:,3), 'filled','y');

    legend(bp,{'vlPFC','TEO','all'},'Location','best')

    ylabel('Decoding Accuracy'); ylim([0.4 1])
    xticks([1 2]); xticklabels({'Real', 'Shuffled'}); xlim([0.25 2.75])
    ax = gca;
    ax.FontSize = 16;

    title(groom_categ_label(groomcat))

end