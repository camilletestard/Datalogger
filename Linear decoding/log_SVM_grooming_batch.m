%% Log_SVM_grooming_batch
%% Run a linear decoder on a the neural activity for different grooming contexts
% This script allows to deocde grooming:
% 1. Start vs. end
% 2. Post-threat or not
% 3. Reciprocated or not
% 4. Initiated or not
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
session_range=[1,2,11,12];

%Set parameters
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
randomsample=0;
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 50; %Number of SVM iterations

%Initialize
mean_hitrate = cell(length(sessions),3);
sd_hitrate = cell(length(sessions),3);
mean_hitrate_shuffled = cell(length(sessions),3);

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];

    chan=1;
    for channel_flag = ["vlPFC", "TEO", "all"]

        %Get data with specified temporal resolution and channels
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final, unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        %[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final, unit_count, groom_labels_all]= log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);

        %filePath is the experimental data path
        %Temp_resolution is the temporal resolution at which we would like to
        %analyze the data
        %Channel_flag specifies with channels to include: only TEO array, only
        %vlPFC array or all channels
        %is_mac is whether a mac or a pc is being used
        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';
        behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject

        %% Add labelling for grooming
        %Beginning vs. end of grooming bout
        %Grooming post-threat vs not.
        %Grooming after reciprocation or not
        %Groom received after grm Prsnt or not

        groom_categ_label = {'Star.vs.end', 'Post-threat','Reciprocated','Initiated'};

        beh = 1;
        for behav = [7,8] %For both groom give and groom receive separately
            groom_behav={'Give','Receive'};

            for groom_categ = 1:4 %For all grooming contexts
                %Note: after threat, almost always groom RECEIVE, not given by
                %subject. Also, grooming after groom present is only for groom
                %RECEIVE.

                groom_labels = groom_labels_all(:,groom_categ+1);


                %% Select behaviors to decode

                % Select behaviors with a minimum # of occurrences
                %behav =[find(matches(behav_categ,'Groom Give'))]; %find(matches(behav_categ,'Groom Give')),


                %Only keep the behaviors of interest
                idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
                Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
                behavior_labels_final = groom_labels(idx,:);%Same as above but in behavior labels

                %If groom label is start vs. end
                if groom_categ==1
                    idx_epoch = find(~ismember(behavior_labels_final,3)); %remove middle chunk
                    Spike_count_raster_final = Spike_count_raster_final(idx_epoch,:);%Only keep timepoints where the behaviors of interest occur in spiking data
                    behavior_labels_final = behavior_labels_final(idx_epoch,:);%Same as above but in behavior labels
                end

                behav_size=tabulate(behavior_labels_final);
                disp('########################')
                tabulate(behavior_labels_final)
                disp('########################')

                channel = char(channel_flag);
                %                 disp('****************************************************************************')
                %                 disp(['Groom categ: ' groom_categ_label{groom_categ} ', Channels: ' channel ', Behavior: Groom ' groom_behav{beh}])
                %                 disp('****************************************************************************')

                %pause(5)
                if all(behav_size(:,2)>=30) && length(behav_size(:,2))>1 %If there are at least 30 occurrence of grooming in context 'b'

                    %Run SVM over multiple iterations

                    disp('Start running SVM...')
                    for iter = 1:num_iter

                        %subsample to match number of neurons across brain areas
                        Labels = behavior_labels_final;
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

                        %disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
                    end

                    mean_hitrate{s,chan}(beh,groom_categ) = mean(hitrate);
                    sd_hitrate{s,chan}(beh,groom_categ) = std(hitrate);
                    mean_hitrate_shuffled{s,chan}(beh,groom_categ) = mean(hitrate_shuffled);
                    sd_hitrate_shuffled = std(hitrate_shuffled);

                else
                    mean_hitrate{s,chan}(beh,groom_categ) = nan;
                    sd_hitrate{s,chan}(beh,groom_categ) = nan;
                end % End of "min number of grooming of category b" clause

                clear labels_id

            end %End of grooming context loop
            beh = beh+1;

        end %End of give vs. receive for loop

        chan = chan +1;
    end %End of channel for loop

    %Plotting results decoding accuracy for grooming context
    figure;  set(gcf,'Position',[150 250 1300 400])
    subplot(1,2,1);hold on; %Groom Give
    for c = 1:3
        y = mean_hitrate{s,c}(1,:);
        std_dev = sd_hitrate{s,c}(1,:);
        scatter(1:4,y,60,'filled', 'MarkerFaceAlpha',0.7)
        %         errorbar(y,std_dev,'s')
    end
    leg = legend(["vlPFC","TEO","All"]);
    title(leg,'Brain Area')
    chance_level = 1/2;
    yline(chance_level,'--','Chance level', 'FontSize',16)
    xticks([1:4]); xlim([0.8 4.2]); ylim([0.4 0.85])
    xticklabels(groom_categ_label)
    ax = gca;
    ax.FontSize = 14;
    ylabel('Mean decoding accuracy','FontSize', 18); xlabel('Grooming Give Context','FontSize', 18)
    title('Decoding accuracy for the context of groom give','FontSize', 14)

    subplot(1,2,2);hold on %Groom Receive
    for c = 1:3
        y = mean_hitrate{s,c}(2,:);
        std_dev = sd_hitrate{s,c}(2,:);
        scatter(1:4,y,60,'filled', 'MarkerFaceAlpha',0.7)
        %         errorbar(y,std_dev,'s')
    end
    leg = legend(["vlPFC","TEO","All"]);
    title(leg,'Brain Area')
    chance_level = 1/2;
    yline(chance_level,'--','Chance level', 'FontSize',16)
    xticks([1:4]); xlim([0.8 4.2]); ylim([0.4 0.85])
    xticklabels(groom_categ_label)
    ax = gca;
    ax.FontSize = 14;
    ylabel('Mean decoding accuracy','FontSize', 18); xlabel('Grooming Receive Context','FontSize', 18)
    title('Decoding accuracy for the context of groom receive','FontSize', 14)

    cd(savePath)
    saveas(gcf,['Decoding grooming given context.png'])
    close all

end %End of session for loop

%Change savePath for all session results folder:
savePath = [home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/'];

%Plotting results decoding accuracy for grooming context
figure;  set(gcf,'Position',[150 250 1300 800])
subplot(2,2,1);hold on; %Groom Give, monkey A
cmap={'b','r','y'};
for s=1:2
    for c = 1:3
        y = mean_hitrate{s,c}(1,:);
        std_dev = sd_hitrate{s,c}(1,:);
        scatter(1:4,y,60,'filled', 'MarkerFaceAlpha',0.7,'MarkerFaceColor',cmap{c})
    end
end
leg = legend(["vlPFC","TEO","All"], 'Location','best');
title(leg,'Brain Area')
chance_level = 1/2;
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([1:4]); xlim([0.8 4.2]); ylim([0.4 0.85])
xticklabels(groom_categ_label)
ax = gca;
ax.FontSize = 14;
ylabel('Mean decoding accuracy','FontSize', 18); xlabel('Grooming Give Context','FontSize', 18)
title('Decoding accuracy for the context of groom give, Monkey A','FontSize', 14)

subplot(2,2,2);hold on; %Groom Receive, monkey A
cmap={'b','r','y'};
for s=1:2
    for c = 1:3
        y = mean_hitrate{s,c}(2,:);
        std_dev = sd_hitrate{s,c}(2,:);
        scatter(1:4,y,60,'filled', 'MarkerFaceAlpha',0.7,'MarkerFaceColor',cmap{c})
    end
end
leg = legend(["vlPFC","TEO","All"], 'Location','best');
title(leg,'Brain Area')
chance_level = 1/2;
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([1:4]); xlim([0.8 4.2]); ylim([0.4 0.85])
xticklabels(groom_categ_label)
ax = gca;
ax.FontSize = 14;
ylabel('Mean decoding accuracy','FontSize', 18); xlabel('Grooming Give Context','FontSize', 18)
title('Decoding accuracy for the context of groom receive, Monkey A','FontSize', 14)

subplot(2,2,3);hold on %Groom Give, monkey H
for s=11:12
    for c = 1:3
    y = mean_hitrate{s,c}(1,:);
    std_dev = sd_hitrate{s,c}(1,:);
    scatter(1:4,y,60,'filled', 'MarkerFaceAlpha',0.7,'MarkerFaceColor',cmap{c})
    end
end
leg = legend(["vlPFC","TEO","All"], 'Location','best');
title(leg,'Brain Area')
chance_level = 1/2;
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([1:4]); xlim([0.8 4.2]); ylim([0.4 0.85])
xticklabels(groom_categ_label)
ax = gca;
ax.FontSize = 14;
ylabel('Mean decoding accuracy','FontSize', 18); xlabel('Grooming Receive Context','FontSize', 18)
title('Decoding accuracy for the context of groom give, Monkey H','FontSize', 14)


subplot(2,2,4);hold on %Groom Receive, monkey A
for s=11:12
    for c = 1:3
    y = mean_hitrate{s,c}(2,:);
    std_dev = sd_hitrate{s,c}(2,:);
    scatter(1:4,y,60,'filled', 'MarkerFaceAlpha',0.7,'MarkerFaceColor',cmap{c})
    %         errorbar(y,std_dev,'s')
    end
end
leg = legend(["vlPFC","TEO","All"], 'Location','best');
title(leg,'Brain Area')
chance_level = 1/2;
yline(chance_level,'--','Chance level', 'FontSize',16)
xticks([1:4]); xlim([0.8 4.2]); ylim([0.4 0.85])
xticklabels(groom_categ_label)
ax = gca;
ax.FontSize = 14;
ylabel('Mean decoding accuracy','FontSize', 18); xlabel('Grooming Receive Context','FontSize', 18)
title('Decoding accuracy for the context of groom receive, Monkey H','FontSize', 14)

cd(savePath)
saveas(gcf,['Decoding grooming given context.png'])
