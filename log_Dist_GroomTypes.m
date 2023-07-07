%% Log_SVM_grooming_batch
% Run a linear decoder on a the neural activity for different grooming contexts
% This script allows to decode grooming:
% 1. Start vs. end
% 2. Post-threat or not
% 3. Reciprocated or not
% 4. Initiated or not
% 5. Neighbor ID
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
session_range_no_partner=[1:6,11:13,15:16,18];
session_range_with_partner=[1:6,11:13,15:16,18];


%Set parameters
with_partner =0;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
randomsample=0;
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 50;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=0;%lump similar behavioral categories together to increase sample size.
threat_precedence =1;
exclude_sq=1;

%Initialize
% mean_hitrate = cell(length(sessions),3);
% sd_hitrate = cell(length(sessions),3);
% mean_hitrate_shuffled = cell(length(sessions),3);

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

    chan=1;
    for channel_flag = ["vlPFC", "TEO", "all"]


        %% Get data with specified temporal resolution and channels
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';
        session_length = size(Spike_count_raster ,1);
        behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        block_labels = cell2mat({labels{:,12}}');

        %% Interpolate short groom bouts
        groomGive = zeros(size(behavior_labels));
        groomGive(behavior_labels== 7)=1;
        groomGive_time = find(groomGive==1);
        time_between_bouts = diff(groomGive_time);
        short_Runs = find(time_between_bouts>1 & time_between_bouts<11);

        for sr =1:length(short_Runs)
            idx_to_fill = groomGive_time(short_Runs(sr))+1:groomGive_time(short_Runs(sr))+time_between_bouts(short_Runs(sr))-1;
            groomGive(idx_to_fill)=1;
        end

        behavior_labels(find(groomGive==1))=7;

        %Groom get
        groomGet = zeros(size(behavior_labels));
        groomGet(behavior_labels== 8)=1;
        groomGet_time = find(groomGet==1);
        time_between_bouts = diff(groomGet_time);
        short_Runs = find(time_between_bouts>1 & time_between_bouts<11);

        for sr =1:length(short_Runs)
            idx_to_fill = groomGet_time(short_Runs(sr))+1:groomGet_time(short_Runs(sr))+time_between_bouts(short_Runs(sr))-1;
            groomGet(idx_to_fill)=1;
        end

        behavior_labels(find(groomGet==1))=8;

        %% Create grooming type labels

        %GROOM GIVE
        groomGive_bout_start = find(diff(groomGive)==1)+1;
        groomGive_bout_end = find(diff(groomGive)==-1);
        if length(groomGive_bout_end)<length(groomGive_bout_start) %can happen if grooming went until very end of session
            groomGive_bout_end(length(groomGive_bout_start))=length(groomGive);
        end
        groomGive_duration = groomGive_bout_end-groomGive_bout_start;
        groomGive_bout=[groomGive_bout_start, groomGive_bout_end, groomGive_duration];

        for bout = 1:length(groomGive_bout_end)
            %Check groom is followed by a threat
            if groomGive_bout_end(bout)+11<session_length
                groomGive_bout(bout,4) = any(behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+11)==9 ...
                    |behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+11)==10);
            else
                groomGive_bout(bout,4) =0;
            end

            %Check if grooming bout was preceded by a threat
            if groomGive_bout_start(bout)-30>0
                groomGive_bout(bout,5) = any(behavior_labels(groomGive_bout_start(bout)-30:groomGive_bout_start(bout)-1)==9 ...
                    |behavior_labels(groomGive_bout_start(bout)-30:groomGive_bout_start(bout)-1)==10);
            else
                groomGive_bout(bout,5)=0;
            end
        end
        groomGive_bout(:,6)=7;
        %4th column: was there a threat right after the groom (which could have
        %cut it short)
        %5th column: was there a threat preceding the grooming bout
        %6th column, grooming receive or give.

        %cut too short duratiin bouts (usually occurs because smth external
        %cuts the bout (aggression)
        %groomGive_bout(groomGive_bout(:,3)<10)

        %GROOM GET
        groomGet_bout_start = find(diff(groomGet)==1)+1;
        groomGet_bout_end = find(diff(groomGet)==-1);
        groomGet_duration = groomGet_bout_end-groomGet_bout_start;
        groomGet_bout=[groomGet_bout_start, groomGet_bout_end, groomGet_duration];

        for bout = 1:length(groomGet_bout_end)


            if groomGet_bout_end(bout)+11<session_length
                groomGet_bout(bout,4) = any(behavior_labels(groomGet_bout_end(bout)+1:groomGet_bout_end(bout)+11)==9 ...
                    |behavior_labels(groomGet_bout_end(bout)+1:groomGet_bout_end(bout)+11)==10);
            else
                groomGet_bout(bout,4)=0;
            end
            if groomGet_bout_start(bout)-30>0
                groomGet_bout(bout,5) = any(behavior_labels(groomGet_bout_start(bout)-30:groomGet_bout_start(bout)-1)==9 ...
                    |behavior_labels(groomGet_bout_start(bout)-30:groomGet_bout_start(bout)-1)==10);
            else
                groomGet_bout(bout,5) =0;
            end
        end
        groomGet_bout(:,6)=8;

        % ALL GROOMS
        allGroomBouts = [groomGet_bout; groomGive_bout];
        [~, idx_sorted] = sort(allGroomBouts(:,1));
        allGroomBouts_sorted = allGroomBouts(idx_sorted,:);
        allGroomBouts_sorted(:,7)=[0;(allGroomBouts_sorted(2:end,1)-allGroomBouts_sorted(1:end-1,2))];
        allGroomBouts_sorted(:,8) = abs([0; diff(allGroomBouts_sorted(:,6))]);
        allGroomBouts_sorted(:,9)=0;
        allGroomBouts_sorted(find(allGroomBouts_sorted(:,8)==1 & allGroomBouts_sorted(:,7)<20),9)=1;


        beh = 1;behav=8;
        for behav = [8,7] %For both groom give and groom receive

            for sim=1:1000
            groombouts = allGroomBouts_sorted(allGroomBouts_sorted(:,6)==behav,:);
            bouts_to_consider = 1:size(groombouts,1);
            idx_all=[];timescale_all=[]; boutid_all=[]; threatid_all=[];  recipid_all=[]; blockid_all=[];
            for b=1:length(bouts_to_consider)
                idx = groombouts(bouts_to_consider(b),1):groombouts(bouts_to_consider(b),2);
                timescale = round(rescale(idx)*10);
                bout_id = ones(size(idx))*b;
                threat_id = ones(size(idx))*groombouts(bouts_to_consider(b),5)+1;
                recip_id = ones(size(idx))*groombouts(bouts_to_consider(b),9)+1;
                block_id = ones(size(idx))*unique(block_labels(idx));

                idx_all = [idx_all, idx];
                timescale_all = [timescale_all, timescale];
                boutid_all = [boutid_all, bout_id];
                threatid_all = [threatid_all, threat_id];
                recipid_all = [recipid_all, recip_id];
                blockid_all = [blockid_all, block_id];

            end
            groom_type_label = [threatid_all; recipid_all; boutid_all; timescale_all; blockid_all];
            groom_categ=1;
            [~, sorted_idx] = sort(groom_type_label(groom_categ,:));
            behavior_labels_interim = groom_type_label(groom_categ,sorted_idx);
            Spike_count_raster_interim = Spike_count_raster(idx_all(sorted_idx),:);%Only keep timepoints where the behaviors of interest occur in spiking data


            %Balance number of trials per class
            Labels = behavior_labels_interim;
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
            minNumTrials = min(num_trials); %use the minimum # of instances

            chosen_trials = [];
            for i = 1:NumOfClasses %for each class
                idx = find(Labels == numericLabels(i)); %find indexes of trials belonging to this class
                rand_i = randsample(length(idx), minNumTrials); %Select a random n number of them
                chosen_trials = [chosen_trials, idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
            end
            Spike_count_raster_final = Spike_count_raster_interim(chosen_trials, :);
            Labels = Labels(chosen_trials);
            Labels_shuffled = Labels(randperm(length(Labels)));

            [~, score]=pca(zscore(Spike_count_raster_final));
            %test=squareform(pdist(zscore(Spike_count_raster_final),'cosine'));
            test=squareform(pdist(score(:,1:50),'cosine'));
%             figure; hold on; imagesc(test); colorbar; %caxis([0 1.5])
%             xline(find(diff(Labels)~=0)-1,'-k')
%             yline(find(diff(Labels)~=0)-1,'-k')

            within_categ(sim) = mean(mean(test(Labels==2,Labels==2)));
            between_categ(sim) = mean(mean(test(Labels==2,Labels==1)));
            end


            %SHUFFLED CONTROL
            for sim=1:1000
                groombouts = allGroomBouts_sorted(allGroomBouts_sorted(:,6)==behav,:);
                groombouts(:,5) = randsample(groombouts(:,5), size(groombouts,1));
                groombouts(:,9) = randsample(groombouts(:,9), size(groombouts,1));
                bouts_to_consider = 1:size(groombouts,1);
                idx_all=[];timescale_all=[]; boutid_all=[]; threatid_all=[];  recipid_all=[]; blockid_all=[];
                for b=1:length(bouts_to_consider)
                    idx = groombouts(bouts_to_consider(b),1):groombouts(bouts_to_consider(b),2);
                    timescale = round(rescale(idx)*10);
                    bout_id = ones(size(idx))*b;
                    threat_id = ones(size(idx))*groombouts(bouts_to_consider(b),5)+1;
                    recip_id = ones(size(idx))*groombouts(bouts_to_consider(b),9)+1;
                    block_id = ones(size(idx))*unique(block_labels(idx));

                    idx_all = [idx_all, idx];
                    timescale_all = [timescale_all, timescale];
                    boutid_all = [boutid_all, bout_id];
                    threatid_all = [threatid_all, threat_id];
                    recipid_all = [recipid_all, recip_id];
                    blockid_all = [blockid_all, block_id];

                end
                groom_type_label = [threatid_all; recipid_all; boutid_all; timescale_all; blockid_all];
                groom_categ=1;
                [~, sorted_idx] = sort(groom_type_label(groom_categ,:));
                behavior_labels_interim = groom_type_label(groom_categ,sorted_idx);
                Spike_count_raster_interim = Spike_count_raster(idx_all(sorted_idx),:);%Only keep timepoints where the behaviors of interest occur in spiking data


                %Balance number of trials per class
                Labels = behavior_labels_interim;
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
                minNumTrials = min(num_trials); %use the minimum # of instances

                chosen_trials = [];
                for i = 1:NumOfClasses %for each class
                    idx = find(Labels == numericLabels(i)); %find indexes of trials belonging to this class
                    rand_i = randsample(length(idx), minNumTrials); %Select a random n number of them
                    chosen_trials = [chosen_trials, idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
                end
                Spike_count_raster_final = Spike_count_raster_interim(chosen_trials, :);
                Labels = Labels(chosen_trials);
                Labels_shuffled = Labels(randperm(length(Labels)));

                [~, score]=pca(zscore(Spike_count_raster_final));
                %test=squareform(pdist(zscore(Spike_count_raster_final),'cosine'));
                test=squareform(pdist(score(:,1:50),'cosine'));
                %             figure; hold on; imagesc(test); colorbar; %caxis([0 1.5])
                %             xline(find(diff(Labels)~=0)-1,'-k')
                %             yline(find(diff(Labels)~=0)-1,'-k')

                within_categ_shuffle(sim) = mean(mean(test(Labels==2,Labels==2)));
                between_categ_shuffle(sim) = mean(mean(test(Labels==2,Labels==1)));
            end

            figure; hold on; histogram(within_categ); histogram(within_categ_shuffle)
            figure; hold on; histogram(between_categ); histogram(between_categ_shuffle)



        end %End of give vs. receive for loop

        chan = chan +1;
    end %End of channel for loop

    close all

    disp('Session done')

end %End of session for loop
