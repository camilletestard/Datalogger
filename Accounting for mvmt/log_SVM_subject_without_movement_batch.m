%% Log_SVM_subject_without_movement_batch
% Run a linear decoder on the neural activity after regressing out movement for the subject's behavior
% (only including behaviors with a minimum # occurrence)
% Batch version of script
% October 2022 - Camille Testard

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
temp = 1; temp_resolution = 30; %set temp resolution at the camera temp resolution (FPS)
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 3;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null 
simplify=0;%lump similar behavioral categories together to increase sample size.
agg_precedence=0;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=15;
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
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma, agg_precedence);
        end

        %Trim neural data and behavioral to align with video data
        camera_start_time = behavior_log{strcmp(behavior_log{:,'Behavior'},"Camera Sync"),"start_time_round"};
        Spike_rasters_trimmed = Spike_rasters(:,camera_start_time:end);
        labels_trimmed = labels(camera_start_time:end,:);

        cd(filePath)

        %Load ME
        load('hooke0819_motion_energy.mat')
        top_view_ME = [0; top_view_ME]; side_view_ME = [0; side_view_ME];

        %Load DLC
        dlc = readtable('hooke0819_dlc_head_alone.csv');% Load DLC key point data
        dlc=dlc(1:end-1,:); %There is an extra datapoint than frame.. for now ignore the first data point

        logger_bottom = table2array(dlc(:,2:4)); logger_bottom(logger_bottom(:,3)<0.8,1:2)=nan;
        logger_top = table2array(dlc(:,5:7)); logger_top(logger_top(:,3)<0.8,1:2)=nan;
        nose = table2array(dlc(:,8:10)); nose(nose(:,3)<0.8,1:2)=nan;

        %Load head derived measures
        head_derived_measures = load('hooke0819_head_direction.mat');
        head_direction = head_derived_measures.head_orientation_dlc; head_direction = head_direction(1:end-1,:);
        quad_position = head_derived_measures.updown_position; quad_position = quad_position(1:end-1,:);

        disp('Data Loaded')

        %% Pool all the data from the alone block

        % Get alone block
        %For behavior labels
        lbls = cell2mat(labels_trimmed(:,3)); 
        lbls=lbls(1:size(dlc,1)); lbls_not_categ = lbls;
        lbls = categorical(lbls); 

        tabulate(lbls)

        %For spike data
        Spike_rasters_final =  zscore(Spike_rasters_trimmed(:,1:size(dlc,1)),0,2)';

        %Combine mvmt predictors
        logger_top_x = logger_top(:,1);
        logger_top_y = logger_top(:,2);
        mvmt_logger_top_x = [0; diff(logger_top(:,1))];
        mvmt_logger_top_y = [0; diff(logger_top(:,2))];
        head_mvmt = [0; diff(head_direction)];

        mvmt = [top_view_ME, side_view_ME,...
            logger_top_x, logger_top_y,...
            mvmt_logger_top_x, mvmt_logger_top_y,...
            head_direction, head_mvmt,...
            quad_position];

        %mvmt = [top_view_ME, side_view_ME,quad_position];


        %Get missing data (from deeplabcut)
        [nanrow, nancol]=find(isnan(mvmt)); length(unique(nanrow))/length(lbls)
        %We get ~70% missing data because Hooke spends a lot of time in a
        %tiny corner.
        idx_to_keep = setdiff(1:length(lbls), unique(nanrow));

        %Remove missing data
        Y = Spike_rasters_final;
        Y_final = Y(idx_to_keep,:);  lbls_not_categ = lbls_not_categ(idx_to_keep,:);
        lbls_final = removecats(lbls(idx_to_keep));
        top_view_ME_final = zscore(top_view_ME(idx_to_keep));
        side_view_ME_final = zscore(side_view_ME(idx_to_keep));
        logger_top_x_final = zscore(logger_top_x(idx_to_keep));
        logger_top_y_final = zscore(logger_top_y(idx_to_keep));
        mvmt_logger_top_x_final = zscore(mvmt_logger_top_x(idx_to_keep));
        mvmt_logger_top_y_final = zscore(mvmt_logger_top_y(idx_to_keep));
        head_direction_final=zscore(head_direction(idx_to_keep));
        head_mvmt_final = zscore(head_mvmt(idx_to_keep));
        quad_position_final = quad_position(idx_to_keep);
        %mvmt_final = mvmt(idx_to_keep,:);


        %% Run regression
        % Use adjusted Rsquared


        for unit = 1:size(Y_final,2) %for now one unit at a time.

            %Set up predictor matrices
            X_mvmt = table(top_view_ME_final, side_view_ME_final,...
                logger_top_x_final,logger_top_y_final,...
                mvmt_logger_top_x_final, mvmt_logger_top_y_final,...
                head_direction_final,head_mvmt_final,quad_position_final,...
                Y_final(:,unit));

            Results.(['unit' num2str(unit)]) = fitlm(X_mvmt); %run linear model with mvmt only for specific unit
            Y_residuals(:, unit) = Results.(['unit' num2str(unit)]).Residuals.Raw;

            X_lbls = table(lbls_final, Y_final(:,unit));
            ResultsBehav.(['unit' num2str(unit)]) = fitlm(X_lbls); %run linear model with behavior lbls only for specific unit
            Y_residuals_behav(:, unit) = ResultsBehav.(['unit' num2str(unit)]).Residuals.Standardized;

            if mod(unit,10)==0
                disp(unit)
            end

        end

        behavior_labels = lbls_not_categ; %Extract unique behavior info for subject
        behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Threat to partner");
        behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Threat to subject");
        behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest


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

        % Remove behaviors that are ill-defined or to not characterized well enough
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

        for res = 1:3
            %1: the raw firing rates
            %2: residuals after regressing out mvmt
            %3: residuals after regressing out behavior

            %Only keep the behaviors of interest
            idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
            
            if res ==1 %firing rate
                Spike_count_raster_final = Y_final(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
            elseif res==2 %residuals after regressing out movement
                Spike_count_raster_final = Y_residuals(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
            elseif res==3 %residuals after regressing out behavior
                Spike_count_raster_final = Y_residuals_behav(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
            end

            behavior_labels_final = behavior_labels(idx,:);%Same as above but in behavior labels
            tabulate(removecats(categorical(behavior_labels_final)));


            %% Run SVM over multiple iterations

            disp('Start running SVM...')
            for iter = 1:num_iter


                %subsample to match number of neurons across brain areas
                Labels = behavior_labels_final;
                if randomsample==1
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
                minNumTrials = 300;%min(num_trials); %30; %find the minimum one %CT change to have 30 of each class
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

                if mod(iter,5)==1
                    disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
                end
            end %End of SVM for loop

            channel = char(channel_flag);
            disp('****************************************************************************')
            disp([num2str(1000/temp_resolution) 'msec resolution, channels: ' channel '. DONE'])
            disp('****************************************************************************')

            mean_hitrate{s}(chan,res) = mean(hitrate)
            sd_hitrate{s}(chan,res) = std(hitrate);
            mean_hitrate_shuffled{s}(chan,res) = mean(hitrate_shuffled)
            sd_hitrate_shuffled{s}(chan,res) = std(hitrate_shuffled);

            C_concat=cat(3,C{:}); %Get confusion matrix
            confusion_mat_avg{s, chan}=round(mean(C_concat,3)*100); %Average over SVM iterations
            rowNames{s} = {labels_id{:,2}}; colNames{s} = {labels_id{:,2}}; %Get behavior names
            C_table{s, chan} = array2table(confusion_mat_avg{s, chan},'RowNames',rowNames{s},'VariableNames',colNames{s});

            
            clear labels_id

        end %end of residuals used for loop
        % chan = chan +1;

    %end %end of channel for loop

    cd(savePath)

    figure; hold on
    y = [mean_hitrate{s} mean(mean_hitrate_shuffled{s})];
    bar(y)
%     std_dev = [sd_hitrate{s} sd_hitrate_shuffled{s}(1)];
%     errorbar(y,std_dev,'s','MarkerSize',10)
    xlim([0.5 4.5]); xticks([1:4])
    yline(1/length(behav),'LineStyle','--')
    xticklabels({'Raw firing rate','Regressing out mvmt & FOV','Regressing out behavior','shuffled (chance)'})
    ylim([0 1]); ylabel('Decoding accuracy')
    ax = gca;
    ax.FontSize = 14;
    saveas(gcf, 'SVM_regressing_outMvmt.pdf')

end %End of session for loop

%% Plot all sessions results

%Change savePath for all session results folder:
cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);
save('SVM_results_subjectBehav.mat', "mean_hitrate","sd_hitrate","mean_hitrate_shuffled","behav","a_sessions","h_sessions","behav_categ")
%load('SVM_results_subjectBehav.mat')
