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
num_iter = 10;%Number of SVM iterations
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
    for channel_flag = ["all"]


        %% Get data with specified temporal resolution and channels
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';
        session_length = size(Spike_count_raster ,1);
        behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject

        %Load cumul groom across all sessions
        cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/'])
        load('ZCumulGroom_AcrossSessions.mat')
        cumul_groom_over_all_sessions = cumul_groom_AcrossSessions{s};


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

        %% Extract grooming stats (bout length and# per session)

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

        %% look whithin a session
        time_integration =300;
        groom_behav_perSession = behavior_labels;
        groom_behav_perSession(groom_behav_perSession~=8 & groom_behav_perSession~=7)=0;
        groom_behav_perSession(groom_behav_perSession==8)=1; groom_behav_perSession(groom_behav_perSession==7)=-1;
        cumul_groom_perSession = cumsum(groom_behav_perSession);
        total_groom_perSession = cumsum(abs(groom_behav_perSession));
        recip_groom_perSession = 1-abs((cumul_groom_perSession./total_groom_perSession));

        cumul_groom_slidingWindow = movsum(groom_behav_perSession,[time_integration 0],"omitnan");
        total_groom_slidingWindow = movsum(abs(groom_behav_perSession),[time_integration 0],"omitnan");
        recip_groom_slidingWindow = (cumul_groom_slidingWindow./total_groom_slidingWindow);


        %% Decode grooming context for groom give and groom receive separately
        %Beginning vs. end of grooming bout
        %Grooming post-threat vs not.
        %Grooming after reciprocation or not
        %Groom received after grm Prsnt or not

        beh = 1;behav=8;


        for behav = [7,8] %For both groom give and groom receive


            for groom_categ = 1:5 %For all grooming contexts
                %Note: after threat, almost always groom RECEIVE, not given by
                %subject. Also, grooming after groom present is only for groom
                %RECEIVE.

                for rand_iter = 1:50
                    allGroomBouts_sorted(:,5) = allGroomBouts_sorted(randsample(size(allGroomBouts_sorted,1),size(allGroomBouts_sorted,1)),5);
                    allGroomBouts_sorted(:,9) = allGroomBouts_sorted(randsample(size(allGroomBouts_sorted,1),size(allGroomBouts_sorted,1)),9);
                    
                    groom_bouts_perSession = allGroomBouts_sorted;
                    n_bouts = size(groom_bouts_perSession ,1);
                    duration_bouts = round(random('Exponential',120, 1,n_bouts));
                    groom_bouts_perSession(:,3)=duration_bouts;
                    groom_bouts_perSession(:,6)=randsample(groom_bouts_perSession(:,6),length(groom_bouts_perSession(:,6)));
                    groom_all = zeros(1,size(Spike_count_raster,1));
                    for bout = 1:size(groom_bouts_perSession,1)

                        idx_bout{bout} = groom_bouts_perSession(bout,1):groom_bouts_perSession(bout,2);
                        if groom_bouts_perSession(bout,6) ==8
                            groom_all(idx_bout{bout})=ones(1,length(idx_bout{bout}));
                        else
                            groom_all(idx_bout{bout}) =ones(1,length(idx_bout{bout}))*(-1);
                        end
                    end
                    idx_all = cell2mat(idx_bout);
                    cumul_groom_perSession = cumsum(groom_all); %figure; plot(cumul_groom_perSession)


                    %% Select behaviors to decode

                    groombouts = allGroomBouts_sorted(allGroomBouts_sorted(:,6)==behav,:);
                    bouts_to_consider = randsample(size(groombouts,1),size(groombouts,1));
                    idx_all=[];timescale_all=[]; boutid_all=[]; threatid_all=[];  recipid_all=[];
                    for b=1:length(bouts_to_consider)
                        idx = groombouts(bouts_to_consider(b),1):groombouts(bouts_to_consider(b),2);
                        timescale = round(rescale(idx)*10);
                        bout_id = ones(size(idx))*b;
                        threat_id = ones(size(idx))*groombouts(bouts_to_consider(b),5)+1;
                        recip_id = ones(size(idx))*groombouts(bouts_to_consider(b),9)+1;

                        idx_all = [idx_all, idx];
                        timescale_all = [timescale_all, timescale];
                        boutid_all = [boutid_all, bout_id];
                        threatid_all = [threatid_all, threat_id];
                        recipid_all = [recipid_all, recip_id];

                    end
                    idx_all=idx_all(1:floor(length(idx_all)./10)*10);
                    timescale_all = timescale_all(1:floor(length(idx_all)./10)*10);
                    recipid_all = recipid_all(1:floor(length(idx_all)./10)*10);
                    threatid_all = threatid_all(1:floor(length(idx_all)./10)*10);
                    boutid_all = boutid_all(1:floor(length(idx_all)./10)*10);
                    cumul_groom = cumul_groom_perSession(idx_all); cumul_groom_categ = round(rescale(cumul_groom)*10);
                    %figure; plot(cumul_groom_categ)
                    %                     recip_groom = recip_groom_perSession(idx_all); recip_groom_categ = round(recip_groom*10);

                    %Shuffle by chunk to keep temporal statistics
                    divisors = alldivisors(length(idx_all));
                    [~, idx_min]=min(abs(divisors-50)); divisor = divisors(idx_min);
                    cumul_groom_categ = reshape(cumul_groom_categ,[],divisor);
                    cumul_groom_categ = reshape(cumul_groom_categ(randperm(size(cumul_groom_categ,1)),:),[],1);
%                     recip_groom_categ = reshape(recip_groom_categ,[],divisor);
%                     recip_groom_categ = reshape(recip_groom_categ(randperm(size(recip_groom_categ,1)),:),[],1);
                    timescale_all = reshape(timescale_all,[],divisor);
                    timescale_all = reshape(timescale_all(randperm(size(timescale_all,1)),:),[],1);


                    %                     groom_type_label = [threatid_all; recipid_all; boutid_all; timescale_all'; recip_groom_categ'; cumul_groom_categ'; cumul_groom_over_all_sessions_categ'];
                    groom_type_label = [threatid_all; recipid_all; boutid_all; timescale_all'; cumul_groom_categ'];
                    behavior_labels_interim = groom_type_label(groom_categ,:);
                    Spike_count_raster_interim = Spike_count_raster(idx_all,:);%Only keep timepoints where the behaviors of interest occur in spiking data


                    behav_size=tabulate(behavior_labels_interim);
                    behav_final=behav_size(behav_size(:,2)>10, 1);
                    behavior_labels_final = behavior_labels_interim(ismember(behavior_labels_interim,behav_final));
                    Spike_count_raster_final = Spike_count_raster_interim(ismember(behavior_labels_interim,behav_final),:);


                    behav_size=tabulate(behavior_labels_final);
                    disp('########################')
                    tabulate(behavior_labels_final)
                    disp('########################')


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
                        end
                        Labels = labels_temp;

                        num_trials = hist(Labels,numericLabels); %number of trials in each class
                        minNumTrials = min(100,min(num_trials)); %use the minimum # of instances

                        chosen_trials = [];
                        for i = 1:NumOfClasses %for each class
                            idx = find(Labels == numericLabels(i)); %find indexes of trials belonging to this class
                            rand_i = randsample(length(idx), minNumTrials); %Select a random n number of them
                            chosen_trials = [chosen_trials, idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
                        end
                        Input_matrix = Input_matrix(chosen_trials, :);
                        Labels = Labels(chosen_trials);
                        Labels_shuffled = Labels(randperm(length(Labels)));

                        % Run svm
                        [hitrate(iter), C{iter}] = log_SVM_basic_function(Input_matrix, Labels', 5, 0, 0);
                        [hitrate_shuffled(iter), C_shuffled{iter}] = log_SVM_basic_function(Input_matrix, Labels_shuffled', 5, 0, 0);

                        if mod(iter,10)==1
                            disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
                        end
                    end %end of SVM iterations

                    hitrate_rand(rand_iter) = mean(hitrate);
                   
                end %end of randomizaion iteration
                mean_hitrate(beh,groom_categ,chan,s) = mean(hitrate_rand);


            end %End of grooming context loop
            beh = beh+1;
        end %End of give vs. receive for loop

        chan = chan +1;
    end %End of channel for loop

    close all

    disp('Session done')

end %End of session for loop
