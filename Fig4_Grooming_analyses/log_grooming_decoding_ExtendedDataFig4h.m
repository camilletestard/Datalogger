%% log_grooming_decoding_ExtendedDataFig4h.m
%Decode grooming reciprocity from neural population activity.

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


%Set parameters
with_partner =0;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = 'all'; %Channels considered
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

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    %     a_sessions = 1:6; h_sessions = [11:13,15:16,18];
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
    for channel_flag = "all" %["vlPFC", "TEO", "all"]

    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')


    Spike_count_raster = Spike_rasters';
    session_length(s) = size(Spike_count_raster,1);
    behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels(behavior_labels==find(behav_categ=='Proximity'))=length(behav_categ);
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

    %     behavior_labels_tosave{1,s} =behavior_labels';
    %     block_labels_tosave{1,s} =block_labels';

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
        if groomGive_bout_end(bout)+11<session_length(s)
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

        % % % %          %Check what the groom is followed by
        % % % %         if groomGive_bout_end(bout)+120<session_length(s)
        % % % %             unique(behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+120),'stable')
        % % % %             groomGive_bout(bout,6) = behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+120);
        % % % %         else
        % % % %             groomGive_bout(bout,6) =0;
        % % % %         end
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
        if groomGet_bout_end(bout)+11<session_length(s)
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
    %Get inter-bout interval
    allGroomBouts_sorted(:,7)=[0;(allGroomBouts_sorted(2:end,1)-allGroomBouts_sorted(1:end-1,2))];
    %Is previous bout grooming in the other direction
    allGroomBouts_sorted(:,8) = abs([0; diff(allGroomBouts_sorted(:,6))]);
    %Mark turn-taking bouts
    allGroomBouts_sorted(:,9)=0;
    allGroomBouts_sorted(find(allGroomBouts_sorted(:,8)==1 & allGroomBouts_sorted(:,7)<20),9)=1;
    allGroomBouts_sorted(:,10)=s;
    %7th column: inter-bout interval
    %8th column: alternation of grooming bout
    %9th column: is grooming bout turn-taking

    obs_groom = allGroomBouts_sorted;
    obs_recip(s) = 1-((sum(obs_groom(obs_groom(:,6)==7,3)) - sum(obs_groom(obs_groom(:,6)==8,3)))./ sum(obs_groom(:,3)));
    obs_recip_bout(s) = 1-((sum(obs_groom(:,6)==7) - sum(obs_groom(:,6)==8))./ size(obs_groom,1));

    %% Remove very short bouts
    short_bouts = find(allGroomBouts_sorted(:,3)<10);
    for b=1:length(short_bouts)
    behavior_labels(allGroomBouts_sorted(short_bouts(b),1):allGroomBouts_sorted(short_bouts(b),2))=length(behav_categ);
    end

    %% Extract grooming variables

    %Initialize variables which will be used to compute the different
    %grooming metrics
    behav_lbls=behavior_labels;
    groom_behav = behavior_labels;
    groomGive = behavior_labels;
    groomGet = behavior_labels;
    total_time=ones(size(behavior_labels));
    cumul_total_time=cumsum(total_time);

    %Groom give only
    groomGive(groom_behav~=7)=0; groomGive(groom_behav==7)=1;
    groomGiveBout=zeros(size(groomGive)); groomGiveBout(find(diff(groomGive)==1)+1)=1;

    %Groom receive only
    groomGet(groom_behav~=8)=0; groomGet(groom_behav==8)=1;
    groomGetBout=zeros(size(groomGive)); groomGetBout(find(diff(groomGet)==1)+1)=1;

    %Time in bout
    time_in_bout = zeros(size(groomGive));
    bout_number = zeros(size(groomGive));
    for b = 1:size(allGroomBouts_sorted,1)
        idx = allGroomBouts_sorted(b,1) : allGroomBouts_sorted(b,2);
        time_in_bout(idx) = 1:length(idx);
        bout_number(idx) = ones(1,length(idx))*b;
    end

    %Groom give or receive
    groom_behav(groom_behav~=8 & groom_behav~=7)=0;
    groom_behav(groom_behav==8)=1; groom_behav(groom_behav==7)=-1;
    groomBout=zeros(size(groomGive)); groomBout(groomGetBout==1)=1; groomBout(groomGiveBout==1)=-1;

    %Cumulative grooming in the session (+1 if you get groomed, -1 if you groom)
    cumul_groomTime = cumsum(groom_behav); %figure; plot(cumul_groom)
    cumul_groomBout= cumsum(groomBout);  %figure; plot(cumul_groomBout)

    %Total amount of grooming thus far in the session (cumulative)
    total_groom = cumsum(abs(groom_behav));  %figure; plot(total_groom)
    total_groomBout = cumsum(abs(groomBout));  %figure; plot(total_groomBout)

    %Total amount of groom give in the session thus far
    cumul_give =cumsum(groomGive); %figure; plot(cumul_give) % Total time grooming partner so far
    %pGive=cumul_give./cumsum(total_time);%figure; plot(pGive) % Total proportion grooming partner so far
    cumul_giveBouts = cumsum(groomGiveBout); %figure; plot(cumul_giveBouts) % Total number of grooming partner bouts so far
    pGive=cumul_give./total_groom;%figure; plot(pGive)


    %Total amount of groom receive in the session thus far
    cumul_GroomGetTime =cumsum(groomGet); %figure; plot(cumul_get)
    %pGet=cumul_get./cumsum(total_time);%figure; plot(pGet)
    cumul_GroomGetBouts = cumsum(groomGetBout); %figure; plot(cumul_getBouts)
    pGet=cumul_GroomGetTime./total_groom;%figure; plot(pGet)


    %Total reciprocity of groom duration
    recip_groom = 1-(cumul_groomTime./total_groom);%figure; plot(recip_groom)
    recip_groomBout = 1-(cumul_groomBout./total_groomBout);%figure; plot(recip_groomBout)

    %Set categorical variables
    cumul_groomTime_categ = discretize(cumul_groomTime,7); %figure; plot(cumul_groomTime_categ)
    cumul_groomBout_categ = discretize(cumul_groomBout,7); %figure; plot(cumul_groomBout_categ)
    time_in_bout_categ = discretize(time_in_bout,7); %figure; plot(time_in_bout_categ)
    bout_number_categ = discretize(bout_number,10);

    %Final set of variables
    groom_variables=[cumul_groomTime_categ, cumul_groomBout_categ, ...
        bout_number_categ];
%     groom_variables=[bout_number_categ];

    groom_var_names={'Net groom time', 'Net groom bout', 'Bout number'};


    %% Decode observed grooming reciprocity from neural data

    %figure(s); hold on
    %     behav=7;
    behav = [7,8];
    for var=1:length(groom_var_names)

        neural_data = zscore(Spike_count_raster);
        groom_var = groom_variables(:,var);
        groom_idx = ismember(behavior_labels,behav);


        behavior_labels_interim = groom_var(groom_idx);
        Spike_count_raster_interim = neural_data(groom_idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data


        behav_size=tabulate(behavior_labels_interim);
        behav_final=behav_size(behav_size(:,2)>30, 1);
        behavior_labels_final = behavior_labels_interim(ismember(behavior_labels_interim,behav_final)); %figure; plot(behavior_labels_final)
        Spike_count_raster_final = Spike_count_raster_interim(ismember(behavior_labels_interim,behav_final),:);

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
                Input_matrix = zscore(Spike_count_raster_final);
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
                chosen_trials = [chosen_trials; idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
            end
            Input_matrix = Input_matrix(sort(chosen_trials), :);
            Labels = Labels(sort(chosen_trials)); %figure; plot(Labels)
            Labels_shuffled = Labels(randperm(length(Labels)));

            % Run svm
            [hitrate(iter), C{iter}, nErr] = log_SVM_basic_function(Input_matrix, Labels, 5, 0, 0);
            %NOTE: we changed the cross-validation to be continuous
            %chunks instead of randomly generated folds
            [hitrate_shuffled(iter), C_shuffled{iter}, nErr] = log_SVM_basic_function(Input_matrix, Labels_shuffled, 5, 0, 0);

            %                     if mod(iter,10)==1
            disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
            %                     end
        end %end of SVM iterations

        mean_hitrate(var,chan,s) = mean(hitrate);
        mean_hitrate_shuffled(var,chan,s) = mean(hitrate_shuffled);

        %clear labels_id
        disp(['Variable:' num2str(var) '/' num2str(length(groom_var_names)) ' done.'])

    end %End of variable loop
    
    chan=chan+1;

    %end %End of behavior loop


    % % %     %% Correlate observed reciprocity to individual firing rate
    % % %     groom_give = nan(size(groom_all)); groom_give(groom_all==1)=2;
    % % %     groom_receive = nan(size(groom_all)); groom_give(groom_all==-1)=2.1;
    % % %     groom_give=groom_give(idx_all); groom_receive=groom_receive(idx_all);
    % % %
    % % %     [~, pca_score]= pca(neural_data(idx_all,:)); X1=score(:,1:10);
    % % %     for n= 1:size(X1,2)%1:10:200
    % % %         raw_activity(:,n) = pca_score(:,n);
    % % %         smooth_activity(:,n) = movmean(raw_activity(:,n),25);
    % % %         smooth_response = movmean(y_obs,50);
    % % %         [correl_per_neuron(n) p(n)] = corr(raw_activity(:,n),y_obs);
    % % %
    % % % %         figure; hold on;
    % % % %         scatter(1:length(groom_give),groom_give,10,'b','filled');
    % % % %         scatter(1:length(groom_receive),groom_receive,10,'c','filled');
    % % % %         plot(smooth_activity(:,n));
    % % %     end
    % % %     figure; histogram(correl_per_neuron)
    % % %     n_of_interest = find(abs(correl_per_neuron)>0.3);
    % % %     for n=1:length(n_of_interest)
    % % %     figure; hold on; plot(smooth_activity(:,n_of_interest(n))); plot(y_obs); title(num2str(correl_per_neuron(n_of_interest(n))))

    end %end of channel loop
    chan=chan+1;

    disp(['Session ' num2str(s) ' done.'])

end

cd('~/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/')
%save('SVM_reciprocity.mat')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bar plot decoding accuracy

figure; hold on
data = squeeze(mean_hitrate(:,1,session_range))';
data_shuffle = squeeze(mean_hitrate_shuffled(:,1,session_range))';
bp = bar([mean(data(:,:)); mean(data_shuffle(:,:))]','FaceAlpha',0.2);

sp1 = scatter(ones(size(data,1))*0.85,data(:,1), 'filled','b');
sp1 = scatter(ones(size(data_shuffle,1))*1.15,data_shuffle(:,1), 'filled','r');

sp1 = scatter(ones(size(data,1))*1.85,data(:,2), 'filled','b');
sp1 = scatter(ones(size(data_shuffle,1))*2.15,data_shuffle(:,2), 'filled','r');

sp1 = scatter(ones(size(data,1))*2.85,data(:,3), 'filled','b');
sp1 = scatter(ones(size(data_shuffle,1))*3.15,data_shuffle(:,3), 'filled','r');

xticks([1 2 3])
xticklabels({"Net duration", "Net number of bouts", "Chronological bout number"})