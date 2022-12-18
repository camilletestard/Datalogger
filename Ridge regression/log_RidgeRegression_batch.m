%% Log_RidgeRegression
% Run a ridge regression (code from Musall et. al)
% Testard C. Feb 2022

%Set session list
is_mac = 0;
is_ron = 0;
if is_mac
    home = '~';
elseif is_ron
    
    home = 'C:\Users\ronwd\Dropbox\Ready to analyze output (1)';
    
else
    home ='C:/Users/GENERAL';
end
%cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
cd(home)
sessions = dir(home); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,15:16,18];
session_range_with_partner=[1:6,11:13,15:16,18];

%Set parameters
with_partner =0;
temp = 1; temp_resolution = 29.97; %frame rate
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1*temp_resolution; %set the smoothing window size (sigma)
simplify=0;
agg_precedence =1;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end


for s =session_range %1:length(sessions)
    %Set path
   
    if is_ron
        
        filePath = [home '\', sessions(s).name];
        savePath = 'C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_results\Mvmt_results';
        
    else
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Mvmt_results'];

    end
    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all"; %2022-12-13: Probably not, but let's check if its a brain region thing/number of neurons thing
    %First trying just TEO, then trying just vlPFC, then just putting in
    %random data as well and seeing if there is still a hold.  Random data
    %works

    %for sig = 1:length(sigma_list)


    %% Get data with specified temporal resolution and channels
    %Get data with specified temporal resolution and channels
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

    cd(filePath) 
% 
%     %2022-12-13: trying totally random data instead.  That runs.  Now
%     %trying non-random behavior data.  This works as well.  So it is
%     %something with Vc that is causing the problem/the neural data in
%     %general.
%     Spike_rasters = poissrnd(.0001,size(Spike_rasters));
%     Spike_rasters(Spike_rasters<0) = 0;

% %Did some simulation work, now setting threshold that cells must be active
% %for at least 0.1% (i.e. 0.001 * session length) of session to be included
% %and seeing if this fixes the issue.
%     criterion = 0.001;
%     low_firing = ((sum(Spike_rasters'>0)/length(Spike_rasters))<criterion);
%     Spike_rasters(low_firing>0,:) = [];
%This doesn't seem to fix the issue.  Might be dropping below activity
%after removing nans, so going to check Vc (neural activity variable)
%before using model
    
    %Load DLC
    dlc = readtable('mvmt_data.csv');% Load DLC key point data
    length(find(sum(isnan(table2array(dlc)),2)==0))/size(dlc,1); %proportion of full data


    %Trim neural and behavioral data to align with video data
    camera_start_time = behavior_log{strcmp(behavior_log{:,'Behavior'},"Camera Sync"),"start_time_round"};
    camera_end_time = camera_start_time + size(dlc,1) -1;

    try
        Spike_rasters_trimmed = Spike_rasters(:,camera_start_time:camera_end_time);
        labels_trimmed = labels(camera_start_time:camera_end_time,:);
    catch %If camera went on past end of recording (occurs in 2 sessions because of an error at the end)
        Spike_rasters_trimmed = Spike_rasters(:,camera_start_time:end);
        labels_trimmed = labels(camera_start_time:end,:);
        dlc = dlc(1:size(labels_trimmed,1),:);
    end

    disp('Data Loaded')


    %% Pool align neural, behavior and video-based mvmt data

    %For behavior labels
    lbls = cell2mat(labels_trimmed(:,3));
    blocks = cell2mat(labels_trimmed(:,12));
    %Note: neural recordings started before camera started and ended after. So we need to trim
    %neural and behavioral data on both ends to match the length of
    %camera-based movement data.

    %Simplify behaviors (lumping categories together)
    lbls(lbls==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

    if simplify
        %Simplify behavioral catagories
        %Lump all aggressive interactions together
        lbls(lbls==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
        lbls(lbls==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");

        %Lump all travel together
        lbls(lbls==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
        lbls(lbls==find(behav_categ=="Leave"))=find(behav_categ=="Travel");
    end

    %tabulate(lbls)

    %Combine mvmt predictors
    mvmt = table2array(dlc);

    %Remove missing data
    lbls_final = lbls;%(idx_to_keep);
    mvmt_final = mvmt;%(idx_to_keep,:);
    block_final = blocks;%(idx_to_keep);


% %     %Extract data table for Seb plotting
% %     toplogger_x = mvmt_final(:,1);
% %     toplogger_y = mvmt_final(:,2);
% %     bottomlogger_x = mvmt_final(:,3);
% %     bottomlogger_y = mvmt_final(:,4);
% %     head_orientation_dlc = mvmt_final(:,5);
% %     dist_traveled = mvmt_final(:,6);
% %     acceleration = mvmt_final(:,7);

    %Format data
    Vc = Spike_rasters_trimmed';

    %% Set options

    %set up sampling rate (not sure this is useful, but it's used in code further down so i'll keep it for now).
    opts.Fs = temp_resolution;

    %Set up windows for time varying kernels. Multiplying factor is in seconds
    %(recall FS is samples/second ),so can adjust window accordingly
    opts.mPreTime = round(1 * opts.Fs); %motor pre time
    opts.mPostTime = round(1 * opts.Fs); %motor post time
    motorIdx = [-(opts.mPreTime: -1 : 1) 0 (1:opts.mPostTime)]; %index for design matrix to cover pre- and post motor action

    opts.folds = 10; %number of folds for cross-validation

    disp('Options set up')

    %% Preprocessing, extract and format behavioral events

    %Extract behavior info:
    behavior_labels_subject = lbls_final;%cell2mat({labels{:,3}}'); %Extract unique behavior info
%     behavior_labels_partner = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
    block_labels = block_final; %cell2mat({labels{:,11}}'); %Extract block info

    %Compute the %of timepoints with unidentified behavior
    %Important note to keep in mind: I still think that the timepoints where
    %there were no behaviors assigned should be excluded...
    no_behavs_timepoints = length(find(behavior_labels_subject == length(behav_categ)))/length(behavior_labels_subject);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('Percent timepoints without a behavior assigned: %s \n', num2str(no_behavs_timepoints));
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    %Create subject behavior regressors:
    subject_behav_reg = zeros(size(Vc,1), length(behav_categ)-1); %initialize
    for b = 1:length(behav_categ)-1 %for all behaviors
        idx = find(behavior_labels_subject == b);
        subject_behav_reg(idx,b) = 1;
    end

%     %Create partner behavior regressors:
%     partner_behav_reg = zeros(size(Vc{s,chan},1), length(behav_categ)-1); %initialize
%     for b = 1:length(behav_categ)-1
%         idx = find(behavior_labels_partner == b);
%         partner_behav_reg(idx,b) = 1;
%     end

    %Create block regressor:
    block_reg = zeros(size(Vc,1), length(unique(block_labels))); %initialize
    for b = 1:length(unique(block_labels))
        idx = find(block_labels == b);
        block_reg(idx,b) = 1;
    end

    % %Create behavior change regressor (exclude for now because we cannot decode behavior shifts irrespective of the shift)
    % subject_behav_change_reg = zeros(size(Vc,1), 1); %initialize
    % subject_behav_change_reg(find(diff(behavior_labels_subject)~=0)+1) = 1;

    %Combine behavioral events
%     behavEvents=[subject_behav_reg, partner_behav_reg, block_reg];
    behavEvents=[subject_behav_reg, block_reg];

    disp('Preprocessing - adding events done')

    %% Set up event types

    %Get variable names
    behavEventNames = [behav_categ(1:end-1),string(block_times{:,"Behavior"})'];%, append('partner.',behav_categ(1:end-1)),string(block_times{:,"Behavior"})'];

    %Set up event type here.
    behavEventTypes=ones(size(behavEventNames))*3;
    % For now I will assign all behavioral events as event type 3
    %Event Type 1 = Whole trial
    %Event Type 2 = from stimulus onset to the rest of the trial
    %Event Type 3 = from before movement onset to after movement onset (pre & post time)
    allEventsInfo = [behavEventNames; behavEventTypes]';

    disp('Preprocessing - Assigning event types done')

    %% Preprocessing - Selecting Analog Tracking

    moveR = mvmt_final;

    %% Setup Design Matrix - Regressor labels

    regLabels = behavEventNames';
    regLabels(length(allEventsInfo)+1) = {'toplogger_x'};
    regLabels(length(allEventsInfo)+2) = {'toplogger_y'};
    regLabels(length(allEventsInfo)+3) = {'bottomlogger_x'};
    regLabels(length(allEventsInfo)+4) = {'bottomlogger_y'};
    regLabels(length(allEventsInfo)+5) = {'head_orientation_dlc'};
    regLabels(length(allEventsInfo)+6) = {'dist_traveled'};
    regLabels(length(allEventsInfo)+7) = {'acceleration'};



    disp('Setup Design Matrix - Regressor labels done')

    %% Preprocessing - Grouping regressors for later cross-validation
    %For now only include behavior and motion energy

    %regGroups = {'BehavioralEvents' 'SubjectBehavior' 'PartnerBehavior' 'block' 'Movements'};
    regGroups = {'BehavioralEvents' 'block' 'Movements'};
    regGroups{2,1} = behavEventNames';
    regGroups{2,2} = behavEventNames(length(behavEventNames)-3+1:end)';
    regGroups{2,3} = regLabels(length(behavEventNames)+1:end);

    disp('Preprocessing - Grouping regressors for later cross-validation done')

    %% Setup Design Matrix - Behavioral Events


    %Creates task regressors with the time varying kernels as described in Churchland
    [behavR, behavIdx] = log_makeDesignMatrix(behavEvents, behavEventTypes, opts);

    %% Setup Design Matrix - Movement
    % Create movement events from analog traces

    %   moveR(:,1:15) = [zeros(1,15); abs(diff(moveR(:,1:15)))]; %Add zeros at the start to keep dimensions the same.

    %Currently subtract 1st percentile to align with Churchland method (instead of min value).
    std_movR = (moveR - prctile(moveR,1))./ nanstd(moveR); %minimum values are at 0, signal in standard deviation units

    [moveMat, traceOut, moveIdx] = log_analogToDesign(std_movR, nanstd(std_movR)*2, motorIdx, behavEventNames);
    % moveMat = "digitized" (or binarized) movement regressors. One cell per movement variable
    % traceOut = original analog traces
    % moveIdx = indices of regressors to keep track of what movement variable a regressor belongs to.

    %moveR contains the movement events (movements above a certain
    %threshold), the original movement variables and ME.
    moveR = [moveMat, traceOut]; %moveR structure: Time X num mvmt regressors
    moveIdx = [moveIdx', (length(behavEventNames)+1:length(behavEventNames)+size(traceOut,2))];

    disp('Setup Design Matrix - Movements done')

    %% Deal with Nans in the behavior data
    %I.e. moments with undefined behaviors.
    % Question: should i remove timepoints where both the partner and the
    % subject have undefined behaviors? only one of the two?

    %Get missing data (from deeplabcut)
    [nanrow, nancol]=find(isnan(mvmt_final)); length(unique(nanrow))/length(lbls)
    %We get ~50-60% missing data.
    all_nan_inds = unique(nanrow);
    percent_lost_to_nans = size(all_nan_inds,1)/size(Vc,1);

    %remove from predictors and neural data
    behavR(all_nan_inds,:) = []; %blank these indices out
    moveR(all_nan_inds,:) = [];
    zero_regressors_behav = find(all(behavR == 0,1));
    zero_regressors_move = find(all(moveR == 0,1));

    Vc(all_nan_inds,:) = [];

    %remove empty regressors
    if ~isempty(zero_regressors_behav)
        behavR(:,zero_regressors_behav) = []; %2022-12-13 RWD: Do you do below here then?  Nvm, done in next cell
        behavIdx(zero_regressors_behav) = []; %CT: Make sure you remove regIdx every time regressors are removed from fullR
    end
    if ~isempty(zero_regressors_move)
        moveR(:,zero_regressors_move) = [];
        moveIdx(zero_regressors_move) = []; %CT: Make sure you remove regIdx every time regressors are removed from fullR
    end

    disp('Nans removed')

    %2020-12-19 CT: Remove regressors with less than 10 events (as in
    %Musall et. al. 2019)
    low_events_idx= find(sum(behavR,1)<10);
    behavR(:,low_events_idx)=[];
    behavIdx(low_events_idx)=[];

    low_events_idx= find(sum(moveR,1)<10);
    moveR(:,low_events_idx)=[];
    moveIdx(low_events_idx)=[];

    %% Combine Design Matrix

    %Combine Design Matrix
    fullR=[behavR , moveR]; %2022-12-13 seeing if one or the other matrix is causing the hang.  Looks like it is moveR...nope it hangs with both of them.
    %Let's try throwing in random, uncorrelated regressors

    %Collect all indecies together.
    regIdx =[ behavIdx; moveIdx']; %see above

% %2022-12-13 RWD trying random regressors to see if these hang.  If they do
% %than issue is probably with neural data
% fullR = normrnd(0,1,size(moveR));
% regIdx = moveIdx';

    disp('Design Matrix Setup Done')


    %% Center and Standardize Continuous Data
    %Center both analog (i.e. movement/tracking data) and neural data.

    % New approach: Find columns that don't just contain zero or one values:
    [~, columns]=find(fullR~=0 & fullR~=1);
    analoginds = unique(columns); clear columns

    %standarize analog regressors
    fullR(:,analoginds) = (fullR(:,analoginds)- mean(fullR(:,analoginds),1))./std(fullR(:,analoginds));

    %Churchland median centered the neuronal data, so we will do the same.
    Vc = (Vc - median(Vc,1));
    
    %2022-12-17 Making sure neurons with sufficient activity to be fit go
    %into the model.  For reasonable tolerance remove any neurons that have
    %less activity than the below criterion.
    criterion = 0.0001; %has to be active for more than .1% of the session
    low_firing = (sum(Vc>0)/length(Vc))<criterion;
    Vc(:,low_firing>0) = []; %remove neurons below this criterion.
    

    disp('Data centered.  Neurons below firing rate criterion removed')

    %% Orthogonalization of movement variables

    %Get all movement variables (i.e. everything except task regressors)
%     allMoveR = fullR(:,length(behavIdx)+1:end); %get instructed + uninstrcuted movements
%     allMove_RegIdx = regIdx(length(behavIdx)+1:end);
    %     disp('Orthogonalize uninstructed variables with respect to instructed movements via QR decomposition')
    %
    %     %Run QR to orthogonalize video variables against all other movement variables
    %     video_reg_id = find(strcmp(regLabels, 'SVD_raw_video')|strcmp(regLabels, 'SVD_ME_video'));
    %     allMoveR_orthog = run_QR_decomp(allMoveR, allMove_RegIdx, video_reg_id);
    %     %Replace relevant columns in fullR
    %     orthog_idx = find(ismember(regIdx, video_reg_id));
    %     fullR(:,orthog_idx) = allMoveR_orthog(:,ismember(allMove_RegIdx, video_reg_id));

    %Replace relevant columns in fullR
    fullR_orthog = fullR;
    regIdx_orthog = regIdx;

    % %     disp('Orthogonalization through QR Done')


    %% run QR and check for rank-defficiency. This will show whether a given regressor is highly collinear with other regressors in the design matrix.

% % %     % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % %     %only consider behavior regressor (not movement)
% % %     fullR_orthog = behavR;
% % %     regIdx_orthog = behavIdx;

    %Run QR for fully orthogonalized uninstructed movement regressors
    disp('Start QR check for colinear regressors')
    rejIdx = log_run_QR(fullR_orthog,0);
    %saveas(gcf, [savepath 'Colinearity values after rejection FULL MODEL.png'])
    %Note: second argument = plotting. If >0 function will output plot. Otherwise no plots.
if is_ron
    
else
    figure;
    histHandle1 = histogram(regIdx_orthog,'BinWidth',1,'FaceAlpha',0.4); hold on
    histHandle1.BinEdges = histHandle1.BinEdges - histHandle1.BinWidth/2;
    histHandle1.BinEdges(length(histHandle1.BinEdges)+1)=max(regIdx_orthog)+0.5;
    
end
    if ~isempty(regIdx_orthog(rejIdx))
        if ~is_ron
        histHandle2 = histogram(regIdx_orthog(rejIdx),'BinWidth',1,'FaceColor',[1 0 0],'FaceAlpha',0.4); %distribution of degenerate regressors
        histHandle2.BinEdges = histHandle2.BinEdges - histHandle2.BinWidth/2;
        histHandle2.BinEdges(length(histHandle2.BinEdges)+1)=max(regIdx_orthog(rejIdx))+0.5;
        title('Rejection of regressors in FullR')
        legend('all','rejected')
        %saveas(gcf,[savepath 'Regressors deletion distribution.png'])
        %What reglabels do they correspond to?
        end
        unique(regLabels(regIdx_orthog(rejIdx)))%which regressors are degenerate

        %Remove rejected regressors
        fullR_orthog(:,rejIdx) = [];
        regIdx_orthog(rejIdx) = []; %2022-12-13 RWD update this was (rejIdx, :) but regIdx is a 1 x regressor vector
        if ~is_ron
        figure;
        histHandle3 = histogram(regIdx_orthog,'BinWidth',1,'FaceAlpha',0.4); title('Regressors after deletion')
        histHandle3.BinEdges = histHandle3.BinEdges - histHandle3.BinWidth/2;
        histHandle3.BinEdges(length(histHandle3.BinEdges)+1)=max(regIdx_orthog)+0.5;
        %saveas(gcf,[savepath 'Regressors after deletion.png'])
        end
    end
    close all

    % %     %Run QR for video regressors orthogonalized only
    % %     disp('Start QR check for colinear regressors')
    % %     rejIdx = run_QR(fullR,1);
    % %     close all
    % %     %2020-12-02 CT: Temporary Check - what are the redundant/degenerate regressors??
    % %     %plot histogram of which regressors are degenerate
    % %     figure;
    % %     histHandle1 = histogram(regIdx,'BinWidth',1,'FaceAlpha',0.4); hold on
    % %     histHandle1.BinEdges = histHandle1.BinEdges - histHandle1.BinWidth/2;
    % %     histHandle1.BinEdges(length(histHandle1.BinEdges)+1)=max(regIdx)+0.5;
    % %     if ~isempty(regIdx(rejIdx))
    % %         histHandle2 = histogram(regIdx(rejIdx),'BinWidth',1,'FaceColor',[1 0 0],'FaceAlpha',0.4); %distribution of degenerate regressors
    % %         histHandle2.BinEdges = histHandle2.BinEdges - histHandle2.BinWidth/2;
    % %         histHandle2.BinEdges(length(histHandle2.BinEdges)+1)=max(regIdx(rejIdx))+0.5;
    % %         title('Rejection of regressors in FullR')
    % %         legend('all','rejected')
    % %         %What reglabels do they correspond to?
    % %         unique(regLabels(regIdx(rejIdx)))%which regressors are degenerate
    % %
    % %         %Remove rejected regressors
    % %         fullR(:,rejIdx) = [];
    % %         regIdx(rejIdx,:) = [];
    % %         figure;
    % %         histHandle3 = histogram(regIdx,'BinWidth',1,'FaceAlpha',0.4); title('Regressors after deletion')
    % %         histHandle3.BinEdges = histHandle3.BinEdges - histHandle3.BinWidth/2;
    % %         histHandle3.BinEdges(length(histHandle3.BinEdges)+1)=max(regIdx)+0.5;
    % %     end
    % %
    % %     close all

    %% Fit full model with cross validation

    disp('Start full model cross validation')
    tic %This is for timing how long it takes to run the entire analysis

    [Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = log_crossValModel(fullR_orthog, Vc, regLabels, regIdx_orthog, regLabels, opts.folds);
    Vfull = Vfull';
    CV_ResultsFull  = modelCorr(Vc,Vfull)
    CV_ResultsFull  = modelCorr(Vc{s,chan},Vfull); %compute model results
    RsqFull(s, chan) = CV_ResultsFull.r_value.^2

    %% Sub_models: (1) regressor groups cvR^2 (shuffle all regressors other than the ones of interest, full contribution)
    %% (2) regressor group dR^2 (shuffle only the regressor of interest, unique contribution)

    % % % % % %         %Save copy of all regressors before shuffling other model fits.
    % % % % % %         fullR_Hold = fullR;
    % % % % % %         groups_of_interest=1:size(regGroups,2);
    % % % % % %
    % % % % % %         for group = groups_of_interest
    % % % % % %
    % % % % % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % %             %%%%% FULL CONTRIBUTION (cvR^2) %%%%%
    % % % % % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % %             %2020-12-27 CT: when computing full contribution use the non-orthogonalized matrix
    % % % % % %             % (except for video variables)
    % % % % % %
    % % % % % %             fullR = fullR_Hold; %reset fullR to the original design matrix with all regressors
    % % % % % %
    % % % % % %             %Get and shuffle everything that is NOT the regressor group of interest
    % % % % % %             shuffleIdx = ismember(regIdx, find(~ismember(regLabels,regGroups{2,group}))); %get index for regressors of interest
    % % % % % %             idx = find(shuffleIdx);
    % % % % % %             shuffling_R = nan(length(fullR), length(idx));
    % % % % % %             for cur_col = 1:length(idx) %for all regressor indices (i.e. columns) that need to be shuffled
    % % % % % %
    % % % % % %                 shuffling_R(:,cur_col) = fullR(randperm(length(shuffling_R)),idx(cur_col)); %shuffle in time each column
    % % % % % %
    % % % % % %             end
    % % % % % %             fullR(:,shuffleIdx) = shuffling_R; %put shuffled columns into full R
    % % % % % %
    % % % % % %             %Fit group model
    % % % % % %             disp(['Start model for group' num2str(group) ', FULL contribution'])
    % % % % % %             [V, Beta, R, Idx, Ridge, Labels] = crossValModel(fullR, Vc, regLabels, regIdx, regLabels, opts.folds);
    % % % % % %
    % % % % % %             CV_Results_full(group) =  modelCorr(Vc,V');
    % % % % % %             Rsq_full(group) =CV_Results_full(group).r_value.^2 %compute explained variance
    % % % % % %             %Rsq_full_adjusted(group) = 1-[(1-Rsq_full(group))*(size(fullR,1)-1)/(size(fullR,1)-size(fullR,2)-1)]
    % % % % % %
    % % % % % %             % % % %         save([savepath 'Only' regGroups{1,group} '_v2.mat'], 'V', 'Beta', 'R', 'Idx', 'Ridge', 'Labels'); %save some results
    % % % % % %
    % % % % % %             disp(['Full contribution model for group' num2str(group) ' cross validation done'])
    % % % % % %             regGroups{3,group}=Rsq_full(group);
    % % % % % %
    % % % % % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % %             %%%%% UNIQUE CONTRIBUTION (dR^2) %%%%%
    % % % % % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % %             %2020-12-27 CT: when computing unique contribution use fully orthogonalized matrix
    % % % % % %
    % % % % % %             fullR = fullR_orthog; %reset fullR to the original design matrix with all regressors
    % % % % % %
    % % % % % %             %Get and shuffle ONLY the regressor group of interest
    % % % % % %             shuffleIdx = ismember(regIdx_orthog, find(ismember(regLabels,regGroups{2,group}))); %get index for regressors of interest
    % % % % % %             idx = find(shuffleIdx);
    % % % % % %             shuffling_R = nan(length(fullR), length(idx));
    % % % % % %             for cur_col = 1:length(idx) %for all regressor indices (i.e. columns) that need to be shuffled
    % % % % % %
    % % % % % %                 shuffling_R(:,cur_col) = fullR(randperm(length(shuffling_R)),idx(cur_col)); %shuffle in time each column
    % % % % % %
    % % % % % %             end
    % % % % % %             fullR(:,shuffleIdx) = shuffling_R; %put shuffled columns into full R
    % % % % % %
    % % % % % %             %Fit group model
    % % % % % %             disp(['Start model for group' num2str(group) ', UNIQUE contribution'])
    % % % % % %             [V, Beta, R, Idx, Ridge, Labels] = crossValModel(fullR, Vc, regLabels, regIdx_orthog, regLabels, opts.folds);
    % % % % % %
    % % % % % %             CV_Results_unique(group) =  modelCorr(Vc,V');
    % % % % % %             Rsq_unique(group) =RsqFull-CV_Results_unique(group).r_value.^2 %compute explained variance
    % % % % % %             Rsq_unique_adjusted(group) = 1-[(1-Rsq_unique(group))*(size(fullR,1)-1)/(size(fullR,1)-size(fullR,2)-1)]
    % % % % % %
    % % % % % %             % % % %         save([savepath 'Except' regGroups{1,group} '_v2.mat'], 'V', 'Beta', 'R', 'Idx', 'Ridge', 'Labels'); %save some results
    % % % % % %
    % % % % % %             disp(['Unique contribution model for group' num2str(group) ' cross validation done'])
    % % % % % %             regGroups{4,group}=Rsq_unique(group);
    % % % % % %
    % % % % % %
    % % % % % %         end
    % % % % % %
    % % % % % %         %%%%%%%%%%%%%%%%%%%
    % % % % % %         %%% Save output %%%
    % % % % % %         %%%%%%%%%%%%%%%%%%%
    % % % % % %
    % % % % % %         output = cell2table(regGroups(3:size(regGroups,1),groups_of_interest),'VariableNames',regGroups(1,groups_of_interest),...
    % % % % % %             'RowNames',{'Full (cvR^2)','Unique (dR^2)'});
    % % % % % %         output.Full = [RsqFull;NaN;];
    % % % % % %         writetable(output,[savePath 'RidgeRegression_results.csv'],'WriteRowNames',true); %save some results

    chan=chwan+1;
end %end of channel for loop


%% Plot some results

for unit =1:size(Vc{s,chan},2)
    correl(unit) = round(corr(Vc{s,chan}(:,unit),Vfull{s,chan}(:,unit)),2);
end

figure; hist(correl)
xlabel('Correlation between predicted and real neural activity')
ylabel('Count')
ax = gca;
ax.FontSize = 14;
xline(mean(correl),'--','LineWidth',2)
cd(savePath)
saveas(gcf,['Histogram_correlation_ridgeReg.png'])

[~,idx_max]=max(correl);

for unit =45;%idx_max

    figure; hold on
    title(['Unit ' num2str(unit) '; r = ' num2str(correl(unit))])
    plot(Vc(:,unit))
    plot(Vfull(:,unit)+6)
    leg = legend('Real', 'Predicted');
    ylabel('Firing rate (Hz)')
    xlabel('Time (in s)')
    ax = gca;
    ax.FontSize = 14;

    cd(savePath)
    saveas(gcf,['Best_example_neuron.png'])

    pause(1)
    close all
end


figure; hold on; set(gcf,'Position',[150 250 1000 500])
scatter(1:6, output{1,:},'b','filled');
scatter(1:6, -output{2,:},'r','filled');
yline(0,'--')
xticks([0.8 1:6 6.2]); xlim([0.8 6.2]); ylim([-0.2 0.2])
xticklabels({'','FULL','All behavior','Subject','Partner','Block ID','Movement',''})
ax = gca;
ax.FontSize = 14;
ylabel('R^2')
leg = legend('Full', 'Unique');
title('Ridge regression Amos 2021-07-29')
cd(savePath)
saveas(gcf,['Ridge_Regression.png'])
close all
