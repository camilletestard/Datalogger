%% Log_RidgeRegression
% Run a ridge regression (code from Musall et. al)
% Testard C. Feb 2022

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
with_partner =1;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units

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


        %% Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        end

        disp('Data Loaded')

        %Format data
        Vc = Spike_rasters';

        %% Set options

        %set up sampling rate (not sure this is useful, but it's used in code further down so i'll keep it for now).
        opts.Fs = 1;

        %Set up windows for time varying kernels. Multiplying factor is in seconds
        %(recall FS is samples/second ),so can adjust window accordingly
        opts.mPreTime = round(1 * opts.Fs); %motor pre time
        opts.mPostTime = round(1 * opts.Fs); %motor post time
        motorIdx = [-(opts.mPreTime: -1 : 1) 0 (1:opts.mPostTime)]; %index for design matrix to cover pre- and post motor action

        opts.folds = 10; %number of folds for cross-validation

        disp('Options set up')

        %% Preprocessing, extract and format behavioral events

        %Extract behavior info:
        behavior_labels_subject = cell2mat({labels{:,3}}'); %Extract unique behavior info
        behavior_labels_partner = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
        block_labels = cell2mat({labels{:,11}}'); %Extract block info

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

        %Create partner behavior regressors:
        partner_behav_reg = zeros(size(Vc,1), length(behav_categ)-1); %initialize
        for b = 1:length(behav_categ)-1
            idx = find(behavior_labels_partner == b);
            partner_behav_reg(idx,b) = 1;
        end

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
        behavEvents=[subject_behav_reg, partner_behav_reg, block_reg];

        disp('Preprocessing - adding events done')

        %% Set up event types

        %Get variable names
        behavEventNames = [behav_categ(1:end-1), append('partner.',behav_categ(1:end-1)),string(block_times{:,"Behavior"})'];

        %Set up event type here.
        behavEventTypes=ones(size(behavEventNames))*3;
        % For now I will assign all behavioral events as event type 3
        %Event Type 1 = Whole trial
        %Event Type 2 = from stimulus onset to the rest of the trial
        %Event Type 3 = from before movement onset to after movement onset (pre & post time)
        allEventsInfo = [behavEventNames; behavEventTypes]';

        disp('Preprocessing - Assigning event types done')

        %% Preprocessing - Selecting Analog Tracking

        moveR = ME_final;

        %% Setup Design Matrix - Regressor labels

        regLabels = behavEventNames';
        regLabels(length(allEventsInfo)+1) = {'Motion Energy L'};
        regLabels(length(allEventsInfo)+2) = {'Motion Energy R'};


        disp('Setup Design Matrix - Regressor labels done')

        %% Preprocessing - Grouping regressors for later cross-validation
        %For now only include behavior and motion energy

        regGroups = {'BehavioralEvents' 'SubjectBehavior' 'PartnerBehavior' 'block' 'Movements'};
        regGroups{2,1} = behavEventNames';
        regGroups{2,2} = behavEventNames(1:length(behav_categ)-1)';
        regGroups{2,3} = behavEventNames(length(behav_categ):end-3)';
        regGroups{2,4} = behavEventNames(end-2:end);
        regGroups{2,5} = {'Motion Energy L' 'Motion Energy R'} ;

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

        all_nan_inds = find(behavior_labels_subject==length(behav_categ));
        percent_lost_to_nans = size(all_nan_inds,1)/size(Vc,1);

        %remove from predictors and neural data
        behavR(all_nan_inds,:) = []; %blank these indices out
        %moveR(all_nan_inds,:) = [];
        zero_regressors_behav = find(all(behavR == 0,1));
        zero_regressors_move = find(all(moveR == 0,1));

        Vc(all_nan_inds,:) = [];

        %remove empty regressors
        if ~isempty(zero_regressors_behav)
            behavR(:,zero_regressors_behav) = [];
            behavIdx(zero_regressors_behav) = []; %CT: Make sure you remove regIdx every time regressors are removed from fullR
        end
        if ~isempty(zero_regressors_move)
            behavR(:,zero_regressors_move) = [];
            behavIdx(zero_regressors_move) = []; %CT: Make sure you remove regIdx every time regressors are removed from fullR
        end

        disp('Nans removed')

        %2020-12-19 CT: Remove regressors with less than 10 events (as in
        %Musall et. al. 2019)
        low_events_idx= find(sum(behavR,1)<10);
        behavR(:,low_events_idx)=[];
        behavIdx(low_events_idx)=[];

        %% Combine Design Matrix

        %Combine Design Matrix
        fullR=[behavR, moveR];

        %Collect all indecies together.
        regIdx =[behavIdx; moveIdx'];

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

        disp('Data centered')

        %% Orthogonalization of movement variables

        %Get all movement variables (i.e. everything except task regressors)
        allMoveR = fullR(:,length(behavIdx)+1:end); %get instructed + uninstrcuted movements
        allMove_RegIdx = regIdx(length(behavIdx)+1:end);
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
        
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        %only consider behavior regressor (not movement)
        fullR_orthog = behavR;
        regIdx_orthog = behavIdx;

        %Run QR for fully orthogonalized uninstructed movement regressors
        disp('Start QR check for colinear regressors')
        rejIdx = log_run_QR(fullR_orthog,0);
        %saveas(gcf, [savepath 'Colinearity values after rejection FULL MODEL.png'])
        %Note: second argument = plotting. If >0 function will output plot. Otherwise no plots.

        figure;
        histHandle1 = histogram(regIdx_orthog,'BinWidth',1,'FaceAlpha',0.4); hold on
        histHandle1.BinEdges = histHandle1.BinEdges - histHandle1.BinWidth/2;
        histHandle1.BinEdges(length(histHandle1.BinEdges)+1)=max(regIdx_orthog)+0.5;
        if ~isempty(regIdx_orthog(rejIdx))
            histHandle2 = histogram(regIdx_orthog(rejIdx),'BinWidth',1,'FaceColor',[1 0 0],'FaceAlpha',0.4); %distribution of degenerate regressors
            histHandle2.BinEdges = histHandle2.BinEdges - histHandle2.BinWidth/2;
            histHandle2.BinEdges(length(histHandle2.BinEdges)+1)=max(regIdx_orthog(rejIdx))+0.5;
            title('Rejection of regressors in FullR')
            legend('all','rejected')
            %saveas(gcf,[savepath 'Regressors deletion distribution.png'])
            %What reglabels do they correspond to?
            unique(regLabels(regIdx_orthog(rejIdx)))%which regressors are degenerate

            %Remove rejected regressors
            fullR_orthog(:,rejIdx) = [];
            regIdx_orthog(rejIdx,:) = [];
            figure;
            histHandle3 = histogram(regIdx_orthog,'BinWidth',1,'FaceAlpha',0.4); title('Regressors after deletion')
            histHandle3.BinEdges = histHandle3.BinEdges - histHandle3.BinWidth/2;
            histHandle3.BinEdges(length(histHandle3.BinEdges)+1)=max(regIdx_orthog)+0.5;
            %saveas(gcf,[savepath 'Regressors after deletion.png'])
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
        CV_ResultsFull  = modelCorr(Vc,Vfull); %compute model results
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

    chan=chan+1;
    end %end of channel for loop
end

%% Plot some results

for unit =1:size(Vc,2)
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
