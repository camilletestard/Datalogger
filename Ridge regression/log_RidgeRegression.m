%% Log_RidgeRegression
% Run a ridge regression (code from Musall et. al)


%% Load data

%Set path
is_mac = 0;
if is_mac
    cd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
else
    cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
end
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)

if is_mac
    cd('~/Dropbox (Penn)/Datalogger/Results/')
else
    cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Results/')
end
savePath = uigetdir('', 'Please select the result directory');

clearvars -except savePath filePath is_mac

%Set temporal resolution, channel
temp = 1; temp_resolution = 1; chan = 1; channel_flag = "all";
%for temp_resolution = [1/30, 1/20,  1/10, 1/5, 1/2, 1, 2, 5, 10]
%1 for second resolution, 10 for 100msec resolution, 100 for 10msec resolution, 1000 for msec resolution. etc.
%0.1 for 10sec resolution, 1/5 for 5sec resolution

%Set channels: 'TEO', 'vlPFC' or 'all'
%for channel_flag = ["vlPFC", "TEO", "all"]

%Get data with specified temporal resolution and channels
[Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac);
%filePath is the experimental data path
%Temp_resolution is the temporal resolution at which we would like to
%analyze the data
%Channel_flag specifies with channels to include: only TEO array, only
%vlPFC array or all channels
%is_mac is whether a mac or a pc is being used
disp('Data Loaded')

%Format data
Spike_count_raster = Spike_rasters';

%% Set options

%set up sampling rate (not sure this is useful, but it's used in code further down so i'll keep it for now).
opts.Fs = 1;

%Set up windows for time varying kernels. Multiplying factor is in seconds
%(recall FS is samples/second ),so can adjust window accordingly
opts.mPreTime = round(1 * opts.Fs); %motor pre time
opts.mPostTime = round(2 * opts.Fs); %motor post time

opts.folds = 10; %number of folds for cross-validation


%% Preprocessing, extract and format behavioral events

%Extract behavior info:
behavior_labels_subject = cell2mat({labels{:,3}}'); %Extract unique behavior info
behavior_labels_partner = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
block_labels = cell2mat({labels{:,10}}'); %Extract block info

%Compute the %of timepoints with unidentified behavior
%Important note to keep in mind: I still think that the timepoints where
%there were no behaviors assigned should be excluded...
no_behavs_timepoints = length(find(behavior_labels_subject == 28))/length(behavior_labels_subject);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
fprintf('Percent timepoints without a behavior assigned: %s \n', num2str(no_behavs_timepoints));
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%Create subject behavior regressors:
subject_behav_reg = zeros(size(Spike_count_raster,1), length(behav_categ)-1); %initialize
for b = 1:length(behav_categ)-1
    idx = find(behavior_labels_subject == b);
    subject_behav_reg(idx,b) = 1;
end

%Create partner behavior regressors:
partner_behav_reg = zeros(size(Spike_count_raster,1), length(behav_categ)-1); %initialize
for b = 1:length(behav_categ)-1
    idx = find(behavior_labels_partner == b);
    partner_behav_reg(idx,b) = 1;
end

%Create block regressor:
block_reg = zeros(size(Spike_count_raster,1), length(unique(block_labels))); %initialize
for b = 1:length(unique(block_labels))
    idx = find(block_labels == b);
    block_reg(idx,b) = 1;
end

% %Create behavior change regressor (exclude for now because we do not seem
% to be able to detect behavior shifts irrespective of the shift)
% subject_behav_change_reg = zeros(size(Spike_count_raster,1), 1); %initialize
% subject_behav_change_reg(find(diff(behavior_labels_subject)~=0)+1) = 1;

%Combine behavioral events
Behavior_events=[subject_behav_reg, partner_behav_reg, block_reg];

disp('Preprocessing - adding events done')

%% Set up event types

%Get variable names
BehaviorEventNames = [behav_categ(1:end-1), append('partner.',behav_categ(1:end-1)),string(block_times{:,"Behavior"})'];

%Set up event type here.
allEventTypes=ones(size(BehaviorEventNames))*3;
% For now I will assign all behavioral events as event type 3
%Event Type 1 = Whole trial
%Event Type 2 = from stimulus onset to the rest of the trial
%Event Type 3 = from before movement onset to after movement onset (pre & post time)
allEventsInfo = [BehaviorEventNames; allEventTypes]';

disp('Preprocessing - Assigning event types done')

%% Preprocessing - Selecting Analog Tracking
%For now we don't need to do anything here...

moveR = ME_final;

%% Setup Design Matrix - Regressor labels

regLabels = BehaviorEventNames';
regLabels(length(allEventsInfo)+1) = {'Motion Energy L'};
regLabels(length(allEventsInfo)+2) = {'Motion Energy R'};


disp('Setup Design Matrix - Regressor labels done')

%% Preprocessing - Grouping regressors for later cross-validation
%For now only include behavior and motion energy

regGroups = {'BehavioralEvents' 'Movements'};
regGroups{2,1} = BehaviorEventNames';
regGroups{2,2} = {'Motion Energy L' 'Motion Energy R'} ;

disp('Preprocessing - Grouping regressors for later cross-validation done')

%% Setup Design Matrix - Behavioral Events


%Creates task regressors with the time varying kernels as described in Churchland
[behavR, behavIdx] = log_makeDesignMatrix(temptaskEvents, behavkEventType, opts);

%% Setup Design Matrix - Movement
% Create movement events from analog traces

[dMat, traceOut, moveIdx] = log_analogToDesign(std_movR, nanstd(std_movR)*2, opts, opts.Fs , opts.Fs , motorIdx, 0, taskEventName);
% dMat = "digitized" (or binarized) movement regressors. One cell per movement variable, per trial
% traceOut = original analog traces
% moveIdx = indices of regressors to keep track of what movement variable a regressor belongs to.

%moveR contains the movement events (movements above a certain
%threshold), the original movement variables and ME.
moveR = [cat(1,temp_moveR{:}), traceOut, video_var_hold]; %moveR structure: Time X num mvmt regressors

disp('Setup Design Matrix - Movements done')

%% Deal with Nans in the behavior data
%I.e. moments with undefined behaviors.
% Question: should i remove timepoints where both the partner and the
% subject have undefined behaviors? only one of the two?

all_nan_inds = [];
percent_lost_to_nans = size(all_nan_inds,1)/size(Spike_count_raster,1)

%remove from predictors and neural data
moveR(all_nan_inds,:) = []; %blank these indices out
behavR(all_nan_inds,:) = []; %blank these indices out
zero_regressors_move = find(all(moveR == 0,1));
zero_regressors_behav = find(all(behavR == 0,1));

%remove empty regressors
if ~isempty(zero_regressors_behav)
    taskR(:,zero_regressors_behav) = [];
    taskIdx(zero_regressors_behav) = []; %CT: Make sure you remove regIdx every time regressors are removed from fullR
end
if ~isempty(zero_regressors_move)
    moveR(:,zero_regressors_move) = [];
    moveIdx(zero_regressors_move) = [];
end

%Remove the nan time points in regressors from the neuronal data
Vc(all_nan_inds,:) = [];


disp('Nans removed')

%2020-12-19 CT: Remove regressors with less than 10 events (as in
%Musall et. al. 2019)
low_events_idx= find(sum(behavR,1)<10);
behavkR(:,low_events_idx)=[];
behavIdx(low_events_idx)=[];

 %% Combine Design Matrix
    
    %Combine Design Matrix
    fullR=[taskR, moveR];
    
    %Collect all indecies together.
    regIdx =[taskIdx; moveIdx'];
        
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
    
    disp('Center and standardize continuous data Done')


%% Orthogonalization
    %2020-11-10 RWD note: I'm not sure if Churchland actually does the QR
    %analysis and throws out redundant regressors or does a QR decomp and
    %takes a specific number of columns of the Q matrix to replace the
    %design matrix...the langauge of their document is slightly unclear.
    %If this is the case we need to simply modifiy the run_QR code to
    %return more than rejIdx and remove the rejIdx portion of this code.
    
    %2020-11-30 RWD Update: Great, turns out I was right.  Below code
    %comments out rejecting regressors based on moveR and instead doing the
    %subtitution
    
    %2020-12-19 CT Update: After reading again the methods, we realized that task regressors
    %also need to be orthoganalized against. Video regressors are always orthogonalized, but for other
    %uninstructed movements that depends on what model we're running.
    
    %2020-12-26 CT Update: After double-checking with Simon Musall, uninstructed movements
    %are ALWAYS orthogonalized against instructed movements when computing explained variance
    %of a group of regressors (not for individual regressors since interpretability is reduced).
    %They also included task regressors that could be observed by the camera (e.g. when spouts
    %or handles were moving an visual stimuli presented). This is not the case in our data since i
    %cut out the visual stimuli presented altogether from the video.
    
    %Get all movement variables (i.e. everything except task regressors)
    allMoveR = fullR(:,length(taskIdx)+1:end); %get instructed + uninstrcuted movements
    allMove_RegIdx = regIdx(length(taskIdx)+1:end);
    disp('Orthogonalize uninstructed variables with respect to instructed movements via QR decomposition')
    
    %Run QR to orthogonalize uninstructed movement against instructed movement variables
    uninstrc_reg_id = find(strcmp(regLabels, 'time'))+1:length(regLabels);
    allMoveR_orthog = run_QR_decomp(allMoveR, allMove_RegIdx, uninstrc_reg_id);
    %Replace relevant columns in fullR
    fullR_orthog = fullR;
    regIdx_orthog = regIdx;
    fullR_orthog(:,length(taskIdx)+1:end) = allMoveR_orthog;
    
    %Run QR to orthogonalize video variables against all other movement variables
    video_reg_id = find(strcmp(regLabels, 'SVD_raw_video')|strcmp(regLabels, 'SVD_ME_video'));
    allMoveR_orthog = run_QR_decomp(allMoveR, allMove_RegIdx, video_reg_id);
    %Replace relevant columns in fullR
    orthog_idx = find(ismember(regIdx, video_reg_id));
    fullR(:,orthog_idx) = allMoveR_orthog(:,ismember(allMove_RegIdx, video_reg_id));
    
    disp('Orthogonalization through QR Done')
    
    
    %% run QR and check for rank-defficiency. This will show whether a given regressor is highly collinear with other regressors in the design matrix.
    
    %Run QR for fully orthogonalized uninstructed movement regressors
    disp('Start QR check for colinear regressors')
    rejIdx = run_QR(fullR_orthog,1);
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
    
    %Run QR for video regressors orthogonalized only
    disp('Start QR check for colinear regressors')
    rejIdx = run_QR(fullR,1);
    close all
    %2020-12-02 CT: Temporary Check - what are the redundant/degenerate regressors??
    %plot histogram of which regressors are degenerate
    figure;
    histHandle1 = histogram(regIdx,'BinWidth',1,'FaceAlpha',0.4); hold on
    histHandle1.BinEdges = histHandle1.BinEdges - histHandle1.BinWidth/2;
    histHandle1.BinEdges(length(histHandle1.BinEdges)+1)=max(regIdx)+0.5;
    if ~isempty(regIdx(rejIdx))
        histHandle2 = histogram(regIdx(rejIdx),'BinWidth',1,'FaceColor',[1 0 0],'FaceAlpha',0.4); %distribution of degenerate regressors
        histHandle2.BinEdges = histHandle2.BinEdges - histHandle2.BinWidth/2;
        histHandle2.BinEdges(length(histHandle2.BinEdges)+1)=max(regIdx(rejIdx))+0.5;
        title('Rejection of regressors in FullR')
        legend('all','rejected')
        %What reglabels do they correspond to?
        unique(regLabels(regIdx(rejIdx)))%which regressors are degenerate
        
        %Remove rejected regressors
        fullR(:,rejIdx) = [];
        regIdx(rejIdx,:) = [];
        figure;
        histHandle3 = histogram(regIdx,'BinWidth',1,'FaceAlpha',0.4); title('Regressors after deletion')
        histHandle3.BinEdges = histHandle3.BinEdges - histHandle3.BinWidth/2;
        histHandle3.BinEdges(length(histHandle3.BinEdges)+1)=max(regIdx)+0.5;
    end
    
    close all
    
    %% Fit model full model with cross validation
    
    disp('Start full model cross validation')
    tic %This is for timing how long it takes to run the entire analysis
    
    [Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR_orthog, Vc, regLabels, regIdx, regLabels, opts.folds);






%end
