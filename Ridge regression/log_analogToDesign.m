function [dMat, traceOut, moveIdx, traceIn_temp] = log_analogToDesign(traceIn, stdThresh, opts, sourceRate, targRate, motorIdx, gaussShift, taskEventName)
% code to create a peri-event design matrix based on an analog trace. Trace
% should be continous and will be reshaped into a trial structure to create
% a design matrix for every trial individually.
% Inputs:   traceIn = analog trace.
%           stdThresh = threshold for event detection in SDUs
%           opts = strucgture containing information about correct trials, duration etc..
%           sourceRate = sampling rate of analog trace
%           targRate = sampling rate of design matrix.
%           motorIdx = index for peri-event matrix.
%           gaussShift = variable for subsampling in case model is used with gaussian convolution.
% Outputs:  dMat cell = structure of size ntrials X nMotorRegressors (e.g.: eyeX, eyeY, HeadX, HeadY...)
%           traceOut = the unmodified analog traces
%           moveIdx = indices of movement variable regressors belong to.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;histogram(traceIn(:,1),50)
traceIn = (traceIn - prctile(traceIn,1))./ nanstd(traceIn); %minimum values are at 0, signal in standard deviation units
% figure;histogram(traceIn(:,1),50)
%NOTE CT: this is done twice (before the function is called and in the function), but this does not
%impact the vector.
traceOut = traceIn; %return normalized analog trace
traceIn = traceIn > stdThresh; %take activity above supplied threshold as indicator for event occurence
traceIn = diff([zeros(1,size(traceIn,2)); traceIn]) == 1; %find event onsets

%Create trial structure. Note: trial length is different for each trial AND we are only considering
%non-ignored trials here (we got rid of the other triials + inter-trial interval before running this
%function)
trialCnt = length(opts.starting_inds);
indices = cumsum(opts.durations);
for ntrial = 1:trialCnt
    if ntrial == 1
        traceIn_temp{ntrial} = traceIn(1:indices(ntrial),:);
    else
        traceIn_temp{ntrial} = traceIn(indices(ntrial-1)+1:indices(ntrial),:);
    end
end
traceIn = traceIn_temp;

%find number of frames per trial. Useful when downsampling hasn't happened yet.
frames = opts.durations / (sourceRate/targRate); %second dimension is trials so first should be frames per trial when taking differences in sampling rate into account

dMat = cell(trialCnt,size(traceIn{1},2));
for iMotorReg = 1:size(traceIn{1},2) %for each analog regressor seperately
    moveIdx(iMotorReg,1:length(motorIdx)) = length(taskEventName)+iMotorReg;
    for iTrials = 1:trialCnt %for each trial
        trace = logical(histcounts(find(traceIn{iTrials}(:,iMotorReg)), 0: sourceRate/targRate : (sourceRate/targRate)*frames(iTrials)))'; %resample to imaging frame rate. This is the zero lag regressor.
        
        % create full design matrix
        cIdx = bsxfun(@plus,find(trace),motorIdx);
        cIdx(cIdx < 1) = 0;
        cIdx(cIdx > frames(iTrials)) = frames(iTrials);
        cIdx = bsxfun(@plus,cIdx,(0:frames(iTrials):frames(iTrials)*length(motorIdx)-1));
        cIdx(cIdx < 1) = frames(iTrials);
        cIdx(cIdx > (frames(iTrials) * length(motorIdx))) = frames(iTrials) * length(motorIdx);
        
        dMat{iTrials,iMotorReg} = false(frames(iTrials), length(motorIdx));
        dMat{iTrials,iMotorReg}(cIdx(:)) = true;
        dMat{iTrials,iMotorReg}(end,:) = false; %don't use last timepoint of design matrix to avoid confusion with indexing.
        dMat{iTrials,iMotorReg}(end,2:end) = dMat{iTrials,iMotorReg}(end-1,1:end-1); %replace with shifted version of previous timepoint
        
        if gaussShift > 1
            dMat{iTrials,iMotorReg} = dMat{iTrials,iMotorReg}(:,1:gaussShift:end);
        end
    end
end
end