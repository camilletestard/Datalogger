function [moveMat, traceOut, moveIdx] = log_analogToDesign(traceIn, stdThresh, motorIdx, behavEventName)
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

%find number of frames in session 
frames = size(traceIn,1);

moveMat = cell(1,size(traceIn,2)); %Initialize event kernels matrix
moveIdx = cell(1,size(traceIn,2)); %Initialize regressor labels matrix

for iMotorReg = 1:size(traceIn,2) %for each analog regressor seperately
    
    numMotorReg = length(behavEventName)+iMotorReg;
        
        trace = logical(traceIn(:,iMotorReg)); %This is the zero lag regressor.

        % create full design matrix
        cIdx = bsxfun(@plus,find(trace),motorIdx);
        cIdx(cIdx < 1) = [];
        cIdx(cIdx > frames) = [];
        
        moveMat{iMotorReg} = false(frames, length(motorIdx)); %initialize matrix
        for c = 1:length(motorIdx)
            moveMat{iMotorReg}(cIdx(:,c),c) = true; %Fill in with time varying kernel; %NOT THE RIGHT INDEXING
        end
        
        cIdx = sum(moveMat{iMotorReg},1) > 0; %don't use empty regressors
        moveMat{iMotorReg} = moveMat{iMotorReg}(:,cIdx);
        moveIdx{iMotorReg} = repmat(numMotorReg,sum(cIdx),1); %keep index on how many regressor were created
        
end

moveMat = cat(2,moveMat{:});
moveIdx = cat(1,moveIdx{:});

end