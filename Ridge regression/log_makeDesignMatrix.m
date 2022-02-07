function [fullMat, eventIdx] = log_makeDesignMatrix(events, eventType, opts)
% function to generate design matrix from a column matrix with binary
% events. eventType defines the type of design matrix that is generated.
% (1 = fullTrial, 2 = post-event, 3 = peri-event)
%CT Note: "events" # cells = # trials considered. Each cell contains Time_trial X #regressors

%2020-05-12 reviewing side by side with tutorial...we end up with a lot
%less empty regressors than they do...so maybe we are thinking about or
%calculating something wrong.

fullMat = cell(1,length(eventType));
eventIdx = cell(1,length(eventType));
trialCnt = size(opts.durations,2); %Durations are for each trial thus this give us the trial number

for iRegs = 1 : length(eventType) %self note: something is going off for regressor 10, not sure why but can check in a bit.
    
    %This one is independent of trial duration
    if eventType(iRegs) == 3
        kernelIdx = [-(opts.mPreTime: -1 : 1) 0 (1:opts.mPostTime)]; %index for design matrix to cover pre- and post event activity
    elseif eventType(iRegs) == 2
        kernelIdx = 0:opts.sPostTime; %this is correct as kernels are same size and below will handle the trials being different lengths. 
    elseif eventType(iRegs) == 1
        kernelIdx = 0:opts.sPostTime; % 2020-06-10 CT: For now leave it like event type 2. Will need to adjust later.
    end
    
    dMat = cell(1,trialCnt);    
    for iTrials = 1 : trialCnt
        
        frames = opts.durations(iTrials);%Duration depend on the trial
        
%         %If eventtype == 1, whole trial.
%         if eventType(iRegs) == 1
%             kernelIdx = 0 : frames-1; #IMPORTANT 2020-06-10 CT: this coding does not work because
              % each trial is of a different length, the number of regressors will then vary per
              % trial and we will not be able to concatenate in the 1st dimension. 
%         end
        
        % Get the zero lag regressor.
        trace = logical(events{iTrials}(:,iRegs));
        
        % create full design matrix (time verying kernel.
        cIdx = bsxfun(@plus,find(trace),kernelIdx); %find time of event and add length of kernel (according to trial type)
        cIdx(cIdx < 1) = 0;% if less than 0, set to 0
        cIdx(cIdx > frames) = frames; %if more than length of trial "frames", set to max length of trial.
        %transform indices to fit matrix of time varying kernel later.
        cIdx = bsxfun(@plus,cIdx,(0:frames:frames*length(kernelIdx)-1));
        cIdx(cIdx < 1) = frames; 
        cIdx(cIdx > (frames * length(kernelIdx))) = frames * length(kernelIdx);
        
        dMat{iTrials} = false(frames, length(kernelIdx)); %initialize matrix
        dMat{iTrials}(cIdx(:)) = true; %Fill in with time varying kernel
        dMat{iTrials}(end,:) = false; %don't use last timepoint of design matrix to avoid confusion with indexing.
        dMat{iTrials}(end,2:end) = dMat{iTrials}(end-1,1:end-1); %replace with shifted version of previous timepoint 
       
    end
    %after doing all off the trials, remove empty regressors.
    
    fullMat{iRegs} = cat(1,dMat{:}); %Output: T_correct_trials X length kernel

    cIdx = sum(fullMat{iRegs},1) > 0; %don't use empty regressors
    fullMat{iRegs} = fullMat{iRegs}(:,cIdx);
    eventIdx{iRegs} = repmat(iRegs,sum(cIdx),1); %keep index on how many regressor were created
    
end

fullMat = cat(2,fullMat{:});
eventIdx = cat(1,eventIdx{:});

end