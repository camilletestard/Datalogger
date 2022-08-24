function [fullMat, eventIdx] = log_makeDesignMatrix(events, eventType, opts)
% function to generate design matrix from a column matrix with binary
% events. eventType defines the type of design matrix that is generated.
% (1 = fullTrial, 2 = post-event, 3 = peri-event)
%CT Note: "events" # cells = # trials considered. Each cell contains Time_trial X #regressors


fullMat = cell(1,length(eventType)); %Initialize event kernels matrix
eventIdx = cell(1,length(eventType)); %Initialize regressor labels matrix

for iRegs = 1 : length(eventType)
    
    %This one is independent of trial duration
    if eventType(iRegs) == 3
        kernelIdx = [-(opts.mPreTime: -1 : 1) 0 (1:opts.mPostTime)]; %index for design matrix to cover pre- and post event activity
    elseif eventType(iRegs) == 2
        kernelIdx = 0:opts.sPostTime; %this is correct as kernels are same size and below will handle the trials being different lengths.
    elseif eventType(iRegs) == 1
        kernelIdx = 0:opts.sPostTime; % 2020-06-10 CT: For now leave it like event type 2. Will need to adjust later.
    elseif eventType(iRegs) == 0 %2022-08-22 RWD: created default case if put in zero where there are no delays added.
        kernelIdx = 0;  %Setting to zero works!
    end
    
    frames = size(events,1);%Duration of session
    
    
    % Get the zero lag regressor.
    trace = logical(events(:,iRegs));
    
    % create full design matrix time verying kernel.
    cIdx = bsxfun(@plus,find(trace),kernelIdx); %find time of original event and add length of kernel (according to trial type)
    cIdx(cIdx < 1) = [];% if less than 0, delete
    cIdx(cIdx > frames) = [];% if more than 0, delete
    
    fullMat{iRegs} = false(frames, length(kernelIdx)); %initialize matrix
    for c = 1:length(kernelIdx)
        fullMat{iRegs}(cIdx(:,c),c) = true; %Fill in with time varying kernel; 
    end
    
    cIdx = sum(fullMat{iRegs},1) > 0; %don't use empty regressors
    fullMat{iRegs} = fullMat{iRegs}(:,cIdx);
    eventIdx{iRegs} = repmat(iRegs,sum(cIdx),1); %keep index on how many regressor were created
    
end

non_empty_cells = find(~cellfun(@isempty,fullMat));
fullMat_nonempty = fullMat(non_empty_cells);
eventIdx_nonempty = eventIdx(non_empty_cells);

fullMat = cat(2,fullMat_nonempty{:});
eventIdx = cat(1,eventIdx_nonempty{:});

end