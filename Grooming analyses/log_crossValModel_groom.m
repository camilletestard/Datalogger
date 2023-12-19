function [Vm, cBeta, cR, subIdx, cRidge, cLabels] =  crossValModel(fullR, Vc, cLabels, regIdx, regLabels, kfolds)
% function to compute cross-validated R^2
%Note: cLabels = labels for regressors of interest. E.g. only task regressors
% regLabels = all regressors
% fullR = all regressor matrix; Vc = neural data; kfolds = #kfolds for cross-validation

cIdx = ismember(regIdx, find(ismember(regLabels,cLabels))); %get index for regressors of interest
cLabels = regLabels(sort(find(ismember(regLabels,cLabels)))); %make sure Labels for regressors of interest are in ascending order

%create new regressor index for labels of interest only ("cLabels")
subIdx = regIdx; %create sub-index vector
subIdx = subIdx(cIdx); %only select regressors of interest
temp = unique(subIdx); %find the unique regressors 
% check regressors included: cLabels(temp)
%Re-label subIdx for the regressors of interest 
for x = 1 : length(temp) %for all regressors of interest
    subIdx(subIdx == temp(x)) = x;
end
cR = fullR(:,cIdx);%sub-select design matrix only including regressors of interest

%Vc = Vc'; %CT change dimensions --> cross validate over TIME not neurons.
Vm = zeros(size(Vc),'single'); %pre-allocate reconstructed V
% rng(1) % for reproducibility

%NEW CROSS_VALIDATION METHOD
%By chunk of time rather than randomly sampling across time
if kfolds>2
    Q = quantile(1:length(Vc),kfolds-1);
    bounds = [1 round(Q) length(Vc)];
elseif kfolds==2
    Q=round(length(Vc)/2);
    bounds = [1 round(Q) length(Vc)];
elseif kfolds==1
    bounds =[1 length(Vc)];
end

for fold=1:kfolds
    indices(bounds(fold):bounds(fold+1))=fold;
end

%OLD CROSS_VALIDATION METHOD
%Create random permutation for cross-validation:
%randIdx = randperm(size(Vc,2)); %generate randomly permuted indices for time

foldCnt = floor(size(Vc,2) / kfolds);
cBeta = cell(1,kfolds); %create nkfolds empty cells

for iFolds = 1:kfolds %for all kfolds
    dataIdx = true(1,size(Vc,2)); %initialize logical vector - considers all time
    
    if kfolds > 1 %if there is more than one fold for cross-validation
          
        %Train with a subset of the neural data
        dataIdx(indices==iFolds) = false; %index for training data, do not include them in testing.
        if iFolds == 1 %ony record ridge penalities for the first fold
            [cRidge, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', cR(dataIdx,:), true); %get beta weights and ridge penalty for model including regressors of interest only
            %disp(['Fold ' num2str(iFolds) 'done.'])

        else %no need to keep track of rigde penalties for other kfolds
            [~, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', cR(dataIdx,:), true, cRidge); %get beta weights. ridge value should be the same as in the first run.
        end
        
        Vm(:,~dataIdx) = (cR(~dataIdx,:) * cBeta{iFolds})'; %predict remaining data
        
%         if rem(iFolds,kfolds/5) == 0
%             fprintf(1, 'Current fold is %d out of %d\n', iFolds, kfolds);
%         end
    else
        [cRidge, cBeta{iFolds}] = ridgeMML(Vc', cR, true); %get beta weights for model including regressors of interest only.
        Vm = (cR * cBeta{iFolds})'; %predict remaining data
        disp('Ridgefold is <= 1, fit to complete dataset instead');
    end
end

% % computed all predicted variance
% Vc = reshape(Vc,size(Vc,1),[]);
% Vm = reshape(Vm,size(Vm,1),[]);
% if length(size(U)) == 3
%     U = arrayShrink(U, squeeze(isnan(U(:,:,1))));
% end
% covVc = cov(Vc');  % S x S
% covVm = cov(Vm');  % S x S
% cCovV = bsxfun(@minus, Vm, mean(Vm,2)) * Vc' / (size(Vc, 2) - 1);  % S x S
% covP = sum((U * cCovV) .* U, 2)';  % 1 x P
% varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
% varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
% stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
% cMap = gather((covP ./ stdPxPy)');
% 
% % movie for predicted variance
% cMovie = zeros(size(U,1),frames, 'single');
% for iFrames = 1:frames
%     
%     frameIdx = iFrames:frames:size(Vc,2); %index for the same frame in each trial
%     cData = bsxfun(@minus, Vc(:,frameIdx), mean(Vc(:,frameIdx),2));
%     cModel = bsxfun(@minus, Vm(:,frameIdx), mean(Vm(:,frameIdx),2));
%     covVc = cov(cData');  % S x S
%     covVm = cov(cModel');  % S x S
%     cCovV = cModel * cData' / (length(frameIdx) - 1);  % S x S
%     covP = sum((U * cCovV) .* U, 2)';  % 1 x P
%     varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
%     varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
%     stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
%     cMovie(:,iFrames) = gather(covP ./ stdPxPy)';
%     clear cData cModel
%     
% end
% fprintf('Run finished. RMSE: %f\n', median(cMovie(:).^2));

end