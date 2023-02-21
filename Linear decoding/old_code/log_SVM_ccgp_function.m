function [hitrate, C] = log_SVM_ccgp_function(trainlbls, traindata, testlbls, testdata, ldaflag, Normalize)
% This script was written by SdT to run SVM on the w8a dataset using the LIBSVM library
%
% Input_matrix = double(Input_matrix); %Augment data to double precision
% Labels = double(Labels);

if ~exist('ldaflag', 'var') %Default option is not to use lda. lda helps with lazy neurons, but overall leads to lower accuracy
    ldaflag = 0;
end

if ~exist('Normalize', 'var') %Default option is not to normalize
    Normalize = 0;
end

if Normalize == 1
    %Z-score data to normalize it. Each neuron FR is normalized to its mean across trials FR
    Input_matrix = zscore(Input_matrix);
end

total_number_trials = size(testdata,1); %Find the total number of trials for which a prediction needs to be made.
Predicted_labels = zeros(size(testlbls)); %Initate matrix of predicted labels for each trial
cumError = 0; %Count the number of prediction errors

if ldaflag == 1 %If you wish to perform LDA
    ldaopts.FisherFace = 1;
    [fdaout, Weights] = fda(trainlbls, ldaopts, traindata);
    testdata = testdata*Weights; % Project testdata to same LDA subspace as traindata:
    traindata = fdaout;
end

% Train/test SVM model:
% !!!!!!!!!! CHANGE SVM KERNEL HERE !!!!!!!!!!!!!!!
model = svmtrain(trainlbls, traindata, '-t, 0, -q'); %train the model using a linear kernel (-t: 0) or a RBF kernel (-t: 2) and default parameters
[svmlbls] = svmpredict(testlbls, testdata, model, '-q'); %get predicted labels given model

nErr= length(find( (testlbls - svmlbls) ~= 0 )); %Find misclassifications
cumError = cumError + nErr; %Count number of errors
Predicted_labels = svmlbls; %Keep track of the predicted labels


%Compute performance of decoder
hitrate = 1 - (cumError/total_number_trials);
%Obtain confusion matrix of predicted against real values.
C = confusionmat(testlbls, Predicted_labels);
C = C ./ repmat(sum(C,2), 1, size(C,2));

end
