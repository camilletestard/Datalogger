function [Y_hat,Betas,CIs,Rsquared_per,LL_per,pvalues_per,problem_neurons] = GLM_CrossVal(X,Y,kfolds,is_smooth)
%GLM_CROSSVAL Custom code to run the cross validation on our glm data
%   Detailed explanation to come.  For now writing this since the built in
%   function for matlab doesn't seem particularly useful and the other
%   function we have from T and T 2022 makes certain assumptions about the
%   format of the data we don't want to deal with.  Designed to work with
%   any of our glm set ups

%INPUTS:
%X - design matrix for these models
%Y - ALL spiking data (will design function to loop over all neurons)
%kfolds - number of folds to use, scalar



%OUTPUTS:


%% Preallocate variables Set up cross validation

%Note: GLM fit also generates a coeff for the intercept term so have to add
%1 to all the size(X,2) things.
Betas = NaN(size(X,2)+1,size(Y,2)); %Average over CV Beta for each regressors (row) for each neuron (column)
CIs = NaN(size(X,2)+1,size(Y,2),2); %Average CI for each regressor, for each neuron, first page lower bound, second upper bound;
Rsquared_per = NaN(1,size(Y,2));%Average R^2 for EACH NEURON MODEL averaged over CV
%For now doing adjusted R^2, but can change this in the CV loop
LL_per = NaN(1,size(Y,2));
pvalues_per = NaN(size(X,2)+1,size(Y,2),kfolds);
Y_hat = NaN(size(Y)); %Prediction for all time points for each neuron, concatenated over CV

rng(1) %For reproducability across runs




%Divide these into approximately kfolds groups (up to a rounding)
foldCnt = floor(size(Y,1)/kfolds);

%% Check each fold of the design matrix is full rank
%If not, keep redoing the permutaion until they are.

if kfolds >1

    good = false;
    counter = 0;
    while ~good
        counter = counter+1;
        display(['Making sure design matrid for all folds is full rank. Try #' num2str(counter)])
        randIdx = randperm(size(Y,1)); %Get randperm of indecies for whole session
        %Need to due permentation instead of just a partition since we have the
        %block design (i.e. don't want fits that are only block 1 or block 2)
        X_check = cell(1,kfolds);
        rank_check = false(1,kfolds);

        for iFolds = 1:kfolds
            dataIdx = true(1,size(Y,1)); %Start with using all of the session (i.e. all true indecies)


            dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false;

            X_check{1,iFolds} = X(dataIdx,:); %Put the design matrix with these random inds into X_check

            rank_check(iFolds) = rank(X_check{iFolds}) == rank(X); %Make sure that the subset matrix is full rank

        end
        rank_check
        sum(rank_check)
        good = all(rank_check); %If all are full rank, set good to true.  Otherwise redo loop

    end
end
%% Loop over each neuron

%Update 2022-08-29 might want to save more stuff than this, but I think
%this is okay for a first pass.
problem_neurons = false(1,size(Y,2));


for n = 1:size(Y,2) %First select neuron to use then do CV proceedure for each neuron separately.
%Update 2022-08-29 I think this will make the loop easier to manage
display(['Working on neuron #' num2str(n)])
Y_cur = Y(:,n); %Select neuron
cBeta = cell(2,kfolds); %One for the beta and one for the confidence interval
cPvalues = NaN(size(X,2)+1,kfolds); %Save the P value for each coefficent for each fold.  Combine later with Fisher's method
cR2 = NaN(1,kfolds); %Save the R^2 value for each neurons fit.
cLL = NaN(1,kfolds); %Save the LL value for each neurons fit.
Y_cur_hat = NaN(size(Y_cur)); %Predicted fits that will be filled in with each fold.
lastwarn('','') %empty warning tracker before each model fit
% Run cross validation loop
    for iFolds = 1:kfolds
        display(['CV fold #' num2str(iFolds)])
        dataIdx = true(1,size(Y_cur,1)); %Start with using all of the session (i.e. all true indecies)
        
        if kfolds >1
        %Set the indecies of the TEST set in the fold to false
        %Only left with the training indecies as true
        dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false;
        end
        
        if is_smooth %Use Gaussian 
        
            mdl = fitglm(X(dataIdx,:),Y_cur(dataIdx),'linear','Distribution','normal'); 
            warning_thrown = lastwarn;
        else %Use Poisson
            
            mdl = fitglm(X(dataIdx,:),Y_cur(dataIdx), 'linear','Distribution','poisson');
            warning_thrown = lastwarn; 
        end
        
        %Put results from each fold in a storage variable
        cBeta{1,iFolds} = mdl.Coefficients.Estimate; %Estimates of coefficents
        cBeta{2,iFolds} = mdl.coefCI; %In the second cell put the confidence intervals
        cPvalues(:,iFolds) = mdl.Coefficients.pValue; %P value for t-test of significance of coefficent
        cR2(1,iFolds) = mdl.Rsquared.Adjusted; %R squared value adjusted for dfe for each neuron and its data.
        cLL(1,iFolds) = mdl.LogLikelihood;% LL for this fit
        
        if kfolds>1
        Y_cur_hat(~dataIdx) = mdl.feval(X(~dataIdx,:)); %Predict from the data that WAS NOT in the training set.
        else
        Y_cur_hat = mdl.feval(X);
        end
    end
    %Check warning at the end, if not empty then indicate cell was a
    %problem cell.
    
    if ~isempty(warning_thrown)
        
        problem_neurons(n) = true;
        
    end
    
    %Collect results for this neuron across folds
    Betas(:,n) = mean(horzcat(cBeta{1,:}),2); %Take mean value across fold
    %Note when we add back actual folds this may not work and we may have
    %to turn this into a for loop or otherwise deal with it.
    CIs(:,n,1) = mean(horzcat(cBeta{2,:}(:,1)),2);
    CIs(:,n,2) = mean(horzcat(cBeta{2,:}(:,2)),2);
    
  
    
    pvalues_per(:,n) = mean(cPvalues,2);   %CHANGE THIS TO FISHERS METHOD!!!!!
    %USING AVERAGE FOR NOW FOR CONVENIENCE
    
    Rsquared_per(n) = mean(cR2,2);
    LL_per(n) = mean(cLL,2);
    Y_hat(:,n) = Y_cur_hat;
    
    
    
    
end
%% Results are collected sufficently above I think, so just put those as the output arguments.
end

