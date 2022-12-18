function [Results] = modelCorr(Vc,Vm)
%Update 2020-06-02: Found code in neuropixel file suggesting that R^2 is
%simply a corr2 between the predicted firing rate by the model and the
%actual firing rate, i.e. just the overall correlation between the two Y
%matrices.   This doesn't really need to be it's own function, but
%just going to leave it like that for sake of consistency in
%implementation/to give flexibilty to add to this if we want to.

%For now going to output results structure.  If this doesn't play nicely
%with other things change.  Doing this for now because it gives us
%flexibility in adding stuff to this function.

%2020-06-03 adding check for higher dimensional arrays so can use this same
%function for the shuffled regressors and such

if numel(size(Vm)) > 2 %check if have higher dimensional model predictions
    
    r_value = nan(1,size(Vm,4));
    
    for iRegs = 1:size(Vm,4) %hard coding using fact that using 4th dimension for regressors in allV
        
        r_value(iRegs) = corr2(Vc, Vm(:,:,1,iRegs)'); %note: added transpose here for multidimensional matrices
        
        
    end
%Corr2 centers the data and divides by the std so don't have to worry about that. 
else
    
r_value = corr2(Vc,Vm);

end

Results.r_value = r_value; %put into results struct

end


%% Posterity code (i.e. code we aren't using anymore)
% test = Vc(:,1); test_mean = mean(test); sum(test-test_mean)
% 
% Vc_squared_error = sum((Vc- mean(Vc,1)).^2);
% Vm_squared_error = sum((Vm- mean(Vc,1)).^2);
% R_squared = Vm_squared_error./Vc_squared_error 
% 
% covVc = cov(Vc');  % S x S
% covVm = cov(Vm');  % S x S
% cCovV = bsxfun(@minus, Vm, mean(Vm,2)) * Vc' / (size(Vc, 2) - 1);  % S x S