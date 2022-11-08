function CV = getCV(x)
% 
% call:
% 
%      CV = getCV(x)
%      
% 
% compute the coefficient of variation (C.V.) of the input vector x.
% The C.V. is defined as std(x)/mean(x), and as such is a measure of 
% the relative variability. The function ignores NaNs.
% 
% INPUT
% 
%      x   :  input vector
%      
% OUTPUT
% 
%      CV  :  coefficient of variation (a scalar)
%      
%      
% R. G. Bettinardi
% -------------------------------------------------------------------

if ~isvector(x)
    error('x has to be a vector!')
end

CV = nanstd(x)/nanmean(x);
     