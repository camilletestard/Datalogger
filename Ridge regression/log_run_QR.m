function [rejIdx median_val min_val max_val] = run_QR(fullR, plotting)
%2020-11-30 RWD Note: this run_QR implements the CSA check.  Use
%run_QR_decomp to do QR decomposition and return the video and video ME
%files orthogonalized.

%Convenience function - task in design matrix, runs QR check for
%mulitcolineraity, and returns indecies of design matrix that can be
%rejected for being colinear.

%   Detailed explanation coming soon!

%plotting is a toggle for having this plot the results of the QR.  For now
%leaving as a required arguments with a >0 toggle.  Will make a parameter
%or setup input parser later.


rejIdx = false(1,size(fullR,2)); %initliaze rejIdx
[~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize normalized design matrix

if plotting > 0 %Plots the colinearity "score" for the regressors in order.
    
    figure; plot(abs(diag(fullQRR)),'linewidth',2); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
    axis square; ylabel('Norm. vector angle'); xlabel('Regressors');
    
end

if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
    temp = ~(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1)));
    fprintf('Design matrix is rank-defficient. Removing %d/%d additional regressors.\n', sum(temp), sum(~rejIdx));
    rejIdx(~rejIdx) = temp; %reject regressors that cause rank-defficint matrix
end

colin_values = abs(diag(fullQRR));
if plotting > 0 %Plots the distirbution of colinearity values
    colin_values(rejIdx(:)) = [];
    figure;
    histogram(colin_values)
    title('distribution of colinearity values after rejection')
end

median_val = median(colin_values);
min_val = min(colin_values);
max_val = max(colin_values);

%% posterity code
% false(1,size(fullR,2));
% [~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize normalized design matrix
% figure; plot(abs(diag(fullQRR)),'linewidth',2); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
% axis square; ylabel('Norm. vector angle'); xlabel('Regressors');
% if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
%     temp = ~(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1)));
%     fprintf('Design matrix is rank-defficient. Removing %d/%d additional regressors.\n', sum(temp), sum(~rejIdx));
%     rejIdx(~rejIdx) = temp; %reject regressors that cause rank-defficint matrix
% end
end

