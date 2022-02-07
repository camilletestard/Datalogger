function [X_r_mod] = run_QR_decomp(X_r, idx, labels)
%2020-11-30 RWD: This function actually performs the QR Decomposition and
%is meant to be used with the motion and video only design matrix (i.e. X_r(educed)) 

%qr Returns Q then R, don't need R so just tildea it out
[X_r_Q, ~] = qr(bsxfun(@rdivide,X_r,sqrt(sum(X_r.^2))),0); %orthogonalize normalized design matrix
%[X_r_Q, ~] = qr(bsxfun(@rdivide,fullR,sqrt(sum(X_r.^2))),0); %orthogonalize normalized design matrix

%There is a 1:1 equivalence of columns, so simply sub the original moveR 
% columns with the orthogonalized ones.  Return the altered matrix.

%2020-12-02: CT added
X_r_mod=X_r;
orthog_idx = find(ismember(idx, labels));
X_r_mod(:,orthog_idx)  = X_r_Q(:,orthog_idx);
%X_r_mod(:,size(X_r,2)-399:size(X_r,2))  = X_r_Q(:,size(X_r,2)-399:size(X_r,2));

end

