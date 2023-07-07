function [beta] = backward(spikes,nStates,dt,A,B,norms)
% Written by Ben Ballintyn (bbal@brandeis.edu) 12/2019.
%
% [beta] = backward(spikes,nStates,dt,A,B,norms)
% Given the norms from the forward algorithm, compute beta(s,t) = P(ot+1...oT | Xt = s)
% Inputs:
%   spikes - N x T matrix of spike counts. Each entry (i,j) holds the # of
%            spikes from neuron i in timebin j
%
%   nStates - # of hidden states predicted to have generated the spikes
%
%   dt - timebin size in seconds (e.g. .001)
%
%   PI - nStates x 1 vector of initial state probabilities
%
%   A - nStates x nStates state transition matrix. Each entry (i,j) gives
%       the probability of transitioning from state i to state j
%
%   B - N x nStates rate matrix. Each entry (i,j) gives the predicted rate
%       of neuron i in state j
% 
% Outputs:
%   beta - nStates x T matrix of backward probabilities P(ot+1...oT | Xt = s)
poiss = @(lambda,n) ((lambda*dt).^n./factorial(n)).*exp(-lambda*dt); % formula for computing probability of each neurons spike count assuming Poisson spiking
nTimeSteps = size(spikes,2);
beta = zeros(nStates,nTimeSteps);
beta(:,end) = 1; % Initialize final beta to be 1 for all states
% Stepping backward in time, compute betas
for t=fliplr(1:(nTimeSteps-1))
    for s=1:nStates
        beta(s,t) = sum(beta(:,t+1).*A(s,:)'.*prod(poiss(B,spikes(:,t+1)),1)'); % making change (08/03)
    end
    beta(:,t) = beta(:,t)./norms(t+1); % scale each beta by norms from forward pass
end
end

