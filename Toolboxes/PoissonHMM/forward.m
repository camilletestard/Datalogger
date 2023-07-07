function [alpha,norms] = forward(spikes,nStates,dt,PI,A,B)
% Written by Ben Ballintyn (bbal@brandeis.edu) 12/2019.
%
% [alpha,norms] = forward(spikes,nStates,dt,PI,A,B)
% Run forward algorithm to compute alpha = P(Xt = i | o1...ot,pi)
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
%   alpha - nStates x T matrix of forward probabilities. Each entry (i,j)
%           gives P(Xt = i | o1...oj,pi)
%
%   norms - 1 x T vector of norms used to normalize alpha to be a
%           probability distribution and also to scale the outputs of the
%           backward algorithm. norms(t) = sum(alpha(:,t))

poiss = @(lambda,n) ((lambda*dt).^n./factorial(n)).*exp(-lambda*dt); % formula for computing probability of each neurons spike count assuming Poisson spiking
nTimeSteps = size(spikes,2);
% For each state, use the initial state distribution and spike counts to
% initialize alpha(:,1)
for i=1:nStates
    alpha(i,1) = PI(i)*prod(poiss(B(:,i),spikes(:,1)));
end
norms(1) = sum(alpha(:,1));
alpha(:,1) = alpha(:,1)./norms(1); % normalize first timestep
% For each timestep, for each state, compute alpha
for t=2:nTimeSteps
    for s=1:nStates
        alpha(s,t) = prod(poiss(B(:,s),spikes(:,t)))*sum(alpha(:,t-1).*A(:,s));
    end
    norms(t) = sum(alpha(:,t));
    alpha(:,t) = alpha(:,t)./norms(t);
end
end

