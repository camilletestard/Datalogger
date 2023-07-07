function [bestPath,maxPathLogProb,T1,T2] = poissViterbi(spikes,PI,A,B,dt)
% Written by Ben Ballintyn (bbal@brandeis.edu) 12/2019.
% 
% [bestPath,maxPathProb,T1,T2] = myViterbi(spikes,PI,A,B,dt)
% Use a lightly modified viterbi algorithm to compute the most likely path
% through a sequence of hidden states given spike counts from N neurons and
% an initial state distribution, transition matrix, and rate matrix
% estimated by the modified Baum-Welch algorithm in poissBaumWelch.m
%
% Inputs:
%   spikes - N x T matrix of spike counts. Each entry (i,j) holds the # of
%            spikes from neuron i in timebin j
%
%   PI - nStates x 1 vector of initial state probabilities
%
%   A - nStates x nStates state transition matrix. Each entry (i,j) gives
%       the probability of transitioning from state i to state j
%
%   B - N x nStates rate matrix. Each entry (i,j) gives the predicted rate
%       of neuron i in state j
%
%   dt - timebin size in seconds (e.g. .001)
%
% Outputs:
%   bestPath - 1 x T sequence of states representing the most likely hidden 
%              state sequence
%
%   maxPathLogProb - log probability of most likely state sequence
%
%   T1 - nStates x T matrix where each entry (i,j) gives the log
%        probability of the most likely path so far ending in state i that
%        generates o1...oj
%
%   T2 - nStates x T matrix of back pointers where each entry (i,j) gives
%        the state xj-1 on the most likely path so far ending in state i
if (size(A,1) ~= size(A,2))
    error('Transition Matrix is not square')
else
    nStates = size(A,1);
end
poiss = @(lambda,n) ((lambda*dt).^n./factorial(n)).*exp(-lambda*dt);
nTimeSteps = size(spikes,2);

T1 = zeros(nStates,nTimeSteps);
T2 = zeros(nStates,nTimeSteps);
for i=1:nStates
    T1(i,1) = log(PI(i)) + log(prod(poiss(B(:,i),spikes(:,1))));
    T2(i,1) = 0;
end

for j=2:nTimeSteps
    for i=1:nStates
        vec1 = T1(:,j-1) + log(A(:,i)) + log(prod(poiss(B(:,i),spikes(:,j))));
        vec2 = T1(:,j-1) + log(A(:,i));
        T1(i,j) = max(vec1);
        [~,ind] = max(vec2);
        T2(i,j) = ind;
    end
end

[maxPathLogProb,bestPathEndState] = max(T1(:,end));
bestPath = zeros(1,nTimeSteps);
bestPath(end) = bestPathEndState;
for t=fliplr(1:(nTimeSteps-1))
    bestPath(t) = T2(bestPath(t+1),t+1);
end
end

