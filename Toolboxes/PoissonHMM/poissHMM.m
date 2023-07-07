function [bestPath,maxPathLogProb,PI,A,B,gamma] = poissHMM(spikes,nStates,dt,maxIter)
% Written by Ben Ballintyn (bbal@brandeis.edu) 12/2019.
%
% [bestPath,maxPathLogProb,PI,A,B,gamma] = poissHMM(spikes,nStates,dt,maxIter)
% Runs a modified Baum-Welch and Viterbi algorithm to compute the most
% likely sequence of hidden states that generated a set of spike trains
%
% Inputs:
% spikes  - the N x T matrix of spikes (N neurons and T timesteps) where
%           spikes are indicated by 1's and all other elements are 0
%
% nStates - The number of hidden states expected
%
% dt      - the size (in seconds) of each timestep
%
% maxIter - maximum number of iterations through the data
%
% Outputs:
%   bestPath - 1 x T sequence of states representing the most likely hidden 
%              state sequence
%
%   maxPathLogProb - log probability of most likely state sequence
%
%   PI - nStates x 1 vector of initial state probabilities
%
%   A - nStates x nStates state transition matrix. Each entry (i,j) gives
%       the estimated probability of transitioning from state i to state j
%
%   B - N x nStates rate matrix. Each entry (i,j) gives the predicted rate
%       of neuron i in state j
%
%   gamma - nStates x T matrix where each entry (i,j) gives
%           P(Xj = i | o1...oT,PI,A,B)

% Estimate state transition matrix and firing rates
[PI,A,B,alpha,beta,gamma,epsilons] = poissBaumWelch(spikes,nStates,dt,maxIter);

% Get most likely state sequence
[bestPath,maxPathLogProb,T1,T2] = poissViterbi(spikes,PI,A,B,dt);
end

