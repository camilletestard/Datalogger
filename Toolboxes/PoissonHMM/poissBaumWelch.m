function [PI,A,B,alpha,beta,gamma,epsilons] = poissBaumWelch(spikes,nStates,dt,maxIter)
% Written by Ben Ballintyn (bbal@brandeis.edu) 12/2019.
%
% [PI,A,B,alpha,beta,gamma,epsilons] = poissBaumWelch(spikes,nStates,dt,maxIter)
% Estimate the initial state probabilities (PI), transition matrix (A), and
% emission matrix (B) using the forward and backward algorithms
% Inputs:
%   spikes - N x T matrix of spike counts. Each entry (i,j) holds the # of
%            spikes from neuron i in timebin j
%
%   nStates - # of hidden states predicted to have generated the spikes
%
%   dt - timebin size in seconds (e.g. .001)
%
%   maxIter - maximum number of updates to perform
%
% Outputs:
%   PI - nStates x 1 vector of initial state probabilities
%
%   A - nStates x nStates state transition matrix. Each entry (i,j) gives
%       the estimated probability of transitioning from state i to state j
%
%   B - N x nStates rate matrix. Each entry (i,j) gives the predicted rate
%       of neuron i in state j
%
%   alpha - output from the forward algorithm
%
%   beta - output from backward algorithm
%
%   gamma - nStates x T matrix where each entry (i,j) gives
%           P(Xj = i | o1...oT,PI,A,B)
%
%   epsilons - nStates x nStates x T tensor where each entry (i,j,k) gives
%              P(Xk = i, Xk+1 = j | o1...oT,PI,A,B)
nNeurons = size(spikes,1);
nTimeSteps = size(spikes,2);
poiss = @(lambda,n) ((lambda*dt).^n./factorial(n)).*exp(-lambda*dt);
minFR = (1/(nTimeSteps*dt));
% Initialize parameter estimates
PI = ones(nStates,1)*(1/nStates); % Initial state distribution
A = zeros(nStates,nStates); % Transition Matrix (columns = post, rows = pre)
for i=1:nStates
    for j=1:nStates
        if (i == j)
            A(i,j) = .99;
        else
            A(i,j) = .01/(nStates-1);
        end
    end
end

B = rand(nNeurons,nStates); % "Emission" (rate) matrix

notConverged = 1;
iterNum = 0;
while (notConverged && (iterNum < maxIter))
    % run forward algorithm
    [alpha,norms] = forward(spikes,nStates,dt,PI,A,B);
    % run backward algorithm
    [beta] = backward(spikes,nStates,dt,A,B,norms);
    % compute temporary variables
    gamma = zeros(nStates,nTimeSteps);
    epsilons = zeros(nStates,nStates,nTimeSteps-1);
    for t=1:nTimeSteps
        if (t < nTimeSteps)
            gamma(:,t) = (alpha(:,t).*beta(:,t))./sum(alpha(:,t).*beta(:,t));
            epsilonNumerator = zeros(nStates,nStates);
            for si=1:nStates
                for sj=1:nStates
                    epsilonNumerator(si,sj) = alpha(si,t)*A(si,sj)*beta(sj,t+1)*prod(poiss(B(:,sj),spikes(:,t+1)));
                end
            end
            epsilons(:,:,t) = epsilonNumerator./sum(epsilonNumerator(:));
        end
    end
    
    % store old paramters for convergence check
    oldPI = PI;
    oldA = A;
    oldB = B;
    % update parameters
    PI = gamma(:,1);
    Anumer = sum(epsilons,3);
    Adenom = sum(gamma,2);
    A = Anumer./Adenom;
    A = A./sum(A,2);
    B = ((spikes*gamma')./sum(gamma,2)')/dt;
    B = max(minFR,B);
    % check for convergence
    notConverged = isConverged(oldPI,oldA,oldB,PI,A,B);
    iterNum=iterNum+1;
    disp(['done with iteration #' num2str(iterNum)])
end
end

function notConverged = isConverged(oldPI,oldA,oldB,PI,A,B)
dPI = sqrt(sum((oldPI - PI).^2));
dA = sqrt(sum((oldA(:) - A(:)).^2));
dB = sqrt(sum((oldB(:) - B(:)).^2));
disp(['dPI = ' num2str(dPI) ' dA = ' num2str(dA) ' dB = ' num2str(dB)])
thresh = 1e-4;
if (dPI < thresh && dA < thresh && dB < thresh)
    notConverged = 0;
else
    notConverged = 1;
end
end