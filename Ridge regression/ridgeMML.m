function [L, betas, convergenceFailures] = ridgeMML(Y, X, recenter, L)
% [lambdas, betas, convergenceFailures] = ridgeMML(Y, X [, recenter] [, lambdas])
% 
% This is an implementation of Ridge regression with the Ridge parameter
% lambda determined using the fast algorithm of Karabatsos 2017 (see
% below). I also made some improvements, described below.
% 
% Inputs are Y (the outcome variables) and X (the design matrix, aka the
% regressors). Y may be a matrix. X is a matrix with as many rows as Y, and
% should *not* include a column of ones.
% 
% A separate value of lambda will be found for each column of Y.
% 
% Outputs are the lambdas (the Ridge parameters, one per column of Y); the
% betas (the regression coefficients, again with columns corresponding to
% columns of Y); and a vector of logicals telling you whether fminbnd
% failed to converge for each column of y (this happens frequently).
% 
% If recenter is 1 (default), the columns of X and Y will be recentered at
% 0. betas will be of size:  size(X, 2) x size(Y, 2)
% To reconstruct the recentered Y, use:
% YRecenteredHat = (X - mean(X, 1)) * betas;
% 
% If recenter is 0, the columns of X and Y will not be recentered. betas
% will be of size:  size(X, 2)+1 x size(Y, 2)
% The first row corresponds to the intercept term.
% To reconstruct Y, you may therefore use either:
% YHat = [ones(size(X, 1), 1), X] * betas;
% or
% YHat = X * betas(2:end, :) + betas(1, :);
% 
% If lambdas is supplied, the optimization step is skipped and the betas
% are computed immediately. This obviously speeds things up a lot.
% 
% 
% TECHNICAL DETAILS:
% 
% To allow for un-centered X and Y, it turns out you can simply avoid
% penalizing the intercept when performing the regression. However, no
% matter what it is important to put the columns of X on the same scale
% before performing the regression (though Matlab's ridge.m does not do
% this if you choose not to recenter). This rescaling step is undone in the
% betas that are returned, so the user does not have to worry about it. But
% this step substantially improves reconstruction quality.
% 
% Improvements to the Karabatsos algorithm: as lambda gets large, local
% optima occur frequently. To combat this, I use two strategies. First,
% once we've passed lambda = 25, we stop using the fixed step size of 1/4
% and start using an adaptive step size: 1% of the current lambda. (This
% also speeds up execution time greatly for large lambdas.) Second, we add
% boxcar smoothing to our likelihood values, with a width of 7. We still
% end up hitting local minima, but the lambdas we find are much bigger and
% closer to the global optimum.
% 
% Source: "Marginal maximum likelihood estimation methods for the tuning
% parameters of ridge, power ridge, and generalized ridge regression" by G
% Karabatsos, Communications in Statistics -- Simulation and Computation,
% 2017. Page 6.
% http://www.tandfonline.com/doi/pdf/10.1080/03610918.2017.1321119
% 
% Written by Matt Kaufman, 2018. mattkaufman@uchicago.edu


%% Optional arguments

if ~exist('recenter', 'var')
  recenter = 1;
end

if ~exist('L', 'var') || isnan(L(1))
  computeL = 1;
else
  computeL = 0;
end


%% Error checking

if size(Y, 1) ~= size(X, 1)
  error('Size mismatch');
end


%% Ensure Y is zero-mean
% This is needed to estimate lambdas, but if recenter = 0, the mean will be
% restored later for the beta estimation

if computeL || recenter
  YMean = mean(Y, 1);
  Y = bsxfun(@minus, Y, YMean);
end

pY = size(Y, 2);


%% Renorm (Z-score)

XStd = std(X, 0, 1);
X = bsxfun(@rdivide, X, XStd);

if computeL || recenter
  XMean = mean(X, 1);
  X = bsxfun(@minus, X, XMean);
  X(isnan(X)) = 0;
end


%% Optimize lambda

if computeL

  
  %% SVD the predictors
  
  [U, S, V] = svd(X, 0);
  
  
  %% Find the valid singular values of X, compute d and alpha
  
  n = size(X, 1);  % Observations
  p = size(V, 2);  % Predictors
  
  d = diag(S);
  
  % Find the number of good singular values. Ensure numerical stability.
  q = sum(d' > eps(U(1)) * (1:p));
  
  d2 = d .^ 2;
  
  % Equation 1
  % Eliminated the diag(1 ./ d2) term: it gets cancelled later and only adds
  % numerical instability (since later elements of d may be tiny).
  % alph = V' * X' * Y;
  alph = S * U' * Y;
  alpha2 = alph .^ 2;
  
  
  %% Compute variance of y
  % In Equation 19, this is shown as y'y
  
  YVar = sum(Y .^ 2, 1);
  
  
  %% Compute the lambdas
  
  L = NaN(1, pY);
  
  convergenceFailures = false(1, pY);
  for i = 1:pY
    [L(i), flag] = ridgeMMLOneY(q, d2, n, YVar(i), alpha2(:, i));
    convergenceFailures(i) = (flag < 1);
  end
  
else
  p = size(X, 2);
end


%% If requested, perform the actual regression

if nargout > 1
  
  if ~recenter
    % Restore the means of X and Y (but don't rescale)
    if computeL
      Y = bsxfun(@plus, Y, YMean);
      X = bsxfun(@plus, X, XMean);
    end
    
    betas = NaN(p + 1, pY);
    
    % Augment X with a column of ones, to allow for a non-zero intercept
    % (offset). This is what we'll use for regression, without a penalty on
    % the intercept column.
    X = [ones(size(X, 1), 1), X];
    
    XTX = X' * X;
    
    % Prep penalty matrix    
    ep = eye(p + 1);
    ep(1, 1) = 0;  % No penalty for intercept column
    
    % For renorming the betas
    % The 1 is so we don't renorm the intercept column.
    % Note that the rescaling doesn't alter the intercept.
    renorm = [1; XStd'];
    
  else
    betas = NaN(p, pY);
    
    % You would think you could compute X'X more efficiently as VSSV', but
    % this is numerically unstable and can alter results slightly. Oh well.
    % XTX = V * bsxfun(@times, V', d2);
    XTX = X' * X;
    
    % Prep penalty matrix
    ep = eye(p);

    % For renorming the betas
    renorm = XStd';
  end
  
 
  % Compute X' * Y all at once, again for speed
  XTY = X' * Y;
  
  % Compute betas for renormed X
  for i = 1:pY
    betas(:, i) = ((XTX + L(i) * ep) \ XTY(:, i));
  end
  
  % Adjust betas to account for renorming.
  betas = bsxfun(@rdivide, betas, renorm);
  betas(isnan(betas)) = 0;
end


%% Display fminbnd failures

if computeL && nargout < 3 && sum(convergenceFailures) > 0
  fprintf('fminbnd failed to converge %d/%d times\n', sum(convergenceFailures), pY);
end



function [L, flag] = ridgeMMLOneY(q, d2, n, YVar, alpha2)
% Compute the lambda for one column of Y

% Width of smoothing kernel to use when dealing with large lambda
smooth = 7;

% Value of lambda at which to switch from step size 1/4 to step size L/stepDenom.
% Value of stepSwitch must be >= smooth/4, and stepSwitch/stepDenom should
% be >= 1/4.
stepSwitch = 25;
stepDenom = 100;


%% Set up smoothing

% These rolling buffers will hold the last few values, to average for
% smoothing
smBuffer = NaN(1, smooth);
testValsL = NaN(1, smooth);

% Initialize index of the buffers where we'll write the next value
smBufferI = 0;


%% Evaluate the log likelihood of the data for increasing values of lambda
% This is step 1 of the two-step algorithm at the bottom of page 6.
% Basically, increment until we pass the peak. Here, I've added trying
% small steps as normal, then switching over to using larger steps and
% smoothing to combat local minima.

% Mint the negative log-likelihood function
NLLFunc = mintNLLFunc(q, d2, n, YVar, alpha2);

% Loop through first few values of k before you apply smoothing.
% Step size 1/4, as recommended by Karabatsos
done = 0;
NLL = Inf;
for k = 0:stepSwitch * 4
  smBufferI = mod(smBufferI, smooth) + 1;
  prevNLL = NLL;
  
  % Compute negative log likelihood of the data for this value of lambda
  NLL = NLLFunc(k / 4);
  
  % Add to smoothing buffer
  smBuffer(smBufferI) = NLL;
  testValsL(smBufferI) = k / 4;
  
  % Check if we've passed the minimum
  if NLL > prevNLL
    % Compute limits for L
    minL = (k - 2) / 4;
    maxL = k / 4;
    done = 1;
    break;
  end
end


% If we haven't already hit the max likelihood, continue increasing lambda,
% but now apply smoothing to try to reduce the impact of local minima that
% occur when lambda is large

% Also increase step size from 1/4 to L/stepDenom, for speed and robustness
% to local minima
if ~done
  L = k / 4;
  NLL = mean(smBuffer);
  while ~done
    L = L + L / stepDenom;
    smBufferI = mod(smBufferI, smooth) + 1;
    prevNLL = NLL;
    
    % Compute negative log likelihood of the data for this value of lambda,
    % overwrite oldest value in the smoothing buffer
    smBuffer(smBufferI) = NLLFunc(L);
    testValsL(smBufferI) = L;
    NLL = mean(smBuffer);
    
    % Check if we've passed the minimum
    if (NLL > prevNLL)
      % Adjust for smoothing kernel (walk back by half the kernel)
      smBufferI = smBufferI - (smooth - 1) / 2;
      smBufferI = smBufferI + smooth * (smBufferI < 1); % wrap around
      maxL = testValsL(smBufferI);
      
      % Walk back by two more steps to find min bound
      smBufferI = smBufferI - 2;
      smBufferI = smBufferI + smooth * (smBufferI < 1); % wrap around
      minL = testValsL(smBufferI);
      
      done = 1;
    end
  end
end


%% Bounded optimization of lambda
% This is step 2 of the two-step algorithm at the bottom of page 6. Note
% that Karabatsos made a mistake when describing the indexing relative to
% k*, which is fixed here (we need to go from k*-2 to k*, not k*-1 to k*+1)

[L, ~, flag] = fminbnd(NLLFunc, max(0, minL), maxL, ...
  optimset('Display', 'off'));


function NLLFunc = mintNLLFunc(q, d2, n, YVar, alpha2)
% Mint an anonymous function with L as the only input parameter, with all
% the other terms determined by the data.

% Equation 19 
% We've modified the math here to eliminate the d^2 term from both alpha
% (Equation 1, in main function) and here (Equation 19), because they
% cancel out and add numerical instability.
NLLFunc = @(L) -(q * log(L) - sum(log(L + d2(1:q))) - ...
  n * log(YVar - sum(alpha2(1:q) ./ (L + d2(1:q)))));
