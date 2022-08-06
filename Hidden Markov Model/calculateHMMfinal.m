function HMM = calculateHMM(spikes, dt, averageStateDuration, verbose, ...
    version, randInit, nStatesFinal, nStatesVector, nRepetitions, ...
    timeSmoothFactor, normalization, timeSmoothFactor2)
%calculateHMM Hidden Markov Model analysis
% From discretized spiking data (nNeurons x nTimeSteps 2D matrix),
% extract symbol = which neuron is spiking at each time step (random choice
% between multiple neurons if > 1 neuron spike simultaneously; if no neuron
% spikes, null symbol). Then, Hidden Markov Model analysis on this 
% 1 x nTimeSteps 1D vector of symbols (= which neurons spike). 
% Initial guess for HMM transition & emission matrix, then optimization.
% Input:
%     - spikes: nNeurons x nTimeSteps 2D binary matrix, containing
%         discretized spike data (did neuron N spike during time bin T, yes
%         or no?)
%     - dt: num scalar, width of time step / time bin
%     - averageStateDuration: num scalar, initial estimation of state period
%         duration (in milliseconds)
%     - verbose: bool scalar, print program working info or not
%     - version: char, what version of data transformation to use.
%         'single spikes' or 'k-means'. Only use 'single spikes' = procedure
%         described above (symbol = identity of 1 neuron spiking at each time
%         step). Inspired by Seidemann / O'Halloran & Cronin.
%         'k-means' = attempt to define symbol through N dimensional
%         network frequency space, where symbol = closest centroid in
%         to vector of network frequency. Inspired by Abeles / Gat / Tishby
%         algorithm. Doesn't seem to work.
%     - randInit: bool scalar, random initialization of transition & emission
%         matrix or not (see detail in code).
%     - nStatesFinal: num scalar, final number of states to choose
%         (irrespective of selection criteria, see code below for more info).
%     - nStatesVector: num vector, multiple HMMs each with a different
%         number of states as defined by this vector.
%     - nRepetitions: num scalar, number of times to repeat HMMs for each
%         number of state.
%     - timeSmoothFactor: num scalar. Only useful for k-means algorithm,
%         specifies milliseconds over which to smooth spikes for frequency
%         estimation. Might not be fully compatible with latest code version.
%     - normalization: char. Only useful for k-means algorithm, specifies
%         frequency normalization method (subtract, divide) between frequency
%         estimated with timeSmoothFactor and timeSmoothFactor2. Might not be
%         fully compatible with latest code version.
%     - timeSmoothFactor2: num scalar. Only useful for k-means algorithm,
%         specifies milliseconds over which to smooth spikes for 2nd longer
%         frequency estimation. Might not be fully compatible with latest
%         code version.
% Output:
%     - HMM: struct scalar, contains many HMM information (the observables,
%         the states, the probability of states, the transition / emission 
%         matrices, etc.)

%% Handle arguments
assert(islogical(spikes) && ~isscalar(spikes), ...
    '1st argument (spikes) must be a logical vector or matrix');

assert(isnumeric(dt) && ~isnan(dt) && ~isinf(dt) && isscalar(dt), ...
    '2nd argument (dt) must be a valid numeric scalar.');

if ~exist('averageStateDuration', 'var') || isempty(averageStateDuration)
    averageStateDuration = 300;
end
assert(isnumeric(averageStateDuration) && ~isnan(averageStateDuration) && ...
    ~isinf(averageStateDuration) && isscalar(averageStateDuration), ...
    '3rd argument (averageStateDuration) must be a valid numeric scalar.');

if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end
assert(islogical(verbose) || ...
    (isnumeric(verbose) && ~isnan(verbose) && ...
    ~isinf(verbose) && isscalar(verbose) && ...
    (verbose == 0 || verbose == 1)), ...
    '4th argument (verbose) must be logical or numeric (1 or 0).');

if ~exist('version', 'var') || isempty(version)
    version = 'single spikes';
end
possibilities = {'single spikes'}; % {'single spikes', 'k-means'}
assert(ischar(version) && any(strcmp(version, possibilities)), ...
    sprintf(['5th argument (version) must be one of the following ' ...
    'chars: ' repmat('''%s'', ', 1, length(possibilities)) '.']));

if ~exist('randInit', 'var') || isempty(randInit)
    randInit = false;
end
assert(islogical(randInit) || ...
    (isnumeric(randInit) && ~isnan(randInit) && ...
    ~isinf(randInit) && isscalar(randInit) && ...
    (randInit == 0 || randInit == 1)), ...
    '6th argument (randInit) must be logical or numeric (1 or 0).');

if ~exist('nStatesFinal', 'var') || isempty(nStatesFinal)
    nStatesFinal = 4;
end
assert(isnumeric(nStatesFinal) && ~isnan(nStatesFinal) && ...
    ~isinf(nStatesFinal) && isscalar(nStatesFinal), ...
    '7th argument (nStatesFinal) must be a valid numeric scalar.');

if ~exist('nStatesVector', 'var') || isempty(nStatesVector)
    nStatesVector = nStatesFinal;
end
assert(isnumeric(nStatesVector) && all(~isnan(nStatesVector)) && ...
    all(~isinf(nStatesVector)) && isvector(nStatesVector), ...
    '8th argument (nStatesVector) must be a valid numeric vector or scalar.');

if ~exist('nRepetitions', 'var') || isempty(nRepetitions)
    nRepetitions = 1;
end
assert(isnumeric(nRepetitions) && ~isnan(nRepetitions) && ...
    ~isinf(nRepetitions) && isscalar(nRepetitions), ...
    '9th argument (nRepetitions) must be a valid numeric scalar.');

if ~exist('timeSmoothFactor', 'var') || isempty(timeSmoothFactor)
    timeSmoothFactor = 500;
end
assert(isnumeric(timeSmoothFactor) && ~isnan(timeSmoothFactor) && ...
    ~isinf(timeSmoothFactor) && isscalar(timeSmoothFactor), ...
    '10th argument (timeSmoothFactor) must be a valid numeric scalar.');

if ~exist('normalization', 'var') || isempty(normalization)
    normalization = 'none';
end
possibilities = {'none', 'subtract'};
assert(ischar(normalization) && any(strcmp(normalization, possibilities)), ...
    sprintf(['11th argument (normalization) must be one of the following ' ...
    'chars: ' repmat('''%s'', ', 1, length(possibilities)) '.']));

if ~exist('timeSmoothFactor2', 'var') || isempty(timeSmoothFactor2)
    timeSmoothFactor2 = 60000;
end
assert(isnumeric(timeSmoothFactor2) && ~isnan(timeSmoothFactor2) && ...
    ~isinf(timeSmoothFactor2) && isscalar(timeSmoothFactor2), ...
    '12th argument (timeSmoothFactor2) must be a valid numeric scalar.');

%% Prepare observations temporal vector
if strcmp(version, 'single spikes')
    %% Average number of spikes per time bin
    nSpikes = 0:4;
    nSpikesPerBin = NaN(1, length(nSpikes));

    for s = 1:length(nSpikes)
        n = nSpikes(s);
        nSpikesPerBin(s) = mean(sum(spikes, 1) >= n & ...
                             sum(spikes, 1) < (n+1));
    end
    pBinWithoutSpikes = nSpikesPerBin(nSpikes == 0);

    % --- Debugging plot ---
    % figure; hold on;
    % plot(nSpikes, nSpikesPerBin);
    % plot(nSpikes, cumsum(nSpikesPerBin));
    
    clear nSpikes nSpikesPerBin s n nSpikesPerBin

    %% Convert multi-variate spike trains to univariate spike idx
    nTimeBins = size(spikes, 2);
    nNeurons = size(spikes, 1);
    nOutputs = nNeurons+1;
    HMM.obs = ones(1, nTimeBins)*nOutputs;
    for kt = 1:nTimeBins
        if any(spikes(:, kt))
            spikeIdx = find(spikes(:, kt));
            if length(spikeIdx) == 1
                HMM.obs(kt) = spikeIdx;
            else
                HMM.obs(kt) = datasample(spikeIdx, 1);
            end
        end
    end
    
    clear nTimeBins nNeurons kt spikeIdx
elseif strcmp(version, 'k-means')
    ti = round(timeSmoothFactor/dt);
    
    % Bins
%     spikes = spikes(:, 1:(round(end/ti)*ti));
%     rates = permute(sum(reshape(spikes, size(spikes, 1), ti, []), 2), ...
%         [1 3 2]);
    
    %% Smooth
    window = ones(1, ti)/ti;
    rates = conv2NoPadding(spikes, window);
    if strcmp(normalization, 'subtract')
        ti2 = round(timeSmoothFactor2/dt);
        window2 = ones(1, ti2)/ti2;
        slow = conv2NoPadding(spikes, window2);
        rates = rates - slow + mean(slow, 2);
        % Other possibilities: high-pass filter, wavelet transform,
        % fourier transform?
        % - Wavelet transform gives same result as conv_fast - conv_slow,
        % but has problems on edges.
        % - High-pass filter is good sometimes, but crazy in other cases.
        % Haven't tried fourier transform, but seems to have other
        % problems (notably, edges).
        % Thus, conv_fast - conv_slow is best.
    end
    
    %% Verification code
%     HMM = struct();
%     nNeurons = size(spikes, 1);
%     rates2 = permute(sum(reshape(spikes(:, 1:floor(end/ti)*ti), ...
%         size(spikes, 1), ti, []), 2), [1 3 2]);
%     
%     t = 0:timeSmoothFactor:1e3;
%     f = @(tau) exp(-t/tau);
%     options = optimoptions('fmincon', 'Display', 'off');
%     taus = NaN(nNeurons, 1);
%     for n = 1:nNeurons
%         temp = autocorr(rates2(n, :), 'NumLags', round(1e3/timeSmoothFactor));
%         error = @(tau) sqrt(sum((f(tau) - temp).^2));
%         taus(n) = fmincon(error, 100, [], [], [], [], [], [], [], ...
%             options);
%     end
%     data = [std(rates, [], 2) ...
%         mean(rates, 2) ...
%         std(rates, [], 2) ./ mean(rates, 2) ...
%         mean(abs(diff(rates, [], 2)), 2) ...
%         taus];
%     data = median(data, 1);
%     fprintf('%g ', data);
%     fprintf('\n');
%     return;

    %% HMM
    nOutputs = size(spikes, 1)*10;
    [HMM.obs, clusters] = kmeans(rates', nOutputs, 'Replicates', 20);
    HMM.obs = HMM.obs';
    
    clear window rates spikes clusters
end

% --- Debugging plot ---
% newSpikes = spikeIdxs;
% newSpikes(newSpikes == nOutputs) = NaN;
% [nAP, tAP] = find(spikes);
% figure; hold on;
% scatter(tAP, nAP, 6, 'r', 'filled');
% scatter(1:nT, newSpikes, 6, 'k', 'filled');

%% Estimate transition & emission matrices
averageTransitionProb = 1/(averageStateDuration/dt);

template = cell(length(nStatesVector), nRepetitions);
transitionMatriceCell = template;
emissionMatriceCell = template;
statesCell = template;
pStatesCell = template;
logProbSeqCell = NaN(length(nStatesVector), nRepetitions);
logLikCell = template; % My extension

clear template

for s = 1:length(nStatesVector)
    nStates = nStatesVector(s);
    fprintf('Number of states: %d. ', nStates);
    for r = 1:nRepetitions
        fprintf('Repetition = %d, ', r);
        %% Transition matrix initialization: Si x Sj (Si --> Sj)
        % Bias transition matrix towards staying in same state (diagonal),
        % with average duration = averageStateDuration
        diagonal = eye(nStates);
        if randInit
            randMat = 2 .* rand(nStates, nStates); % uniform [0 2]
        else
            randMat = ones(nStates, nStates); % all ones
        end
        transitionMatrixInit = ...
            diagonal .* (1-averageTransitionProb) + ... Diagonal non-random
            (1 - diagonal) .* ...
            (averageTransitionProb / (nStates-1) .* randMat);
        if ~verLessThan('matlab', '9.0')
             transitionMatrixInit = transitionMatrixInit ./ ...
                 sum(transitionMatrixInit, 2);
        else
             transitionMatrixInit = transitionMatrixInit ./ ...
                 repmat(sum(transitionMatrixInit, 2), 1, nStates);
        end

        %% Emission matrix initialization: S x O
        % Biais emission matrix with high probability for no spike
        emissionMatrixInit = NaN(nStates, nOutputs);
        
        if randInit
            randMat = 2 .* rand(nStates, nOutputs); % uniform [0 2]
        else
            randMat = ones(nStates, nOutputs); % all ones
        end
        
        if strcmp(version, 'single spikes')
            emissionMatrixInit(:, end) = pBinWithoutSpikes;
            randMat(:, end) = 1; % Spike absence non-random
            emissionMatrixInit(:, 1:end-1) = ...
                (1-pBinWithoutSpikes) / (nOutputs-1);
        elseif strcmp(version, 'k-means')
            emissionMatrixInit(:) = 1/nOutputs;
        end
        emissionMatrixInit = emissionMatrixInit .* randMat;
        if ~verLessThan('matlab', '9.0')
            emissionMatrixInit = emissionMatrixInit ./ ...
                sum(emissionMatrixInit, 2);
        else
            emissionMatrixInit = emissionMatrixInit ./ ...
                repmat(sum(emissionMatrixInit, 2), 1, nOutputs);
        end

        %% Estimate transition & emission matrices
        [transitionMatrix, emissionMatrix, loglik] = ...
            hmmtrain(HMM.obs, transitionMatrixInit, emissionMatrixInit, ...
            'Verbose', verbose); % I added loglik
        states = hmmviterbi(HMM.obs, transitionMatrix, emissionMatrix);
        [pStates, logProbSeq] = ...
            hmmdecode(HMM.obs, transitionMatrix, emissionMatrix);

        %% Compare values for different number of states
        transitionMatriceCell{s, r} = transitionMatrix;
        emissionMatriceCell{s, r} = emissionMatrix;
        statesCell{s, r} = states;
        pStatesCell{s, r} = pStates;
        logProbSeqCell(s, r) = logProbSeq;
        logLikCell{s, r} = loglik; % my addition
    end
    fprintf('\n');
end

clear s nStates r diagonal randMat transitionMatrixInit emissionMatrixInit
clear transitionMatrix emissionMatrix states pStates logProbSeq

%% Choose best HMM for each number of states, then best number of states
% [AIC, BIC, p] = compareHMMNumberOfStates(nStatesVector, HMM.obs, ...
%     logProbSeqCell, transitionMatriceCell, emissionMatriceCell, ...
%     pStatesCell, statesCell, nOutputs, size(spikes, 2));

% Get best HMM (highest log-likelihood) for each number of states?
% Or just choose number of states?
% [~, repMax] = max(logProbSeqCell, [], 2);
% maxIdx = sub2ind(size(logProbSeqCell), 1:length(nStatesVector), repMax')';
maxIdx = find(nStatesFinal == nStatesVector);

transitionMatriceCell = transitionMatriceCell(maxIdx);
emissionMatriceCell = emissionMatriceCell(maxIdx);
statesCell = statesCell(maxIdx);
pStatesCell = pStatesCell(maxIdx);
logProbSeqCell = logProbSeqCell(maxIdx);
logLikCell = logLikCell(maxIdx);

% Choose best number of states
stateIdx = nStatesVector == nStatesFinal;
HMM.transitionMatrix = transitionMatriceCell{stateIdx};
HMM.emissionMatrix = emissionMatriceCell{stateIdx};
HMM.states = statesCell{stateIdx};
HMM.pStates = pStatesCell{stateIdx};
HMM.logProbSeq = logProbSeqCell(stateIdx);
HMM.logLik = logLikCell{stateIdx};

clear repMax logProbSeqCell maxIdx nStatesVector transitionMatriceCell
clear emissionMatriceCell statesCell pStatesCell logProbSeqCell stateIdx

%% Add null state
nullStateThreshold = 0.8;
existingStates = HMM.pStates > nullStateThreshold;
[~, HMM.newStates] = max(existingStates, [], 1);
nullState = ~any(existingStates, 1);
HMM.newStates(nullState) = 0;

clear nullStateThreshold existingStates nullState

%% Concatenate states when null state in between
paddedStates = [NaN HMM.newStates NaN];
stateDiffs = diff(paddedStates);
newStateChangesIdx = find(stateDiffs);
newStateChangesIdx(end) = newStateChangesIdx(end)-1;
% stateDiffsAtChanges = stateDiffs(newStateChangesIdx);
newStateAtChanges = HMM.newStates(newStateChangesIdx);

for s = 1:nStatesFinal
    for t = 1:length(newStateAtChanges)-2
        if all(newStateAtChanges(t:t+2) == [s 0 s])
            HMM.newStates(newStateChangesIdx(t+1):newStateChangesIdx(t+2)-1) = s;
        end
    end
end

clear paddedStates stateDiffs s t newStateChangesIdx newStateAtChanges

%% Calculate state periods
% Prepare data
paddedStates = [NaN HMM.newStates NaN];
stateDiffs = diff(paddedStates);
HMM.newStateChangesIdx = find(stateDiffs);
HMM.newStateChangesIdx(end) = HMM.newStateChangesIdx(end)-1;
HMM.newStateAtChanges = HMM.newStates(HMM.newStateChangesIdx);

% Calculate all durations
allNewDurations = diff(HMM.newStateChangesIdx)*dt;

clear paddedStates stateDiffs

%% Calculate durations of states
% Calculate durations of all states
HMM.newDurations = allNewDurations(HMM.newStateAtChanges(1:end-1) > 0);
HMM.newAverageDuration = mean(HMM.newDurations);
HMM.newMedianDuration = median(HMM.newDurations);

% Calculate durations for each state
HMM.newDurationsInState = cell(1, nStatesFinal);
HMM.newAverageDurationInState = cell(1, nStatesFinal);
HMM.newMedianDurationInState = cell(1, nStatesFinal);
for s = 1:nStatesFinal
    stateNewDurations = allNewDurations(HMM.newStateAtChanges(1:end-1) == s);
    HMM.newDurationsInState{s} = stateNewDurations;
    HMM.newAverageDurationInState{s} = mean(HMM.newDurationsInState{s});
    HMM.newMedianDurationInState{s} = median(HMM.newDurationsInState{s});
end

% Calculate null durations
HMM.newDurationsInNull = allNewDurations(HMM.newStateAtChanges(1:end-1) == 0);
HMM.newAverageDurationInNull = mean(HMM.newDurationsInNull);
HMM.newMedianDurationInNull = median(HMM.newDurationsInNull);

clear allNewDurations s stateNewDurations

%% Calculate time spent in states
% Calculate time spent of all states
data = HMM.newDurations;
edges = (min(data)-dt/2) : dt : (max(data)+dt/2);
temp = reverseHist(data, edges);
HMM.newTimeSpent = temp{1};
HMM.newAverageTimeSpent = mean(HMM.newTimeSpent);
HMM.newMedianTimeSpent = median(HMM.newTimeSpent);

% Calculate time spent for each states
HMM.newTimeSpentInState = cell(1, nStatesFinal);
HMM.newAverageTimeSpentInState = cell(1, nStatesFinal);
HMM.newMedianTimeSpentInState = cell(1, nStatesFinal);
for s = 1:nStatesFinal
    data = HMM.newDurationsInState{s};
    edges = (min(data)-dt/2) : dt : (max(data)+dt/2);
    temp = reverseHist(data, edges);
    HMM.newTimeSpentInState{s} = temp{1};
    HMM.newAverageTimeSpentInState{s} = mean(HMM.newTimeSpentInState{s});
    HMM.newMedianTimeSpentInState{s} = median(HMM.newTimeSpentInState{s});
end

% Calculate time spent of in null state
data = HMM.newDurationsInNull;
edges = (min(data)-dt/2) : dt : (max(data)+dt/2);
temp = reverseHist(data, edges);
HMM.newTimeSpentInNull = temp{1};
HMM.newAverageTimeSpentInNull = mean(HMM.newTimeSpentInNull);
HMM.newMedianTimeSpentInNull = median(HMM.newTimeSpentInNull);

clear data edges temp s

%% Last calculations
HMM.averageStateDurations = dt ./ (1 - diag(HMM.transitionMatrix));
% if ~verLessThan('matlab', '9.0')
%     HMM.neuronStateProb = emissionMatrix ./ sum(emissionMatrix, 1);
% else
%     HMM.neuronStateProb = emissionMatrix ./ ...
%         repmat(sum(emissionMatrix, 1), nStatesFinal, 1);
% end
% HMM.neuronStateProb = HMM.neuronStateProb(:, 1:end-1);
% [~, HMM.neuronsInState] = max(HMM.neuronStateProb, [], 1);

%% Storing the Viterbi-based results

% Computing the original Markovian states
fprintf('\tcomputing state chnages...\n');
paddedStates = [NaN HMM.states NaN];
stateDiffs = diff(paddedStates);
HMM.stateChangesIdx = find(stateDiffs);
HMM.stateChangesIdx(end) = HMM.stateChangesIdx(end)-1;
HMM.stateAtChanges = HMM.states(HMM.stateChangesIdx);

% Calculate all durations
allDurations = diff(HMM.stateChangesIdx)*dt;
clear paddedStates stateDiffs

% Calculate durations of states
fprintf('\tcomputing durations...\n');
% Calculate durations of all states
HMM.durations = allDurations(HMM.stateAtChanges(1:end-1) > 0);
HMM.averageDuration = mean(HMM.durations);
HMM.medianDuration = median(HMM.durations);

% Calculate durations for each state
HMM.durationsInState = cell(1, nStatesFinal);
HMM.averageDurationInState = cell(1, nStatesFinal);
HMM.medianDurationInState = cell(1, nStatesFinal);
for s = 1:nStatesFinal
    stateDurations = allDurations(HMM.stateAtChanges(1:end-1) == s);
    HMM.durationsInState{s} = stateDurations;
    HMM.averageDurationInState{s} = mean(HMM.durationsInState{s});
    HMM.medianDurationInState{s} = median(HMM.durationsInState{s});
end

clear allDurations s stateDurations

% Calculate time spent in state
fprintf('\tcomputing time spent...\n');
% Calculate time spent of all states
data = HMM.durations;
edges = (min(data)-dt/2) : dt : (max(data)+dt/2);
temp = reverseHist(data, edges);
HMM.timeSpent = temp{1};
HMM.averageTimeSpent = mean(HMM.timeSpent);
HMM.medianTimeSpent = median(HMM.timeSpent);

% Calculate time spent for each states
HMM.timeSpentInState = cell(1, nStatesFinal);
HMM.averageTimeSpentInState = cell(1, nStatesFinal);
HMM.medianTimeSpentInState = cell(1, nStatesFinal);
for s = 1:nStatesFinal
    data = HMM.durationsInState{s};
    edges = (min(data)-dt/2) : dt : (max(data)+dt/2);
    temp = reverseHist(data, edges);
    HMM.timeSpentInState{s} = temp{1};
    HMM.averageTimeSpentInState{s} = mean(HMM.timeSpentInState{s});
    HMM.medianTimeSpentInState{s} = median(HMM.timeSpentInState{s});
end

clear data edges temp s

%% Keep data
HMM.version = version;
HMM.dt = dt;
HMM.nStates = nStatesFinal;

end

% To re-calculate everything else quickly, most important variables:
%     HMM.obs, HMM.transitionMatrix, HMM.emissionMatrix
%     HMM.stateChangesIdx, HMM.stateAtChanges

% Automatically detect number of states
    % * Repetitions?
        % * Without random initialization, multiple repetitions give exactly 
        % same transition / emission matrix / probability & chosen state
        % sequence for N states.
        % * With random initialization, probably different.
    % * Increasing number of states?
        % * Increases log-prob, partial structure kept.
        % * transition times ~< 300ms.
        % * To choose, BIC/AIC, cross-validation, (duration of states), 
        % Bayesian estimation of the parameters. No fit-all procedure.
% Smoothed probabilities vs real probabilities with threshold
% Test multiple random initializations