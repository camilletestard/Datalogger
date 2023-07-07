% HMM vs. poissHMM
clear all;
nNeurons = 100;
nActualStates = 3;
nPredictedStates = 3;
maxFR = 15;
stateSeqLength = 10;
dt = .01;
maxIterations = 100;
stateDurMeans = normrnd(1,.1,1,nActualStates);
stateDurStds = .2*normrnd(1,.2,1,nActualStates);
stateSeq = ceil(rand(1,stateSeqLength)*nActualStates);
rates = ceil(rand(nNeurons,nActualStates)*maxFR);
[spikes,stateSeqDurs,stateSeqVec] = genPoissSpikes(rates,stateSeq,dt,stateDurMeans,stateDurStds);

%%
for i=1:size(spikes,2)
    spikingNeurons = find(spikes(:,i) > 0);
    if (length(spikingNeurons) > 1)
        ind = ceil(rand*length(spikingNeurons));
        v = zeros(nNeurons,1); v(ind) = 1;
        soloSpikes(:,i) = v;
        symbols(i) = spikingNeurons(ind)+1;
    elseif (length(spikingNeurons) == 1)
        v=zeros(nNeurons,1); v(spikingNeurons) = 1;
        soloSpikes(:,i) = v;
        symbols(i) = spikingNeurons+1;
    else
        soloSpikes(:,i) = zeros(nNeurons,1);
        symbols(i) = 1;
    end
end
%%
% poissHMM
[bestPath,maxPathLogProb,PI,A,B,gamma] = poissHMM(spikes,nPredictedStates,dt,maxIterations);

[correspondingStates] = getCorrespondingStates(rates,B);
for i=1:length(bestPath)
    bestPath(i) = find(bestPath(i) == correspondingStates);
end
tvec=linspace(dt,sum(stateSeqDurs),length(bestPath));
figure;
subplot(2,1,1)
plot(tvec,bestPath);
hold on;
plot(tvec,stateSeqVec);
legend({'Predicted','Actual'})

subplot(2,1,2);
plot(tvec,gamma'); legend()

% normal HMM
for i=1:nPredictedStates
    for j=1:nPredictedStates
        if (i == j)
            TRGUESS(i,j) = .99;
        else
            TRGUESS(i,j) = .01/(nPredictedStates-1);
        end
    end
end
EMITGUESS=rand(nPredictedStates,nNeurons+1);
[ESTTR,ESTEMIT] = hmmtrain(symbols,TRGUESS,EMITGUESS);
STATES = hmmviterbi(symbols,ESTTR,ESTEMIT);
PSTATES = hmmdecode(symbols,ESTTR,ESTEMIT);
figure;
subplot(2,1,1)
plot(tvec,STATES); hold on; plot(tvec,stateSeqVec)
subplot(2,1,2)
plot(PSTATES')