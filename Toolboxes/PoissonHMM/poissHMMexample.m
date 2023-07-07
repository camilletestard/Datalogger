% HMMExample
nNeurons = 20;
nActualStates = 3;
nPredictedStates = 2;
maxFR = 100;
stateSeqLength = 10;
dt = .001;
maxIterations = 20;
stateSeq = ceil(rand(1,stateSeqLength)*nActualStates);
rates = ceil(rand(nNeurons,nActualStates)*maxFR);
[spikes] = genSpikes(rates,stateSeq,dt);

[bestPath,maxPathLogProb,PI,A,B,gamma] = poissHMM(spikes,nPredictedStates,dt,maxIterations);

tvec=dt:dt:stateSeqLength;
figure;
subplot(2,1,1)
plot(tvec,bestPath);
hold on;
plot(0:stateSeqLength-1,stateSeq);
legend({'Predicted','Actual'})

subplot(2,1,2);
plot(tvec,gamma'); legend()