function [spikes,stateSeqDurs,stateSeqVec] = genSpikes(rates,stateSeq,dt,varargin)
if (nargin > 3)
    stateDurMeans = varargin{1};
    stateDurStds = varargin{2};
else
    nStates = size(rates,2);
    stateDurMeans = ones(1,nStates);
    stateDurStds = zeros(1,nStates);
end
for i=1:length(stateSeq)
    stateSeqDurs(i) = max(.05,normrnd(stateDurMeans(stateSeq(i)),stateDurStds(stateSeq(i))));
end
spikes = zeros(size(rates,1),ceil(sum(stateSeqDurs)/dt));
stateSeqVec = [];
for i=1:size(spikes,1)
    startInd=1;
    for j=1:length(stateSeq)
        curDur = stateSeqDurs(j);
        curRate = rates(i,stateSeq(j));
        v = rand(1,ceil(curDur/dt)); v = v < curRate*dt;
        endInd = startInd + length(v) - 1;
        spikes(i,startInd:endInd) = v;
        startInd=startInd+length(v);
        if (i == 1)
            stateSeqVec = [stateSeqVec repelem(stateSeq(j),length(v))];
        end
    end
end
end

