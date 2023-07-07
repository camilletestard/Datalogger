function [states,spikes] = poissHMMgenerate(seqLength,T,RATES,PI,dt)
nNeurons = size(RATES,1);
states = zeros(1,seqLength);
spikes = zeros(nNeurons,seqLength);
states(1) = randomDraw(PI,1);
spikes(:,1) = rand(nNeurons,1) < RATES(:,states(1))*dt;

for t=2:seqLength
    newState = randomDraw(T(states(t-1),:),1);
    newSpikes = rand(nNeurons,1) < RATES(:,newState)*dt;
    states(t) = newState;
    spikes(:,t) = newSpikes;
end
end