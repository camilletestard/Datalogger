function [correspondingStates] = getCorrespondingStates(realR,estR)
M = [realR estR];
D = dist(M);
nstates = size(realR,2);
correspondingStates=zeros(1,nstates);
C = D(1:nstates,(nstates+1):end);
usedStates=[];
for i=1:nstates
    [~,inds] = sort(C(i,:));
    for j=1:length(inds)
        if (ismember(inds(j),usedStates))
            continue;
        else
            correspondingStates(i) = inds(j);
            usedStates = [usedStates inds(j)];
            break;
        end
    end
end
end

