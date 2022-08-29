function [cls, pnode] = knncluster(X, K)
%{
Variation of the clustering algorithm based on directed graphs described
in "A graph-theoretic approach to nonparametric cluster analysis"
by Koontz et al. in IEEE Trans. on Computers, September 1976. Instead
of using a threshold distance to determine if two objects are neighbors,
we determine up to a fixed number of nearest neighbors

Inputs:
    X: Matrix of data to be clustered. Each row corresponds
       to an object in the data set and each column is an
       attribute of the object.
    K: The maximum number of nearest neighbors
Outputs:
    cls: A vector of the cluster indices of the objects
    pnode: A vector of the parent object of each object
%}
[N,~] = size(X);
numNeighbors = zeros(N,1);
maxDistance = zeros(N,1);
neighborID = zeros(N,K);
neighborDist = zeros(N,K);
for i=2:N
    for j=1:i-1 % for all i, j such that j < i
        dij = sqrt(sum((X(i,:)-X(j,:)).^2));
        % See if j can be a neighbor of i
        if numNeighbors(i) < K % neighborhood not yet full
            n = numNeighbors(i) + 1;
            neighborID(i,n) = j;
            neighborDist(i,n) = dij;
            numNeighbors(i) = n;
        elseif dij < maxDistance(i) % is this object nearer?
            % replace farthest neighbor
            [~,nmax] = max(neighborDist(i,:));
            neighborID(i,nmax) = j;
            neighborDist(i,nmax) = dij;
        end
        % update maximum distance
        maxDistance(i) = max(neighborDist(i,:));
        % Repeat, swapping i and j
        if numNeighbors(j) < K
            n = numNeighbors(j) + 1;
            neighborID(j,n) = i;
            neighborDist(j,n) = dij;
            numNeighbors(j) = n;
        elseif dij < maxDistance(j)
            [~,nmax] = max(neighborDist(j,:));
            neighborID(j,nmax) = i;
            neighborDist(j,nmax) = dij;
        end
        maxDistance(j) = max(neighborDist(j,:));
    end
end
% Now find parents and orphans based on minimum maximum distance
M=0;
pnode = 1:N;
cls=zeros(N,1);
for n=1:N
    dmin = maxDistance(n);
    for m = 1:numNeighbors(n)
        k = neighborID(n,m);
        if maxDistance(k) < dmin
            dmin = maxDistance(k);
            pnode(n) = k;
        end
    end
    if pnode(n)==n % orphan?
        % Create new cluster
        M=M+1;
        cls(n)=M;
    end
end
%
% Cluster the rest recursively
%
for n=1:N
    pn=n;
    while cls(pn)==0
        pn=pnode(pn);
    end
    cls(n)=cls(pn);
end