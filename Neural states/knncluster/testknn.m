load iris.dat
X=iris(:,1:4);
[N,~]=size(X);
[cls,pnode]=knncluster(X,25);
% Color-coded scatter plot of first two attributes
scatter(X(:,1),X(:,2),[],cls,'filled')
xlabel('Sepal Length')
ylabel('Sepal Width')
% Generate and plot directed graph
G=digraph(1:N,pnode,'OmitSelfLoops');
figure(2)
% Plot in data space
plot(G,'XData',X(:,1),'YData',X(:,2),'MarkerSize',4)
xlabel('Sepal Length')
ylabel('Sepal Width')
figure(3)
% Plot cluster hierarchy
plot(G,'Layout','Layered')
% Compare clusters with actual classifications
M=max(cls);
conf=zeros(M,3);
for n=1:N
    i=cls(n);
    j=iris(n,5);
    conf(i,j)=conf(i,j)+1;
end
fprintf('Cluster     Setosa Versicolor  Virginica\n')
for i=1:M
    fprintf('%7d ',i)
    for j=1:3
        fprintf('%10d ',conf(i,j))
    end
    fprintf('\n')
end