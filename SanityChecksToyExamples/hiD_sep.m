%% Running Notes
%Developing set of tests for checking intuition about higher dimensional
%spaces and for developing potential null model checks for datalogger
%project as well as other high neuronal data, population level analyses.


%% Linear Separability test
%Curious if nearest neighbor failure (i.e. ratio between nearest and
%furthest converges to 1 as dim -> inf) has consequences for linear
%separability.

%Fairly satisified with the fact that it doesn't unless there is some
%difference in means for random gaussian data

dims = [20,50,100,500,1000];

cv_perform = zeros(length(dims),1);

num_cat = 2; %note: as expected takes much longer if you increase categories and make the problem harder.

dat_cat = cell(length(dims),num_cat);%just going to start with two categories

labels_cat = cell(size(dat_cat));

dat_tog = cell(length(dims),1);

samples = 25000; %may also want explore this as a function of sample size

means = [0 0.001 .4]; %note have to manually set

mdls = cell(1,length(dims));

%Essentially creating high-d spheres by making each category random pulls
%from a gaussian of dimensionality dims with identity correlation matrix.
%Keep variance tight, but move average to see separability

for i = 1:length(dims)
    
    count = 0;
    
    for cat = 1:num_cat
        
        dat_cat{i,cat} = normrnd(means(cat),1, samples, dims(i)); %keep observation x source structure matlab prefers
        
        labels_cat{i,cat} = count*ones(size(dat_cat{i,cat},1),1); %create vectors of labels
        
        count = count+1;
        
    end
    
    dat_tog{i} = [vertcat(labels_cat{i,:}) vertcat(dat_cat{i,:})];  %Put everything together with labels as leading column
    
    %default learner for these models is an SVM so don't have to set it in
    %the code.
    
    if num_cat<3
    
       mdls{i} = fitclinear(dat_tog{i}(:,2:end),dat_tog{i}(:,1),'KFold',5);
    
    else
     
        mdls{i} = fitcecoc(dat_tog{i}(:,2:end),dat_tog{i}(:,1),'KFold',5);
        
    end
    
    cv_perform(i) = 100 - kfoldLoss(mdls{i})*100  %Loss as percent misclassified, so subtract to get percent correct
    
    
end

figure()
plot(dims,cv_perform,'.k')
title(['Difference in means: ' num2str(means(1:num_cat))])
ylabel('Classification performance')
xlabel('Dimensionality of data')
ylim([(100/num_cat)-5 , 105])
xlim([min(dims)-5 max(dims)+5])




