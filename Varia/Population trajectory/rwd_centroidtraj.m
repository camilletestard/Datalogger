%% Running Notes
%As first pass analysis for trajectories/manifold idea, finding population
%average (centroid) for each behavior and using this to predict moment to
%moment behavioral state.

%May need to reduce the number of components used so that MLE has an easier
%time for the multinomial regression (i.e. manual number of components)

%% Load in data and preprocess


is_ron = 1;

if is_ron >0
    load_folder = 'C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_data\';
else
    load_folder = '~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/';
end

session = 'Amos_2021-07-29';

%Pull need variables, push struct into separate variables in case need info
%later

Label_struct = load([load_folder session '/Labels_per_sec.mat']);
Data_struct = load([load_folder session '/Neural_data.mat']);

Labels = cell2mat(Label_struct.labels(:,3)); %Get numerical labels based on prioritization
Data = Data_struct.Unit_rasters; %Get neural data

Z_data = zscore(Data,[],2)';%Data';
%note: Z-score increases dim relative to not z-scoring as neurons that fire
%more can't dominate the overall signal.  Not an issue but may be worth
%keeping in mind
%% Sort and summarize data

%Quick histogram of behaviors

C = categorical(Labels,1:20,Label_struct.behav_categ);
figure()
histogram(C,Label_struct.behav_categ(1:end-1)) %removing rest since it dominates data


%Take advantage of numeric labels to get neural activity associated with
%each time a behavior was present.  %Not sure the sorting really helps but
%makes it easier to check things.


LD_holding = [Labels Z_data];%[Labels Data_group];% %This keeps matrix for all behaviors, now first column is labels

boi = [1:17]; %manually set behaviors of interest

index_use = LD_holding(:,1)==boi; %creates numel(boi) vectors of logical indecies for every time point
index_use = sum(index_use,2)>0; %has ones know where one of the behaviors of interest happened

LD_tog = LD_holding(index_use,:);

%undoing the sorting is not working as predicted...need to do a differen
%procedure it seems fix that here and then fix it in the main
%code...removing the sorting is now causing some kind of weird response and
%tanking performance.  Not really clear on what is happening but it seems like
%some offset between the data and behaviors is being introduced as
%performance drops to the same level as when didn't sort LD_tog but s_inds
%was still defined.  Going to try shuffling behaviors than sorting to make
%sure this isn't bs.

%Sanity check only

%LD_tog(:,1) = LD_tog(randperm(length(LD_tog)),1); %shuffle labels while sorting to see if the metric is false as there is now no relation between states and neurons

%Look like metric is inflated by this since now it just need to get the
%number of point right rather than the actual timing.  Removing all
%sorting.


%% Dim reduction, calculate centriod
% Only doing this to check linear embedding dimension and to aid visualization. 
%2021-12-3 copying Camille's method to select behaviors to analyze instead
%of using whole data set at once
%Just set boi to different values if you want different numbers of
%behaviors

man_num = 1; %user chooses number of components to use

[coeff, score,latent,~,explained] = pca(LD_tog(:,2:end), 'Centered', true);


centriods = cell(length(Label_struct.behav_categ),1);
%radii = zeros(length(Label_struct.behav_categ),1); %drawing spheres at some point?

var_sum = cumsum(explained);
        figure; plot(var_sum(1:end),'.k'); xlabel('PCs used'); ylabel('var explained');
        
if man_num
    
    num_component = input('Select number of PCs to use: ');
    
else
    
   num_component = find(cumsum(explained)>85,1,'first'); %use number of components that explains more than 85% of the variance



end

if num_component > 20 %want to use 10 or less, 20 is limit where nearest neighbor falls apart for high dim data

    disp(['number of components = ' num2str(num_component)])
    warning('Using more than 20 components in dim reduced space.  Changing from L2 to L1 measure')

end

%note, don't need to do reconstruction to see things in pca space, just use
%score variable

DR_data = score(:,1:num_component); 


figure; 
scatter3(DR_data(:,1), DR_data(:,2),DR_data(:,3),12); title('Data in PCA space');
figure; hold on
color = hsv(length(boi));



        
        for b = 1:length(boi)
            
            scatter3(DR_data(LD_tog(:,1)==boi(b),1), DR_data(LD_tog(:,1)==boi(b),2),DR_data(LD_tog(:,1)==boi(b),3),...
                12,color(b,:),'filled');
            
            centriods{boi(b)} = mean(DR_data(LD_tog(:,1)==boi(b),:),1);
            
           
            
        end
       
%         figure;hold on
%         for b = 1:length(boi)
%             
%             scatter3(centriods{boi(b)}(1),centriods{boi(b)}(2),centriods{boi(b)}(3),12,color(b,:),'filled')
%             
%         end


%end


%% Predict state from closest centriod
%removing this for now and switching to the multinomial regression
% %Try training and testing
% %Try moving centriod closer
% 
dis_mat = nan(length(DR_data),length(boi)); %initialize distance matrix for distance between 

%As a check add noise which should tank performance
%Check passed...still surprised at how well this works.


use_noise = 0; %Note: check average distance and make sure to change noise parameters accordingly

%Trying other noise implementation too where use it on DR_data directly;
%DR_data=DR_data_hold;

noise_dr = 0;

DR_data_hold = DR_data;
DR_data = DR_data+noise_dr*normrnd(3,1,size(DR_data));

for b=1:length(boi) 
    
    if num_component > 20 %use L1 norm
    
    dis_mat(:,b) = vecnorm(DR_data'-centriods{boi(b)}',1)+use_noise*normrnd(200,100,1,length(dis_mat)); %works more easily if everything is column vector
    
    else %use L2 norm
        
    dis_mat(:,b) = vecnorm(DR_data'-centriods{boi(b)}',2)+use_noise*normrnd(200,100,1,length(dis_mat));
        
    end
    
    

end

[~,cs]=  min(dis_mat,[],2); %Find which centriod the activity was closer to

preds = boi(cs); %Predict behavior that was that centriod



%plot step function
figure
plot(LD_tog(:,1), 'LineWidth',2)
hold on
plot(preds, 'LineWidth',2)
yticks([0:18]);
ticklabs = Label_struct.behav_categ;
yticklabels({'',ticklabs{boi},''})
ylim([0 max(boi)+2])
legend('Real','Predicted State', 'FontSize',16)
ylabel('Behavior', 'FontSize',18)
xlabel('Time', 'FontSize',18)
%xlabel('index')
ax = gca;
ax.FontSize = 14; 
title('Predicted behavioral state based on neural data vs. Real behavioral state')


per_cor_centroid = sum(LD_tog(:,1)==preds')/length(preds)*100


%%  Trying DBSCAN for unsupervised clustering
%2021-12-06 update need to mess with paramters for this to work
%This requires a lot of massaging to work so I'm not sure it's worth using
%Update feel like this nees to be used on heavily reduced data (i.e. with
%UMAP or something)
% [idxs, corepts] = dbscan(DR_data,250,5);
% unique(idxs)
% gscatter(DR_data(:,1),DR_data(:,2),idx)
% 


%% Trying multinomial regression instead of centriod analysis
%Note: should implement a k-fold cross validation loop for this for greater
%confidence in the results
%Further note: keep an eye for when MLE does and doesn't converge

behavs = categorical(LD_tog(:,1),boi,{Label_struct.behav_categ{boi}}); %create categorical array for use with mnrfit

%Matlab has a built in k-fold function but it isn't clear if it interfaces
%with mnrfit, so just going to code my own.  Going with 10 portions instead
%of 10 grabs

folds = 10;

sample_size = floor(length(LD_tog(:,1))/folds);

sam_if = ones(folds,1)*sample_size; %samples in each fold


if mod(length(LD_tog(:,1)),folds)>0 %check if not divisible by 10
    disp('Samples do not divide evenly accross folds')
    %if it does not divide evenly, add remainder randomly to one of the
    %folds
    ind = randsample(folds,1);
    sam_if(ind) = sam_if(ind) + mod(length(LD_tog(:,1)),folds);
    
end


%inds_groups = [1 cumsum(sam_if)'];  %Doing the other implementation of shuffling all of the inds and grabbing the number of inds set above

shuffledind = randperm(length(behavs));

ind_if = cell(folds,1);

ind_if{1} = shuffledind(1:sam_if(1)); %set the first group out of the loop

groups = cumsum(sam_if);


if sum(sam_if)~=length(LD_tog(:,1))
    
    error('Missing inds from shuffle')
    
end


for i=2:folds
    
    ind_if{i} = shuffledind(groups(i-1)+1:groups(i));
    
    
end

%2021-12-07 update: non-trival to increase the number of iterations mnrfit
%does...so will need to work around this or take the time to make something more custom...or switch to python.

cv_per_cor = nan(folds,1);
cv_preds = cell(folds,1);
cv_behavs = cell(folds,1);

clear i %just incase messed up index somewhere

for k = 1:folds
   
  disp(['iteration: ' num2str(k)])
    
 states = behavs(ind_if{k},1);
 
 cv_behavs{k} = states;
 
 foldsidx = 1:folds;
 
 trainingidx = horzcat(ind_if{foldsidx~=k})'; %train on all indecies that aren't in the current fold

[Betas,~,stats] = mnrfit(DR_data(trainingidx,:),behavs(trainingidx,1)); %note always get B is dim predictors+1 for the intercept term x length(boi)-1  as one behavior is selected as reference

[pihat, ~,~] = mnrval(Betas,DR_data(ind_if{k},:),stats); %give probabilities for each behavior

%take max of each predicted probability as the predicted behavioral state
%like above


[~,cs]= max(pihat,[],2); %Find which centriod the activity was closer to

preds = boi(cs); %Predict behavior that was that centriod

cv_preds{k} = preds;

%preds = categorical(preds,boi,{Label_struct.behav_categ{boi}}); %use same categorical trick above so can do string compare



% figure() %Have to wait for end to concatenate everything and unscramble
% plot(LD_tog(s_inds,1))
% hold on
% plot(preds(s_inds))
% title('Multinomial Regression Predictions')



per_cor_mnr = sum(LD_tog(ind_if{k},1)==preds')/length(preds)*100; %strcmp(behavs(ind_if{k}),preds)

cv_per_cor(k) = per_cor_mnr;



end

%% Plot CV results

%% Run and plot without cv just for convenience for now

behavs = categorical(LD_tog(:,1),boi,{Label_struct.behav_categ{boi}}); %create categorical array for use with mnrfit

[Betas,dev,stats] = mnrfit(DR_data,behavs); %note always get B is dim predictors+1 for the intercept term x length(boi)-1  as one behavior is selected as reference

[pihat, dlow,dhi] = mnrval(Betas,DR_data,stats); %give probabilities for each behavior

%take max of each predicted probability as the predicted behavioral state
%like above

[~,cs]= max(pihat,[],2);

preds = boi(cs); %Predict behavior that was that centriod



per_cor_nocv = sum(LD_tog(:,1)==preds')/length(preds)*100

figure() %Have to wait for end to concatenate everything and unscramble
plot(LD_tog(:,1))
hold on
plot(preds)
title('Multinomial Regression Predictions')

%% posterity code

