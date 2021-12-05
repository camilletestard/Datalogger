%% Running Notes
%As first pass analysis for trajectories/manifold idea, finding population
%average (centroid) for each behavior and using this to predict moment to
%moment behavioral state.

%% Load in data and preprocess


is_ron = 0;

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

Z_data = zscore(Data,[],2)';
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

[~, s_inds] = sort(LD_tog(:,1));

LD_tog = LD_tog(s_inds,:);

%% Dim reduction, calculate centriod
% Only doing this to check linear embedding dimension and to aid visualization. 
%2021-12-3 copying Camille's method to select behaviors to analyze instead
%of using whole data set at once
%Just set boi to different values if you want different numbers of
%behaviors

% use_boi = 1;
% 
% if use_boi <1
% 
% [coeff, score,~,~,explained] = pca(Z_data, 'Centered', false); %center data in algorithm
% 
% num_component = find(cumsum(explained)>85,1,'first'); %use number of components that explains more than 85% of the variance
% 
%     if num_component > 20 %want to use 10 or less, 20 is limit where nearest neighbor falls apart for high dim data
% 
%         disp(['number of components = ' num2str(num_component)])
%         warning('Using more than 20 components in dim reduced space.  Changing from L2 to L1 measure')
% 
%     end
% 
% 
% DR_data = score(:,1:num_component); 
% 
% figure; plot(cumsum(explained)); xlabel('PCs used'); ylabel('var explained');
%         figure; 
%         scatter3(DR_data(:,1), DR_data(:,2),DR_data(:,3),12); title('Data in PCA space');
%         figure; hold on
%         color = hsv(length(1:19));
%         for b = 1:19 %plot everything but rest.
%             scatter3(DR_data(LD_holding(:,1)==b,1), DR_data(LD_holding(:,1)==b,2),DR_data(LD_holding(:,1)==b,3),...
%                 12,color(b,:),'filled');
%         end
% 
% 
% else %only consider behaviors of interest
    
[coeff, score,latent,~,explained] = pca(LD_tog(:,2:end), 'Centered', false);

num_component = find(cumsum(explained)>85,1,'first'); %use number of components that explains more than 85% of the variance

    if num_component > 20 %want to use 10 or less, 20 is limit where nearest neighbor falls apart for high dim data

        disp(['number of components = ' num2str(num_component)])
        warning('Using more than 20 components in dim reduced space.  Changing from L2 to L1 measure')

    end


%note, don't need to do reconstruction to see things in pca space, just use
%score variable

DR_data = score(:,1:num_component); 



centriods = cell(length(Label_struct.behav_categ),1);
%radii = zeros(length(Label_struct.behav_categ),1); %drawing spheres at some point?

var_sum = cumsum(explained);
        figure; plot(1:10:length(var_sum),var_sum(1:10:end),'.k'); xlabel('PCs used'); ylabel('var explained');
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

%Try training and testing
%Try moving centriod closer

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

[rs,cs]= find(dis_mat == min(dis_mat,[],2)); %Find which centriod the activity was closer to

preds = boi(cs); %Predict behavior that was that centriod

%plot step function
figure
plot(LD_tog(s_inds,1), 'LineWidth',2)
hold on
plot(preds(s_inds), 'LineWidth',2)
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


per_cor = sum(LD_tog(:,1)==preds')/length(preds)*100


