%% Running Notes
%As first pass analysis for trajectories/manifold idea, finding population
%average (centroid) for each behavior and using this to predict moment to
%moment behavioral state.

%% Load in data

is_ron = 1;

if is_ron >0
    
    load_folder = 'C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_data\';
    
    
else
    
    'Put your data folder path here'
    
end

session = 'Amos_2021-07-29';

%Pull need variables, push struct into separate variables in case need info
%later

Label_struct = load([load_folder session '\Labels_per_sec.mat']);
Data_struct = load([load_folder session '\Neural_data.mat']);

Labels = cell2mat(Label_struct.labels(:,3)); %Get numerical labels based on prioritization
Data = Data_struct.Unit_rasters; %Get neural data



%% Linear embedding dim and dim reduction
% Only doing this to check linear embedding dimension and to aid visualization. 


[coeff, score,~,~,explained,mu] = pca(Data', 'Centered', true); %center data in algorithm

num_component = find(cumsum(explained)>85,1,'first'); %use number of components that explains more than 85% of the variance

if num_component > 20 %want to use 10 or less, 20 is limit where nearest neighbor falls apart for high dim data
    
    disp(['number of components = ' num2str(num_component)])
    warning('Using more than 20 components in dim reduced space.  Changing from L2 to L1 measure')
    
end


DR_data = mu + score*coeff;
DR_data = DR_data(:,1:num_component)'; %back to our normal form of neurons


%% Sort and summarize data, calculation centriod

%Quick histogram of behaviors

C = categorical(Labels,1:20,Label_struct.behav_categ);
figure()
histogram(C,Label_struct.behav_categ(1:end-1)) %removing rest since it dominates data


%Take advantage of numeric labels to get neural activity associated with
%each time a behavior was present

[sorted,s_inds] = sort(Labels);

Data_group = DR_data(:,s_inds);

LD_tog = [sorted'; Data_group];

centroids = cell(length(Label_struct.behav_categ),1); %Will be same size in principle but making cell in case this gets violated later on

%pausing for now as not seeing clustering that Camille saw
figure()
for i = [6 17] %1:length(centroids) %starting with just groom give and receive and self groom since they happen the most
    
    %plot clusters and their centroid on first 3 pcs
    
    plot3(LD_tog(2,LD_tog(1,:) == i),LD_tog(3,LD_tog(1,:) == i),LD_tog(4,LD_tog(1,:) == i),'.') %recal row 1 is the label
    hold on
    centroids{i} = mean(LD_tog(2:end,LD_tog(1,:) == i),2); %get average population response for each behavior (i.e. centroid of the cluster)
    
    
    
end