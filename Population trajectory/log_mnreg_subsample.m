%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,15:16,18];
session_range_with_partner=[1:6,11:13,15:16,18];

%Set parameters
with_partner =0;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1;%set the smoothing window size (sigma)
var_explained_threshold=90;
num_iter = 500; num_units = 100;
min_occurrences = 60;
behav=5;%[1,5,18,29];

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

dim=nan([length(session_range),2,length(behav),2,num_iter]);

s=1;

for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Dimensionality_results/'];

    chan = 1;

    %pre-choose number of features to use or use 85% variance threshold for PCA
    choose_numcom = 1; man_num = 20; %update 2021-12-06 this doesn't seem to effect the trend of vlPFC being worse prediction wise than TEO for the centroid analysis.
    randomize = 0;

    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all";
    for channel_flag = ["vlPFC", "TEO", "all"]
        
        channel = char(channel_flag);%for later saving
        
        %Get data with specified temporal resolution and channels
        [Spike_rasters, labels, behav_categ, block_times, monkey]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac);
        %filePath is the experimental data path
        %Temp_resolution is the temporal resolution at which we would like to
        %analyze the dat
        %Channel_flag specifies with channels to include: only TEO array, only
        %vlPFC array or all channels
        disp('Data Loaded')
        
        Labels = cell2mat(labels(:,3)); %Get numerical labels based on prioritization
        tabulate(Labels)
        Data = Spike_rasters; %Get neural data
        
        Z_data = zscore(Data,[],2)';
        %% Sort and summarize data
        
        %Quick histogram of behaviors
        
        C = categorical(Labels,1:length(behav_categ),behav_categ);
        figure()
        histogram(C,behav_categ(1:end-1)) %removing "rest" since it dominates data
        pause(1); close all
        
        
        %Take advantage of numeric labels to get neural activity associated with
        %each time a behavior was present.  %Not sure the sorting really helps but
        %makes it easier to check things.
        
        
        LD_holding = [Labels Z_data];%[Labels Data_group];% %This keeps matrix for all behaviors, now first column is labels
        
        boi = [5,7:10]; %[4:8, 17] %manually set behaviors of interest
        not_boi = setdiff(1:length(behav_categ), boi);
        
        index_use = LD_holding(:,1)==boi; %creates numel(boi) vectors of logical indecies for every time point
        index_use = sum(index_use,2)>0; %has ones know where one of the behaviors of interest happened
        
        LD_tog_temp = LD_holding(index_use,:);
        
        %% Run MNRegression over multiple iterations
        num_iter = 100;
        disp('Start running MNRegression...')
        
        for iter = 1:num_iter
            
            Labels = LD_tog_temp(:,1);
            unqLabels = unique(Labels); NumOfClasses = length(unqLabels);
            minNumTrials = 90;%min(num_trials); %find the minimum one %CT change to have 200 of each class
            chosen_trials = [];
            for i = 1:NumOfClasses %for each class
                idx = find(Labels == unqLabels(i)); %find indexes of trials belonging to this class
                rand_i = randsample(length(idx), minNumTrials); %Select a random n number of them
                chosen_trials = [chosen_trials; idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
            end
            
            LD_tog = LD_tog_temp(chosen_trials,:);
            
            % Run PCA
            [coeff, score,latent,~,explained] = pca(LD_tog(:,2:end), 'Centered', false);
            
            if choose_numcom
                num_component = man_num;% set desired number of components to use instead of using what explains 85% of the variance
            else
                num_component = find(cumsum(explained)>85,1,'first'); %use number of components that explains more than 85% of the variance
            end
            
            DR_data = score(:,1:num_component);
        
            % Run k-fold cross-validation
            
            behavs = categorical(LD_tog(:,1),boi,{behav_categ{boi}}); %create categorical array for use with mnrfit
            if randomize
                behavs = behavs(randperm(length(behavs),length(behavs)));
            end
            
            % Set up k-fold cross-validation
            %Matlab has a built in k-fold function but it isn't clear if it interfaces
            %with mnrfit, so just going to code my own.  Going with 10 portions instead
            %of 10 grabs
            
            folds = 5;
            
            sample_size = floor(length(LD_tog(:,1))/folds);
            
            sam_if = ones(folds,1)*sample_size; %samples in each fold
            
            
            if mod(length(LD_tog(:,1)),folds)>0 %check if not divisible by folds
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
            
            %Run multinomial regression
            %2021-12-07 update: non-trival to increase the number of iterations mnrfit
            %does...so will need to work around this or take the time to make something more custom...or switch to python.
            
            cv_per_cor = nan(folds,1);
            cv_preds = cell(folds,1);
            cv_behavs = cell(folds,1);
            
            for k = 1:folds
                
                disp(['iteration: ' num2str(k)])   
                states = behavs(ind_if{k},1);
                cv_behavs{k} = states;
                foldsidx = 1:folds;
                trainingidx = horzcat(ind_if{foldsidx~=k})'; %train on all indices that aren't in the current fold
                
                [Betas,~,stats] = mnrfit(DR_data(trainingidx,:),behavs(trainingidx,1)); %note always get B is dim predictors+1 for the intercept term x length(boi)-1  as one behavior is selected as reference
                
                [pihat, ~,~] = mnrval(Betas,DR_data(ind_if{k},:),stats); %give probabilities for each behavior
                
                %take max of each predicted probability as the predicted behavioral state
                %like above

                [~,cs]= max(pihat,[],2); %Find which centriod the activity was closer to
                
                preds = boi(cs); %Predict behavior that was that centriod
                
                cv_preds{k} = preds;
                
                %preds = categorical(preds,boi,{Label_struct.behav_categ{boi}}); %use same categorical trick above so can do string compare
                per_cor_mnr = sum(LD_tog(ind_if{k},1)==preds')/length(preds)*100; %strcmp(behavs(ind_if{k}),preds)
                
                cv_per_cor(k) = per_cor_mnr;

            end
            
            per_cor_cvmnr(iter) = mean(cv_per_cor);
            
        
            
            disp(['MNReg run' num2str(iter) '/' num2str(num_iter)])
        end
        
        per_cor_cvmnr_mean(temp, chan) = mean(per_cor_cvmnr)
        per_cor_cvmnr_sd(temp, chan) = std(per_cor_cvmnr)
        
        close all
        
        clearvars -except randomize temp chan channel_flag temp_resolution per_cor savePath filePath boi per_cor_cvmnr_sd choose_numcom man_num is_mac per_cor_cvmnr_mean
        chan = chan +1;
    end
    temp = temp+1;
end

rowNames = ["1sec", "500msec", "200msec", "100msec"]; colNames = ["vlPFC","TEO","all"];
%rowNames = ["1sec"]; colNames = ["vlPFC","TEO","all"];
result_hitrate = array2table(per_cor_cvmnr_mean,'RowNames',rowNames,'VariableNames',colNames)

save([savePath '\MNRegression_results_allBehav.mat'], 'per_cor_cvmnr_mean', 'per_cor_cvmnr_sd', 'boi')
writetable(result_hitrate,[savePath '\MNRegression_results_allBehav.csv'],'WriteRowNames',true,'WriteVariableNames',true);

%Plot multinomial regression results
%load([savePath '\MNRegression_results_allBehav.mat'])
figure; hold on; set(gcf,'Position',[150 250 1000 500])
cmap = cool(size(per_cor_cvmnr_mean,1));
for b = 1:size(per_cor_cvmnr_mean,1)
    y = per_cor_cvmnr_mean(b,:);
    std_dev = per_cor_cvmnr_sd(b,:);
    errorbar(y,std_dev,'s','MarkerSize',10,...
    'MarkerEdgeColor',cmap(b,:),'MarkerFaceColor',cmap(b,:))
end
chance_level = 1/length(boi)*100; yline(chance_level,'--','Chance level', 'FontSize',16)
leg = legend("1sec","500msec","200msec","100msec","chance", 'Location','southwest');
title(leg,'Window size')
xticks([0.8 1 2 3 3.2]); xlim([0.8 3.2]); ylim([0 100])
xticklabels({'','vlPFC','TEO','all',''})
ax = gca;
ax.FontSize = 14;
ylabel('%Accuracy','FontSize', 18); xlabel('Brain area','FontSize', 18)
title('Accuracy of behavioral states prediction based on neural data','FontSize', 20)
cd(savePath)
saveas(gcf,['MNRegression_results_6Behav.png'])

%
% is_still_ron = 1;
%
% if is_still_ron
%
%     cd('C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger\Population trajectory')
%
% end