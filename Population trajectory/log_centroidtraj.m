%% Running Notes
%As first pass analysis for trajectories/manifold idea, finding population
%average (centroid) for each behavior and using this to predict moment to
%moment behavioral state.

%% Load in data and preprocess

%Set path
is_mac = 0; is_camille = 0;

if is_camille
    if is_mac
        cd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
    else
        cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
    end
    filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)

    if is_mac
        cd('~/Dropbox (Penn)/Datalogger/Results/')
    else
        cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Results/')
    end
    savePath = uigetdir('', 'Please select the result directory');

else

    %"[RON ADD PATHS]"    
    
     addpath('C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger\Behavior')
     addpath('C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger\Neural preprocessing')
     cd('C:\Users\ronwd\Dropbox\Ready to analyze output'); 
   
     filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
     
     cd('C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_results')
     
     savePath = uigetdir('', 'Please select the results directory');
     
     
     
     
end

clearvars -except savePath filePath temp_resolution channel_flag is_mac

%Set temporal resolution
temp = 1; temp_resolution = 1;
run_cvmnr = 1; %toggle running the cross validated multinomial regression

%pre-choose number of features to use or use 85% variance threshold for PCA
choose_numcom = 1; man_num = 20; %update 2021-12-06 this doesn't seem to effect the trend of vlPFC being worse prediction wise than TEO for the centroid analysis.
randomize = 1;
for temp_resolution = [1, 2, 5, 10] %1sec, 500msec, 200msec, 100msec
    %temp_resolution = [1/5, 1/2, 1, 5, 10] %5sec, 2sec, 1sec,500msec, 100msec
    %1 for second resolution, 10 for 100msec resolution, 100 for 10msec resolution, 1000 for msec resolution. etc.
    %0.1 for 10sec resolution, 1/5 for 5sec resolution

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
        Data = Spike_rasters; %Get neural data

        Z_data = zscore(Data,[],2)';
        %% Sort and summarize data

        %Quick histogram of behaviors

        C = categorical(Labels,1:20,behav_categ);
        figure()
        histogram(C,behav_categ(1:end-1)) %removing "rest" since it dominates data
        pause(1); close all


        %Take advantage of numeric labels to get neural activity associated with
        %each time a behavior was present.  %Not sure the sorting really helps but
        %makes it easier to check things.


        LD_holding = [Labels Z_data];%[Labels Data_group];% %This keeps matrix for all behaviors, now first column is labels

        boi = [4:8 17]; %manually set behaviors of interest
        not_boi = setdiff(1:length(behav_categ), boi);

        index_use = LD_holding(:,1)==boi; %creates numel(boi) vectors of logical indecies for every time point
        index_use = sum(index_use,2)>0; %has ones know where one of the behaviors of interest happened

        LD_tog = LD_holding(index_use,:);

        %% Dim reduction, calculate centriod
        % Only doing this to check linear embedding dimension and to aid visualization.
        %2021-12-3 copying Camille's method to select behaviors to analyze instead
        %of using whole data set at once
        %Just set boi to different values if you want different numbers of
        %behaviors
        
        [coeff, score,latent,~,explained] = pca(LD_tog(:,2:end), 'Centered', false);
        
        if choose_numcom
            
            num_component = man_num;% set desired number of components to use instead of using what explains 85% of the variance
            
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


        centriods = cell(length(behav_categ),1);
        %radii = zeros(length(behav_categ),1); %drawing spheres at some point?

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

        [~,cs] = min(dis_mat,[],2); %Find which centriod the activity was closer to

        preds = boi(cs); %Predict behavior that was that centriod

        %plot centroid results
        close all; figure; set(gcf,'Position',[150 250 1700 800])
        plot(LD_tog(:,1), 'LineWidth',2)
        hold on
        plot(preds, 'LineWidth',2)
        yticks([0:length(behav_categ)+1]);
        ticklabs = behav_categ; ticklabs(not_boi)={' '};
        yticklabels({'',ticklabs{:},''})
        ylim([0 length(behav_categ)+1])
        legend('Real','Predicted State', 'FontSize',16)
        ylabel('Behavior', 'FontSize',18)
        xlabel('Time', 'FontSize',18)
        %xlabel('index')
        ax = gca;
        ax.FontSize = 14;
        per_cor(temp, chan) = sum(LD_tog(:,1)==preds')/length(preds)*100
        title(['Predicted based on neural data vs. Real behavioral state. Resolution: ' num2str(1000/temp_resolution) 'msec. Area: ' channel '. Accuracy: ' num2str(round(per_cor(temp, chan))) '%'])
        

%         cd(savePath)
%         saveas(gcf,[savePath '/Centroid_' num2str(1000/temp_resolution) 'msec_' channel '.png'])
         pause(1); close all
        
        %% Predict behavioral state from neural data using multinomial regression
        

if run_cvmnr        
       

 behavs = categorical(LD_tog(:,1),boi,{behav_categ{boi}}); %create categorical array for use with mnrfit
    if randomize %As it should be, get same result shuffling either way.
        behavs = behavs(randperm(length(behavs))); %trying shuffling the data instead
       % DR_data = DR_data(randperm(length(DR_data)),:);
    end
        
        % Set up k-fold cross-validation
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
        
        per_cor_cvmnr(temp,chan) = mean(cv_per_cor)
        
end        
        
        %Think about how to make plot of CV tomorrow
        
        

        close all


        clearvars -except randomize temp chan channel_flag temp_resolution per_cor savePath filePath boi choose_numcom man_num is_mac per_cor_cvmnr run_cvmnr

        chan = chan +1;
    end
    temp = temp+1;
end

rowNames = ["1sec", "500msec", "200msec", "100msec"]; colNames = ["vlPFC","TEO","all"];
%rowNames = ["1sec"]; colNames = ["vlPFC","TEO","all"];
result_hitrate = array2table(per_cor_cvmnr,'RowNames',rowNames,'VariableNames',colNames)

save([savePath '\MNRegression_results_allBehav.mat'], 'per_cor', 'per_cor_cvmnr', 'boi')
writetable(result_hitrate,[savePath '\MNRegression_results_allBehav.csv'],'WriteRowNames',true,'WriteVariableNames',true); 

%Plot centroid result
figure; hold on; set(gcf,'Position',[150 250 1000 500])
cmap = cool(size(per_cor,1));
for b = 1:size(per_cor,1)
    y = per_cor(b,:);
    x = 1:size(per_cor,2);
    plot(x,y,'s','MarkerSize',15,...
    'MarkerEdgeColor',cmap(b,:),'MarkerFaceColor',cmap(b,:))
    %plot(x,y,'Color','k')
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
saveas(gcf,['Centroid_results_allBehav.png'])
close all

%Plot multinomial regression results
%load([savePath '\MNRegression_results_allBehav.mat'])
figure; hold on; set(gcf,'Position',[150 250 1000 500])
cmap = cool(size(per_cor,1));
for b = 1:size(per_cor,1)
    y = per_cor_cvmnr(b,:);
    x = 1:size(per_cor,2);
    plot(x,y,'s','MarkerSize',15,...
    'MarkerEdgeColor',cmap(b,:),'MarkerFaceColor',cmap(b,:))
    %plot(x,y,'Color','k')
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
saveas(gcf,['MNRegression_results_allBehav.png'])

% 
% is_still_ron = 1;
% 
% if is_still_ron
%     
%     cd('C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger\Population trajectory')
%     
% end