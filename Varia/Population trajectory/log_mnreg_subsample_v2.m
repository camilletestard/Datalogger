%% Multinomial regressionn
% Get probability of a behavior based on neural data

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
num_iter = 500; 

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=1;

for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/'];

    chan = 1;

    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "vlPFC";

    for channel_flag = ["vlPFC", "TEO", "all"]
        

        channel = char(channel_flag);%for later saving
        
        %% Get data with specified temporal resolution and channels
        %Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        end

        disp('Data Loaded')
        
        %% Extract labels

        Labels = cell2mat(labels(:,3)); %Get numerical labels based on prioritization

        %Adjust labels
        Labels(Labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Threat to partner");
        Labels(Labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Threat to subject");
        Labels(Labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

        %Lump all aggressive interactions together
        Labels(Labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
        Labels(Labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");
  
        %Lump all travel together
        Labels(Labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
        Labels(Labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

        tabulate(Labels)
        Data = Spike_rasters; %Get neural data


        Z_data = zscore(Data,[],2)';


        %% Sort and summarize data
        
        %Quick histogram of behaviors
        
        C = categorical(Labels,1:length(behav_categ),behav_categ);
        figure()
        histogram(C,behav_categ) %removing "rest" since it dominates data
        pause(1); close all
        
        
        %Take advantage of numeric labels to get neural activity associated with
        %each time a behavior was present.  %Not sure the sorting really helps but
        %makes it easier to check things.
        
        
        LD_holding = [Labels Z_data];%[Labels Data_group];% %This keeps matrix for all behaviors, now first column is labels
        
        boi = [5,7,24,29]; %[4:8, 17] %manually set behaviors of interest
        not_boi = setdiff(1:length(behav_categ), boi);
        
        index_use = LD_holding(:,1)==boi; %creates numel(boi) vectors of logical indecies for every time point
        index_use = sum(index_use,2)>0; %has ones know where one of the behaviors of interest happened
        
        LD_tog_temp = LD_holding(index_use,:);
        
        labels_temp = LD_tog_temp(:,1);
        classes = unique(labels_temp);
        for i = 1:length(classes) %for each class
            idx = find(labels_temp == classes(i)); %find indexes of trials belonging to this class
            labels_temp(idx)=i;
        end
        LD_tog_temp(:,1) = labels_temp;

% % % %         %%CAMILLE TEST
% % % %         % Run PCA
% % % %         
% % % %         [coeff, score,latent,~,explained] = pca(LD_tog_temp(:,2:end), 'Centered', false);
% % % %         num_component = find(cumsum(explained)>85,1,'first');
% % % %         figure; plot(cumsum(explained))
% % % % 
% % % %         data = score(:,200);
% % % %         [B,~,stats] = mnrfit(data,LD_tog_temp(:,1));
% % % %         pihat = mnrval(B,data,stats);
% % % %         figure; hist(pihat(:,4))
% % % % 
% % % %         behav = 4;
% % % %         figure; hold on
% % % %         plot(pihat)
% % % %         labels_plot = nan(size(Labels)); labels_plot(Labels==behav)=1;
% % % %         plot(labels_plot)
        %%%%%%%%%%%%%%%%%%%%
        
        %% Run MNRegression over multiple iterations
        num_iter = 100;
        disp('Start running MNRegression...')
        
        for iter = 1:num_iter

            Lbls = LD_tog_temp(:,1);
            data_details = tabulate(Lbls);
            unqLabels = unique(Lbls); NumOfClasses = length(unqLabels);
            minNumTrials = min(data_details(:,2)); %find the minimum one %CT change to have 200 of each class
            chosen_trials = [];
            for i = 1:NumOfClasses %for each class
                idx = find(Lbls == unqLabels(i)); %find indexes of trials belonging to this class
                rand_i = randsample(length(idx), minNumTrials); %Select a random n number of them
                chosen_trials = [chosen_trials; idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
            end

            LD_tog = LD_tog_temp(chosen_trials,:);

            % Run PCA
            [coeff, score,latent,~,explained] = pca(LD_tog_temp(:,2:end), 'Centered', false);
            [~, score_allD,~,~,explained_allD] = pca(LD_holding(:,2:end), 'Centered', false);

            %             if choose_numcom
            %                 num_component = man_num;% set desired number of components to use instead of using what explains 85% of the variance
            %             else
            num_component = find(cumsum(explained)>90,1,'first'); %use number of components that explains more than 85% of the variance
            num_component_allD = find(cumsum(explained_allD)>90,1,'first');
            %             end

            %DR_data = score(chosen_trials,1:num_component);
            DR_data = LD_tog_temp(chosen_trials,2:end);

            %CT test code
            [B,~,stats] = mnrfit(DR_data,LD_tog(:,1));
            %pihat = mnrval(B,score(:,1:num_component),stats);
            pihat = mnrval(B,LD_tog_temp(:,2:end),stats);
            pihat_allD = mnrval(B,score_allD(:,1:num_component),stats);
            [~,cs]= max(pihat,[],2); %get max predictions
            [~,cs_allD]= max(pihat_allD,[],2);

            %check the distribution is bimodal
            behav = 7;
            figure; hist(pihat(:,find(boi==behav)))
            figure; hist(pihat_allD(:,find(boi==behav)))

            %Only considering behaviors of interest
            figure; hold on
            plot(pihat(:,find(boi==behav)))
            pred_plot = nan(size(cs)); pred_plot(cs==find(boi==behav))=1;
            plot(pred_plot, 'LineWidth',6)
            labels_plot = nan(size(Lbls)); labels_plot(Lbls==find(boi==behav))=1.05;
            plot(labels_plot, 'LineWidth',6)
            ylim([0 1.1])

            figure; hold on
            for b=1:length(boi)
                pred_plot = nan(size(cs)); pred_plot(cs==b)=b;
                plot(pred_plot, 'LineWidth',3)
                labels_plot = nan(size(Lbls)); labels_plot(Lbls==b)=b-0.1;
                plot(labels_plot, 'LineWidth',3)
            end
            ylim([0 length(boi)+1])

            %Plotting pihat before and after start of event
            behav_times = zeros(size(Lbls)); behav_times(Lbls==find(boi==behav))=1;
            behav_starts = find(diff(behav_times)>0)+1;
            behav_ends = find(diff(behav_times)<0);
            bout_length = behav_ends-behav_starts;
            behav_starts_final = behav_starts(bout_length>=10);

%             %Check the ends and start fall at the right times
%             figure; hold on
%             plot(behav_times)
%             for i = 1:length(behav_starts)
%                 xline(behav_starts(i), 'Color','r','LineWidth',2)
%                 xline(behav_ends(i), 'Color','g','LineWidth',2)
%             end

            window = 5;

            pihat_around_behav = nan(length(behav_starts_final),window*2+1);
            for i = 1:length(behav_starts_final)
                pihat_around_behav(i,:) = pihat_allD(behav_starts(i)-window: behav_starts(i)+window,find(boi==behav));    
            end

            figure; hold on
            ylim([0 1.05]); xlim([0 window*2+1]);
            xline(window+1, 'Color','r','LineWidth',3)
            plot(mean(pihat_around_behav), 'LineWidth',3)
            plot(pihat_around_behav', 'Color',[0.5 0.5 0.5])

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Considering the whole session

            %Plot pihat, predicted times for behavior 'behav' and actual
            %times.
            figure; hold on
            plot(pihat_allD(:,find(boi==behav)))
            pred_plot = nan(size(cs_allD)); pred_plot(cs_allD==find(boi==behav))=1;
            plot(pred_plot, 'LineWidth',6)
            labels_plot = nan(size(LD_holding(:,1))); labels_plot(LD_holding(:,1)==behav)=1.05;
            plot(labels_plot, 'LineWidth',6)
            ylim([0 1.1])

            %Plot predicted vs actual behaviors
            figure; hold on
            for b=1:length(boi)
                pred_plot = nan(size(cs_allD)); pred_plot(cs_allD==b)=b;
                plot(pred_plot, 'LineWidth',3)
                labels_plot = nan(size(LD_holding(:,1))); labels_plot(LD_holding(:,1)==boi(b))=b-0.1;
                plot(labels_plot, 'LineWidth',3)
            end
            ylim([0 length(boi)+1])

            %Plotting pihat before and after start of event
            idx_behav = find(~cellfun(@isempty,(strfind(table2cell(behavior_log(:,'Behavior')), behav_categ_original(behav)))));
            idx_dur = find(table2array(behavior_log(:,'duration_round'))>3);
            idx = intersect (idx_behav, idx_dur);
            behav_starts = table2array(behavior_log(idx,"start_time_round")); %start times in msec
            
            window = 5;

            pihat_around_behav = nan(length(behav_starts),window*2+1);
            for i = 1:length(behav_starts)
                pihat_around_behav(i,:) = pihat_allD(behav_starts(i)-window: behav_starts(i)+window,find(boi==behav));    
            end

            figure; hold on
            ylim([0 1.05]); xlim([0 window*2+1]);
            xline(window+1, 'Color','r','LineWidth',3)
            plot(mean(pihat_around_behav), 'LineWidth',3)
            plot(pihat_around_behav', 'Color',[0.5 0.5 0.5])



% %             %preds = categorical(preds,boi,{Label_struct.behav_categ{boi}}); %use same categorical trick above so can do string compare
% %             per_cor_mnr = sum(LD_tog(ind_if{k},1)==preds')/length(preds)*100; %strcmp(behavs(ind_if{k}),preds)
% % 
% %             cv_per_cor(k) = per_cor_mnr;
% % 
% %             per_cor_cvmnr(iter) = mean(cv_per_cor);

        
            
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