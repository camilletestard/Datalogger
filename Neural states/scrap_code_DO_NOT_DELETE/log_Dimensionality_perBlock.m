%%log_Dimensionality_perBehav
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
min_occurrences = 200;
agg_precedence=0;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

% dim=nan([length(session_range),2,length(behav),2,num_iter]);

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Dimensionality_results/'];

    chan = 1;

    for channel_flag = ["vlPFC", "TEO", "all"]
        %channel_flag = "vlPFC";

        %% Get data with specified temporal resolution and channels
        %Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma, agg_precedence);
        end

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';

        block_labels = cell2mat({labels{:,13}}');


        %% Select behaviors to decode
        block=[0,1];

        % For all behaviors
        for bl=1:length(block)

            %Only keep the behaviors of interest
            idx = find(ismember(block_labels,block(bl))); %find the indices of the behaviors considered

            Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
            behavior_labels_final = block_labels(idx,:);%Same as above but in behavior labels

                    for iter = 1:num_iter

                        %Select time points to run PCA
                        idx = randsample(1:length(Spike_count_raster_final),min_occurrences);

                        %Select unit to run PCA
                        Input_matrix = zscore(Spike_count_raster(idx,randsample(size(Spike_count_raster,2), num_units)));

                        %figure; hold on; hist(corr(Input_matrix))

                        %PCA
                        [coeff,score,~,~,explained] = pca(Input_matrix);

                        %Get dimensionality
                        var_explained = cumsum(explained);
                        idxl = find(var_explained>=var_explained_threshold);
                        dim(s,chan,bl,iter) = min(idxl);
                    end


        end % end of block loop

        chan=chan+1;
        disp([channel_flag ' done'])
    end %end of channel loop

    %% Plot results for the session

    Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0.1 0.8 0.9];[0 0.7 0];[1 0 1];[0 1 1];...
        [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0.7 0 1];...
        [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0.8 0 0];[1 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
        [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];
% 
%     blck_lbls = {'Paired','Alone'};
% 
%     figure; hold on; set(gcf,'Position',[150 250 1000 400])
% 
%     for b=length(behav)
%         subplot(1,length(behav),b); hold on
%         dim_vlpfc = squeeze(squeeze(dim(s,1,b,:,:)))';
%         mean_dim_vlpfc = squeeze(nanmean(dim(s,1,b,:,:),5)); [~, orderIdx] = sort(mean_dim_vlpfc);
%         %sd_dim_vlpfc = squeeze(nanstd(dim(s,1,:,:),0,4));
%         violin(dim_vlpfc(:,orderIdx), 'facecolor',Cmap(orderIdx,:))
%         xticks([1:length(blocks)]); xlim([0.5 length(blocks)+0.5]);
%         xticklabels(blck_lbls(orderIdx)); ylim([0 35])
%         ax = gca;
%         ax.FontSize = 14;
%         ylabel(['Dimensionality'],'FontSize', 14);
%         title(behav_categ(behav(b)))
%     end
% sgtitle('vlPFC')
% 
%     figure; hold on; set(gcf,'Position',[150 250 1000 400])
%     for b=1:length(behav)
%         subplot(1,length(behav),b); hold on
%         dim_vlpfc = squeeze(squeeze(dim(s,2,b,:,:)))';
%         mean_dim_vlpfc = squeeze(nanmean(dim(s,2,b,:,:),5)); [~, orderIdx] = sort(mean_dim_vlpfc);
%         %sd_dim_vlpfc = squeeze(nanstd(dim(s,1,:,:),0,4));
%         violin(dim_vlpfc(:,orderIdx), 'facecolor',Cmap(orderIdx,:))
%         xticks([1:length(blocks)]); xlim([0.5 length(blocks)+0.5]);
%         xticklabels(blck_lbls(orderIdx)); ylim([0 35])
%         ax = gca;
%         ax.FontSize = 14;
%         ylabel(['Dimensionality'],'FontSize', 14);
%     end
%     title('vlPFC')
%    
%         close all

end%end of session loop

cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/Dimensionality_results/']);
save('Dimensionality_block.mat', "dim","block","a_sessions","h_sessions","behav_categ","Cmap")

%% Plot results across sessions

load('Dimensionality_block.mat')

dim_all = squeeze(nanmean(dim(:,3,:,:),1));
blck_lbls = {'Alone','Paired'};

%Pooling monkeys and areas
figure; hold on; set(gcf,'Position',[150 250 1000 500]); lowlimit=24; uplimit=29;
mean_dim = mean(dim_all,2); [~, orderIdx] = sort(mean_dim);
violin(dim_all(orderIdx,:)', 'facecolor',Cmap(orderIdx,:))
xticks([1:length(blck_lbls)]); xlim([0.5 length(blck_lbls)+0.5]);
xticklabels(blck_lbls(orderIdx)); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
[h, p]=ttest(dim_all(1,:), dim_all(2,:)); abs(computeCohen_d(dim_all(1,:), dim_all(2,:), 'paired'))
saveas(gcf,'DimensionalityPerBlock.pdf')


%Pooling both monkeys
figure; hold on; set(gcf,'Position',[150 250 1000 500]); lowlimit=15; uplimit=30;
subplot(1,2,1); hold on
mean_dim_vlpfc = squeeze(mean(dim_all(1,:,:),3)); [~, orderIdx] = sort(mean_dim_vlpfc);
violin(squeeze(dim_all(1,orderIdx,:))', 'facecolor',Cmap(orderIdx,:))
xticks([1:length(blck_lbls)]); xlim([0.5 length(blck_lbls)+0.5]);
xticklabels(blck_lbls(orderIdx)); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('vlPFC')

subplot(1,2,2); hold on
mean_dim_teo = squeeze(mean(dim_all(2,:,:),3)); [~, orderIdx] = sort(mean_dim_teo);
violin(squeeze(dim_all(2,orderIdx,:))', 'facecolor',Cmap(orderIdx,:))
xticks([1:length(blck_lbls)]); xlim([0.5 length(blck_lbls)+0.5]);
xticklabels(blck_lbls(orderIdx)); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('TEO')
saveas(gcf,'DimensionalityPerBehav_GroomCateg.pdf')



%Separated by monkey
figure; hold on; set(gcf,'Position',[150 250 1000 800]); lowlimit=25; uplimit=45;

%Amos
subplot(2,2,1); hold on
mean_dim_vlpfc = squeeze(mean(dim_amos(1,:,:),3)); [~, orderIdx] = sort(mean_dim_vlpfc);
violin(squeeze(dim_amos(1,orderIdx,:))', 'facecolor',Cmap(block(orderIdx),:))
xticks([1:length(block)]); xlim([0.5 length(block)+0.5]);
xticklabels(behav_categ(block(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Amos, vlPFC')

subplot(2,2,2); hold on
mean_dim_teo = squeeze(mean(dim_amos(2,:,:),3)); [~, orderIdx] = sort(mean_dim_teo);
violin(squeeze(dim_amos(2,orderIdx,:))', 'facecolor',Cmap(block(orderIdx),:))
xticks([1:length(block)]); xlim([0.5 length(block)+0.5]);
xticklabels(behav_categ(block(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Amos, TEO')

%Hooke
subplot(2,2,3); hold on
mean_dim_vlpfc = squeeze(mean(dim_hooke(1,:,:),3)); [~, orderIdx] = sort(mean_dim_vlpfc);
violin(squeeze(dim_hooke(1,orderIdx,:))', 'facecolor',Cmap(block(orderIdx),:))
xticks([1:length(block)]); xlim([0.5 length(block)+0.5]);
xticklabels(behav_categ(block(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Hooke, vlPFC')

subplot(2,2,4); hold on
mean_dim_teo = squeeze(mean(dim_hooke(2,:,:),3)); [~, orderIdx] = sort(mean_dim_teo);
violin(squeeze(dim_hooke(2,orderIdx,:))', 'facecolor',Cmap(block(orderIdx),:))
xticks([1:length(block)]); xlim([0.5 length(block)+0.5]);
xticklabels(behav_categ(block(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Hooke, TEO')

saveas(gcf,'DimensionalityPerBehav_GroomCateg.pdf')