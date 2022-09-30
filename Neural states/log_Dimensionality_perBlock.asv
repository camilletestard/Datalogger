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

    for channel_flag = ["vlPFC", "TEO"]
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
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        end

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';

        behavior_labels = cell2mat({labels{:,3}}');
        behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Threat to partner");
        behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Threat to subject");
        behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest
        block_labels = cell2mat({labels{:,12}}');


        %% Select behaviors to decode
        %Simplify behavioral catagories
        %Lump all aggressive interactions together
        behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
        behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");
        behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Aggression");
        behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Aggression");

        %Lump all travel together
        behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
        behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

        %Lump Drinking and foraging
        behavior_labels(behavior_labels==find(behav_categ=="Drinking"))=find(behav_categ=="Foraging");

        % For all behaviors
        for b=1:length(behav)

            disp(['Behavior: ' behav_categ(behav(b))])

            %Only keep the behaviors of interest
            idx = find(ismember(behavior_labels,behav(b))); %find the indices of the behaviors considered

            Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
            behavior_labels_final = block_labels(idx,:);%Same as above but in behavior labels
            behavior_labels_final(behavior_labels_final==2)=1; %pool the paired blocks together
            behavfreq= tabulate(behavior_labels_final);
            blocks = unique(behavior_labels_final);

            if all(behavfreq([1,3],2)>=min_occurrences)

                for bl = 1:length(blocks)

                    for iter = 1:num_iter

                        %Select time points to run PCA
                        idx= find(ismember(behavior_labels_final,blocks(bl)));
                        idx_beh = idx(randsample(1:length(idx),min_occurrences));

                        %Select unit to run PCA
                        Input_matrix = zscore(Spike_count_raster(idx_beh,randsample(size(Spike_count_raster,2), num_units)));

                        %figure; hold on; hist(corr(Input_matrix))

                        %PCA
                        [coeff,score,~,~,explained] = pca(Input_matrix);

                        %Get dimensionality
                        var_explained = cumsum(explained);
                        idxl = find(var_explained>=var_explained_threshold);
                        dim(s,chan,b,bl,iter) = min(idxl);

                    end % end of interation loop

                end %end of block loop

            else 

                disp(['Not enough occurrences across blocks for: ' behav_categ(behav(b))])

            end % end of if clause

        end % end of behavior loop

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
save('Dimensionality_GroomCateg.mat', "dim","behav","a_sessions","h_sessions","behav_categ","Cmap")

%% Plot results across sessions

load('Dimensionality_GroomCateg.mat')

dim_amos = squeeze(nanmean(dim(a_sessions,:,1,:,:),1));
dim_hooke = squeeze(nanmean(dim(h_sessions,:,1,:,:),1));
dim_all = squeeze(nanmean(dim(:,:,1,:,:),1));
blck_lbls = {'Paired','Alone'};

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
violin(squeeze(dim_amos(1,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Amos, vlPFC')

subplot(2,2,2); hold on
mean_dim_teo = squeeze(mean(dim_amos(2,:,:),3)); [~, orderIdx] = sort(mean_dim_teo);
violin(squeeze(dim_amos(2,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Amos, TEO')

%Hooke
subplot(2,2,3); hold on
mean_dim_vlpfc = squeeze(mean(dim_hooke(1,:,:),3)); [~, orderIdx] = sort(mean_dim_vlpfc);
violin(squeeze(dim_hooke(1,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Hooke, vlPFC')

subplot(2,2,4); hold on
mean_dim_teo = squeeze(mean(dim_hooke(2,:,:),3)); [~, orderIdx] = sort(mean_dim_teo);
violin(squeeze(dim_hooke(2,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Hooke, TEO')

saveas(gcf,'DimensionalityPerBehav_GroomCateg.pdf')