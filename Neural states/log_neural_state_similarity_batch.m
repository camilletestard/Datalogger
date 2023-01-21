%% log_neural_state_similarity_batch
% This script applied unsupervised umap on smoothed firing
% Ref: Connor Meehan, Jonathan Ebrahimian, Wayne Moore, and Stephen Meehan (2022). Uniform Manifold Approximation and Projection (UMAP) (https://www.mathworks.com/matlabcentral/fileexchange/71902), MATLAB Central File Exchange.

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
temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 0; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=1; %lump similar behavioral categories together
agg_precedence=0; % 1: aggression takes precedence; 0: Threat to partner and subject states take precedence

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1; chan=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];


%     for channel_flag = ["vlPFC", "TEO", "all"]
        %channel_flag = "TEO";

        %% Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma, agg_precedence);
        end

        disp('Data Loaded')

        %Raw data
        Spike_count_raster = Spike_rasters';

        %Low-pass filter
        %Spike_count_raster = lowpass(Spike_rasters',0.005,1);
        %PCA
        % [coeff,score,latent,tsquared,explained] = pca(Spike_rasters');
        % Spike_count_raster = score(:,1:15);


        %% Select behaviors to visualize

        %Extract behavior labels and frequency
        behavior_labels = cell2mat({labels{:,3}}');
        behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

        %Lump all travel together
        behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
        behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

        %Extract block labels
        block_labels = cell2mat({labels{:,12}}');

        % Select behaviors

        %Compute freq of behavior for the session
        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

        % Select behaviors with a minimum # of occurrences
        min_occurrences = 30;
        behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences

        %Remove behaviors we're not interested in for now
        behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rowdy Room')));%excluding Rowdy Room which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Other monkeys vocalize')));

        % OR select behaviors manually
        behav=[behav;29];
        %behav =[7] ;%unique(behavior_labels); %[4,5,7,8,9,10,24];% [4:10, 23]; %[4:8,17]; %manually select behaviors of interest

        %Print behaviors selected
        behavs_eval = behav_categ(behav);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        fprintf('Behaviors evaluated are: %s \n', behavs_eval);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        %Only consider indices with behavior of interest
        idx= find(ismember(behavior_labels,behav));
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx);%Same as above but in behavior labels
        block_labels_final =  block_labels(idx);
        behavior_labels_final_rand = randsample(behavior_labels_final, length(behavior_labels_final));

        if null
            %Simulate fake labels
            [sim_behav] = GenSimBehavior(behavior_labels_final,behav_categ, temp_resolution);
            behavior_labels_final = sim_behav;
        end


        %Set colormap
        Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0.1 0.8 0.9];[0 0.7 0];[1 0 1];[0 1 1];...
            [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0.7 0 1];...
            [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0.8 0 0];[1 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
            [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];

        Cmap_block = [[0.9 0.7 0.12];[0 0.6 0.8];[0.5 0 0]];



        %% Visualize distances between behavior states
        %Within a block

        plot_behav_labels = behavior_labels_final(block_labels_final==1);
        plot_neural_states = Spike_count_raster_final(block_labels_final==1,:);

        idx_groom = find(ismember(plot_behav_labels,7));
        idx_travel = find(ismember(plot_behav_labels,18));
        idx_getgroom = find(ismember(plot_behav_labels,8));
        %idx_selfgroom = find(ismember(plot_behav_labels,24));
        idx_forage = find(ismember(plot_behav_labels,5));
        idx_hip = find(ismember(plot_behav_labels,9));
        idx_his = find(ismember(plot_behav_labels,10));
        idx_rest = find(ismember(plot_behav_labels,29));
        idx_length = 15;

        for iter=1:500
            idx_ordered = [randsample(idx_rest,idx_length);...
                randsample(idx_forage,idx_length);...
                %idx_selfgroom(randsample(length(idx_selfgroom),idx_length));...
                randsample(idx_groom,idx_length);...
                randsample(idx_getgroom,idx_length);...
                randsample(idx_travel,idx_length);...
                randsample(idx_hip,idx_length);...
                randsample(idx_his,idx_length)];

            D_subsample(iter,:,:)=corrcoef(plot_neural_states(idx_ordered,:)');
        end
        %corrcoef(D_subsample_pca'); correlation across iterations is much
        %lower for pca than umap. Umap is more repeatable than pca.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Check for whole session
        test = reshape(corrcoef(Spike_count_raster(randsample(size(Spike_count_raster,1),30000),:)'),1,[]);
        figure; hist(test,500); xlabel('Correlation coefficient across time points'); title('Datalogger (bin size = 33msec)')
        figure; scatter(Spike_count_raster(2039,:), Spike_count_raster(6039,:)); corr([Spike_count_raster(2039,:); Spike_count_raster(6039,:)]')

        for r=1:size(Spike_count_raster,1)
        shuffledData(r,:) = Spike_count_raster(r,randsample(size(Spike_count_raster,2),size(Spike_count_raster,2)));
        end

        test2 = reshape(corrcoef(shuffledData'),1,[]);
        figure; hist(test2,500); xlabel('Correlation coefficient across time points'); title('Datalogger (Shuffled)')

        test2 = reshape(corrcoef(Spike_count_raster),1,[]);
        figure; hist(test2,20); xlabel('Correlation coefficient across neurons')

        [acf,lags] = autocorr(Spike_count_raster(:,1));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        D = squeeze(mean(D_subsample,1)); 
        figure; set(gcf,'Position',[150 250 800 600]); heatmap(D,'Colormap',flipud(hot))
        caxis([0.7 1])

        D_lowRes= sepblockfun(D,[idx_length,idx_length],@nanmean);
        figure; set(gcf,'Position',[150 250 800 600]); heatmap(D_lowRes,'Colormap',flipud(hot))
        caxis([0.7 1])



        %% Visualize distances paired vs. alone

        plot_neural_states = Spike_count_raster_final;

        idx_block1= find(block_labels_final==1);
        idx_block2= find(block_labels_final==2);
        idx_block3= find(block_labels_final==3);

        idx_length = 100;

        for iter=1:500
            idx_ordered = [idx_block1(randsample(length(idx_block1),idx_length));...
                idx_block2(randsample(length(idx_block2),idx_length));...
                idx_block3(randsample(length(idx_block3),idx_length))];
            D_subsample_block(iter,:,:)=corrcoef(plot_neural_states(idx_ordered,:)');
        end

       D = squeeze(mean(D_subsample_block,1)); 
        figure; set(gcf,'Position',[150 250 800 600]); heatmap(D,'Colormap',flipud(hot))
        caxis([0.8 1])

        D_lowRes= sepblockfun(D,[idx_length,idx_length],@nanmean);
        figure; set(gcf,'Position',[150 250 800 600]); heatmap(D_lowRes,'Colormap',flipud(hot))
        caxis([0.8 1])


        %% Visualize distances between social context given a particular behavior
        beh=29;

        if ismember(beh,1)
            %Pull threats to partner and subject together
            behavior_labels_final(behavior_labels_final==9)=1;
            behavior_labels_final(behavior_labels_final==10)=1;
        end

        plot_behav_labels = behavior_labels_final(behavior_labels_final==beh);
        plot_block_labels = block_labels_final(behavior_labels_final==beh);
        plot_neural_states = Spike_count_raster_final(behavior_labels_final==beh,:);

        idx_block1= find(plot_block_labels==1);
        idx_block2= find(plot_block_labels==2);
        idx_block3= find(plot_block_labels==3);

        for iter=1:500
            if ismember(beh, [7,8])
                idx_length = 100;
                idx_ordered = [idx_block1(randsample(length(idx_block1),idx_length));...
                    idx_block2(randsample(length(idx_block2),idx_length))];
            else
                idx_length = min([length(idx_block1), length(idx_block2), length(idx_block3)]);
                idx_ordered = [idx_block1(randsample(length(idx_block1),idx_length));...
                    idx_block2(randsample(length(idx_block2),idx_length));...
                    idx_block3(randsample(length(idx_block3),idx_length))];
            end
            D_subsample_blockBehav(iter,:)=corrcoef(plot_neural_states(idx_ordered,:)');
        end

        D = squeeze(mean(D_subsample_blockBehav,1)); 
        figure; set(gcf,'Position',[150 250 800 600]); heatmap(D,'Colormap',flipud(hot))
        caxis([0.8 1])

        D_lowRes= sepblockfun(D,[idx_length,idx_length],@nanmean);
        figure; set(gcf,'Position',[150 250 800 600]); heatmap(D_lowRes,'Colormap',flipud(hot))
        caxis([0.8 1])

         close all

%     end %end of channel loop

    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    disp(s)
    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

end %end of session for loop

%Save data
cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/UMAP_results/'])
% save('State_distances.mat', "behav_categ","Z_behav_umap","Z_behav_pca",...
%     "Z_behav_umap_lowRes", "Z_behav_pca_lowRes","Z_umap_block","Z_pca_block",...
%     "Z_umap_block_lowRes","Z_pca_block_lowRes","Z_umap_blockBehav","Z_pca_blockBehav",...
%     "Z_umap_blockBehav_lowRes","Z_pca_blockBehav_lowRes")
load('State_similarity.mat')

%Plot distance to threat to subject
behav=[5,7,8,18,9];
data = cat(3,Z_behav_umap_lowRes{:}); data = squeeze(data(6,1:5,:));
[~, idx_sorted]=sort(nanmean(data,2),1,'descend');
figure; hold on
bp = bar([nanmean(data(idx_sorted,:),2)],'FaceAlpha',0.2);
xticks([1:5]); xticklabels(behav_categ(behav(idx_sorted)))
ylabel('Eucledian distance to Threat to Subject')
scatter(ones(size(data,2))*1,data(idx_sorted(1),:), 'filled','r');
scatter(ones(size(data,2))*2,data(idx_sorted(2),:), 'filled','y');
scatter(ones(size(data,2))*3,data(idx_sorted(3),:), 'filled','g');
scatter(ones(size(data,2))*4,data(idx_sorted(4),:), 'filled','b');
scatter(ones(size(data,2))*5,data(idx_sorted(5),:), 'filled','c');
saveas(gcf,['DistanceBetweenBehavStates_allSessions.pdf'])

data_similarity = 1-data./(max(max(data)));
[~, idx_sorted]=sort(nanmean(data_similarity,2),1,'descend');
figure; hold on
bp = bar([nanmean(data_similarity(idx_sorted,:),2)],'FaceAlpha',0.2);
xticks([1:5]); xticklabels(behav_categ(behav(idx_sorted)))
ylabel('Eucledian distance to Threat to Subject')
scatter(ones(size(data_similarity,2))*1,data_similarity(idx_sorted(1),:), 'filled','r');
scatter(ones(size(data_similarity,2))*2,data_similarity(idx_sorted(2),:), 'filled','y');
scatter(ones(size(data_similarity,2))*3,data_similarity(idx_sorted(3),:), 'filled','g');
scatter(ones(size(data_similarity,2))*4,data_similarity(idx_sorted(4),:), 'filled','b');
scatter(ones(size(data_similarity,2))*5,data_similarity(idx_sorted(5),:), 'filled','c');

%UMAP
D_behav= mean(cat(3,Z_behav_umap_lowRes{:}),3);
AxesLabels = {'Forage','Groom','Get groomed','Travel','Threat partner','Threat subject'};
idx = tril(D_behav);
D_behav(~idx) = nan;
figure; set(gcf,'Position',[150 250 800 600]); hp= heatmap(D_behav,'Colormap',hot, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
hp.XDisplayLabels = AxesLabels; hp.YDisplayLabels = AxesLabels;
saveas(gcf,['DistanceBetweenBehav_allSessions.pdf'])

%PCA
D_behav= mean(cat(3,Z_behav_pca_lowRes{:}),3);
AxesLabels = {'Forage','Groom','Get groomed','Travel','Threat partner','Threat subject'};
idx = tril(D_behav);
D_behav(~idx) = nan;
figure; set(gcf,'Position',[150 250 800 600]); hp= heatmap(D_behav,'Colormap',cool, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
hp.XDisplayLabels = AxesLabels; hp.YDisplayLabels = AxesLabels;

%Behavior, separated by monkey
figure; set(gcf,'Position',[150 250 1500 300]); subplot(1,2,1);
behav_order = [2,3,4,1,5,6];
D_behav= mean(cat(3,Z_behav_lowRes{a_sessions}),3);
D_behav= D_behav(behav_order,behav_order);
AxesLabels = {'Forage','Groom','Get groomed','Travel','Threat partner','Threat subject'};
AxesLabels_ordered = AxesLabels(behav_order);
idx = tril(D_behav);
D_behav(~idx) = nan;
hp= heatmap(D_behav,'Colormap',hot, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
caxis([0 10]); title('Amos')
hp.XDisplayLabels = AxesLabels_ordered; hp.YDisplayLabels = AxesLabels_ordered;

subplot(1,2,2);
D_behav= mean(cat(3,Z_behav_lowRes{h_sessions(1:end-1)}),3);
D_behav= D_behav(behav_order,behav_order);
idx = triu(D_behav);
D_behav(~idx) = nan;
hp= heatmap(D_behav,'Colormap',hot, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
caxis([0 10]); title('Hooke')
hp.XDisplayLabels = AxesLabels_ordered; hp.YDisplayLabels = AxesLabels_ordered;
saveas(gcf,['DistanceBetweenBehav_allSessions.pdf'])

%%%%%%%%%%%%%%%%%
% By block
D_block= mean(cat(3,Z_block_lowRes{:}),3);
AxesLabels = {'Female Neighbor','Male Neighbor','Alone'};
figure; set(gcf,'Position',[150 250 800 600]); hp= heatmap(D_block,'Colormap',hot);
hp.XDisplayLabels = AxesLabels; hp.YDisplayLabels = AxesLabels;
saveas(gcf,['DistanceBetweenBlocks_allSessions.pdf'])

%Block, separated by monkey
figure; set(gcf,'Position',[150 250 1500 300]); subplot(1,2,1);
behav_order = [1,2,3];
D_block= mean(cat(3,Z_block_lowRes{a_sessions}),3);
D_block= D_block(behav_order,behav_order);
AxesLabels = {'Female Neighbor','Male Neighbor','Alone'};
AxesLabels_ordered = AxesLabels(behav_order);
idx = tril(D_block);
D_block(~idx) = nan;
hp= heatmap(D_block,'Colormap',hot, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
caxis([5 12]); title('Amos')
hp.XDisplayLabels = AxesLabels_ordered; hp.YDisplayLabels = AxesLabels_ordered;

subplot(1,2,2);
D_block= mean(cat(3,Z_block_lowRes{h_sessions}),3);
D_block= D_block(behav_order,behav_order);
idx = triu(D_block);
D_block(~idx) = nan;
hp= heatmap(D_block,'Colormap',hot, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
caxis([5 12]); title('Hooke')
hp.XDisplayLabels = AxesLabels_ordered; hp.YDisplayLabels = AxesLabels_ordered;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Behav and block

%Block by behav, separated by monkey
figure; set(gcf,'Position',[150 250 1500 300]); subplot(1,2,1);
behav_order = [1,2,3];
D_blockBehav= mean(cat(3,Z_blockBehav_lowRes{a_sessions}),3);
D_blockBehav= D_blockBehav(behav_order,behav_order);
AxesLabels = {'Female Neighbor','Male Neighbor','Alone'};
AxesLabels_ordered = AxesLabels(behav_order);
idx = tril(D_blockBehav);
D_blockBehav(~idx) = nan;
hp= heatmap(D_blockBehav,'Colormap',hot, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
caxis([1 6]); title('Amos')
hp.XDisplayLabels = AxesLabels_ordered; hp.YDisplayLabels = AxesLabels_ordered;

subplot(1,2,2);
D_blockBehav= mean(cat(3,Z_blockBehav_lowRes{h_sessions(1:end-1)}),3);
D_blockBehav= D_blockBehav(behav_order,behav_order);
idx = triu(D_blockBehav);
D_blockBehav(~idx) = nan;
hp= heatmap(D_blockBehav,'Colormap',hot, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
caxis([1 6]); title('Hooke')
hp.XDisplayLabels = AxesLabels_ordered; hp.YDisplayLabels = AxesLabels_ordered;
