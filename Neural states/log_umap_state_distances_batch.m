%% log_umap_visualization_batch
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
session_range_no_partner=[1:6,11:13,15:16];
session_range_with_partner=[1:6,11:13,15:16];

%Set parameters
with_partner =0;
temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
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

s=2; chan=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];


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

        %Pull threats to partner and subject together
        behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Threat to partner");
        behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Threat to subject");
        
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

        %% Run umap

        %Supervised
%         data = [Spike_count_raster_final, behavior_labels_final];
%         [umap_result{s,chan}]=run_umap(data, 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3,'label_column', 'end'); %Run umap to get 2d embedded states

        %Unsupervised
        [umap_result]=run_umap(Spike_count_raster_final, 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 20); %Run umap to get 2d embedded states
        close

        channel = char(channel_flag);

        %Set colormap
        Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0.1 0.8 0.9];[0 0.7 0];[1 0 1];[0 1 1];...
            [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0.7 0 1];...
            [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0.8 0 0];[1 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
            [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];

        Cmap_block = [[0.9 0.7 0.12];[0 0.6 0.8];[0.5 0 0]];


        %% Plot UMAP projection in 3D space

        figure; hold on; set(gcf,'Position',[150 250 1200 500])

        %Plot UMAP results color-coded by behavior
        ax1=subplot(1,2,1);
        scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),12,Cmap(behavior_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior')
        set(gca,'FontSize',12);
        %saveas(gcf,[savePath '/umap_supervised_ColorCodedByBehav_' channel 'Units.png'])
        %pause(5)

        %Color-coded by block
        ax2=subplot(1,2,2);
        scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),12,Cmap_block(block_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Block')
        set(gca,'FontSize',12);

        sgtitle([channel ' units, UMAP, ' sessions(s).name])

        hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'});
        rotate3d on


        %% Visualize distances between beahvior states

        plot_behav_labels = behavior_labels_final(block_labels_final==1);
        plot_neural_states = umap_result(block_labels_final==1,:);
        idx_groom = find(ismember(plot_behav_labels,7));
        idx_travel = find(ismember(plot_behav_labels,18));
        idx_getgroom = find(ismember(plot_behav_labels,8));
        %idx_selfgroom = find(ismember(plot_behav_labels,24));
        idx_forage = find(ismember(plot_behav_labels,5));
        idx_hip = find(ismember(plot_behav_labels,9));
        idx_his = find(ismember(plot_behav_labels,10));
        idx_length = 15;

        for iter=1:500
        idx_ordered = [idx_forage(randsample(length(idx_forage),idx_length));...
                       %idx_selfgroom(randsample(length(idx_selfgroom),idx_length));...
                       idx_groom(randsample(length(idx_groom),idx_length));...
                       idx_getgroom(randsample(length(idx_getgroom),idx_length));...
                       idx_travel(randsample(length(idx_travel),idx_length));...
                       idx_hip(randsample(length(idx_hip),idx_length));...
                       idx_his(randsample(length(idx_his),idx_length))];
        D_subsample(iter,:)=pdist(plot_neural_states(idx_ordered,:));
        end

        D = D_subsample(randsample(500,1),:);%mean(D_subsample);
        Z_behav{s} = squareform(D); Z_behav{s}(eye(size(Z_behav{s}))==1) = nan; 
        figure; set(gcf,'Position',[150 250 800 600]); heatmap(Z_behav{s},'Colormap',flipud(hot))
%         saveas(gcf,[savePath '/DistanceBetweenStates_behavior.pdf'])
        
        Z_behav_lowRes{s} = sepblockfun(Z_behav{s},[idx_length,idx_length],@nanmean);
        %figure; set(gcf,'Position',[150 250 800 600]); heatmap(Z_behav_lowRes{s},'Colormap',flipud(jet))

        %Visualize distances paired vs. alone
        plot_neural_states = umap_result;

        idx_block1= find(block_labels_final==1);
        idx_block2= find(block_labels_final==2);
        idx_block3= find(block_labels_final==3);
        
        idx_length = 50;

        for iter=1:500
        idx_ordered = [idx_block1(randsample(length(idx_block1),idx_length));...
                       idx_block2(randsample(length(idx_block2),idx_length));...
                       idx_block3(randsample(length(idx_block3),idx_length))];
        D_subsample_block(iter,:)=pdist(plot_neural_states(idx_ordered,:), 'cityblock');
        end

        D = D_subsample_block(randsample(500,1),:);%mean(D_subsample_block);
        Z_block{s} = squareform(D); Z_block{s}(eye(size(Z_block{s}))==1) = nan; 
        figure; set(gcf,'Position',[150 250 800 600]); heatmap(Z_block{s},'Colormap',flipud(cool))
        %saveas(gcf,[savePath '/DistanceBetweenStates_PairedAloneContext.pdf'])

        Z_block_lowRes{s} = sepblockfun(Z_block{s},[idx_length,idx_length],@nanmean);
        %figure; set(gcf,'Position',[150 250 800 600]); heatmap(Z_block_lowRes{s},'Colormap',flipud(jet))

        %Visualize distances between social context given a particular
        %behavior
        beh=8;

        if ismember(beh,1)
            %Pull threats to partner and subject together
            behavior_labels_final(behavior_labels_final==9)=1;
            behavior_labels_final(behavior_labels_final==10)=1;
        end

        plot_behav_labels = behavior_labels_final(behavior_labels_final==beh);
        plot_block_labels = block_labels_final(behavior_labels_final==beh);
        plot_neural_states = umap_result(behavior_labels_final==beh,:);

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
            D_subsample_blockBehav(iter,:)=pdist(plot_neural_states(idx_ordered,:), 'cityblock');
        end

        D = D_subsample_blockBehav(randsample(500,1),:); % mean(D_subsample_blockBehav);
        Z_blockBehav{s} = squareform(D); Z_blockBehav{s}(eye(size(Z_blockBehav{s}))==1) = nan; 
        figure; set(gcf,'Position',[150 250 800 600]); heatmap(Z_blockBehav{s},'Colormap',flipud(cool))
        %saveas(gcf,[savePath '/DistanceBetweenStates_NeighborContext.pdf'])

        Z_blockBehav_lowRes{s} = sepblockfun(Z_blockBehav{s},[idx_length,idx_length],@nanmean);
        %figure; set(gcf,'Position',[150 250 800 600]); heatmap(Z_blockBehav_lowRes{s},'Colormap',flipud(jet))
        close all


end %end of session for loop

savePath = [home '/Dropbox (Penn)/Datalogger/Results/All_sessions/UMAP_results/'];

D_behav= mean(cat(3,Z_behav_lowRes{:}),3);
AxesLabels = {'Forage','Groom','Get groomed','Travel','Threat partner','Threat subject'};
figure; set(gcf,'Position',[150 250 800 600]); hp= heatmap(D_behav,'Colormap',flipud(jet));
hp.XDisplayLabels = AxesLabels; hp.YDisplayLabels = AxesLabels;
saveas(gcf,[savePath '/DistanceBetweenBehav_allSessions.pdf'])
        
D_block= mean(cat(3,Z_block_lowRes{:}),3);
AxesLabels = {'Female Neighbor','Male Neighbor','Alone'};
figure; set(gcf,'Position',[150 250 800 600]); hp= heatmap(D_block,'Colormap',flipud(jet));
hp.XDisplayLabels = AxesLabels; hp.YDisplayLabels = AxesLabels;
saveas(gcf,[savePath '/DistanceBetweenBlocks_allSessions.pdf'])
        
D_blockBehav= mean(cat(3,Z_blockBehav_lowRes{:}),3);
AxesLabels = {'Female Neighbor','Male Neighbor'};
figure; set(gcf,'Position',[150 250 800 600]); hp= heatmap(D_blockBehav,'Colormap',flipud(jet));
hp.XDisplayLabels = AxesLabels; hp.YDisplayLabels = AxesLabels;