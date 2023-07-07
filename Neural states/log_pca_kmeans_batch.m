%% log_pca_kmeans_batch
% This script applied unsupervised umap on smoothed firing 
% Ref: Connor Meehan, Jonathan Ebrahimian, Wayne Moore, and Stephen Meehan (2022). Uniform Manifold Approximation and Projection (PCA) (https://www.mathworks.com/matlabcentral/fileexchange/71902), MATLAB Central File Exchange.

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
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null 
simplify=1; %lump similar behavioral categories together
threat_precedence=0; % 1: aggression takes precedence; 0: Threat to partner and subject states take precedence
exclude_sq=1;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/PCA_results/'];

    chan = 1;

    for channel_flag = ["vlPFC", "TEO"]
        %channel_flag = "TEO";


        %% Get data with specified temporal resolution and channels
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

        disp('Data Loaded')


        %Raw data
        Spike_count_raster = zscore(Spike_rasters');


        %% Select behaviors to visualize

        %Extract behavior labels and frequency
        behavior_labels = cell2mat({labels{:,3}}');
       
        %Extract block labels
        block_labels = cell2mat({labels{:,12}}');
        block_categ={"F","M","Alone"};
        
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
        %behav =[29] ;%unique(behavior_labels); %[4,5,7,8,9,10,24];% [4:10, 23]; %[4:8,17]; %manually select behaviors of interest

        %Print behaviors selected
        behavs_eval = behav_categ(behav);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        fprintf('Behaviors evaluated are: %s \n', behavs_eval);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        %Only consider indices with behavior of interest
        idx= find(ismember(behavior_labels,behav));% & ismember(block_labels,3));
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx);%Same as above but in behavior labels
        block_labels_final =  block_labels(idx);
 
        if null
            %Simulate fake labels
            [sim_behav] = GenSimBehavior(behavior_labels_final,behav_categ, temp_resolution);
            behavior_labels_final = sim_behav;
        end

        %% Run pca

        [umap_result]=run_umap(Spike_count_raster_final, 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3); %Run umap to get 2d embedded states
        close

        [~, pca_result]=pca(Spike_count_raster_final); %Run umap to get 2d embedded states
        pca_result=pca_result(:,1:50);
        close

        channel = char(channel_flag);


        %% Run Kmeans clustering

%         eva = evalclusters(pca_result,'kmeans','CalinskiHarabasz','KList',3:20); eva.OptimalK
        n_clusters = length(behav);%eva.OptimalK; %30
        idx_cluster_pca = kmeans(pca_result,n_clusters);
        idx_cluster_umap = kmeans(umap_result,n_clusters);


        %% Set colormap
        Cmap = [[1 0 0];...%Aggression; red
            [1 0.4 0.1];...%Approach; dark orange
            [0 0 0];...%But sniff; NA
            [0.3 0.7 1];...%Drinking; light blue
            [0 0.7 0];...%Foraging; dark green
            [1 0 1];...%Groom sollicitation; magenta
            [0 1 1];...%Groom partner; cyan
            [0 0 1];...%Getting groomed; dark blue
            [0.8 0 0];...%Threat to partner; dark red
            [1 0 0];...%Threat to subject; red
            [0.9 0.9 0];...%leave; dark yellow
            [0 0 0];...%Lipsmack
            [0.2 0.9 0.76];...%Masturbating; turquoise
            [0.7 0 1];...%Mounting; light purple
            [0.9 0.5 0];...%Other monkeys vocalize; orange
            [1 0.8 0.1];...%Travel; yellow orange
            [0 0 0];...%Proximity; NA
            [0 0 0];...%Rowdy room; NA
            [1 0 0];...%SP; red
            [1 0 0];...%SS; red
            [0.6314 0.5059 0.0118];...%Scratch; maroon
            [0.5 0.2 0.5];...%Self-groom; dark purple
            [ 1 0.07 0.65];...%Submission; dark pink
            [0 0.4 0.5];...%Vocalzation; blue green
            [0 0 0];...%Yawning; NA
            [0.8 0.8 0.8]];%Rest; grey

        Cmap_block = [[0.9 0.7 0.12];[0 0.6 0.8];[0.5 0 0]];

        Cmap_time = copper(length(idx));

        Cmap_kmeans = hsv(n_clusters);

        mean_activity = int32(round(rescale(mean(Spike_count_raster_final,2)),2)*100); mean_activity(mean_activity==0)=1;
        std_activity = int32(round(rescale(std(Spike_count_raster_final,[],2)),2)*100);std_activity(std_activity==0)=1;
        Cmap_firingRate = jet(double(max(mean_activity)));
        Cmap_firingRateStd = jet(double(max(std_activity)));
        %Note: mean activity correlates with variation in activity.


        %% Plot PCA projection in 3D space

        figure; hold on; set(gcf,'Position',[150 250 1500 500])

%         for i=1:n_clusters
%             scatter3(umap_result(idx==i,1), umap_result(idx==i,2),umap_result(idx==i,3),8,Cmap_kmeans(i,:),'filled')
%             xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
%             %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
%             title('Block')
%             set(gca,'FontSize',12);
%             pause(2)
%         end

        %Plot PCA results color-coded by behavior
        ax1=subplot(2,2,1);
        scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),8,Cmap(behavior_labels_final,:),'filled')
%         scatter3(pca_result(:,1), pca_result(:,2),pca_result(:,3),8,Cmap(behavior_labels_final,:),'filled')
        xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior')
        set(gca,'FontSize',12);
        %saveas(gcf,[savePath '/umap_supervised_ColorCodedByBehav_' channel 'Units.png'])
        %pause(5)


        %Color-coded by block
        ax2=subplot(2,2,2);
         scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),8,Cmap_block(block_labels_final,:),'filled')
%         scatter3(pca_result(:,1), pca_result(:,2),pca_result(:,3),8,Cmap_block(block_labels_final,:),'filled')
        xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Block')
        set(gca,'FontSize',12);


        %Color-coded by kmeans cluster in pca space
        ax3=subplot(2,2,3);
        scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),8,Cmap_kmeans(idx_cluster_pca,:),'filled')
%         scatter3(pca_result(:,1), pca_result(:,2),pca_result(:,3),8,Cmap_kmeans(idx_cluster,:),'filled')
        xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Kmeans cluster in pca space')
        set(gca,'FontSize',12);

        %Color-coded by kmeans cluster in umap space
        ax4=subplot(2,2,4);
        scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),8,Cmap_kmeans(idx_cluster_umap,:),'filled')
%         scatter3(pca_result(:,1), pca_result(:,2),pca_result(:,3),8,Cmap_kmeans(idx_cluster,:),'filled')
        xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Kmeans cluster in umap space')
        set(gca,'FontSize',12);

%         %Color-coded by time
%         ax2=subplot(1,2,2);
%         scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),8,Cmap_time,'filled')
%         xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
%         %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
%         title('Time')
%         set(gca,'FontSize',12);

%         %Color-coded by mean activity
%         ax2=subplot(1,2,2);
%         scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),8,Cmap_firingRate(mean_activity,:),'filled')
%         xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
%         %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
%         title('Mean activity')
%         set(gca,'FontSize',12);

%         %Color-coded by std activity
%         ax3=subplot(1,3,3);
%         scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),8,Cmap_firingRateStd(std_activity,:),'filled')
%         xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
%         %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
%         title('Variation in of activity')
%         set(gca,'FontSize',12);
% 
%         %sgtitle([channel ' units, PCA, ' sessions(s).name])

        hlink = linkprop([ax1,ax2,ax3,ax4],{'CameraPosition','CameraUpVector'});
        rotate3d on

        %savefig([savePath 'Umap_3Dprojection_' channel '.fig'])
        %saveas(gcf,[savePath '/umap_ColorCodedByBlock_' channel 'Units.png'])


        % % % %         %plot by cluster
        % % % %         [cls, pnode] = knncluster(umap_result, 200); unique(cls)
        % % % %         Cmap_cls = jet(length(unique(cls)));
        % % % %         figure;scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),8,Cmap_cls(cls,:),'filled')
        % % % %         xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        % % % %         %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        % % % %         title('Behavior')
        % % % %         set(gca,'FontSize',12);
        % % % %
        % % % %         %Find correlation between low-D embedding and high-D embedding
        % % % %         corr(s,chan) = pdist_plot(Spike_count_raster_final);
        % % % %         corr_supervised(s,chan) = pdist_plot(data,'label_column',307);

        %         cd(savePath)
        %         OptionZ.FrameRate=30;OptionZ.Duration=15;OptionZ.Periodic=true;
        %         CaptureFigVid([-20,10;-380,190],['Umap_3Dprojection_' channel],OptionZ)

        %behavior_labels=categorical(behavior_labels);

        behaviors = unique(behavior_labels_final);
        behav_mat = zeros(n_clusters, length(behaviors));
        for cl = 1:n_clusters
            
            cluster_behav = behavior_labels_final(idx_cluster_pca==cl);
            behav_tabulation = tabulate(cluster_behav); 
            behav_mat(cl,behav_tabulation(:,1))= behav_tabulation(:,3);
%             cluster_behav_distr{1,cl}=behav_tabulation{cl}(behav_tabulation{cl}(:,3)~=0,:);
%             [max_behav, idx_max] =  max(cluster_behav_distr{1,cl}(:,3));
%             cluster_behav_distr{2,cl} = behav_categ(cluster_behav_distr{1,cl}(idx_max,1));
%             cluster_behav_distr{3,cl} = max_behav;
% 
            cluster_block = block_labels_final(idx_cluster_pca==cl);
            block_tabulation = tabulate(cluster_block); 
            block_mat(cl,block_tabulation(:,1))= block_tabulation(:,3);

%             cluster_block_distr{1,cl}=behav_tabulation(block_tabulation(:,3)~=0,:);

        end

        behav_mat_final = behav_mat(:,behaviors);

        %Get list of behaviors which actually occur
        figure;
        hm=heatmap(behav_mat_final);
        hm.XDisplayLabels =behav_categ(behaviors);
        colormap jet

        figure;
        hm2=heatmap(block_mat);
        hm2.XDisplayLabels =block_categ;

% % %         %For session 1, with 30 clusters self-groom separates into two
% % %         %clusters. What may be different between the two?
% % %         coi = find(behav_mat(:,find(strcmp(behav_categ,"Self-groom"))) >50);
% % %         times=find(idx_cluster == coi(1));
% % %         time_plot = zeros(1,length(behavior_labels));
% % %         time_plot(times)=1; 
% % % 
% % %         times2=find(idx_cluster == coi(2));
% % %         time_plot2 = zeros(1,length(behavior_labels));
% % %         time_plot2(times2)=1; 
% % % 
% % %         figure; hold on; plot(time_plot); plot(time_plot2); ylim([-0.5 1.5])
% % %         figure; hold on; histogram(times); histogram(times2)
% % %         hrs=floor(median(times)/3600)
% % %         mins=floor((median(times)/3600-hrs)*60)
% % %         secs=((median(times)/3600-hrs)*60-mins)*60
% % %         
% % %         hrs=floor(8380/3600)
% % %         mins=floor((8380/3600-hrs)*60)
% % %         secs=((8380/3600-hrs)*60-mins)*60
        


     end %end of channel for loop

     figure; hold on; set(gcf,'Position',[150 250 1000 800])

        %Plot PCA results color-coded by behavior vlPFC
        ax1=subplot(2,1,1); chan =1;
        scatter3(pca_result(:,1), pca_result(:,2),pca_result(:,3),8,Cmap(behavior_labels_final,:),'filled')
        xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior, vlPFC')
        set(gca,'FontSize',12);

        %Color-coded by block vlPFC
        ax2=subplot(2,1,2);
        scatter3(pca_result(:,1), pca_result(:,2),pca_result(:,3),8,Cmap_block(block_labels_final,:),'filled')
        xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Social context, vlPFC')
        set(gca,'FontSize',12);

        sgtitle([sessions(s).name])
        hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'});
        rotate3d on


        figure; hold on; set(gcf,'Position',[150 250 800 800])
        %Plot PCA results color-coded by behavior TEO
        ax3=subplot(2,1,1); chan =2;
        scatter3(pca_result(:,1), pca_result(:,2),pca_result(:,3),8,Cmap(behavior_labels_final,:),'filled')
        xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior, TEO')
        set(gca,'FontSize',12);

        ax4=subplot(2,1,2);
        scatter3(pca_result(:,1), pca_result(:,2),pca_result(:,3),8,Cmap_block(block_labels_final,:),'filled')
        xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Social context, TEO')
        set(gca,'FontSize',12);

        sgtitle([sessions(s).name])
        hlink = linkprop([ax3,ax4],{'CameraPosition','CameraUpVector'});
        rotate3d on

      
        saveas(gcf,[savePath '/umap_ColorCodedBy_Units.pdf'])
% % % %         savefig([savePath 'Umap_3Dprojection_bothAreas.fig'])

%         cd(savePath)
%         OptionZ.FrameRate=30;OptionZ.Duration=15;OptionZ.Periodic=true;
%         CaptureFigVid([-20,10;-380,190],['Umap_3Dprojection_bothAreas'],OptionZ)


end %end of session for loop

figure
for s = 1:length(a_sessions)
    subplot(length(a_sessions),1,s)
    gscatter(pca_result{a_sessions(s),3}(:,1), pca_result{a_sessions(s),3}(:,2), labels_plot{a_sessions(s)},[],[],10)
    xlabel('PCA 1'); ylabel('PCA 2');
    set(gca,'xtick',[]); set(gca,'ytick',[])
end

for s = 1:length(h_sessions)
    subplot(length(h_sessions),1,s)
    gscatter(pca_result{a_sessions(s),3}(:,1), pca_result{a_sessions(s),3}(:,2), labels_plot{a_sessions(s)},[],[],10)
    xlabel('PCA 1'); ylabel('PCA 2');
    set(gca,'xtick',[]); set(gca,'ytick',[])
end

