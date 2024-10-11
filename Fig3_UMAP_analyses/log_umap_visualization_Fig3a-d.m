%% log_umap_visualization_Fig3a-d
% This script applies unsupervised umap on smoothed firing rate 
% to visualize neural states in a dimensionally reduced space.
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
session_range=[1:6,11:13,15:16,18];
a_sessions = 1:6; h_sessions = [11:13,15:16,18];

%Set parameters
with_partner =0;
temp_resolution = 1;
channel_flag = "TEO";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_MU =1;%0: MU cluster is excluded; 1:MU cluster is included; 2:ONLY multi-unit cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null 
simplify=1; %lump similar behavioral categories together
threat_precedence=0; % 1: aggression takes precedence; 0: Threat to partner and subject states take precedence
exclude_sq=1;

s=2;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];

    chan = 1;

    for channel_flag = ["vlPFC", "TEO"]
        %channel_flag = "TEO";


        %% Get data with specified temporal resolution and channels
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_MU, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

        disp('Data Loaded')

        %Raw data
        Spike_count_raster = zscore(Spike_rasters');


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

        %% Run umap

        %Supervised
%         data = [Spike_count_raster_final, behavior_labels_final];
%         [umap_result{s,chan}]=run_umap(data, 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3,'label_column', 'end'); %Run umap to get 2d embedded states

        %Unsupervised
        [umap_result{s,chan}]=run_umap(Spike_count_raster_final, 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3); %Run umap to get 2d embedded states
        close

        channel = char(channel_flag);

        %Set colormap
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
            [1 0 0];...%SP; NA
            [1 0 0];...%SS; NA
            [0.6314 0.5059 0.0118];...%Scratch; maroon
            [0.5 0.2 0.5];...%Self-groom; dark purple
            [ 1 0.07 0.65];...%Submission; dark pink
            [0 0.4 0.5];...%Vocalzation; blue green
            [0 0 0];...%Yawning; NA
            [0.8 0.8 0.8]];%Rest; grey

        Cmap_block = [[0.9 0.7 0.12];[0 0.6 0.8];[0.5 0 0]];

        Cmap_time = copper(length(idx));

        mean_activity = int32(round(rescale(mean(Spike_count_raster_final,2)),2)*100); mean_activity(mean_activity==0)=1;
        std_activity = int32(round(rescale(std(Spike_count_raster_final,[],2)),2)*100);std_activity(std_activity==0)=1;
        Cmap_firingRate = jet(double(max(mean_activity)));
        Cmap_firingRateStd = jet(double(max(std_activity)));
        %Note: mean activity correlates with variation in activity.


        %% Plot UMAP projection in 3D space

        figure; hold on; set(gcf,'Position',[150 250 1000 500])

        %Plot UMAP results color-coded by behavior
        ax1=subplot(1,2,1);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap(behavior_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior')
        set(gca,'FontSize',12);
        %saveas(gcf,[savePath '/umap_supervised_ColorCodedByBehav_' channel 'Units.png'])
        %pause(5)


        %Color-coded by block
        ax2=subplot(1,2,2);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_block(block_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Block')
        set(gca,'FontSize',12);


%         %Color-coded by time
%         ax3=subplot(1,3,3);
%         scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_time,'filled')
%         xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
%         %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
%         title('Block')
%         set(gca,'FontSize',12);

%         %Color-coded by mean activity
%         ax2=subplot(1,2,2);
%         scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_firingRate(mean_activity,:),'filled')
%         xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
%         %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
%         title('Mean activity')
%         set(gca,'FontSize',12);

%         %Color-coded by std activity
%         ax3=subplot(1,3,3);
%         scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_firingRateStd(std_activity,:),'filled')
%         xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
%         %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
%         title('Variation in of activity')
%         set(gca,'FontSize',12);
% 
%         %sgtitle([channel ' units, UMAP, ' sessions(s).name])

        hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'});
        rotate3d on
        grid on

        %savefig([savePath 'Umap_3Dprojection_' channel '.fig'])
        %saveas(gcf,[savePath '/umap_ColorCodedByBlock_' channel 'Units.png'])


        % % % %         %plot by cluster
        % % % %         [cls, pnode] = knncluster(umap_result{s,chan}, 200); unique(cls)
        % % % %         Cmap_cls = jet(length(unique(cls)));
        % % % %         figure;scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_cls(cls,:),'filled')
        % % % %         xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
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

        chan=chan+1;

     end %end of channel for loop

     figure; hold on; set(gcf,'Position',[150 250 1000 800])

        %Plot UMAP results color-coded by behavior vlPFC
        ax1=subplot(2,1,1); chan =1;
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap(behavior_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior, vlPFC')
        set(gca,'FontSize',12);

        %Color-coded by block vlPFC
        ax2=subplot(2,1,2);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_block(block_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Social context, vlPFC')
        set(gca,'FontSize',12);

        sgtitle([sessions(s).name])
        hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'});
        rotate3d on


        figure; hold on; set(gcf,'Position',[150 250 800 800])
        %Plot UMAP results color-coded by behavior TEO
        ax3=subplot(2,1,1); chan =2;
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap(behavior_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior, TEO')
        set(gca,'FontSize',12);

        ax4=subplot(2,1,2);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_block(block_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
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
    gscatter(umap_result{a_sessions(s),3}(:,1), umap_result{a_sessions(s),3}(:,2), labels_plot{a_sessions(s)},[],[],10)
    xlabel('UMAP 1'); ylabel('UMAP 2');
    set(gca,'xtick',[]); set(gca,'ytick',[])
end

for s = 1:length(h_sessions)
    subplot(length(h_sessions),1,s)
    gscatter(umap_result{a_sessions(s),3}(:,1), umap_result{a_sessions(s),3}(:,2), labels_plot{a_sessions(s)},[],[],10)
    xlabel('UMAP 1'); ylabel('UMAP 2');
    set(gca,'xtick',[]); set(gca,'ytick',[])
end

