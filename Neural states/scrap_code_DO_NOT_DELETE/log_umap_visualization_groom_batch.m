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
session_range_with_partner=[1:3,11:13];

%Set parameters
with_partner =0;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)
null=0;%Set whether we want the null 
simplify=0; %lump similar behavioral categories together

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=1;
for s =session_range(2:end) %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];

    chan = 1;

    for channel_flag = ["vlPFC", "TEO"]
        %channel_flag = "vlPFC";


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
                is_mac, with_NC, isolatedOnly, smooth, sigma);
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
%         behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Threat to partner");
%         behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Threat to subject");
%         behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

        if simplify
            %Simplify behavioral catagories
            %Lump all aggressive interactions together
            behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
            behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");
            behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Aggression");
            behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Aggression");

            %Lump all travel together
            behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
            behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");
        end

        if null
            %Simulate fake labels
            [sim_behav] = GenSimBehavior(behavior_labels,behav_categ, temp_resolution);
            behavior_labels = sim_behav;
        end
        
        %Extract block labels
        block_labels = cell2mat({labels{:,12}}');
        block_categ = string(block_times{:,1})';

        % Select behaviors

        %Compute freq of behavior for the session
        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

        % Select grooming manually
        behav =[7,8] ;%unique(behavior_labels); %[4,5,7,8,9,10,24];% [4:10, 23]; %[4:8,17]; %manually select behaviors of interest

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

        %Choose grooming category to represent
        groom_categ= {'Star.vs.end', 'Post-threat','Reciprocated','Initiated'}; %label grooming categ
        sum(groom_labels_all==2)
        categ =3;
        groom_labels_final = groom_labels_all(idx,categ+1);
        groom_categ_label = {'Not Reciprocated', 'Reciprocated'};

        %% Run umap

        %Supervised
%         data = [Spike_count_raster_final, behavior_labels_final];
%         [umap_result{s,chan}]=run_umap(data, 'n_neighbors', 50, 'min_dist', 0.5, 'n_components', 3,'label_column', 'end'); %Run umap to get 2d embedded states

        %Unsupervised
        [umap_result{s,chan}]=run_umap(Spike_count_raster_final, 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3); %Run umap to get 2d embedded states
        close

        channel = char(channel_flag);

        %Set colormap
        Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0.1 0.8 0.9];[0 0.7 0];[1 0 1];[0 1 1];...
            [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0.7 0 1];...
            [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0.8 0 0];[1 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
            [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];

        Cmap_block = [[0.9 0.7 0.12];[0 0.6 0.8];[0.5 0 0]];

        %Cmap_time = copper(size(Spike_count_raster,1));
        Cmap_time = copper(length(idx));

        Cmap_groom = [[0 1 0];[1 0 0]];

        %% Plot UMAP projection in 3D space

        figure; hold on; set(gcf,'Position',[150 250 1500 500])

        %Plot UMAP results color-coded by behavior
        ax1=subplot(1,4,1);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap(behavior_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior')
        set(gca,'FontSize',12);
        %saveas(gcf,[savePath '/umap_supervised_ColorCodedByBehav_' channel 'Units.png'])
        %pause(5)

        %Plot UMAP results color-coded by time
        ax2=subplot(1,4,2);
        %scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_time(idx,:),'filled')
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_time,'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Time')
        set(gca,'FontSize',12);

        %Color-coded by groom-context
        ax3=subplot(1,4,3);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_groom(groom_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Block')
        set(gca,'FontSize',12);

        %Color-coded by block
        ax4=subplot(1,4,4);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_block(block_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Block')
        set(gca,'FontSize',12);
      
        sgtitle([channel ' units, UMAP, ' sessions(s).name])

        hlink = linkprop([ax1,ax2,ax3,ax4],{'CameraPosition','CameraUpVector'});
        rotate3d on

        savefig([savePath 'Umap_3Dprojection_' channel '.fig'])
        %saveas(gcf,[savePath '/umap_ColorCodedByBlock_' channel 'Units.png'])

%         cd(savePath)
%         OptionZ.FrameRate=30;OptionZ.Duration=15;OptionZ.Periodic=true;
%         CaptureFigVid([-20,10;-380,190],['Umap_3Dprojection_' channel],OptionZ)
        
        chan=chan+1;

    end %end of channel for loop

     figure; hold on; set(gcf,'Position',[150 250 1000 800])

        %Plot UMAP results color-coded by behavior vlPFC
        ax1=subplot(2,2,1); chan =1;
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap(behavior_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior, vlPFC')
        set(gca,'FontSize',12);

        %Color-coded by block vlPFC
        ax2=subplot(2,2,2);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_block(block_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Social context, vlPFC')
        set(gca,'FontSize',12);

        %Plot UMAP results color-coded by behavior TEO
        ax3=subplot(2,2,3); chan =2;
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap(behavior_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior, TEO')
        set(gca,'FontSize',12);

        ax4=subplot(2,2,4);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_block(block_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Social context, TEO')
        set(gca,'FontSize',12);

        sgtitle([sessions(s).name])
        hlink = linkprop([ax1,ax2,ax3,ax4],{'CameraPosition','CameraUpVector'});
        rotate3d on

        savefig([savePath 'Umap_3Dprojection_bothAreas.fig'])

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

