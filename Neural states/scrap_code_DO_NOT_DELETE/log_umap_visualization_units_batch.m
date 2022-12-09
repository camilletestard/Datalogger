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
smooth= 0; % 1: smooth the data; 0: do not smooth
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
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= ...
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
        
        %Unsupervised
        [umap_result{s,chan}]=run_umap(Spike_rasters, 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3); %Run umap to get 2d embedded states
        close

        channel = char(channel_flag);

        %Set colormap
        for n=1:length(brain_label)
            if brain_label(n)=="vlPFC"
                brain_label_num(n)=1;
            else
                brain_label_num(n)=2;
            end
        end

        Cmap=[1 0 0; 0 0 1];
        %% Plot UMAP projection in 3D space

        figure; hold on;% set(gcf,'Position',[150 250 1500 500])

        %Plot UMAP results color-coded by behavior
        figure; scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),15,Cmap(brain_label_num,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')

        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior')
        set(gca,'FontSize',12);
        %saveas(gcf,[savePath '/umap_supervised_ColorCodedByBehav_' channel 'Units.png'])
        %pause(5)

        %Plot UMAP results color-coded by time
        ax2=subplot(1,3,2);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_time,'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Time')
        set(gca,'FontSize',12);

        %Color-coded by block
        ax3=subplot(1,3,3);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_block(block_labels_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Block')
        set(gca,'FontSize',12);
      
        sgtitle([channel ' units, UMAP, ' sessions(s).name])

        hlink = linkprop([ax1,ax2,ax3],{'CameraPosition','CameraUpVector'});
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

