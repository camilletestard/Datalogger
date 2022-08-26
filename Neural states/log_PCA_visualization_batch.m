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
sigma = 5;%set the smoothing window size (sigma)

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
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];

    chan = 1;

    for channel_flag = ["vlPFC", "TEO"]
        %channel_flag = "vlPFC";

        %% Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        end

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';


        %PCA
        [coeff,score,latent,tsquared,explained] = pca(Spike_count_raster);
        var_explained = cumsum(explained); idxl = find(var_explained>=80); dim = min(idxl);
        Spike_count_raster = score(:,1:15);

        %         %FA
        %         [lambda,psi,T,stats,F] = factoran(Spike_count_raster,3)
        %         Spike_count_raster = F;


        %% Select behaviors to decode
        %Compute freq of behavior for the session
        behavior_labels = cell2mat({labels{:,3}}');
        behavior_labels(behavior_labels==21)=9;
        behavior_labels(behavior_labels==22)=10;
        behavior_labels(behavior_labels==19)=29;
        block_labels = cell2mat({labels{:,11}}');

        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=0,:); % Discard 0 (non-defined behaviors)
        block_categ = string(block_times{:,1})';

        % Select behaviors

        %Compute freq of behavior for the session
        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

        % Select behaviors with a minimum # of occurrences
        min_occurrences = 30;
        behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
        behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rowdy Room')));%excluding Rowdy Room which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Other monkeys vocalize')));

        % OR select behaviors manually
        %behav = unique(behavior_labels); %[4,5,7,8,9,10,24];% [4:10, 23]; %[4:8,17]; %manually select behaviors of interest

        %Print behaviors selected
        behavs_eval = behav_categ(behav);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        fprintf('Behaviors evaluated are: %s \n', behavs_eval);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        %Subset data to only consider epochs of interest
        idx= find(ismember(behavior_labels,behav));

        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx);%Same as above but in behavior labels
        block_labels_final =  block_labels(idx);
        behavior_labels_final_rand = randsample(behavior_labels_final, length(behavior_labels_final));

        %% Plot PCA

        %Set colormap
         Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0 0.4 0.5];[0 0.7 0];[1 0 1];[0 1 1];...
            [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0 0 0];...
            [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
            [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];

        Cmap_block = [[0.9 0.7 0.12];[0 0.6 0.8];[0.5 0 0]];

        Cmap_time = copper(length(idx));
        

        figure; hold on; set(gcf,'Position',[150 250 1500 500])

        %Plot PCA results color-coded by behavior
        ax1=subplot(1,3,1);
        scatter3(Spike_count_raster_final(:,1), Spike_count_raster_final(:,2),Spike_count_raster_final(:,3),8,Cmap(behavior_labels_final,:),'filled')
        sum(explained(1:3))
        xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior')
        set(gca,'FontSize',12);
        %saveas(gcf,[savePath '/umap_supervised_ColorCodedByBehav_' channel 'Units.png'])
        %pause(5)

        %Plot UMAP results color-coded by time
        ax2=subplot(1,3,2);
        scatter3(Spike_count_raster_final(:,1), Spike_count_raster_final(:,2),Spike_count_raster_final(:,3),8,Cmap_time,'filled')
        xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Time')
        set(gca,'FontSize',12);

        %Color-coded by block
        ax3=subplot(1,3,3);
        scatter3(Spike_count_raster_final(:,1), Spike_count_raster_final(:,2),Spike_count_raster_final(:,3),8,Cmap_block(block_labels_final,:),'filled')
        xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Block')
        set(gca,'FontSize',12);
        %saveas(gcf,[savePath '/umap_ColorCodedByBlock_' channel 'Units.png'])

        sgtitle([channel_flag ' units, PCA'])

        hlink = linkprop([ax1,ax2,ax3],{'CameraPosition','CameraUpVector'});
        rotate3d on

        %clearvars -except temp chan savePath filePath temp_resolution is_mac
        chan=chan+1;

    end %end of channel for loop

end %end of session for loop


