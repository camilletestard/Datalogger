%Set path
is_mac = 0;
if is_mac
    cd('~/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
else
    cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/')
end
filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)

if is_mac
    cd('~/Dropbox (Penn)/Datalogger/Results/')
else
    cd('C:/Users/GENERAL/Dropbox (Penn)/Datalogger/Results/')
end
savePath = uigetdir('', 'Please select the result directory');

clearvars -except savePath filePath temp_resolution channel_flag is_mac

%Set temporal resolution
temp = 1; temp_resolution = 1;
for temp_resolution = [1 2, 5, 10] %1sec, 500msec, 100msec
    %temp_resolution = [1/5, 1/2, 1, 5, 10] %5sec, 2sec, 1sec,500msec, 100msec
    %1 for second resolution, 10 for 100msec resolution, 100 for 10msec resolution, 1000 for msec resolution. etc.
    %0.1 for 10sec resolution, 1/5 for 5sec resolution

    %Set channels: 'TEO', 'vlPFC' or 'all'
    chan = 1; channel_flag = "all";
    for channel_flag = ["vlPFC", "TEO", "all"]

       %Get data with specified temporal resolution and channels
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac);
        %filePath is the experimental data path
        %Temp_resolution is the temporal resolution at which we would like to
        %analyze the dat
        %Channel_flag specifies with channels to include: only TEO array, only
        %vlPFC array or all channels
        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';

        %% Select behaviors to decode
        %Compute freq of behavior for the session
        behavior_labels = cell2mat({labels{:,3}}');
        behavior_labels_context = cell2mat({labels{:,10}}');
        ME_label = cell2mat({labels{:,11}}');
        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=0,:); % Discard 0 (non-defined behaviors)
        block_categ = string(block_times{:,1})';

        % Select behaviors
        behav = [1:6,9:17, 19:27];% [4:10, 23]; %[4:8,17];%[1:6,9:11,16,17]; %manually select behaviors of interest
        behavs_eval = behav_categ(behav);

        idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx);%Same as above but in behavior labels
        behavior_labels_context_final =  behavior_labels_context(idx);
        ME_label_final = ME_label(idx);
        tabulate(behavior_labels_final);

        %% Run umap
        [umap_result]=run_umap(Spike_count_raster_final, 'n_neighbors', 15, 'min_dist', 0.05); %Run umap to get 2d embedded states
        channel = char(channel_flag);
        saveas(gcf,[savePath '/umap_unsupervised_' num2str(1000/temp_resolution) 'msec_' channel '.png'])
        
        %Order for later plotting
        labels_plot = categorical(behav_categ(behavior_labels_final));
        labels_plot_context = categorical(block_categ(behavior_labels_context_final));
        MElabels_plot = categorical(ME_label_final)';
%         if strcmp(monkey, 'Amos')
%             %labels_order = reordercats(labels_plot,{'Foraging','Threat to partner', 'Threat to subject','Groom Receive','Groom Give','Self-groom',});
        labels_order = reordercats(labels_plot,{'Foraging','Self-groom','Threat to partner', 'Threat to subject','Groom Receive','Groom Give'});
        %labels_order = reordercats(labels_plot_context,{'Alone.block','Pair2.block','Pair1.block'});
%         else
%             labels_order = reordercats(labels_plot,{'Foraging','Groom Give','Groom Receive','Threat to partner', 'Threat to subject'});
%         end
        
        %Plot results color coded by behavior
        figure
        gscatter(umap_result(:,1), umap_result(:,2), labels_plot,[],[],10)
        xlabel('UMAP 1'); ylabel('UMAP 2');
        set(gca,'xtick',[]); set(gca,'ytick',[])
        saveas(gcf,[savePath '/umap_ColorCoded_' num2str(1000/temp_resolution) 'msec_' channel '.png'])

        %Color-coded by amount of ME
        figure
        gscatter(umap_result(:,1), umap_result(:,2), MElabels_plot,[],[],10)
        xlabel('UMAP 1'); ylabel('UMAP 2');
        set(gca,'xtick',[]); set(gca,'ytick',[])


        disp('****************************************************************************')
        disp([num2str(1000/temp_resolution) 'msec resolution, channels: ' channel '. DONE'])
        disp('****************************************************************************')

        %pause(2)
        close all

        clearvars -except temp chan savePath filePath temp_resolution is_mac
    end
end
