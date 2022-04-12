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
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units

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

    for channel_flag = ["vlPFC", "TEO", "all"]

        
        %% Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        end

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';

        %% Select behaviors to decode
        %Compute freq of behavior for the session
        behavior_labels = cell2mat({labels{:,3}}');
        block_labels = cell2mat({labels{:,11}}');

        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=0,:); % Discard 0 (non-defined behaviors)
        block_categ = string(block_times{:,1})';

        % Select behaviors
        behav = [5,7,8,9,10,24];% [4:10, 23]; %[4:8,17];%[1:6,9:11,16,17]; %manually select behaviors of interest
        behavs_eval = behav_categ(behav);

        idx = find(ismember(behavior_labels,behav)); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx);%Same as above but in behavior labels
        block_labels_final =  block_labels(idx);
        %ME_label_final = ME_label(idx);
        tabulate(behavior_labels_final);

        %% Run umap
        [umap_result{s,chan}]=run_umap(Spike_count_raster_final, 'n_neighbors', 15, 'min_dist', 0.05); %Run umap to get 2d embedded states
        channel = char(channel_flag);
        saveas(gcf,[savePath '/umap_unsupervised_' num2str(1000/temp_resolution) 'msec_' channel '.png'])
        
        %Order for later plotting
        labels_plot{s} = categorical(behav_categ(behavior_labels_final));
        labels_plot_context{s} = categorical(block_categ(block_labels_final));
        %MElabels_plot = categorical(ME_label_final)';
%         if strcmp(monkey, 'Amos')
%         labels_order = reordercats(labels_plot{s},{'Foraging','Threat to partner', 'Threat to subject','Groom Receive','Groom Give','Self-groom'});
%         labels_order = reordercats(labels_plot{s},{'Foraging','Self-groom','Threat to partner', 'Threat to subject','Groom Receive','Groom Give'});
%         labels_order = reordercats(labels_plot_context{s},{'Alone.block','Pair2.block','Pair1.block'});
%         else
%             labels_order = reordercats(labels_plot,{'Foraging','Groom Give','Groom Receive','Threat to partner', 'Threat to subject'});
%         end
        
        %Plot results color coded by behavior
        figure
        gscatter(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2), labels_plot{s},[],[],10)
        xlabel('UMAP 1'); ylabel('UMAP 2');
        set(gca,'xtick',[]); set(gca,'ytick',[])
        saveas(gcf,[savePath '/umap_ColorCodedByBehav_' num2str(1000/temp_resolution) 'msec_' channel '.png'])

%         %Color-coded by amount of ME
%         figure
%         gscatter(umap_result(:,1), umap_result(:,2), MElabels_plot,[],[],10)
%         xlabel('UMAP 1'); ylabel('UMAP 2');
%         set(gca,'xtick',[]); set(gca,'ytick',[])

        %Color-coded by block
        figure
        gscatter(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2), labels_plot_context{s},[],[],10)
        xlabel('UMAP 1'); ylabel('UMAP 2');
        set(gca,'xtick',[]); set(gca,'ytick',[])
        saveas(gcf,[savePath '/umap_ColorCodedByBlock_' num2str(1000/temp_resolution) 'msec_' channel '.png'])

        disp('****************************************************************************')
        disp([num2str(1000/temp_resolution) 'msec resolution, channels: ' channel '. DONE'])
        disp('****************************************************************************')

        pause(2)
        close all

        %clearvars -except temp chan savePath filePath temp_resolution is_mac
        chan=chan+1;
    end %end of channel for loop

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

