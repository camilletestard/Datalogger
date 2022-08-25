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
with_partner =1;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1;%set the smoothing window size (sigma)


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

        Spike_count_raster = Spike_rasters';

        %Extract subejct and partner video
        behavior_labels_subject_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        behavior_labels_partner_init = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
        
        %Extract block labels
        block_labels = cell2mat({labels{:,12}}');
        block_categ = string(block_times{:,1})';

        %% Select behaviors to decode

        %Compute freq of behavior for the session
        behav_freq_table = tabulate(behavior_labels_partner_init);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

        % Select behaviors with a minimum # of occurrences
        behav = behav_freq_table(behav_freq_table(:,2)>=30,1);%Get behaviors with a min number of occurrences
        behav = behav(behav~=find(matches(behav_categ,'Proximity')));%excluding proximity which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.

        % Then select non-reciprocal behaviors
        behav = setdiff(behav, reciprocal_set);

        % OR select behaviors manually
        %behav = [4,5,17,23,25];%manually select behaviors of interest
        %Select behaviors manually to ensure that the same
        %behaviors are considered for the partner and subject comparisons.
        %This list could change from session to session.. I'll have to
        %think of a way to automatize this.

        %Only keep the behaviors of interest
        idx = find(ismember(behavior_labels_partner_init,behav)); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
       
        %check the amount of labels that differ in the partner vs. subject
        %labels after selecting the behaviors of interest.
        subject_behav_after_selection = behavior_labels_subject_init(idx);
        partner_behav_after_selection = behavior_labels_partner_init(idx);
        block_after_selection = block_labels(idx);

        overlap_partner_subject_after_selection = length(find(partner_behav_after_selection == subject_behav_after_selection))/length(idx);
        fprintf('Percent overlap between subject and partner labels AFTER selecting behaviors: %s \n', num2str(overlap_partner_subject_after_selection))
        alone_block_obs = length(find(block_after_selection==3))/length(idx);
        fprintf('Percent observations in "alone" block AFTER selecting behaviors: %s \n', num2str(alone_block_obs))

        %Only consider windows where the behaviors of subject and
        %partner do not overlap
        diff_idx = find(partner_behav_after_selection ~= subject_behav_after_selection); %find the indices where subject and partner behavior do not overlap
        Spike_count_raster_final = Spike_count_raster_final(diff_idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_partner_final = partner_behav_after_selection(diff_idx);%Same as above but in behavior labels
        behavior_labels_subject_final = subject_behav_after_selection(diff_idx);
        block_after_selection_final = block_after_selection(diff_idx);

   
        %% Run umap

        %Supervised
%         data = [Spike_count_raster_final, behavior_labels_final];
%         [umap_result{s,chan}]=run_umap(data, 'n_neighbors', 50, 'min_dist', 0.5, 'n_components', 3,'label_column', 'end'); %Run umap to get 2d embedded states

        %Unsupervised
        [umap_result{s,chan}]=run_umap(Spike_count_raster_final, 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3); %Run umap to get 2d embedded states
        close

        channel = char(channel_flag);

        %Set colormap
        Cmap = [[0 0 0];[1 0.4 0.1];[0 0 0];[0 0.6 0.8];[0 0.7 0];[1 0 1];[0 1 1];...
            [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0 0 0];...
            [0 0 0];[0 0 0];[0.8 0.8 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
            [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];

        Cmap_block = [[1 0 1];[0 0.6 0.8];[0.5 0 0]];

        Cmap_time = copper(length(idx));

        %% Plot UMAP projection in 3D space

        figure; hold on; set(gcf,'Position',[150 250 1500 500]); pointsize = 10;

        %Plot UMAP results color-coded by subject behavior
        ax1=subplot(1,3,1);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),pointsize,Cmap(behavior_labels_subject_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior subject')
        set(gca,'FontSize',12);
        
        %Plot UMAP results color-coded by partner behavior
        ax2=subplot(1,3,2);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),pointsize,Cmap(behavior_labels_partner_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Behavior partner')
        set(gca,'FontSize',12);

%         %Plot UMAP results color-coded by time
%         ax2=subplot(1,3,2);
%         scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),8,Cmap_time,'filled')
%         xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
%         %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
%         title('Time')
%         set(gca,'FontSize',12);

        %Color-coded by block
        ax3=subplot(1,3,3);
        scatter3(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2),umap_result{s,chan}(:,3),pointsize,Cmap_block(block_after_selection_final,:),'filled')
        xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
        %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
        title('Block')
        set(gca,'FontSize',12);
      
        sgtitle([channel ' units, UMAP, ' sessions(s).name])

        hlink = linkprop([ax1,ax2,ax3],{'CameraPosition','CameraUpVector'});
        rotate3d on

        savefig([savePath 'Umap_3Dprojection_' channel '.fig'])

         chan=chan+1;

    end %end of channel for loop

end %end of session for loop

figure
for s = 1:length(a_sessions)
    subplot(length(a_sessions),1,s)
    gscatter(umap_result{a_sessions(s),1}(:,1), umap_result{a_sessions(s),1}(:,2), labels_plot{a_sessions(s)},[],[],5)
    xlabel('UMAP 1'); ylabel('UMAP 2');
    set(gca,'xtick',[]); set(gca,'ytick',[])
end

for s = 1:length(h_sessions)
    subplot(length(h_sessions),1,s)
    gscatter(umap_result{a_sessions(s),1}(:,1), umap_result{a_sessions(s),1}(:,2), labels_plot{a_sessions(s)},[],[],10)
    xlabel('UMAP 1'); ylabel('UMAP 2');
    set(gca,'xtick',[]); set(gca,'ytick',[])
end

