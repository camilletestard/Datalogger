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
with_partner =2;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0;%Only consider isolated units. 0=all units; 1=only well isolated units
min_occurrences=30;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

s=11;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];

    chan = 1;

    %for channel_flag = ["vlPFC", "TEO", "all"]

        
        %% Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        elseif with_partner ==2
            [Spike_rasters, labels, labels_partner, labels_neighbor, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function_neighbor(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        end

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';
        behavior_labels_subject_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        behavior_labels_partner_init = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner
        behavior_labels_neighbor_init = cell2mat({labels_neighbor{:,3}}');
        block_labels = cell2mat({labels{:,11}}'); %Extract block info
        alone_block_num = find(strcmp(block_times{:,1},"Alone.block"));

        %SANITY CHECK: Compute overlap between partner behavior and subject behavior
        perc_overlap_subject_partner = length(find(behavior_labels_subject_init == behavior_labels_partner_init))/length(behavior_labels_subject_init);
        perc_overlap_subject_neighbor = length(find(behavior_labels_subject_init == behavior_labels_neighbor_init))/length(behavior_labels_subject_init);
        perc_overlap_partner_neighbor = length(find(behavior_labels_partner_init == behavior_labels_neighbor_init))/length(behavior_labels_subject_init);   
    


        behavior_labels = behavior_labels_neighbor_init;

        %% Select behaviors to decode

        %Compute freq of behavior for the session
        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

        % Select behaviors with a minimum # of occurrences
        behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences
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

        %Print behaviors selected
        behavs_eval = behav_categ(behav);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        fprintf('Behaviors evaluated are: %s \n', behavs_eval);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        %Only keep the behaviors of interest
        idx = find(ismember(behavior_labels,behav) & block_labels ==2);%~=alone_block_num); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster(idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels(idx,:);%Same as above but in behavior labels
        tabulate(removecats(categorical(behavior_labels_final)));

        %check the amount of labels that differ in the partner vs. subject
        %labels after selecting the behaviors of interest.
        subject_behav_after_selection = behavior_labels_subject_init(idx);
        partner_behav_after_selection = behavior_labels_partner_init(idx);
        neighbor_behav_after_selection = behavior_labels_neighbor_init(idx);
        block_after_selection = block_labels(idx);


        overlap_neighbor_subject_after_selection = length(find(neighbor_behav_after_selection == subject_behav_after_selection))/length(idx);
        fprintf('Percent overlap between subject and neighbor labels AFTER selecting behaviors: %s \n', num2str(overlap_neighbor_subject_after_selection))
        
        overlap_neighbor_partner_after_selection = length(find(neighbor_behav_after_selection == partner_behav_after_selection))/length(idx);
        fprintf('Percent overlap between partner and neighbor labels AFTER selecting behaviors: %s \n', num2str(overlap_neighbor_partner_after_selection))

        alone_block_obs = length(find(block_after_selection==3))/length(idx);
        fprintf('Percent observations in "alone" block AFTER selecting behaviors: %s \n', num2str(alone_block_obs))

        %Only consider windows where the behaviors of subject and
        %partner do not overlap
        diff_idx = find(neighbor_behav_after_selection ~= subject_behav_after_selection); %find the indices where subject and partner behavior do not overlap
        Spike_count_raster_final = Spike_count_raster_final(diff_idx,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels_final(diff_idx,:);%Same as above but in behavior labels
        block_after_selection_final = block_after_selection(diff_idx);

        %Make sure that only the behaviors with min occurrences are
        %considered
        behav_freq_table2= tabulate(behavior_labels_final);
        behav = behav_freq_table2(behav_freq_table2(:,2)>=min_occurrences,1);
        idxfinal = find(ismember(behavior_labels_final,behav)); %find the indices of the behaviors considered
        Spike_count_raster_final = Spike_count_raster_final(idxfinal,:);%Only keep timepoints where the behaviors of interest occur in spiking data
        behavior_labels_final = behavior_labels_final(idxfinal,:);%Same as above but in behavior labels
        block_after_selection_final = block_after_selection(idxfinal);

        tabulate(removecats(categorical(behavior_labels_final)));
        tabulate(block_after_selection_final)
        crosstab(removecats(categorical(behavior_labels_final)), block_after_selection_final)
        %Note: not all behaviors are equally happening across blocks. It's
        %possible that what we decode is actually block ID and not the
        %behavior itself...

        %Display which behaviors will be decoded
        behavs_eval = behav_categ(behav);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        fprintf('Behaviors evaluated are: %s \n', behavs_eval);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

   
        %% Run umap

        %Supervised
        data = [Spike_count_raster_final, behavior_labels_final];
        [umap_result{s,chan}]=run_umap(data, 'n_neighbors', 50, 'min_dist', 1, 'label_column', 'end'); %Run umap to get 2d embedded states

        %Unsupervised
        %[umap_result{s,chan}]=run_umap(Spike_count_raster_final, 'n_neighbors', 10, 'min_dist', 0.5); %Run umap to get 2d embedded states
        
        channel = char(channel_flag);
        %saveas(gcf,[savePath '/umap_unsupervised_' num2str(1000/temp_resolution) 'msec_' channel '.png'])
        
        %Order for later plotting
        labels_plot{s} = categorical(behav_categ(behavior_labels_final));
        ordered_labels = string(unique(labels_plot{s}));
        labels_order = reordercats(labels_plot{s},ordered_labels);
  
        %Plot results color-coded by behavior
        figure
        gscatter(umap_result{s,chan}(:,1), umap_result{s,chan}(:,2), labels_plot{s} ,[],[],10, 'MarkerFaceAlpha',0.7)
        legend(ordered_labels,'Location','best')
        xlabel('UMAP 1'); ylabel('UMAP 2');
        set(gca,'xtick',[]); set(gca,'ytick',[])
        title('Neighbor behavior')
        saveas(gcf,[savePath '/umap_supervised_ColorCodedNeighborBehv_' channel 'Units.png'])
        
        close all

        %clearvars -except temp chan savePath filePath temp_resolution is_mac
        chan=chan+1;

    %end %end of channel for loop

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

