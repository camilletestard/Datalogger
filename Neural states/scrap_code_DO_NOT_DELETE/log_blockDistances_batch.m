%% log_blockDistanc_batch
% This script applied unsupervised umap on smoothed firing
% Ref: Connor Meehan, Jonathan Ebrahimian, Wayne Moore, and Stephen Meehan (2022). Uniform Manifold Approximation and Projection (UMAP) (https://www.mathworks.com/matlabcentral/fileexchange/71902), MATLAB Central File Exchange.
% Jan 2023, C. Testard

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
agg_precedence=0; % 1: aggression takes precedence; 0: Threat to partner and subject states take precedence

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=2; 
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/UMAP_results/'];


    chan=1;
    for channel_flag = ["vlPFC", "TEO", "all"]
        %channel_flag = "TEO";

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



        %% Select behaviors to visualize

        %Extract behavior labels and frequency
        behavior_labels = cell2mat({labels{:,3}}');
        behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

        %Lump all travel together
        behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
        behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

        %Extract block labels
        block_labels = cell2mat({labels{:,12}}');
        alone_label = cell2mat({labels{:,13}}');

        % Select behaviors

        %Compute freq of behavior for the session
        behav_freq_table = tabulate(behavior_labels);
        behav_freq_table = behav_freq_table(behav_freq_table(:,1)~=length(behav_categ),:); % Discard 0 (non-defined behaviors)

        % Select behaviors with a minimum # of occurrences
        min_occurrences = 30;
        behav = behav_freq_table(behav_freq_table(:,2)>=min_occurrences,1);%Get behaviors with a min number of occurrences

        %Remove behaviors we're not interested in for now
        behav = behav(behav~=find(matches(behav_categ,'Scratch')));%excluding scratch which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rowdy Room')));%excluding Rowdy Room which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Rest')));%excluding rest which is a source of confusion.
        behav = behav(behav~=find(matches(behav_categ,'Other monkeys vocalize')));

        % OR select behaviors manually
         behav =[5,7,8,9,10,16] ;%unique(behavior_labels); %[4,5,7,8,9,10,24];% [4:10, 23]; %[4:8,17]; %manually select behaviors of interest

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
        alone_labels_final = alone_label(idx);


        %% Run umap

        %Supervised
        %         data = [Spike_count_raster_final, behavior_labels_final];
        %         [umap_result{s,chan}]=run_umap(data, 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3,'label_column', 'end'); %Run umap to get 2d embedded states

        %Unsupervised
        [umap_result]=run_umap(Spike_count_raster_final, 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3); %Run umap to get 2d embedded states
        close

        %% Run PCA
        [~, pca_results] = pca(zscore(Spike_count_raster_final));



        %% Plot UMAP projection in 3D space

        %Set colormap
        Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0.1 0.8 0.9];[0 0.7 0];[1 0 1];[0 1 1];...
            [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0.7 0 1];...
            [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0.8 0 0];[1 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
            [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];

        Cmap_block = [[0.9 0.7 0.12];[0 0.6 0.8];[0.5 0 0]];

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
        
                %sgtitle([channel ' units, UMAP, ' sessions(s).name])
        
                hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'});
                rotate3d on


        

        %% Visualize distances between social context given a particular behavior

        boi =[5,7,8,9,10,16];
        paired_presence = crosstab(behavior_labels_final, block_labels_final); paired_presence=paired_presence(:,1:2);
        paired_alone_presence = crosstab(behavior_labels_final, alone_labels_final);
        distances{s,chan}=nan(length(behav),3);

        for beh= 1:length(boi) %for each behavior

            distances{s,chan}(beh,1)=boi(beh);

            if all(paired_alone_presence(beh,:)>15) %if the behavior occurs both when monkey is paired and alone
            %Measure the distance between the centres of mass neural state for this behavior in paired and alone blocks

            paired_activity_beh_com = mean(umap_result(behavior_labels_final==boi(beh)...
                & alone_labels_final==1,:));

            alone_activity_beh_com = mean(umap_result(behavior_labels_final==boi(beh)...
                & alone_labels_final==0,:));

            distances{s,chan}(beh,2)=pdist([paired_activity_beh_com;alone_activity_beh_com], 'cityblock');

            end

            if all(paired_presence(beh,:)>15) %if the behavior occurs both when monkey is paired and alone
            %Measure the distance between the centres of mass neural state for this behavior in paired and alone blocks

            paired_activity_beh_com = mean(umap_result(behavior_labels_final==boi(beh)...
                & block_labels_final==1,:));

            alone_activity_beh_com = mean(umap_result(behavior_labels_final==boi(beh)...
                & block_labels_final==2,:));

            distances{s,chan}(beh,3)=pdist([paired_activity_beh_com;alone_activity_beh_com], 'cityblock');

            end

        end

        chan=chan+1;

    end %end of channel loop

    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    disp(s)
    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

end %end of session for loop

%Save data
cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/UMAP_results/'])
%save("blockDist_behav.mat",'distances','behav_categ','boi')
load("blockDist_behav.mat")

i=2; %set channel
zscored_distances = distances;
for s=a_sessions
    x=distances{s,i}(:,2:3);
    zscored_distances{s,i}(:,2:3)= (x - nanmean(reshape(x,1,[])))/nanstd(reshape(x,1,[]));
end

rescaled_distances = distances;
for s=a_sessions
    x=distances{s,i}(:,2:3);
    rescaled_distances{s,i}(:,2:3)= rescale(x,0,1);
end

distances_allSession=cat(3,rescaled_distances{a_sessions,i});
distances_allSession_mean = nanmean(distances_allSession,3);
figure; hold on; bar(distances_allSession_mean(:,2:3))
scatter(0.85*ones(1,12),squeeze(distances_allSession(1,2,:)),'b','filled')
scatter(1.15*ones(1,12),squeeze(distances_allSession(1,3,:)),'r','filled')
scatter(1.85*ones(1,12),squeeze(distances_allSession(2,2,:)),'b','filled')
scatter(2.15*ones(1,12),squeeze(distances_allSession(2,3,:)),'r','filled')
scatter(2.85*ones(1,12),squeeze(distances_allSession(3,2,:)),'b','filled')
scatter(3.15*ones(1,12),squeeze(distances_allSession(3,3,:)),'r','filled')
scatter(3.85*ones(1,12),squeeze(distances_allSession(4,2,:)),'b','filled')
scatter(4.15*ones(1,12),squeeze(distances_allSession(4,3,:)),'r','filled')
scatter(4.85*ones(1,12),squeeze(distances_allSession(5,2,:)),'b','filled')
scatter(5.15*ones(1,12),squeeze(distances_allSession(5,3,:)),'r','filled')
scatter(5.85*ones(1,12),squeeze(distances_allSession(6,2,:)),'b','filled')
scatter(6.15*ones(1,12),squeeze(distances_allSession(6,3,:)),'r','filled')
xticks([1:6])
xticklabels(behav_categ(boi))
ylabel('Scaled eucledian distance in umap')
legend({'Alone vs. Paired', 'Female vs. Male neighbor'})
set(gca,'fontsize', 16)

%Need to show raw datapoints.
