%% Log_grooming_behavior.m

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
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
randomsample=0;
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 50;%Number of SVM iterations
smooth= 0; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=0;%lump similar behavioral categories together to increase sample size.
threat_precedence =1;
exclude_sq=1;

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
Cmap_block = [[0.9 0.7 0.12];[0 0.6 0.8];[0.5 0 0];[0.8 0.8 0.8]];

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
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];


    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')


    Spike_count_raster{s} = Spike_rasters';
    session_length(s) = size(Spike_count_raster{s},1);
    behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ);
    block_labels = cell2mat({labels{:,12}}');


    %% Extract rest bouts
    RestBouts = zeros(size(behavior_labels));
    RestBouts(behavior_labels== length(behav_categ))=1;

    rest_bout_start = [1; find(diff(RestBouts)==1)+1];
    rest_bout_end = find(diff(RestBouts)==-1);

    if length(rest_bout_end)<length(rest_bout_start) %can happen if grooming went until very end of session
        rest_bout_end(length(rest_bout_start))=length(RestBouts);
    end
    rest_duration = rest_bout_end-rest_bout_start;
    rest_bout=[rest_bout_start, rest_bout_end, rest_duration];
    rest_bout_final = rest_bout(rest_bout(:,3)>9,:);


    %% Run umap
    neural_data = zscore(Spike_count_raster{s});
    time = 1:length(behavior_labels);

    bouts_to_consider = 1:size(rest_bout_final,1);
    rest_idx=[]; timescale_all=[]; boutid_all=[];
    for b=1:length(bouts_to_consider)
        idx = rest_bout_final(bouts_to_consider(b),1):rest_bout_final(bouts_to_consider(b),2);
        timescale = rescale(idx);
        bout_id = ones(size(idx))*b;

        rest_idx = [rest_idx, idx];
        timescale_all = [timescale_all, timescale];
        boutid_all = [boutid_all, bout_id];

    end

    [umap_result]=run_umap(neural_data(rest_idx,:), 'n_neighbors', 15, 'min_dist', 0.1, 'n_components', 3); %Run umap to get 2d embedded states
    close

    figure;
    ax1= subplot(1,2,1);
    scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),10,Cmap(behavior_labels(rest_idx),:),'filled')
    xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
    colorbar
    %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
    set(gca,'FontSize',12);
    title('Behavior')

%     figure;
%     scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),10,Cmap_block(block_labels(rest_idx),:),'filled')
%     xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
%     %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
%     set(gca,'FontSize',12);
%     title('Neighbor ID')

    ax2= subplot(1,2,2);
    scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),10,boutid_all,'filled')%Cmap_recip(recip,:),'filled')
    xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
    colormap(jet)
    colorbar
    %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
    set(gca,'FontSize',12);
    title('Bout Number')

%     figure
%     scatter3(umap_result(:,1), umap_result(:,2),umap_result(:,3),10,time(rest_idx),'filled')%Cmap_recip(recip,:),'filled')
%     xlabel('UMAP 1'); ylabel('UMAP 2'); zlabel('UMAP 3')
%     colormap(hot)
%     %colorbar
%     %set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[])
%     set(gca,'FontSize',12);
%     title('Time in session')

    hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'});
        rotate3d on

    %% Decoding

    neural_data = zscore(Spike_count_raster{s});
    behav=length(behav_categ);
    rest_idx = find(ismember(behavior_labels,behav));

    for categ = 1:2
        for subsample_iter = 1:10
            bouts_to_consider = randsample(size(rest_bout_final,1),16);
            idx_all=[]; timescale_all=[]; boutid_all=[];

            for b=1:length(bouts_to_consider)
                idx = rest_bout_final(bouts_to_consider(b),1):rest_bout_final(bouts_to_consider(b),2);
                timescale = round(rescale(idx)*10);
                bout_id = ones(size(idx)).*b;

                idx_all = [idx_all, idx];
                timescale_all = [timescale_all, timescale];
                boutid_all = [boutid_all, bout_id];
            end
            rest_categ = [timescale_all;boutid_all];
            neural_data_final=neural_data(idx_all,:);
            behavior_labels_final = rest_categ(categ,:);

            disp('Start running SVM...')
            for iter = 1:num_iter

                %subsample to match number of neurons across brain areas
                Labels = behavior_labels_final;
                Input_matrix = neural_data_final;

                %Balance number of trials per class
                uniqueLabels = unique(Labels); %IDentify unique labels (useful when not numbers)
                NumOfClasses = length(uniqueLabels); % Total number of classes
                numericLabels = 1:NumOfClasses; %Numeric name of labels

                labels_temp = Labels;
                for i=1:NumOfClasses
                    idx = Labels == uniqueLabels(i);
                    labels_temp(idx) = numericLabels(i);
                end
                Labels = labels_temp;

                num_trials = hist(Labels,numericLabels); %number of trials in each class
                minNumTrials = min(num_trials); %use the minimum # of instances

                chosen_trials = [];
                for i = 1:NumOfClasses %for each class
                    idx = find(Labels == numericLabels(i)); %find indexes of trials belonging to this class
                    rand_i = randsample(length(idx), minNumTrials); %Select a random n number of them
                    chosen_trials = [chosen_trials, idx(rand_i)]; %Put the selected trials in a matrix, ordered by class
                end
                Input_matrix = Input_matrix(chosen_trials, :);
                Labels = Labels(chosen_trials);
                Labels_shuffled = Labels(randperm(length(Labels)));

                % Run svm
                [hitrate(iter), C{iter}] = log_SVM_basic_function(Input_matrix, Labels', 5, 0, 0);
                [hitrate_shuffled(iter), C_shuffled{iter}] = log_SVM_basic_function(Input_matrix, Labels_shuffled', 5, 0, 0);

                if mod(iter,10)==1
                    disp(['SVM run' num2str(iter) '/' num2str(num_iter)])
                end
            end %end of SVM iterations

            mean_hitrate(categ, subsample_iter) = mean(hitrate);
            mean_hitrate_shuffled(categ, subsample_iter) = mean(hitrate_shuffled);
        end

    end

end
