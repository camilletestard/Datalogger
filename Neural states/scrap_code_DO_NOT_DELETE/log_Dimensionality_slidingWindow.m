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
channel_flag = "vlPFC";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1;%set the smoothing window size (sigma)
var_explained_threshold=90; num_units = 100;
num_iter = 50;
window_size = 50; step = 5;
simplify =1;

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
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey,...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, ...
                with_NC, isolatedOnly, smooth, sigma);
        end

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';
        behavior_labels = cell2mat({labels{:,3}}');

        behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Threat to partner");
        behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Threat to subject");
        behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

        if simplify == 1 %Compare across 5 behavioral categories that are common across sessions
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


        %% Compute dimensionality for 100 units, over multiple iterations,over the whole session

        timepoints = window_size/2+1:step:size(Spike_count_raster,1)-window_size/2;

        for tp=1:length(timepoints)

            dim = nan(1,num_iter);
            for iter = 1:num_iter

                %Select unit to run SVM
                Input_matrix = Spike_count_raster(timepoints(tp)-window_size/2:timepoints(tp)+window_size/2-1,randsample(size(Spike_count_raster,2), num_units));

                %PCA
                [coeff,score,~,~,explained] = pca(Input_matrix);

                %Get dimensionality
                var_explained = cumsum(explained);
                idxl = find(var_explained>=var_explained_threshold);
                dim(iter) = min(idxl);
                

            end % end of interation loop
   
            disp([num2str(tp) '/' num2str(length(timepoints))])
            mean_dim(chan, tp) = mean(dim);
            std_dim(chan, tp) = std(dim);
        end

        chan=chan+1;
        disp([channel_flag ' done'])
    end %end of channel loop

    %% Plot results for the session
    
    vlpfc_dim = interp(mean_dim(1,:), step); vlpfc_dim = [ones(1,size(Spike_count_raster,1)-length(vlpfc_dim))*vlpfc_dim(1), vlpfc_dim];
    teo_dim = interp(mean_dim(2,:), step); teo_dim = [ones(1,size(Spike_count_raster,1)-length(teo_dim))*teo_dim(1), teo_dim];
    corrcoef(vlpfc_dim,teo_dim)

    figure; hold on; set(gcf,'Position',[150 250 1500 300])

    %Get behavior indices
    idx_rest = find(behavior_labels == length(behav_categ)); 
    idx_groom = find(behavior_labels == 7); 
    idx_getgroom = find(behavior_labels == 8); 
    idx_selfgroom = find(behavior_labels == 24);
    idx_forage = find(behavior_labels == 5 | behavior_labels == 4);
    idx_agg = find(behavior_labels == 9 | behavior_labels == 10 |behavior_labels == 1 |behavior_labels == 21 |behavior_labels == 22);
    idx_travel = find(behavior_labels == 2 | behavior_labels == 12 |behavior_labels == 18);
    idx_scratch = find(behavior_labels == 23);

    xline([block_times.end_time_round(1), block_times.end_time_round(2)], "-",["Block 1 end", "Block 2 end"], "LineWidth",2);
    
    %Full dimensionality trace
    session_length = size(Spike_count_raster,1);
    trace = vlpfc_dim;
    plot(trace)

    %Rest
    to_plot=nan(1,session_length); to_plot(idx_rest)=trace(idx_rest);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0.5 0.5 0.5])

    %Getting groomed
    to_plot=nan(1,session_length); to_plot(idx_getgroom)=trace(idx_getgroom);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",'b')

    %Grooming
    to_plot=nan(1,session_length); to_plot(idx_groom)=trace(idx_groom);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",'c')

    %Self-groom
    to_plot=nan(1,session_length); to_plot(idx_selfgroom)=trace(idx_selfgroom);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0.8 0 0.8])

    %Foraging
    to_plot=nan(1,session_length); to_plot(idx_forage)=trace(idx_forage);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0 0.7 0])

    %Aggression
    to_plot=nan(1,session_length); to_plot(idx_agg)=trace(idx_agg);
    plot(1:session_length, to_plot, "LineWidth",2, "Color","r")

    %Travel
    to_plot=nan(1,session_length); to_plot(idx_travel)=trace(idx_travel);
    plot(1:session_length, to_plot, "LineWidth",2, "Color","y")

    %Scratch
    to_plot=nan(1,session_length); to_plot(idx_scratch)=trace(idx_scratch);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0.8 0.6 0])

    legend({'','','', 'Rest','Getting groomed','Grooming','Self-groom','Foraging','Aggression','Travel','Scratch'},...
        'Location','best')

    set(gca,'FontSize',15);
    ylabel(['Dimensionality (' num2str(var_explained_threshold) '% of variation)']); xlabel('Time (s)');
    title('Dimensionality through time')
    ax = gca;
    ax.FontSize = 14;


end%end of session loop
