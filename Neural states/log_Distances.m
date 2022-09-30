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
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1;%set the smoothing window size (sigma)
var_explained_threshold=90;
num_iter = 500; num_units = 100;
simplify=0;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end

if simplify ==1
    dim = nan(max(session_range),2,5,num_iter);
    min_occurrences = 50;
elseif simplify ==2
    dim = nan(max(session_range),2,2,num_iter);
    min_occurrences = 230;
elseif simplify ==3
    dim = nan(max(session_range),2,3,num_iter);
    min_occurrences = 50;
else
    dim = nan(max(session_range),2,8,num_iter);
    min_occurrences = 30;
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Dimensionality_results/'];

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

            %Lump Drinking and foraging
            behavior_labels(behavior_labels==find(behav_categ=="Drinking"))=find(behav_categ=="Foraging");

            %Lump all grooming together
            behavior_labels(behavior_labels==find(behav_categ=="Getting groomed"))=find(behav_categ=="Groom partner");
            behavior_labels(behavior_labels==find(behav_categ=="Groom sollicitation"))=find(behav_categ=="Groom partner");
            behavior_labels(behavior_labels==find(behav_categ=="Self-groom"))=find(behav_categ=="Groom partner");

            behav = [1,5,7,18,29];

        elseif simplify == 2 % Compare one behavior across all others

            %Lump all that is not grooming together
            behavior_labels(ismember(behavior_labels,find(behav_categ~="Groom partner" & behav_categ~="Getting groomed")))=find(behav_categ=="Rest");

            behav = [find(behav_categ=="Getting groomed"),find(behav_categ=="Rest")];

            %             %Lump all that is not foraging
            %             behavior_labels(ismember(behavior_labels,find(behav_categ~="Foraging" )))=find(behav_categ=="Rest");
            %
            %             behav = [find(behav_categ=="Foraging"),find(behav_categ=="Rest")];

        elseif simplify == 3 %Comapre grooming behaviors between each other

            behav = [find(behav_categ=="Groom partner"),find(behav_categ=="Getting groomed"),find(behav_categ=="Self-groom")];

        else %Compare all behaviors separately (without pooling across)
            %behav = [4,5,7,8,9,10,24,29];
            behav = [7,8];

        end %End of simplifying loop

        behav_freq_table = tabulate(behavior_labels);

        %% Compute dimensionality over increasing numbers of units, over multiple iterations

        idx_all_beh= find(ismember(behavior_labels,behav));
        behav_freq=tabulate(behavior_labels(idx_all_beh));
        n_per_behav{s} = behav_freq(behav_freq(:,2)>0,2);
        if all(behav_freq(behav_freq(:,2)>0,2)>min_occurrences) && length(find(behav_freq(:,2)>0))>=length(behav)

            for iter = 1:num_iter

                for b = 1:length(behav)

                    disp(['Behavior: ' behav_categ(behav(b))])

                    %Select time points to run PCA
                    idx= find(ismember(behavior_labels,behav(b)));
                    idx_beh = idx(randsample(1:length(idx),min_occurrences));

                    %Select unit to run PCA
                    Input_matrix{b} = zscore(Spike_count_raster(idx_beh,randsample(size(Spike_count_raster,2), num_units)))';

                end
                Input_matrix_full = cell2mat(Input_matrix);
                %figure; hold on; hist(corr(Input_matrix))

                %PCA
                [coeff,score,~,~,explained] = pca(Input_matrix_full');
                %[coeff,score,~,~,explained] = pca(Spike_count_raster);

                %Get dimensionality
                var_explained = cumsum(explained);
                idxl = find(var_explained>=var_explained_threshold);

                %Get distances
                D=pdist(score(:,1:20), 'cityblock');
                Z = squareform(D);
                heatmap(Z,'Colormap',jet)

            end % end of interation loop

        end % end of if clause

        chan=chan+1;
        disp([channel_flag ' done'])
    end %end of channel loop

    %% Plot results for the session

    Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0.1 0.8 0.9];[0 0.7 0];[1 0 1];[0 1 1];...
        [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0.7 0 1];...
        [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0.8 0 0];[1 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
        [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.5 0.5 0.5]];

    %     figure; hold on; set(gcf,'Position',[150 250 1000 400])
    %
    %     subplot(1,2,1); hold on
    %     dim_vlpfc = squeeze(squeeze(dim(s,1,:,:)))';
    %     mean_dim_vlpfc = squeeze(mean(dim(s,1,:,:),4)); [~, orderIdx] = sort(mean_dim_vlpfc);
    %     sd_dim_vlpfc = squeeze(std(dim(s,1,:,:),0,4));
    %     violin(dim_vlpfc(:,orderIdx), 'facecolor',Cmap(behav(orderIdx),:))
    %     xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
    %     xticklabels(behav_categ(behav(orderIdx))); ylim([0 35])
    %     ax = gca;
    %     ax.FontSize = 14;
    %     ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
    %     title('vlPFC')
    %
    %     subplot(1,2,2); hold on
    %     dim_teo = squeeze(squeeze(dim(s,2,:,:)))';
    %     mean_dim_teo = squeeze(mean(dim(s,2,:,:),4)); [~, orderIdx] = sort(mean_dim_teo);
    % %     sd_dim_vlpfc = squeeze(std(dim(s,2,:,:),0,4));
    %     violin(dim_teo(:,orderIdx), 'facecolor',Cmap(behav(orderIdx),:))
    % %     scatter(1:length(behav), mean_dim_vlpfc(orderIdx),40,Cmap(behav(orderIdx),:),'filled')
    % %     errorbar(1:length(behav), mean_dim_vlpfc(orderIdx), sd_dim_vlpfc,'-s','MarkerSize',1, ...
    % %         'MarkerEdgeColor','k')
    %     xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
    %     xticklabels(behav_categ(behav(orderIdx))); ylim([0 35])
    %     ax = gca;
    %     ax.FontSize = 14;
    %     ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
    %     title('TEO')
    %     sgtitle(string(sessions(s).name))
    %
    %     saveas(gcf,[savePath '/DimensionalityPerBehav_2categ.pdf'])
    %
    %     close all

end%end of session loop

cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/Dimensionality_results/']);
save('Dimensionality_2categ.mat', "dim","behav","a_sessions","h_sessions","behav_categ","Cmap")

%% Plot results across sessions

load('Dimensionality_2categ.mat')

dim_amos = squeeze(nanmean(dim(a_sessions,:,:,:),1));
dim_hooke = squeeze(mean(dim(h_sessions,:,:,:),1));


figure; hold on; set(gcf,'Position',[150 250 1000 800]); lowlimit=25; uplimit=55;

%Amos
subplot(2,2,1); hold on
mean_dim_vlpfc = squeeze(mean(dim_amos(1,:,:),3)); [~, orderIdx] = sort(mean_dim_vlpfc);
violin(squeeze(dim_amos(1,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Amos, vlPFC')

subplot(2,2,2); hold on
mean_dim_teo = squeeze(mean(dim_amos(2,:,:),3)); [~, orderIdx] = sort(mean_dim_teo);
violin(squeeze(dim_amos(2,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Amos, TEO')

%Hooke
subplot(2,2,3); hold on
mean_dim_vlpfc = squeeze(mean(dim_hooke(1,:,:),3)); [~, orderIdx] = sort(mean_dim_vlpfc);
violin(squeeze(dim_hooke(1,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Hooke, vlPFC')

subplot(2,2,4); hold on
mean_dim_teo = squeeze(mean(dim_hooke(2,:,:),3)); [~, orderIdx] = sort(mean_dim_teo);
violin(squeeze(dim_hooke(2,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
title('Hooke, TEO')

saveas(gcf,'DimensionalityPerBehav_2categ.pdf')