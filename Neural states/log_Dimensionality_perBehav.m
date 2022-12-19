%%log_Dimensionality_perBehav
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
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 1;%set the smoothing window size (sigma)
var_explained_threshold=90;
num_iter = 500; num_units = 200;
simplify=4;
agg_precedence =1;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16];
end


Cmap = [[1 0 0];[1 0.4 0.1];[0 0 0];[0.1 0.8 0.9];[0 0.7 0];[1 0 1];[0 1 1];...
    [0 0 1];[0.8 0 0];[1 0 0];[0 0 0];[0.2 0.9 0.76];[0 0 0];[0 0 0];[0.7 0 1];...
    [0 0 0];[0 0 0];[0.9 0.5 0];[0 0 0];[0 0 0];[0.8 0 0];[1 0 0];[0.9 0.7 0.12];[0.5 0.2 0.5];...
    [0 0 0];[0 0 0];[0.8 0.4 0.4];[0 0 0];[0.8 0.8 0.8]];

if simplify ==1
    dim = nan(max(session_range),2,5,num_iter);
    min_occurrences = 50;
    Cmap = [[1 0 0]; [0 0.7 0]; [0 1 1]; [0.9 0.5 0]; [0.5 0.5 0.5]];
elseif simplify ==2
    dim = nan(max(session_range),2,2,num_iter);
    min_occurrences = 230;
    Cmap = [[0 1 1]; [0.5 0.5 0.5]];
elseif simplify ==3
    dim = nan(max(session_range),2,4,num_iter);
    min_occurrences = 100;
else
    dim = nan(max(session_range),2,8,num_iter);
    min_occurrences = 50;
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/Dimensionality_results/'];

    chan = 1;

    %for channel_flag = ["vlPFC", "TEO"]
        %channel_flag = "vlPFC";

        %% Get data with specified temporal resolution and channels
        %Get data with specified temporal resolution and channels
        if with_partner ==1
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= ...
                log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma);
        else
            [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
                reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= ...
                log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
                is_mac, with_NC, isolatedOnly, smooth, sigma, agg_precedence );
        end

        disp('Data Loaded')

        Spike_count_raster = Spike_rasters';

        behavior_labels = cell2mat({labels{:,3}}');
        behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %Make proximity equal to rest

        if simplify == 1 %Compare across 5 behavioral categories that are common across sessions
            %Simplify behavioral catagories
            %Lump all aggressive interactions together
            behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
            behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");

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

        elseif simplify == 2 % Compare one behavior to all others

            %Lump all that is not grooming together
            behavior_labels(ismember(behavior_labels,find(behav_categ~="Groom partner" & behav_categ~="Getting groomed")))=find(behav_categ=="Rest");

            behav = [find(behav_categ=="Getting groomed"),find(behav_categ=="Rest")];

% %             %Lump all that is not foraging
% %             behavior_labels(ismember(behavior_labels,find(behav_categ~="Foraging" )))=find(behav_categ=="Rest");
% % 
% %             behav = [find(behav_categ=="Foraging"),find(behav_categ=="Rest")];

        elseif simplify == 3 %Comapre grooming behaviors between each other

            %Lump all that is not grooming together
            behavior_labels(ismember(behavior_labels,find(behav_categ~="Groom partner" & behav_categ~="Getting groomed" & behav_categ~="Self-groom")))=find(behav_categ=="Rest");

            behav = [find(behav_categ=="Groom partner"),find(behav_categ=="Getting groomed"),find(behav_categ=="Self-groom",find(behav_categ=="Rest"))];

        else %Compare all behaviors separately (without pooling across)
            behav = [4,5,7,8,9,10,24,29];
        end

        behav_freq_table = tabulate(behavior_labels);

        %% Compute dimensionality over increasing numbers of units, over multiple iterations

        idx_all_beh= find(ismember(behavior_labels,behav));
        behav_freq=tabulate(behavior_labels(idx_all_beh));
        n_per_behav{s} = behav_freq(behav_freq(:,2)>0,2);

        if all(behav_freq(behav_freq(:,2)>0,2)>min_occurrences) && length(find(behav_freq(:,2)>0))>=length(behav)

            for b = 1:length(behav)

                disp(['Behavior: ' behav_categ(behav(b))])


                for iter = 1:num_iter

                    %Select time points
                    idx= find(ismember(behavior_labels,behav(b)));
                    idx_beh = idx(randsample(1:length(idx),min_occurrences)); %subsample number of indices


                   %Select neural data
                   neural_data = (Spike_count_raster(idx_beh,randsample(size(Spike_count_raster,2), num_units)));

                   %Compute mean firing rate, variance and correlation in response
                    mean_Hz{s, b}(iter,:) = mean(neural_data);
                    std_Hz{s, b}(iter,:) = std(neural_data);
                    cv_Hz{s, b}(iter,:) = std(neural_data)./mean(neural_data);
                    corr_data=reshape(triu(corrcoef(neural_data)), 1,[]);
                    correl_Hz{s, b, iter} = corr_data(corr_data~=0 & corr_data~=1);

                    %PCA
                    Input_matrix = zscore(neural_data);
                    [coeff{b,iter},score,~,~,explained] = pca(Input_matrix);

                    %Get dimensionality
                    var_explained{b}(:,iter) = cumsum(explained);
                    idxl = find(var_explained{b}(:,iter)>=var_explained_threshold);
                    dim(s,chan,b,iter) = min(idxl);
                    var_expl_dim3(s,chan,b,iter) = var_explained{b}(3,iter);

                    %get volume occupied
                    P = score(:,1:3);
                    [k,vol(s,chan,b,iter)] = convhulln(P);
%                     trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','cyan')


                end % end of interation loop

                var_explained_mean{s}(:,b)= mean(var_explained{b},2);

            end %end of behavior loop

        

        chan=chan+1;
        disp([channel_flag ' done'])
    %end %end of channel loop

    %% Plot results for the session


% % % %     figure; hold on
% % % %     for b=1:size(var_explained_mean{s},2)
% % % %         plot(var_explained_mean{s}(:,b),'LineWidth',3, 'Color', Cmap(behav(b),:))
% % % %     end
% % % %     yline(90,'LineStyle','--')
% % % %     xline(3,'LineStyle',':')
% % % %     xlabel('Dimensions')
% % % %     ylabel('Var. explained')
% % % %     legend(behav_categ(behav))
% % % %     ax = gca;
% % % %     ax.FontSize = 14;
% % % %     %saveas(gcf,[savePath '/DimensionalityPerBehav_allExplainedVar.pdf'])
% % % % 
% % % %         figure; hold on; set(gcf,'Position',[150 250 600 400])
% % % %     
% % % %         dim_plot = squeeze(squeeze(dim(s,1,:,:)))';
% % % %         mean_dim = squeeze(mean(dim(s,1,:,:),4)); [~, orderIdx] = sort(mean_dim);
% % % %         violin(dim_plot(:,orderIdx), 'facecolor',Cmap(behav(orderIdx),:))
% % % %         xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
% % % %         xticklabels(behav_categ(behav(orderIdx))); %ylim([0 35])
% % % %         ax = gca;
% % % %         ax.FontSize = 14;
% % % %         ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
% % % %         title('Dimensionality across behaviors')
    %
    %     saveas(gcf,[savePath '/DimensionalityPerBehav_2categ.pdf'])


% % %         figure; hold on; set(gcf,'Position',[150 250 600 400])
% % %     
% % %         dim_plot = squeeze(squeeze(var_expl_dim3(s,1,:,:)))';
% % %         mean_dim = squeeze(mean(var_expl_dim3(s,1,:,:),4)); [~, orderIdx] = sort(mean_dim);
% % %         violin(dim_plot(:,orderIdx), 'facecolor',Cmap(behav(orderIdx),:))
% % %         xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
% % %         xticklabels(behav_categ(behav(orderIdx))); %ylim([0 35])
% % %         ax = gca;
% % %         ax.FontSize = 14;
% % %         ylabel(['Variance explained at 3D'],'FontSize', 14);
% % %         title('Variance explained at 3D')


% % %         figure; hold on; set(gcf,'Position',[150 250 600 400])
% % %     
% % %         vol_plot = squeeze(squeeze(vol(s,1,:,:)))';
% % %         mean_vol = squeeze(mean(vol(s,1,:,:),4)); [~, orderIdx] = sort(mean_vol);
% % %         violin(vol_plot(:,orderIdx), 'facecolor',Cmap(behav(orderIdx),:))
% % %         xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
% % %         xticklabels(behav_categ(behav(orderIdx))); %ylim([0 35])
% % %         ax = gca;
% % %         ax.FontSize = 14;
% % %         ylabel(['Volume in 3D space'],'FontSize', 14);
% % %         title('Volume across behaviors')
        %
        %     close all

        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(s)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        end % end of if clause

end%end of session loop

cd([home '/Dropbox (Penn)/Datalogger/Results/All_sessions/Dimensionality_results/']);
save('Dimensionality_allBehav.mat', "dim","behav","a_sessions","h_sessions","behav_categ","Cmap")

%% Plot results across sessions

load('Dimensionality_GroomCateg.mat')

dim_amos = squeeze(nanmean(dim(a_sessions,:,:,:),1));
dim_hooke = squeeze(nanmean(dim(h_sessions,:,:,:),1));
dim_all = squeeze(nanmean(dim,1));

%All Pooled
figure; hold on; lowlimit=17; uplimit=35;
mean_dim = squeeze(mean(dim_all(1,:,:),3)); [~, orderIdx] = sort(mean_dim);
violin(squeeze(dim_all(1,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
ax = gca;
ax.FontSize = 14;
ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);



mean_alldims = mean(cat(3,(var_explained_mean{:})),3);
sd_alldims =  std(cat(3,(var_explained_mean{:})),[],3);
    figure; hold on
    for b=1:size(mean_alldims,2)
        y=mean_alldims(:,b);
        sd=sd_alldims(:,b);
        upper_lim=y+sd;
        lower_lim=y-sd;
%         p = fill([1:length(y) length(y):-1:1],[upper_lim; flip(lower_lim)],'red');
%         p.FaceColor = Cmap(behav(b),:);
%         p.FaceAlpha = 0.3;
%         p.EdgeColor = 'none';

        plot(y,'LineWidth',3, 'Color', Cmap(behav(b),:))
    end
    yline(90,'LineStyle','--')
    xline(3,'LineStyle',':')
    xlabel('Dimensions')
    ylabel('Var. explained')
    legend(behav_categ(behav))
    ax = gca;
    ax.FontSize = 14;
    %saveas(gcf,[savePath '/DimensionalityPerBehav_allExplainedVar.pdf'])

% % % % %Separate per brain area, pooling both monkeys
% % % % figure; hold on; set(gcf,'Position',[150 250 1000 500]); lowlimit=20; uplimit=35;
% % % % subplot(1,2,1); hold on
% % % % mean_dim = squeeze(mean(dim_all(1,:,:),3)); [~, orderIdx] = sort(mean_dim);
% % % % violin(squeeze(dim_all(1,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
% % % % xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
% % % % xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
% % % % ax = gca;
% % % % ax.FontSize = 14;
% % % % ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
% % % % title('vlPFC')
% % % % 
% % % % subplot(1,2,2); hold on
% % % % mean_dim_teo = squeeze(mean(dim_all(2,:,:),3)); [~, orderIdx] = sort(mean_dim_teo);
% % % % violin(squeeze(dim_all(2,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
% % % % xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
% % % % xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
% % % % ax = gca;
% % % % ax.FontSize = 14;
% % % % ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
% % % % title('TEO')
% % % % saveas(gcf,'DimensionalityPerBehav_GroomCateg.pdf')
% % % % 
% % % % 
% % % % %Separated by monkey and brain area
% % % % figure; hold on; set(gcf,'Position',[150 250 1000 800]); lowlimit=25; uplimit=45;
% % % % 
% % % % %Amos
% % % % subplot(2,2,1); hold on
% % % % mean_dim = squeeze(mean(dim_amos(1,:,:),3)); [~, orderIdx] = sort(mean_dim);
% % % % violin(squeeze(dim_amos(1,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
% % % % xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
% % % % xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
% % % % ax = gca;
% % % % ax.FontSize = 14;
% % % % ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
% % % % title('Amos, vlPFC')
% % % % 
% % % % subplot(2,2,2); hold on
% % % % mean_dim_teo = squeeze(mean(dim_amos(2,:,:),3)); [~, orderIdx] = sort(mean_dim_teo);
% % % % violin(squeeze(dim_amos(2,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
% % % % xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
% % % % xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
% % % % ax = gca;
% % % % ax.FontSize = 14;
% % % % ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
% % % % title('Amos, TEO')
% % % % 
% % % % %Hooke
% % % % subplot(2,2,3); hold on
% % % % mean_dim = squeeze(mean(dim_hooke(1,:,:),3)); [~, orderIdx] = sort(mean_dim);
% % % % violin(squeeze(dim_hooke(1,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
% % % % xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
% % % % xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
% % % % ax = gca;
% % % % ax.FontSize = 14;
% % % % ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
% % % % title('Hooke, vlPFC')
% % % % 
% % % % subplot(2,2,4); hold on
% % % % mean_dim_teo = squeeze(mean(dim_hooke(2,:,:),3)); [~, orderIdx] = sort(mean_dim_teo);
% % % % violin(squeeze(dim_hooke(2,orderIdx,:))', 'facecolor',Cmap(behav(orderIdx),:))
% % % % xticks([1:length(behav)]); xlim([0.5 length(behav)+0.5]);
% % % % xticklabels(behav_categ(behav(orderIdx))); ylim([lowlimit uplimit])
% % % % ax = gca;
% % % % ax.FontSize = 14;
% % % % ylabel(['Dims needed to explain ' num2str(var_explained_threshold) '% of variation'],'FontSize', 14);
% % % % title('Hooke, TEO')
% % % % 
% % % % saveas(gcf,'DimensionalityPerBehav_GroomCateg.pdf')