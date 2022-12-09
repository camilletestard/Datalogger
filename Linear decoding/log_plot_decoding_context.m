%% log_plot_decoding_context
%Script which plots Fig 3 Panel H.
% poolig results from decoding paired vs. alone and neighbor ID
% C. Testard Nov. 2022

%% load data
cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/SVM_results/']);
load('SVM_allBehavs_AloneVSPaired.mat')
mean_hitrate_AlonePaired = mean_hitrate;
mean_hitrate_AlonePaired_shuffled = mean_hitrate_shuffled;

load('SVM_allBehavs_NeighborID.mat')
mean_hitrate_NeighborID = mean_hitrate;
mean_hitrate_NeighborID_shuffled = mean_hitrate_shuffled;


% Bar plot decoding accuracy
figure; hold on; set(gcf,'Position',[150 250 800 600]);

for b=1:length(behav)

    subplot(2,2,b); hold on

    data_AlonePaired = mean_hitrate_AlonePaired{b}; data_AlonePaired(data_AlonePaired==0)=nan;
    data_NeighborID = mean_hitrate_NeighborID{b}; data_NeighborID(data_NeighborID==0)=nan;
    data_shuffled = mean_hitrate_shuffled{b}; data_shuffled(data_shuffled==0)=nan;

    bp = bar([nanmean(data_AlonePaired(:,:)); nanmean(data_NeighborID(:,:)); nanmean(data_shuffled(:,:))],'FaceAlpha',0.2);

    sp1 = scatter(ones(size(data_AlonePaired,1))*0.77,data_AlonePaired(:,1),8, 'filled','b');
    sp1 = scatter(ones(size(data_AlonePaired,1)),data_AlonePaired(:,2),8, 'filled','r');
    sp1 = scatter(ones(size(data_AlonePaired,1))*1.22,data_AlonePaired(:,3),8, 'filled','y');

    sp1 = scatter(ones(size(data_NeighborID,1))*1.77,data_NeighborID(:,1),8, 'filled','b');
    sp1 = scatter(ones(size(data_NeighborID,1))*2,data_NeighborID(:,2),8, 'filled','r');
    sp1 = scatter(ones(size(data_NeighborID,1))*2.22,data_NeighborID(:,3),8, 'filled','y');

    sp1 = scatter(ones(size(data_NeighborID,1))*2.77,data_shuffled(:,1),8, 'filled','b');
    sp1 = scatter(ones(size(data_NeighborID,1))*3,data_shuffled(:,2),8, 'filled','r');
    sp1 = scatter(ones(size(data_NeighborID,1))*3.22,data_shuffled(:,3),8, 'filled','y');

    %legend(bp,{'vlPFC','TEO','all'},'Location','best')
    title(behav_categ(behav(b)))

    ylabel('Decoding Accuracy'); ylim([0.45 1])
    xticks([1 2 3]); xticklabels({'Alone vs. Paired', 'Neighbor ID', 'Shuffled'}); xlim([0.25 3.75])
    ax = gca;
    ax.FontSize = 16;

end
saveas(gcf,['Decoding_SocialContext.pdf'])
