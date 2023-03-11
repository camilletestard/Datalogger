%% log_plot5F.m

cd(['~/Dropbox (Penn)/Datalogger/Results/All_sessions/Mvmt_results/']);
load(['SVM_results_mvmtControlled_rawRes_c0.2.mat'])

data = cell2mat(mean_hitrate');
data_shuffle = cell2mat(mean_hitrate_shuffled');
data_relative(:,1) = data(:,2)./data(:,1); % Decoding behavior from neural activity without mvmt
data_relative(:,3) = data(:,3)./data(:,1); % Decoding behavior from neural activity without behavior
data_relative(:,4) = data_shuffle(:,1)./data(:,1); %Chance level for decoding behavior

load(['SVM_from_mvmts_c0.2.mat'])
data = cell2mat(mean_hitrate');
data_relative(:,2) = data(:,2)./data(:,1);% Decoding behavior from mvmt

%Relative to full
figure; hold on
bp = bar(mean(data_relative),'FaceAlpha',0.2);
sp1 = scatter(ones(size(data_relative,1))*1,data_relative(:,1), 'filled','b');
sp1 = scatter(ones(size(data_relative,1))*2,data_relative(:,2), 'filled','r');
sp1 = scatter(ones(size(data_relative,1))*3,data_relative(:,3), 'filled','g');
sp1 = scatter(ones(size(data_relative,1))*4,data_relative(:,4), 'filled','g');
ylabel('Accuracy relative to full model'); ylim([0 1])
xticks([1:4]); xticklabels({'Neural (-) mvmt', 'Mvmt', 'Neural (-) behavior','Chance'}); 
xlim([0.25 4.75])
ax = gca;
ax.FontSize = 16;
saveas(gcf, ['SVM_RegressingOutMvmt_RelativeToFull_c' num2str(c_cutoff) '.pdf'])
