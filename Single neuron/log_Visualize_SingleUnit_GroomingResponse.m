%% Log_Visualize_SingleUnit_GroomingResponse
%  This script produces two plots for a given neuron to illustrate its
%  response profile. First, it computes firing rate of individual neuron 
%  under different behavioral conditions. 
% 1st plot: a proportional histogram of firing rates for a given behavior, compared to rest. 
% 2nd plot: smooth spike trace with labeled time points according to a set of behaviors.
% C. Testard, July 2022

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range=[1:6,11:13,15:16,18];

%Set parameters
with_partner =0;
temp = 1; temp_resolution = 1;
channel_flag = "all";
randomsample=0; %subsample neurons to match between brain areas
unq_behav=0; %If only consider epochs where only 1 behavior happens
with_NC =1;%0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly= 0;%Only consider isolated units. 0=all units; 1=only well isolated units
smooth= 1; %smooth the data
sigma = 10*temp_resolution;%set the smoothing window size (sigma)
threat_precedence =0;
exclude_sq=1;
plot_toggle=0;

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/Example_units'];

    %% Load data

    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label{s}, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')

    Spike_count_raster = zscore(Spike_rasters');

    %Extract behavior labels
    behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure

    %Set colormap
    %uisetcolor([0.6 0.8 1])
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
        [0 0 0];...%SP; NA
        [0 0 0];...%SS; NA
        [0.6314 0.5059 0.0118];...%Scratch; maroon
        [0.5 0.2 0.5];...%Self-groom; dark purple
        [ 1 0.07 0.65];...%Submission; dark pink
        [0 0.4 0.5];...%Vocalzation; blue green
        [0 0 0];...%Yawning; NA
        [0.8 0.8 0.8]];%Rest; grey

    %% Pre-process behavior

    %Simplify behavioral categories
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

    %Lump all aggressive interactions together
    behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
    behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");

    %Lump all travel together
    behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
    behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");


    %% Plot actvity for one session in relation to behavior

    %Exclude axctivity at the very beginning and very end of the
    %session
    session_with_buffer = [250*temp_resolution:size(Spike_count_raster,1)-250*temp_resolution];
    
    %% Plot firing rate profile across macaques behavioral repertoire

    %Set parameters
    unqLabels = [7,8];%1:max(behavior_labels)-1; %Get unique behavior labels (exclude rest)
    n_neurons(s) = size(Spike_rasters,1); %Get number of neurons
    n_behav = length(unqLabels); %Get number of unique behavior labels

    neurons = [1:7, 77, 99, 130, 135, 171, 187, 252, 266];%randsample(n_neurons(s),1)
    %135: Figure 2 paper. To add in supplements: 77, 252; 204; 187; 130
    %266: alone block table top
    %99: super activated during grooming
    %171: example grooming neuron. Also 6
    %4: Super sensitive to aggression
    %1: Sensitive to grooming

    n=49;
    for n=1:size(Spike_rasters,1)
    
        idx_gg = find(behavior_labels == 7); %get idx where behavior b occurred
        idx_gr = find(behavior_labels == 8);

        mean_gg{s}(n) = mean(Spike_count_raster(idx_gg,n));
        mean_gr{s}(n) = mean(Spike_count_raster(idx_gr,n));

        if plot_toggle
        %plot firing rate distribution for example unit
        figure; hold on;
        histogram(Spike_rasters(n, idx_gg),'BinWidth',0.5,'Normalization','pdf', 'FaceColor',Cmap(7,:))
        histogram(Spike_rasters(n, idx_gr),'BinWidth',0.5,'Normalization','pdf', 'FaceColor',Cmap(8,:))
        legend({'Groom give','Groom receive'})

        set(gca,'FontSize',15);
        xlabel('Firing rate (Hz)'); ylabel('Normalized Frequency');

        sgtitle(['Unit #' num2str(n)])
        %saveas(gcf, [savePath '/ExampleUnitSelectivity.eps'])


        %Plot smoothed firing rate over the length of the session and highlight
        %behavior tipes

        %Get behavior indices
        idx_groom = find(behavior_labels == 7);
        idx_getgroom = find(behavior_labels == 8);
        idx_rest = find(behavior_labels == length(behav_categ));
      
        unit=n;
        session_length = size(Spike_rasters,2);
        figure; hold on; set(gcf,'Position',[150 250 1500 500]);
        xline([block_times.end_time_round(1), block_times.end_time_round(2)], "-",["Block 1 end", "Block 2 end"], "LineWidth",2);
        session_with_buffer = [120:session_length-120];

        %Full spike rate trace
        plot(1:session_length, Spike_rasters(unit,1:session_length));

        %Rest
        to_plot=nan(1,session_length); to_plot(idx_rest)=Spike_rasters(unit,idx_rest);
        plot(1:session_length, to_plot, "LineWidth",2, "Color",[0.5 0.5 0.5])

        %Getting groomed
        to_plot=nan(1,session_length); to_plot(idx_getgroom)=Spike_rasters(unit,idx_getgroom);
        plot(1:session_length, to_plot, "LineWidth",2, "Color",'b')

        %Grooming
        to_plot=nan(1,session_length); to_plot(idx_groom)=Spike_rasters(unit,idx_groom);
        plot(1:session_length, to_plot, "LineWidth",2, "Color",'c')

% %         legend({'','','Smooth spike trace', 'Rest','Getting groomed','Grooming','Self-groom','Foraging','Aggression','Travel','Scratch'},...
% %             'Location','bestoutside')

        set(gca,'FontSize',15);
        ylabel('Firing rate (Hz)'); xlabel('Time (s)');
        title(['Example Unit #' num2str(n) ' from ' brain_label(n)])
        %saveas(gcf, [savePath '/ExampleUnitLabeledSmoothTrace.eps'])
        end
    
    end

%     figure; hold on
%     scatter(mean_gg{s},mean_gr{s})
%     xlabel("Zscore activity guring groom give")
%     ylabel("Zscore activity durin groom receive")
%     hL=plot([-1.5 1.5],[-1.5 1.5],'DisplayName','Diagonal'); 
    
    
    close all

end

mean_gg_all = cat(2,mean_gg{:});
mean_gr_all = cat(2,mean_gr{:});
brain_label_all = cat(2,brain_label{:});


figure; hold on
scatter(mean_gg_all(strcmp(brain_label_all,"TEO")),mean_gr_all(strcmp(brain_label_all,"TEO")), 'Color','g')
scatter(mean_gg_all(strcmp(brain_label_all,"vlPFC")),mean_gr_all(strcmp(brain_label_all,"vlPFC")), 'Color', 'b')
xlabel("Zscore activity guring groom give")
ylabel("Zscore activity durin groom receive")
hL=plot([-2 2],[-2 2],'DisplayName','Diagonal');

% [p, h, stats ]= signrank(mean_gg_all,mean_gr_all);
% [h,p,ci,stats]=ttest(mean_gg_all,mean_gr_all);
% 
% figure
% histogram(mean_gg_all-mean_gr_all, 'BinWidth',0.1)
% length(find(mean_gg_all-mean_gr_all>0))/length(mean_gg_all)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Old code
%     %Plot example spike trains
%     %figure; hold on
%     spiketrains = 50;
%     length_bout = 10;
%     %xlim([0 length_bout]); ylim([0 spiketrains])
%     spike_times_behav = []; spike_times_rest = []; trials_behav = []; trials_rest = [];
%     for st=1:spiketrains
%         Hz_behav =Spike_rasters(n,randsample(idx,length_bout));
%         Hz_rest =Spike_rasters(n,randsample(idx_rest,length_bout));
% 
%         total_spikes_behav = 0; total_spikes_rest = 0;
%         for lb = 1:length_bout
%             ticks_behav = rand(Hz_behav(lb),1);
%             spike_times_behav = [spike_times_behav; lb-1+ticks_behav];
%             total_spikes_behav = total_spikes_behav + Hz_behav(lb);
% 
%             ticks_rest = rand(Hz_rest(lb),1);
%             spike_times_rest=[spike_times_rest; lb-1+ticks_rest];
%             total_spikes_rest = total_spikes_rest+Hz_rest(lb);
%         end
% 
%         trials_rest = [trials_rest, ones(1,total_spikes_rest)*st];
%         trials_behav = [trials_behav, ones(1,total_spikes_behav)*st];
%     end
%     figure; subplot(2,1,1); spikeRasterPlot(seconds(spike_times_behav), trials_behav,'ColorOrder',[0 0 1])
%     subplot(2,1,2); spikeRasterPlot(seconds(spike_times_rest), trials_rest,'ColorOrder',[1 0 0])