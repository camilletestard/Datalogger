%% Log_Visualize_SingleUnit_Selectivity
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
session_range_no_partner=[1:6,11:13,15:16,18];
session_range_with_partner=[1:6,11:13,15:16,18];

%Set parameters
with_partner = 0; %need partner info? 0: No; 1:yes
temp_resolution = 1; %Temporal resolution of firing rate. 1: 1sec; 10:100msec; 0.1: 10sec
channel_flag = "all"; %Channels considered. vlPFC, TEO or all
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
min_occurrence =30; %Minimum number of occurrences in the session needed to be considered for this analysis.
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 10*temp_resolution;%set the smoothing window size in sec (sigma)


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
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/Example_units'];

    %% Load data

    %Get data with specified temporal resolution and channels
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= ...
            log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
            reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= ...
            log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, ...
            is_mac, with_NC, isolatedOnly, smooth, sigma);
    end

    session_length = size(Spike_rasters,2); % get session length

    %% Pre-process behavior
    %Extract behavior labels
    behavior_labels = cell2mat({labels{:,3}}');%Get behavior label from labels structure

    %Simplify behavioral categories
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").
    
    %Lump all aggressive interactions together
    behavior_labels(behavior_labels==find(behav_categ=="Threat to partner"))=find(behav_categ=="Aggression");
    behavior_labels(behavior_labels==find(behav_categ=="Threat to subject"))=find(behav_categ=="Aggression");
    behavior_labels(behavior_labels==find(behav_categ=="Squeeze partner"))=find(behav_categ=="Aggression");
    behavior_labels(behavior_labels==find(behav_categ=="Squeeze Subject"))=find(behav_categ=="Aggression");

    %Lump all travel together
    behavior_labels(behavior_labels==find(behav_categ=="Approach"))=find(behav_categ=="Travel");
    behavior_labels(behavior_labels==find(behav_categ=="Leave"))=find(behav_categ=="Travel");

    %% Set parameters
    unqLabels = [1,5,7,8,24,18];%1:max(behavior_labels)-1; %Get unique behavior labels (exclude rest)
    n_neurons(s) = size(Spike_rasters,1); %Get number of neurons
    n_behav = length(unqLabels); %Get number of unique behavior labels

    %Estimate "baseline" neural firing distribution.
    if with_partner ==0
        idx_rest=find(behavior_labels==length(behav_categ));%Get idx of "rest" epochs.
    else
        behavior_labels_partner = cell2mat({labels_partner{:,3}}');
        behavior_labels_partner(behavior_labels_partner==find(behav_categ=="Proximity"))=length(behav_categ);
        idx_rest = intersect(find(behavior_labels ==length(behav_categ)), find(behavior_labels_partner ==length(behav_categ)));
        percent_corest = length(intersect(find(behavior_labels ==length(behav_categ)),...
            find(behavior_labels_partner ==length(behav_categ))))/length(find(behavior_labels ==length(behav_categ)));
    end
    mean_baseline = mean(Spike_rasters(:,idx_rest),2);
    std_baseline = std(Spike_rasters(:,idx_rest),0,2);


    %% Compute cohen's d

    cohend = nan(n_neurons(s),n_behav);
    cohend_shuffle = nan(n_neurons(s),n_behav);
    mean_beh = nan(n_neurons(s), n_behav);
    mean_beh_shuffle = nan(n_neurons(s), n_behav);
    std_beh = nan(n_neurons(s), n_behav);
    std_beh_shuffle = nan(n_neurons(s), n_behav);
    p = nan(n_neurons(s), n_behav);
    p_rand = nan(n_neurons(s), n_behav);

    for n = 1:n_neurons(s) %for all neurons

        for b = 1:n_behav %for all behaviors
            idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
            n_per_behav(s,b)=length(idx);%get the number of time points for behavior b

            if n_per_behav(s,b)>min_occurrence % if behavior occurs at least during 'min_occurrence' time points

                if length(idx)<length(idx_rest)% sample the same number of rest time points 
                    idx_rand = randsample(idx_rest,length(idx));
                else
                    idx_rand = randsample(idx_rest,length(idx),true);
                end

                mean_beh(n,b)=mean(Spike_rasters(n, idx),2); %get the mean firing rate during behavior b
                std_beh(n,b)=std(Spike_rasters(n, idx),0,2); %get the standard deviation firing rate during behavior b

                %extract a "shuffled" mean and standard deviation firing
                %rate. This is to control for stochastic differences in
                %firing rate compared to baseline.
                mean_beh_shuffle(n,b)=mean(Spike_rasters(n, idx_rand),2); %get the mean firing rate during rest, randomly sampling the same #data points as behavior b
                std_beh_shuffle(n,b)=std(Spike_rasters(n, idx_rand),0,2); %get the standard deviation firing rate during rest, randomly sampling the same #data points as behavior b

                %Compute a cohen d between the distribution of firing rate
                %during behavior b and a baseline state (rest)
                cohend(n,b) = (mean_beh(n,b)-mean_baseline(n)) ./ sqrt( ((n_per_behav(s,b)-1)*(std_beh(n,b).^2) + (length(idx_rest)-1)*(std_baseline(n).^2)) / (n_per_behav(s,b)+length(idx_rest)-2) ); %Compute cohen d
                cohend_shuffle(n,b) = (mean_beh_shuffle(n,b)-mean_baseline(n)) ./ sqrt( ((n_per_behav(s,b)-1)*(std_beh_shuffle(n,b).^2) + (length(idx_rest)-1)*(std_baseline(n).^2)) / (n_per_behav(s,b)+length(idx_rest)-2) );%Compute cohen d

                %get p-value from ttest comparing the distribution of
                %firing rate during behavior b and rest.
                [~, p(n,b)] = ttest2(Spike_rasters(n, idx), Spike_rasters(n,idx_rest));
                [~, p_rand(n,b)] = ttest2(Spike_rasters(n, idx_rand), Spike_rasters(n,idx_rest));

            end

        end
    end

    %% Plot firing rate profile across macaques behavioral repertoire

    %Set color scheme
   color_scheme = {'r','','','',[0 0.7 0],'','c','b','','','','','','','','','','y','','','','',[0.8 0.6 0],[0.8 0 0.8]};
%     Cmap = [[0 0 0];[1 0.4 0.1];[0 0 0];[0 0.6 0.8];[0 0.7 0];[1 0 0.65];[0 1 1];...
%             [0 0 1];[0.5 0 0];[1 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];...
%             [0 0 0];[0 0 0];[0.9 0.7 0.12];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0.8 0.6 0];[0.8 0 0.8];...
%             [0 0 0];[0 0 0];[0.2 0.2 1]];

    n= 138;%randsample(unit_count(3),1);%42;138;218;222
    cohend_per_neuron = cohend(n,:);
    [~, orderIdx] = sort(cohend_per_neuron, 'ascend');

    %Plot distribution of firing rate across behaviors for example neuron n
    figure; hold on; set(gcf,'Position',[150 250 1500 300]); i=1;
    for b = unqLabels(orderIdx)
        subplot(1, length(orderIdx), i); hold on
        idx = find(behavior_labels == b); %get idx where behavior b occurred

        %plot firing rate distribution for example unit
        histogram(Spike_rasters(n, idx),'BinWidth',1,'Normalization','pdf', 'FaceColor',color_scheme{b})
        histogram(Spike_rasters(n, idx_rest),'BinWidth',1,'Normalization','pdf', 'FaceColor',[0.5 0.5 0.5])
        legend({behav_categ(b),'Rest'})
        title([behav_categ(b)])

        set(gca,'FontSize',15);
        xlabel('Firing rate (Hz)'); ylabel('Normalized Frequency');
        i=i+1;
    end
    %saveas(gcf, [savePath '/ExampleUnitSelectivity.eps'])


    %Plot smoothed firing rate over the length of the session and highlight
    %behavior tipes

    %Get behavior indices
    idx_groom = find(behavior_labels == 7); 
    idx_getgroom = find(behavior_labels == 8); 
    idx_selfgroom = find(behavior_labels == 24);
    idx_forage = find(behavior_labels == 5);
    idx_agg = find(behavior_labels == 1);
    idx_travel = find(behavior_labels == 2 | behavior_labels == 12 |behavior_labels == 18);

    unit=n;
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

    %Self-groom
    to_plot=nan(1,session_length); to_plot(idx_selfgroom)=Spike_rasters(unit,idx_selfgroom);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0.8 0 0.8])

    %Foraging
    to_plot=nan(1,session_length); to_plot(idx_forage)=Spike_rasters(unit,idx_forage);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0 0.7 0])

    %Aggression
    to_plot=nan(1,session_length); to_plot(idx_agg)=Spike_rasters(unit,idx_agg);
    plot(1:session_length, to_plot, "LineWidth",2, "Color","r")

    %Travel
    to_plot=nan(1,session_length); to_plot(idx_travel)=Spike_rasters(unit,idx_travel);
    plot(1:session_length, to_plot, "LineWidth",2, "Color","y")

    legend({'','','Smooth spike trace', 'Rest','Getting groomed','Grooming','Self-groom','Foraging','Aggression','Travel','Scratch'},...
        'Location','bestoutside')

    set(gca,'FontSize',15);
    ylabel('Firing rate (Hz)'); xlabel('Time (s)');
    title(['Example Unit from ' brain_label(n)])
    %saveas(gcf, [savePath '/ExampleUnitLabeledSmoothTrace.eps'])
    
    close all

end

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