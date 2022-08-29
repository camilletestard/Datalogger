%% Log_Visualize_SingleUnit_Selectivity
%  This script computes firing rate of individual neuron under different
%  behavioral conditions. Then, it computes a cohen's d (or effect size)
%  difference between the distribution of firing rates during behavior X
%  with a baseline firing rate.

%Set session list
is_mac = 1;
if is_mac
    home = '~';
else
    home ='C:/Users/GENERAL';
end
cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
session_range_no_partner=[1:6,11:13,16];
session_range_with_partner=[1:3,11:13];

%Set parameters
plot_toggle = 0;
with_partner = 0;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
min_occurrence =100;
cohend_cutoff=0; p_cutoff=0.01;%Set thresholds

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    a_sessions = 1:3; h_sessions = 11:13;
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,16];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/Example_units'];

    %% Load data

    %Get data with specified temporal resolution and channels
    if with_partner ==1
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
    else
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all, brain_label]= log_GenerateDataToRes_function_temp(filePath, temp_resolution, channel_flag, is_mac, with_NC, isolatedOnly);
    end

    session_length = size(Spike_rasters,2); % get session length

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
    unqLabels = 1:max(behavior_labels)-1; %Get unique behavior labels (exclude rest)
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

    for n = 1:n_neurons(s)

        for b = 1:n_behav
            idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
            n_per_behav(s,b)=length(idx);

            if n_per_behav(s,b)>min_occurrence

                if length(idx)<length(idx_rest)
                    idx_rand = randsample(idx_rest,length(idx));
                else
                    idx_rand = randsample(idx_rest,length(idx),true);
                end

                mean_beh(n,b)=mean(Spike_rasters(n, idx),2);
                std_beh(n,b)=std(Spike_rasters(n, idx),0,2);

                mean_beh_shuffle(n,b)=mean(Spike_rasters(n, idx_rand),2);
                std_beh_shuffle(n,b)=std(Spike_rasters(n, idx_rand),0,2);

                cohend(n,b) = (mean_beh(n,b)-mean_baseline(n)) ./ sqrt( ((n_per_behav(s,b)-1)*(std_beh(n,b).^2) + (length(idx_rest)-1)*(std_baseline(n).^2)) / (n_per_behav(s,b)+length(idx_rest)-2) );
                cohend_shuffle(n,b) = (mean_beh_shuffle(n,b)-mean_baseline(n)) ./ sqrt( ((n_per_behav(s,b)-1)*(std_beh_shuffle(n,b).^2) + (length(idx_rest)-1)*(std_baseline(n).^2)) / (n_per_behav(s,b)+length(idx_rest)-2) );

                [~, p(n,b)] = ttest2(Spike_rasters(n, idx), Spike_rasters(n,idx_rest));
                [~, p_rand(n,b)] = ttest2(Spike_rasters(n, idx_rand), Spike_rasters(n,idx_rest));

            end

        end
    end

    %Smooth data
    sigma = 10;% 0.045 * opts.Fs; %og SDF window was 45 ms, so simply mutiply 45 ms by fs
    gauss_range = -3*sigma:3*sigma; %calculate 3 stds out, use same resolution for convenience
    smoothing_kernel = normpdf(gauss_range,0,sigma); %Set up Gaussian kernel
    smoothing_kernel = smoothing_kernel/sum(smoothing_kernel);
    smoothing_kernel = smoothing_kernel * 1; %Rescale to get correct firing rate
    Spike_rasters_smooth = conv2(Spike_rasters, smoothing_kernel,'same');

    %Low-pass filter
    Spike_rasters_LowPass = lowpass(Spike_rasters',0.00005,1);
%     figure; hold on; plot(Spike_rasters_LowPass(:,138)); plot(Spike_rasters_smooth(138,:))
    Spike_rasters_smooth = Spike_rasters_LowPass';

    color_scheme = {'r','','','',[0 0.7 0],'','c','b','','','','','','','','','','y','','','','',[0.8 0.6 0],[0.8 0 0.8]};
%     Cmap = [[0 0 0];[1 0.4 0.1];[0 0 0];[0 0.6 0.8];[0 0.7 0];[1 0 0.65];[0 1 1];...
%             [0 0 1];[0.5 0 0];[1 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];...
%             [0 0 0];[0 0 0];[0.9 0.7 0.12];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0.8 0.6 0];[0.8 0 0.8];...
%             [0 0 0];[0 0 0];[0.2 0.2 1]];

    n= 138;%randsample(unit_count(3),1);%42;138;218;222
    cohend_per_neuron = cohend(n,:);
    [~, orderIdx] = sort(cohend_per_neuron, 'ascend');
    cohend_per_neuron_sorted = cohend_per_neuron(:,orderIdx); 
    behav_sorted = behav_categ(orderIdx);
    is_not_nan = ~isnan(cohend_per_neuron_sorted);

    figure; hold on; set(gcf,'Position',[150 250 1500 500]); i=1
    for b = orderIdx(is_not_nan)
        subplot(1, length(orderIdx(is_not_nan)), i); hold on
        %plot firing rate distribution for example unit
        %Specifically comparing groom give to rest
        %b=7 ; [~, n]=max(abs(cohend(:,b))); %Groom give (suppressed activity) 209 is the min separability considered. (For session s=1)
        %b=10 ;[max_val, n]=max(abs(cohend(:,b))); %Threat to subject (increased activity)
        idx = find(behavior_labels == unqLabels(b)); %get idx where behavior b occurred
        %idx_rand = randsample(idx_rest,length(idx));
        %figure; hold on
        histogram(Spike_rasters_smooth(n, idx),'BinWidth',1,'Normalization','pdf', 'FaceColor',color_scheme{b})
        histogram(Spike_rasters_smooth(n, idx_rest),'BinWidth',1,'Normalization','pdf', 'FaceColor',[0.5 0.5 0.5])
        legend({behav_categ(b),'Rest'})
        title([behav_categ(b)])
        %     histogram(Spike_rasters(n, idx),20, 'FaceColor','r')
        %     histogram(Spike_rasters(n, idx_rand),20, 'FaceColor',[0.5 0.5 0.5])
        %     legend({'Aggression','Rest'})
        %     title('Firing rate during aggression vs. rest for example unit')
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
    idx_forage = find(behavior_labels == 5 | behavior_labels == 4);
    idx_agg = find(behavior_labels == 9 | behavior_labels == 10 |behavior_labels == 1 |behavior_labels == 21 |behavior_labels == 22);
    idx_travel = find(behavior_labels == 2 | behavior_labels == 12 |behavior_labels == 18);
    idx_scratch = find(behavior_labels == 23);

    unit=n;
    figure; hold on; set(gcf,'Position',[150 250 1500 500]);
    xline([block_times.end_time_round(1), block_times.end_time_round(2)], "-",["Block 1 end", "Block 2 end"], "LineWidth",2);
    session_with_buffer = [120:session_length-120];
    %Full spike rate trace
    plot(1:session_length, Spike_rasters_smooth(unit,1:session_length));

    %Rest
    to_plot=nan(1,session_length); to_plot(idx_rest)=Spike_rasters_smooth(unit,idx_rest);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0.5 0.5 0.5])

    %Getting groomed
    to_plot=nan(1,session_length); to_plot(idx_getgroom)=Spike_rasters_smooth(unit,idx_getgroom);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",'b')

    %Grooming
    to_plot=nan(1,session_length); to_plot(idx_groom)=Spike_rasters_smooth(unit,idx_groom);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",'c')

    %Self-groom
    to_plot=nan(1,session_length); to_plot(idx_selfgroom)=Spike_rasters_smooth(unit,idx_selfgroom);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0.8 0 0.8])

    %Foraging
    to_plot=nan(1,session_length); to_plot(idx_forage)=Spike_rasters_smooth(unit,idx_forage);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0 0.7 0])

    %Aggression
    to_plot=nan(1,session_length); to_plot(idx_agg)=Spike_rasters_smooth(unit,idx_agg);
    plot(1:session_length, to_plot, "LineWidth",2, "Color","r")

    %Travel
    to_plot=nan(1,session_length); to_plot(idx_travel)=Spike_rasters_smooth(unit,idx_travel);
    plot(1:session_length, to_plot, "LineWidth",2, "Color","y")

    %Scratch
    to_plot=nan(1,session_length); to_plot(idx_scratch)=Spike_rasters_smooth(unit,idx_scratch);
    plot(1:session_length, to_plot, "LineWidth",2, "Color",[0.8 0.6 0])

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