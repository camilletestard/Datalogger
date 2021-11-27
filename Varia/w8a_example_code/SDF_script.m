%% w8a_SpikeDensityFunction_perTarget_batch
% This is a batch implementation of w8a_SpikeDensityFunction_perTarget
% This script produces spike density function plots that pool trials according to target identity (star or circle).
% The script produces the raster matrix for selected blocks and convolves the raster to produce the SDF.
% The SDF is aligned to three events of interest (Cue onset, Target onset, and Response) and concatenated into a single SDF.
% Created by Testard and Tremblay 02/2017

    
    %Identify hit trials
        %Identify the trial numbers part of each block
        Trials_block1 = find(BHV.BlockNumber == block(1));
        Trials_block2 = find(BHV.BlockNumber == block(2));
        %Identify the hit trials in each block
        hit_block1 = hits(ismember(hits, Trials_block1));
        hit_block2 = hits(ismember(hits, Trials_block2));
        %Pool the hits from both blocks (Corrupts chronology of trials. Shouldn't matter for SDF)
    
    %Get the timestamp of cue and target onset in the neural data for each hit trial (Not in chronological order)
    for i = 1:length(hits)
        Cue_onset_neural(i) = words_per_trials(hits(i)).timestamps(words_per_trials(hits(i)).words == 12);
        Targets_onset(i) = words_per_trials(hits(i)).timestamps(words_per_trials(hits(i)).words == 25);
        Response(i) = words_per_trials(hits(i)).timestamps(words_per_trials(hits(i)).words == 40);
    end
    Cue_onset_neural = Cue_onset_neural';
    Targets_onset = Targets_onset';
    Response = Response';
    
    %% Create matrix for rasters aligned to three events
    %Define your analysis periods (in seconds)
    Analysis_period_lenght = 5000; %msec
    Cue_window_length = 3000; %msec
    Target_window_length = 750;
    Response_window_length = 1250;
    
    %Define the timestamps limits for the analysis period for each trial
    Analysis_period_cue = [(Cue_onset_neural - 1.5) (Cue_onset_neural + 1.5)];
    Analysis_period_target = [(Targets_onset - .5) (Targets_onset + .25)];
    Analysis_period_response = [(Response - .25) (Response + 1)];
    
    %Create structure with rasters per trials for all units
    hWaitbar = waitbar(0, 'Creating Raster Matrix');
    unit=1;
    for i = 1:length(fields(SpikeData)) %For all channels
        waitbar(i/length(fields(SpikeData)), hWaitbar)
        
        for j = 2:length(SpikeData.(['Channel_' num2str(i)])) %For all units
            
            for k = 1:length(hits) %For all hit trials
                
                Unit_rasters{unit}(k,:) = zeros(1,Analysis_period_lenght); %Fill the line with zeros to initiate raster for that trial
                temp_cue = find(SpikeData.(['Channel_' num2str(i)]){j} > Analysis_period_cue(k,1) & SpikeData.(['Channel_' num2str(i)]){j} < Analysis_period_cue(k,2)); %Find indicies of spikes during this analysis period
                tick_cue = SpikeData.(['Channel_' num2str(i)]){j}(temp_cue) - Analysis_period_cue(k,1); %Align those spikes to cue onset -1500ms for that trial
                tick_cue = ceil(tick_cue*1000); %Convert spikes timings (in raster time) to miliseconds
                Unit_rasters{unit}(k, tick_cue) = 1; %Fill in spikes in the raster
                temp_target = find(SpikeData.(['Channel_' num2str(i)]){j} > Analysis_period_target(k,1) & SpikeData.(['Channel_' num2str(i)]){j} < Analysis_period_target(k,2)); %Find indicies of spikes during this analysis period
                tick_target = SpikeData.(['Channel_' num2str(i)]){j}(temp_target) - Analysis_period_target(k,1); %Align those spikes to target - 500ms for that trial
                tick_target = ceil(tick_target*1000) + Cue_window_length; %Convert spikes timings (in raster time) to miliseconds, add the cue analysis period length so as to concatenate both matrices.
                Unit_rasters{unit}(k, tick_target) = 1; %Fill in spikes in the raster
                temp_response = find(SpikeData.(['Channel_' num2str(i)]){j} > Analysis_period_response(k,1) & SpikeData.(['Channel_' num2str(i)]){j} < Analysis_period_response(k,2)); %Find indicies of spikes during this analysis period
                tick_response = SpikeData.(['Channel_' num2str(i)]){j}(temp_response) - Analysis_period_response(k,1); %Align those spikes to response - 250ms for that trial
                tick_response = ceil(tick_response*1000) + Cue_window_length + Target_window_length; %Convert spikes timings (in raster time) to miliseconds, add the cue analysis period length so as to concatenate both matrices.
                Unit_rasters{unit}(k, tick_response) = 1; %Fill in spikes in the raster
                clear tick_cue temp_cue tick_target temp_target temp_response tick_response
            end
            Electrode_ID(unit) = i;
            unit = unit+1;
        end
        
    end
    
    close(hWaitbar)
    
    %% Compute SDF
    %Define kernel
    sigma = .045; %Define width of kernel (in sec)
    edges = -3*sigma:.001:3*sigma;
    kernel = normpdf(edges,0,sigma); %Use a gaussian function
    kernel = kernel*.001; %Time 1/1000 so the total area under the gaussian is 1
    
    % Compute Spike Density Function for all hit trials
    hWaitbar = waitbar(0, 'Convolving rasters');
    for i = 1:length(Unit_rasters) % for all units
        waitbar(i/length(Unit_rasters), hWaitbar)
        for j = 1:size(Unit_rasters{1,i},1) % for all hit trials
            
            sdf = conv(Unit_rasters{1,i}(j,:),kernel); %Convolve the gaussian
            center = ceil(length(edges)/2);
            SDF{1,i}(j,:) = sdf(center:length(sdf)-(center-1)).*1000; %Substract the irrelevant edges due to convolution operation
            clear sdf
            
        end
    end
    close(hWaitbar)
    clear Unit_rasters

    
    %Identify trials where the correct target was a star or a circle among the hit trials
    star = strfind(BHV.TimingFileByCond, 'visual_star.m'); star_conditions = find(not(cellfun('isempty', star)));
    circle = strfind(BHV.TimingFileByCond, 'visual_circle.m'); circle_conditions = find(not(cellfun('isempty', circle)));
    star_hit_trials = find(ismember(BHV.ConditionNumber(hits), star_conditions)); % indices of star trials in the hit trials vector (hits)
    circle_hit_trials = find(ismember(BHV.ConditionNumber(hits), circle_conditions));
    
    % Compute Mean SDF for correct target star and correct target circle
    for i = 1:length(SDF) %For all units
        Mean_SDF{1,i} = mean(SDF{1,i}(star_hit_trials,:));
        Error_SDF{1,i} = std(SDF{1,i}(star_hit_trials,:))/(size(SDF{1,i}(star_hit_trials,:),1))^(1/2);
        Mean_SDF{2,i} = mean(SDF{1,i}(circle_hit_trials,:));
        Error_SDF{2,i} = std(SDF{1,i}(circle_hit_trials,:))/(size(SDF{1,i}(circle_hit_trials,:),1))^(1/2);
    end
    
    
    %% Plotting
    hWaitbar = waitbar(0, ['Creating SDF plots for each ' num2str(length(SDF)) ' unit']);
    x = 1:length(Mean_SDF{1,1}); % Define x axis for plotting
    C = {'g', 'b'}; % Define colors for SDF, green for star, blue for circle
    
    for unit = 1:length(SDF) % For all units
        waitbar(unit/length(SDF), hWaitbar);
        
        for j = 1:size(Mean_SDF,1) % For all task conditions
            
            figure(unit);
            hold on
            
            %Reference timings
                line([1500 1500], [0 max(max(Mean_SDF{j,unit}))+max(max(Mean_SDF{j,unit}))/10], 'Color', 'g') % mark cue onset
                line([2500 2500], [0 max(max(Mean_SDF{j,unit}))+max(max(Mean_SDF{j,unit}))/10], 'Color', 'g') % mark cue offset
                line([3000 3000], [0 max(max(Mean_SDF{j,unit}))+max(max(Mean_SDF{j,unit}))/10], 'Color', 'k', 'LineStyle', '--') % mark delay axis break
                line([3500 3500], [0 max(max(Mean_SDF{j,unit}))+max(max(Mean_SDF{j,unit}))/10], 'Color', 'r') % mark target onset
                line([3750 3750], [0 max(max(Mean_SDF{j,unit}))+max(max(Mean_SDF{j,unit}))/10], 'Color', 'k', 'LineStyle', '--') % mark reaction time axis break
                line([4000 4000], [0 max(max(Mean_SDF{j,unit}))+max(max(Mean_SDF{j,unit}))/10], 'Color', 'c') % mark response onset
            
            %Plot SDF
            fill([x fliplr(x)], [Mean_SDF{j,unit}(1,:)+Error_SDF{j,unit}(1,:), fliplr(Mean_SDF{j,unit}(1,:)-Error_SDF{j,unit}(1,:))], C{j})
            plot(Mean_SDF{j,unit}(1,:), 'r')
            
        end
            
        title(['Star (green) and Circle (blue) hit trials, Neuron # '  num2str(unit) ' , Channel # ' num2str(Electrode_ID(unit))])
        xlabel('Time (msec)')
        ylabel('Firing rate (Hz)')
        ylim([-Inf Inf])
        xlim([500 Analysis_period_lenght-200])
        set(gca,'TickDir','out') % draw the tick marks on the outside

        pause(1)
        close all
        
    end %End of loop over all units in this session
    
    close(hWaitbar)
    
    disp(['SDF plots done in ', num2str(toc), ' seconds for session ' session_name])
    
    clearvars -except saving_folder Session PathName neural_dir correct_drift

