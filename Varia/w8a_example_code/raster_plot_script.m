
    
    %% Plot the rasters aligned to Cue onset
    close all
    tic
    
    %Identify hit trials
    hits = find(BHV.TrialError == 0);
    
    %Identify timing of relevant events (all timings are relative to cue onset on a trial by trial basis)
    cue_ticks = 1500; %Identify time of significant cue onset on raster (as specified above)
    for y = 1:length(hits) %identify other timings that are trial specific
        targets_onset(y) = (words_per_trials(hits(y)).timestamps(words_per_trials(hits(y)).words == 25) - words_per_trials(hits(y)).timestamps(words_per_trials(hits(y)).words == 12))*1000 +1500;
        target_selection(y) = (words_per_trials(hits(y)).timestamps(words_per_trials(hits(y)).words == 40)- words_per_trials(hits(y)).timestamps(words_per_trials(hits(y)).words == 12))*1000 +1500;
    end
    
    %Identify trials where the correct target was a star or a circle among the hit trials
    star = strfind(BHV.TimingFileByCond, 'visual_star.m'); star_conditions = find(not(cellfun('isempty', star)));
    circle = strfind(BHV.TimingFileByCond, 'visual_circle.m'); circle_conditions = find(not(cellfun('isempty', circle)));
    star_hit_trials = find(ismember(BHV.ConditionNumber(hits), star_conditions)); % indices of star trial among hit trials vector (hits)
    circle_hit_trials = find(ismember(BHV.ConditionNumber(hits), circle_conditions));
    
    %Find the indices of hit trials when targets onset time are organized in ascending order
    [~,index_star] = sort(targets_onset(star_hit_trials));
    [~,index_circle] = sort(targets_onset(circle_hit_trials));
    
    clear BHV

    hWaitbar = waitbar(0, ['Creating raster plots for each unit in ' session_name]);
    
    % Draw the raster plot figures for each neuron
    for i = 1:length(Unit_rasters) %For all units
        
        waitbar(i/length(Unit_rasters), hWaitbar);
        
        figure
        hold on
        
        subplot(2,1,1) %Star trials
        for j = 1:length(star_hit_trials) %For all hit trials with star as a target
            
            line([cue_ticks cue_ticks], [j-.3 j+.3], 'Color', 'g', 'LineWidth', 1.5) %Cue onset
            line([cue_ticks+1000 cue_ticks+1000], [j-.3 j+.3], 'Color', 'r', 'LineWidth', 1.5) %Cue Offset
            line([targets_onset(star_hit_trials(index_star(j))) targets_onset(star_hit_trials(index_star(j)))], [j-.3 j+.3], 'Color', 'g', 'LineWidth', 1.5)
            line([target_selection(star_hit_trials(index_star(j))) target_selection(star_hit_trials(index_star(j)))], [j-.3 j+.3], 'Color', 'r', 'LineWidth', 1.5)
            
            ticks = find(Unit_rasters{i}(star_hit_trials(index_star(j)),:)); %Find the index of this spike
            hold on; scatter(ticks, ones(1,length(ticks))*j, .5, 'b', 'filled')
            
        end
        title(['Star target, Neuron # '  num2str(i) ' , Channel # ' num2str(Electrode_ID(i))])
        xlabel('Bin # (1 msec)')
        ylabel('Trial #')
        ylim([0 length(star_hit_trials)+1])
        xlim([0 Analysis_period_lenght*1000])
        set(gca,'TickDir','out') % draw the tick marks on the outside
        
        subplot(2,1,2) %Circle trials
        for j = 1:length(circle_hit_trials) %For all hit trials with circle as a target
            
            line([cue_ticks cue_ticks], [j-.3 j+.3], 'Color', 'g', 'LineWidth', 1.5) %Cue onset
            line([cue_ticks+1000 cue_ticks+1000], [j-.3 j+.3], 'Color', 'r', 'LineWidth', 1.5) %Cue Offset
            line([targets_onset(circle_hit_trials(index_circle(j))) targets_onset(circle_hit_trials(index_circle(j)))], [j-.3 j+.3], 'Color', 'g', 'LineWidth', 1.5)
            line([target_selection(circle_hit_trials(index_circle(j))) target_selection(circle_hit_trials(index_circle(j)))], [j-.3 j+.3], 'Color', 'r', 'LineWidth', 1.5)
            
            ticks = find(Unit_rasters{i}(circle_hit_trials(index_circle(j)),:));  %Find the index of this spike
            hold on; scatter(ticks, ones(1,length(ticks))*j, .5, 'b', 'filled')
            
        end
        title(['Circle target, Neuron # '  num2str(i) ' , Channel # ' num2str(Electrode_ID(i))])
        xlabel('Bin # (1 msec)')
        ylabel('Trial #')
        ylim([0 length(circle_hit_trials)+1])
        xlim([0 Analysis_period_lenght*1000])
        set(gca,'TickDir','out') % draw the tick marks on the outside
        
        
        % Saving the plots
        obj = findobj('type', 'figure');
        
        %Adjust paths
        if ~exist([saving_folder '/' session_name '/Neuron_'  num2str(i)], 'dir');
            mkdir([saving_folder '/' session_name '/Neuron_'  num2str(i)])
        end
        cd([saving_folder '/' session_name '/Neuron_'  num2str(i)])
        saveas(obj, ['Neuron_'  num2str(i) '_Raster_aligned_RT'], 'tif')
        close all
        
    end %End of loop over all units in this session
    
    close(hWaitbar)
    
    disp(['Raster plots done in ', num2str(toc), ' seconds for session ' session_name])
    
    clearvars -except saving_folder Session PathName neural_dir

end %end of iteration across sessions
