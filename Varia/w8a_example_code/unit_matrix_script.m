%Practice script for making unit matrix.
% Info to know: word 12 = cue onset

%Load data

%Define your analysis period after target onset (in seconds)
Analysis_period_lenght = 5.5;

% Identify indices of hit trials
hits = find(BHV.TrialError == 0);

%Define your event of interest to which neural activity is aligned
for i = 1:length(hits) %for all hit trials
    Cue_onset_neural(i) = words_per_trials(hits(i)).timestamps(words_per_trials(hits(i)).words == 12); %Get the timestamp of the cue
end
Cue_onset_neural = Cue_onset_neural';

%Define the timestamps limits for the analysis period for each trial
Analysis_period = [Cue_onset_neural - 1.5 (Cue_onset_neural + 4)];

hWaitbar = waitbar(0, 'Creating raster matrix');

%Create structure with rasters per hit trials for all units
unit=1;
for i = 1:length(fields(SpikeData)) %For all channels
    waitbar(i/length(fields(SpikeData)), hWaitbar);
    
    for j = 2:length(SpikeData.(['Channel_' num2str(i)])) %For all units, except the unsorted cluster
        
        for k = 1:length(Cue_onset_neural) %For all hit trials
            
            Unit_rasters{unit}(k,:) = zeros(1,Analysis_period_lenght*1000); %Fill the line with zeros to initiate raster for that trial
            temp = find(SpikeData.(['Channel_' num2str(i)]){j} > Analysis_period(k,1) & SpikeData.(['Channel_' num2str(i)]){j} < Analysis_period(k,2)); %Find indicies of spikes during this analysis period
            ticks = SpikeData.(['Channel_' num2str(i)]){j}(temp) - Analysis_period(k,1); %Align those spikes to cue onset - 1500 msec for that trial
            ticks = ceil(ticks*1000); %Convert spikes timings (in raster time) to miliseconds
            Unit_rasters{unit}(k, ticks) = 1; %Fill in spikes in the raster
            clear ticks temp
        end
        Electrode_ID(unit) = i;
        unit = unit+1;
    end
    
end

close(hWaitbar)

clearvars -except BHV words_per_trials Unit_rasters Electrode_ID Analysis_period_lenght saving_folder session_name Session PathName neural_dir

disp(['Creation of the raster matrix done in ', num2str(toc), ' seconds.'])
