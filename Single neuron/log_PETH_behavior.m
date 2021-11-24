%% log_PETH_behavior
% Firing rate aligned to onset or offset of behavior

filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)

session = filePath(end-9:end);
monkey = filePath(end-14:end-10);

%behav_log = readtable('behavioral_log_session1.csv');
behavior_log = readtable(['EVENTLOG_restructured_',monkey,session,'.csv']);% Behavioral data
load('Neural_data.mat') % Neural data; array1 is in TEO and array2 is in vlPFC
session_length = size(Unit_rasters,2);
n_neurons = size(Unit_rasters,1);

%Preprocessing: round times in behavioral log
behavior_log{:,'start_time_round'}=round(behavior_log{:,'start_time'});
behavior_log{:,'end_time_round'}=round(behavior_log{:,'end_time'});
behavior_log{:,'duration_round'}=behavior_log{:,'end_time_round'}-behavior_log{:,'start_time_round'};

for unit = 1:n_neurons

    events = find(strcmp('Groom Give', behavior_log{:,'Behavior'}));
    event_times = behavior_log{events,'start_time_round'};
    intervals = [event_times-4 event_times+4];
    figure; hold on
    for e = 1:length(event_times)
        PETH(e,:) = Unit_rasters(unit,intervals(e,1):intervals(e,2));
        plot(PETH(e,:), 'Color',[0, 0, 1, 0.1])
    end
    plot(mean(PETH,1), LineWidth=5), 'Color',[0, 0.8, 1, 0.2]; ylabel('Firing rate (Hz)'); title(['Event onset Unit #' num2str(unit)])
    line([5 5],[0, max(PETH(:))], 'LineStyle','--', 'Color','red')

    pause(1)

    close all
    
end