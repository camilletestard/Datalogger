
%% Load data

filePath = uigetdir('', 'Please select the experiment directory'); % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
cd(filePath)

session = filePath(end-9:end);
monkey = filePath(end-14:end-10);

%behav_log = readtable('behavioral_log_session1.csv');
behavior_log = readtable(['EVENTLOG_restructured_',monkey,session,'.csv']);% Behavioral data
behav_log = readtable('behavioral_log_session1.csv');
behavior_log = readtable('behav_log_restructured_session1.csv');% Behavioral data
behavior_log{:,'start_time_round'}=round(behavior_log{:,'start_time'});
behavior_log{:,'stop_time_round'}=round(behavior_log{:,'stop_time'});
behavior_log{:,'duration_round'}=behavior_log{:,'stop_time_round'}-behavior_log{:,'start_time_round'};

load('Raster_per_second_Allchannels.mat') % Neural data
% load('SpikeData_timestamps.mat')


%% Plot SDF (intergrated signal at the second resolution)
% Check what the neural signal looks like

figure; hold on
for neuron = 1:size(Unit_rasters,1)
    plot(Unit_rasters(neuron,:)+ones(1,size(Unit_rasters,2))*300*neuron)
end
pause(1); close all

%For now exclude unit 4 and 37 which have a weird activity
included_units = [1:3 5:36 38:59];
Unit_rasters = Unit_rasters(included_units,:);

%% Neural trajectory by wakefulness state

%Get wakefulness events:
wakingup = round(behav_log{strcmp(behav_log{:,'Behavior'},'wakingup'),'Time'}/1000); %time in sec
awake = round(behav_log{strcmp(behav_log{:,'Behavior'},'awake'),'Time'}/1000); %time in sec
fullyawake = round(behav_log{strcmp(behav_log{:,'Behavior'},'fullyawake'),'Time'}/1000); %time in sec

%Create input structure for Data High:
D =struct();
D(1).type = 'traj';
D(1).data = Unit_rasters;
D(1).epochStarts=[1 wakingup awake ];%fullyawake];
D(1).epochColors=[1,0,0;0,1,0;0,0,1];%0,1,1];

cd('C:\Users\Camille Testard\Desktop\Grants_Fellowships\R37_Grant_renewal\DataHigh1.3')
DataHigh(D,'DimReduce')

%% Create trajectory video

%I manually extracted the trajectory to have total plotting control.
%Load 3D trajectory coordinates:
load('trajectory_allUnits.mat')

%get new epoch starts
load('trajectory_allUnits_info.mat')
D.epochStarts = handles.D.epochStarts; close all
%set frame rate
window_size = 5; step = 1;
no_of_windows = max(1 + floor((size(p, 2) - window_size) / step), 0);
plotting_windows=buffer([1:size(p,2)], window_size, step);
plotting_windows(plotting_windows==0)=1;
plotting_windows(find(plotting_windows(:,end)==1),end)=max(max(plotting_windows));
C=[0.8 0.8 0.8; 0.7 0 0.7;0.9 0.6 0.1];

%Get new times for behaviors
ratio = size(p,2)/size(Unit_rasters,2);
behavior_log.new_start_time = round(behavior_log.start_time*ratio); behavior_log.new_start_time(1)=1;
behavior_log.new_stop_time = round(behavior_log.stop_time*ratio);
behavior_log.new_duration=behavior_log.new_stop_time - behavior_log.new_start_time;
all_Colors = [1 0 0; 1 0 1; 0.8 1 0.6; 0 0 1; 0 1 1; 0.8 0.8 1; 0 0.6 0.4; 0.6 0.6 0.6; 1 0.5 0; 0.6 0.8 0];

for row_num = 1:length(behavior_log.new_start_time)
    if strcmp(behavior_log.Behavior(row_num), 'Aggression')
        behavior_log.Color(row_num) = {[1 0 0]};
    elseif strcmp(behavior_log.Behavior(row_num),'Drinking')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Feeding')
        behavior_log.Color(row_num) = {[0.8 1 0.6]};
    elseif strcmp(behavior_log.Behavior(row_num),'Grooming')
        behavior_log.Color(row_num) = {[0 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Proximity')
        behavior_log.Color(row_num) = {[0 1 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Resting')
        behavior_log.Color(row_num) = {[0.8 0.8 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Self-directed behavior')
        behavior_log.Color(row_num) = {[0 0.6 0.4]};
    elseif strcmp(behavior_log.Behavior(row_num),'Sleeping')
        behavior_log.Color(row_num) = {[0.6 0.6 0.6]};
    elseif strcmp(behavior_log.Behavior(row_num),'Submission')
        behavior_log.Color(row_num) = {[1 0.5 0]};
    else
        behavior_log.Color(row_num) = {[0.6 0.8 0]};
    end
end

% % figure; hold on; ylim([1 2]); xlim([0 455])
% % set(gca,'Xticklabel',[]); set(gca, 'XTick', []);
% % set(gca,'Yticklabel',[]); set(gca, 'YTick', []);
% % set(gca,'Zticklabel',[]); set(gca, 'ZTick', []);
% % for b = 1:length(behavior_log.new_start_time)
% %     rectangle('Position',[behavior_log.new_start_time(b),1,...
% %         behavior_log.new_stop_time(b)-behavior_log.new_start_time(b),1], 'FaceColor',behavior_log.Color{b},...
% %         'EdgeColor',behavior_log.Color{b})
% % end
% % %Just to make a legend...
% % for c=1:size(all_Colors,1)
% %     scatter(0,380,1,all_Colors(c,:),'filled')
% % end
% % legend({'Aggression','Drinking','Feeding','Grooming','Proximity',...
% %     'Resting','Self-groooming','Sleeping','Submission','Walking'},'Location','East')
% % set(gcf, 'Position',  [100, 100, 1000, 200])

vidfile = VideoWriter('TrajectoryMovie.mp4','MPEG-4');
open(vidfile);
for k = 1:size(p,2)%size(plotting_windows,2)
    
    disp([num2str(k) '/' num2str(size(p,2))])
    figure; hold on; set(gcf, 'Position',  [100, 100, 1000, 600])
    
    ethogramplot = subplot(2,1,1); hold on; ylim([1 2]); xlim([0 372])
    set(gca,'Xticklabel',[]); set(gca, 'XTick', []);
    set(gca,'Yticklabel',[]); set(gca, 'YTick', []);
    set(gca,'Zticklabel',[]); set(gca, 'ZTick', []);
    ethogramplot.Position = ethogramplot.Position + [0 0.08 0 -0.05];
    
    behav2plot=find(behavior_log.new_start_time<=k);
    
    for i = 1:length(behav2plot)
        b=behav2plot(i);
        rectangle('Position',[behavior_log.new_start_time(b),1,...
            k-behavior_log.new_start_time(b)+1,1], 'FaceColor',behavior_log.Color{b},...
            'EdgeColor',behavior_log.Color{b})
    end
    %Just to make a legend...
    xline(D.epochStarts(2),'--',{'WakingUp'},'LineWidth',2,'Color',C(2,:));
    xline(D.epochStarts(3),'--',{'Awake'},'LineWidth',2,'Color',C(3,:));
    for c=1:size(all_Colors,1)
        scatter(0,380,0.1,all_Colors(c,:),'filled')
    end
    legend({'Waking up','Awake','Aggression','Drinking','Feeding','Grooming','Proximity',...
        'Resting','Self-groooming','Sleeping','Submission','Walking'},'Location','eastoutside')
    t=title('Monkey Behavior','FontSize', 16);
    set(t,'position',get(t,'position')-[130 0 0])
    
    % Neural trajectory plot
    trajplot=subplot(2,1,2); hold on;
    trajplot.Position = trajplot.Position + [0 -0.1 0 0.3];
    xlim([min(p(1,:)), max(p(1,:))]); set(gca,'Xticklabel',[])
    ylim([min(p(2,:)), max(p(2,:))]); set(gca,'Yticklabel',[])
    zlim([min(p(3,:)), max(p(3,:))]); set(gca,'Zticklabel',[])
    view(40,30)
    grid on
    
    if k<=D.epochStarts(2)
        line(p(1,1:k), p(2,1:k),...
            p(3,1:k), 'LineWidth', 2, 'Color', C(1,:));
        
        legend('Sedated','Location','SouthEast')
    end
    if k>D.epochStarts(2)
        line(p(1,1:D.epochStarts(2)), p(2,1:D.epochStarts(2)),...
            p(3,1:D.epochStarts(2)), 'LineWidth', 2, 'Color', C(1,:));
        line(p(1,D.epochStarts(2):k), p(2,D.epochStarts(2):k),...
            p(3,D.epochStarts(2):k), 'LineWidth', 2, 'Color', C(2,:));
        legend('Sedated', 'Waking up','Location','SouthEast')
    end
    if k>D.epochStarts(3)
        line(p(1,D.epochStarts(3)+1:k), p(2,D.epochStarts(3)+1:k),...
            p(3,D.epochStarts(3)+1:k), 'LineWidth', 2, 'Color', C(3,:));
        legend('Sedated', 'Waking up', 'Awake','Location','SouthEast')
    end
    t=title('Neural Trajectory','FontSize', 16);
    set(t,'position',get(t,'position')-[650 0 800])
    
    frame = getframe(gcf);
    writeVideo(vidfile,frame);
    
    close all
end

close(vidfile);

%% Neural states during "awake" state. 1- Compare active behaviors states

% Start with "long behaviors":
min_duration = 5;%in sec
max_duration = 100;%in sec % exclude long behaviors as they may not be "clean". I.e. they co-occur with a number of other behaviors.
idx = find(behavior_log{:,'duration_round'}>=min_duration &...
    behavior_log{:,'duration_round'}<=max_duration &...
    ~strcmp(behavior_log{:,'Behavior'},'Resting')); %exclude resting for now
new_log=behavior_log(idx,:);

% Create structure with chunked neural data per event types.
time_around_epoch = 1;%in sec
event_window=min_duration;%in sec
B = struct();
pseudo_trial = 1;
for behav = 1:size(new_log,1)
    
    total_time = new_log{behav,'start_time_round'}-time_around_epoch ...
        :new_log{behav,'stop_time_round'}+time_around_epoch;
    
    divide=floor(length(total_time)./event_window);
    mode=mod(length(total_time),event_window);
    new_chunks = reshape(total_time(1:end-mode),divide,[]);
    for ii = 1:divide
        
        B(pseudo_trial).data = Unit_rasters(:,new_chunks(ii,:));
        
        B(pseudo_trial).condition = char(new_log{behav,'Behavior'});
        %     if strcmp(new_log{behav,'Behavior'},'Resting')
        %         B(behav).condition = 'Non-Resting';
        %     else
        %         B(behav).condition = 'Resting';
        %     end
        
        B(pseudo_trial).epochStarts=1;
        if strcmp(new_log{behav,'Behavior'},'Aggression')
            B(pseudo_trial).epochColors=[1,0,0];
        elseif strcmp(new_log{behav,'Behavior'},'Grooming')
            B(pseudo_trial).epochColors=[0,0,1];
        elseif strcmp(new_log{behav,'Behavior'},'Proximity')
            B(pseudo_trial).epochColors=[0,0.75,1];
        elseif strcmp(new_log{behav,'Behavior'},'Self-directed behavior')
            B(pseudo_trial).epochColors=[0,1,0];
        elseif strcmp(new_log{behav,'Behavior'},'Submission')
            B(pseudo_trial).epochColors=[1,0.75,0];
        elseif strcmp(new_log{behav,'Behavior'},'Walking')
            B(pseudo_trial).epochColors=[1,0,1];
        elseif strcmp(new_log{behav,'Behavior'},'Resting')
            B(pseudo_trial).epochColors=[0,0.5,1];s
        elseif strcmp(new_log{behav,'Behavior'},'Feeding')
            B(pseudo_trial).epochColors=[0.5,0.5,0.5];
        elseif strcmp(new_log{behav,'Behavior'},'Drinking')
            B(pseudo_trial).epochColors=[0.5,0,1];
        end
        B(pseudo_trial).type='state';
        
        pseudo_trial=pseudo_trial+1;
    end
end

% Average "chunked" trials (for better visualization):
table_B = struct2table(B); avg_B=struct();
tabulate(table_B{:,'condition'});
behav_list = unique(table_B{:,'condition'});
num_events=4; avg_state=1;
for b = 1:length(behav_list)
    
    idx=find(ismember(table_B{:,'condition'},behav_list{b}));
    divide=floor(length(idx)./num_events);
    mode=mod(length(idx),num_events);
    average_over = reshape(idx(1:end-mode),divide,[]);
    
    for avg = 1:size(average_over,1)
        avg_B(avg_state).data = mean(cat(3,B(average_over(avg,1)).data,B(average_over(avg,2)).data,...
            B(average_over(avg,3)).data,B(average_over(avg,4)).data),3);
        avg_B(avg_state).condition = behav_list{b};
        avg_B(avg_state).epochStarts = 1;
        avg_B(avg_state).epochColors = B(average_over(avg,1)).epochColors;
        avg_B(avg_state).type='state';
        
        avg_state=avg_state+1;
    end
end

% Call DataHigh
DataHigh(avg_B,'DimReduce')

%% Scrap code
for row_num = 1:length(labels)
    if strcmp(behavior_log.Behavior(row_num), 'Aggression')
        behavior_log.Color(row_num) = {[1 0 0]};
    elseif strcmp(behavior_log.Behavior(row_num),'Approach')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Drinking')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Foraging')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Groom Give')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Groom Receive')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'HIP')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'HIS')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Leave')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Other monkeys vocalize')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Pacing/Travel')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Proximity')
        behavior_log.Color(row_num) = {[0.8 1 0.6]};
    elseif strcmp(behavior_log.Behavior(row_num),'RR')
        behavior_log.Color(row_num) = {[0 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'SP')
        behavior_log.Color(row_num) = {[1 0 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'SS')
        behavior_log.Color(row_num) = {[0 1 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'SP')
        behavior_log.Color(row_num) = {[0.8 0.8 1]};
    elseif strcmp(behavior_log.Behavior(row_num),'Scratch')
        behavior_log.Color(row_num) = {[0 0.6 0.4]};
    elseif strcmp(behavior_log.Behavior(row_num),'Self-groom')
        behavior_log.Color(row_num) = {[0.6 0.6 0.6]};
    elseif strcmp(behavior_log.Behavior(row_num),'Vocalization')
        behavior_log.Color(row_num) = {[1 0.5 0]};
    else
        behavior_log.Color(row_num) = {[0.6 0.8 0]};
    end
end