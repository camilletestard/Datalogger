%% Log_CreateTrajectoryVideo

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