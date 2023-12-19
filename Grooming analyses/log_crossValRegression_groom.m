%% log_crossValReg_groom
% Run cross-validated regression to predict reciprocity from neural data

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


%Set parameters
with_partner =0;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
channel_flag = 'all'; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
randomsample=0;
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 100;%Number of SVM iterations
smooth= 1; % 1: smooth the data; 0: do not smooth
sigma = 1*temp_resolution;%set the smoothing window size (sigma)
null=0;%Set whether we want the null
simplify=0;%lump similar behavioral categories together to increase sample size.
threat_precedence =1;
exclude_sq=1;

%Select session range:
if with_partner ==1
    session_range = session_range_with_partner;
    %     a_sessions = 1:6; h_sessions = [11:13,15:16,18];
else
    session_range = session_range_no_partner;
    a_sessions = 1:6; h_sessions = [11:13,15:16,18];
end

s=1;
for s =session_range %1:length(sessions)

    %Set path
    filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
    savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SVM_results/'];


    chan=1;
    for channel_flag = "all"%["vlPFC", "TEO", "all"]

    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')


    Spike_count_raster = Spike_rasters';
    session_length(s) = size(Spike_count_raster,1);
    behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels(behavior_labels==find(behav_categ=='Proximity'))=length(behav_categ);
    block_labels = cell2mat({labels{:,12}}');


    %% Interpolate short groom bouts
    groomGive = zeros(size(behavior_labels));
    groomGive(behavior_labels== 7)=1;
    groomGive_time = find(groomGive==1);
    time_between_bouts = diff(groomGive_time);
    short_Runs = find(time_between_bouts>1 & time_between_bouts<11);

    for sr =1:length(short_Runs)
        idx_to_fill = groomGive_time(short_Runs(sr))+1:groomGive_time(short_Runs(sr))+time_between_bouts(short_Runs(sr))-1;
        groomGive(idx_to_fill)=1;
    end

    behavior_labels(find(groomGive==1))=7;

    %Groom get
    groomGet = zeros(size(behavior_labels));
    groomGet(behavior_labels== 8)=1;
    groomGet_time = find(groomGet==1);
    time_between_bouts = diff(groomGet_time);
    short_Runs = find(time_between_bouts>1 & time_between_bouts<11);

    for sr =1:length(short_Runs)
        idx_to_fill = groomGet_time(short_Runs(sr))+1:groomGet_time(short_Runs(sr))+time_between_bouts(short_Runs(sr))-1;
        groomGet(idx_to_fill)=1;
    end

    behavior_labels(find(groomGet==1))=8;

    %     behavior_labels_tosave{1,s} =behavior_labels';
    %     block_labels_tosave{1,s} =block_labels';

    %% Extract grooming stats (bout length and# per session)

    %GROOM GIVE
    groomGive_bout_start = find(diff(groomGive)==1)+1;
    groomGive_bout_end = find(diff(groomGive)==-1);
    if length(groomGive_bout_end)<length(groomGive_bout_start) %can happen if grooming went until very end of session
        groomGive_bout_end(length(groomGive_bout_start))=length(groomGive);
    end
    groomGive_duration = groomGive_bout_end-groomGive_bout_start;
    groomGive_bout=[groomGive_bout_start, groomGive_bout_end, groomGive_duration];

    for bout = 1:length(groomGive_bout_end)
        %Check groom is followed by a threat
        if groomGive_bout_end(bout)+11<session_length(s)
            groomGive_bout(bout,4) = any(behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+11)==9 ...
                |behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+11)==10);
        else
            groomGive_bout(bout,4) =0;
        end

        %Check if grooming bout was preceded by a threat
        if groomGive_bout_start(bout)-30>0
            groomGive_bout(bout,5) = any(behavior_labels(groomGive_bout_start(bout)-30:groomGive_bout_start(bout)-1)==9 ...
                |behavior_labels(groomGive_bout_start(bout)-30:groomGive_bout_start(bout)-1)==10);
        else
            groomGive_bout(bout,5)=0;
        end

        % % % %          %Check what the groom is followed by
        % % % %         if groomGive_bout_end(bout)+120<session_length(s)
        % % % %             unique(behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+120),'stable')
        % % % %             groomGive_bout(bout,6) = behavior_labels(groomGive_bout_end(bout)+1:groomGive_bout_end(bout)+120);
        % % % %         else
        % % % %             groomGive_bout(bout,6) =0;
        % % % %         end
    end
    groomGive_bout(:,6)=7;
    %4th column: was there a threat right after the groom (which could have
    %cut it short)
    %5th column: was there a threat preceding the grooming bout
    %6th column, grooming receive or give.

    %cut too short duratiin bouts (usually occurs because smth external
    %cuts the bout (aggression)
    %groomGive_bout(groomGive_bout(:,3)<10)

    %GROOM GET
    groomGet_bout_start = find(diff(groomGet)==1)+1;
    groomGet_bout_end = find(diff(groomGet)==-1);
    groomGet_duration = groomGet_bout_end-groomGet_bout_start;
    groomGet_bout=[groomGet_bout_start, groomGet_bout_end, groomGet_duration];

    for bout = 1:length(groomGet_bout_end)
        if groomGet_bout_end(bout)+11<session_length(s)
            groomGet_bout(bout,4) = any(behavior_labels(groomGet_bout_end(bout)+1:groomGet_bout_end(bout)+11)==9 ...
                |behavior_labels(groomGet_bout_end(bout)+1:groomGet_bout_end(bout)+11)==10);
        else
            groomGet_bout(bout,4)=0;
        end
        if groomGet_bout_start(bout)-30>0
            groomGet_bout(bout,5) = any(behavior_labels(groomGet_bout_start(bout)-30:groomGet_bout_start(bout)-1)==9 ...
                |behavior_labels(groomGet_bout_start(bout)-30:groomGet_bout_start(bout)-1)==10);
        else
            groomGet_bout(bout,5) =0;
        end
    end
    groomGet_bout(:,6)=8;

    % ALL GROOMS
    allGroomBouts = [groomGet_bout; groomGive_bout];
    [~, idx_sorted] = sort(allGroomBouts(:,1));
    allGroomBouts_sorted = allGroomBouts(idx_sorted,:);
    %Get inter-bout interval
    allGroomBouts_sorted(:,7)=[0;(allGroomBouts_sorted(2:end,1)-allGroomBouts_sorted(1:end-1,2))];
    %Is previous bout grooming in the other direction
    allGroomBouts_sorted(:,8) = abs([0; diff(allGroomBouts_sorted(:,6))]);
    %Mark turn-taking bouts
    allGroomBouts_sorted(:,9)=0;
    allGroomBouts_sorted(find(allGroomBouts_sorted(:,8)==1 & allGroomBouts_sorted(:,7)<20),9)=1;
    allGroomBouts_sorted(:,10)=s;
    %7th column: inter-bout interval
    %8th column: alternation of grooming bout
    %9th column: is grooming bout turn-taking

    obs_groom = allGroomBouts_sorted;
    obs_recip(s) = 1-((sum(obs_groom(obs_groom(:,6)==7,3)) - sum(obs_groom(obs_groom(:,6)==8,3)))./ sum(obs_groom(:,3)));
    obs_recip_bout(s) = 1-((sum(obs_groom(:,6)==7) - sum(obs_groom(:,6)==8))./ size(obs_groom,1));


    %% Remove very short bouts
    short_bouts = find(allGroomBouts_sorted(:,3)<10);
    for b=1:length(short_bouts)
        behavior_labels(allGroomBouts_sorted(short_bouts(b),1):allGroomBouts_sorted(short_bouts(b),2))=length(behav_categ);
    end

    %% Extract grooming variables

    %Initialize variables which will be used to compute the different
    %grooming metrics
    behav_lbls=behavior_labels;
    groom_behav = behavior_labels;
    groomGive = behavior_labels;
    groomGet = behavior_labels;
    total_time=ones(size(behavior_labels));
    cumul_total_time=cumsum(total_time);

    %Groom give only
    groomGive(groom_behav~=7)=0; groomGive(groom_behav==7)=1;
    groomGiveBout=zeros(size(groomGive)); groomGiveBout(find(diff(groomGive)==1)+1)=1;

    %Groom receive only
    groomGet(groom_behav~=8)=0; groomGet(groom_behav==8)=1;
    groomGetBout=zeros(size(groomGive)); groomGetBout(find(diff(groomGet)==1)+1)=1;

    %Time in bout
    time_in_bout = zeros(size(groomGive));
    bout_number = zeros(size(groomGive));
    for b = 1:size(allGroomBouts_sorted,1)
        idx = allGroomBouts_sorted(b,1) : allGroomBouts_sorted(b,2);
        time_in_bout(idx) = 1:length(idx);
        bout_number(idx) = ones(1,length(idx))*b;
    end

    %Groom give or receive
    groom_behav(groom_behav~=8 & groom_behav~=7)=0;
    groom_behav(groom_behav==8)=1; groom_behav(groom_behav==7)=-1;
    groomBout=zeros(size(groomGive)); groomBout(groomGetBout==1)=1; groomBout(groomGiveBout==1)=-1;

    %Cumulative grooming in the session (+1 if you get groomed, -1 if you groom)
    cumul_groomTime = cumsum(groom_behav)*0.5; %figure; plot(cumul_groom)
    cumul_groomBout= cumsum(groomBout);  %figure; plot(cumul_groomBout)

    %Total amount of grooming thus far in the session (cumulative)
    total_groom = cumsum(abs(groom_behav));  %figure; plot(total_groom)
    total_groomBout = cumsum(abs(groomBout));  %figure; plot(total_groomBout)

    %Total amount of groom give in the session thus far
    cumul_give =cumsum(groomGive); %figure; plot(cumul_give) % Total time grooming partner so far
    %pGive=cumul_give./cumsum(total_time);%figure; plot(pGive) % Total proportion grooming partner so far
    cumul_giveBouts = cumsum(groomGiveBout); %figure; plot(cumul_giveBouts) % Total number of grooming partner bouts so far
    pGive=cumul_give./total_groom;%figure; plot(pGive)


    %Total amount of groom receive in the session thus far
    cumul_GroomGetTime =cumsum(groomGet); %figure; plot(cumul_get)
    %pGet=cumul_get./cumsum(total_time);%figure; plot(pGet)
    cumul_GroomGetBouts = cumsum(groomGetBout); %figure; plot(cumul_getBouts)
    pGet=cumul_GroomGetTime./total_groom;%figure; plot(pGet)


    %Total reciprocity of groom duration
    recip_groom = 1-(cumul_groomTime./total_groom);%figure; plot(recip_groom)
    recip_groomBout = 1-(cumul_groomBout./total_groomBout);%figure; plot(recip_groomBout)

    %Set categorical variables
    cumul_groomTime_categ = discretize(cumul_groomTime,5); %figure; plot(cumul_groomTime_categ)
    cumul_groomBout_categ = discretize(cumul_groomBout,5); %figure; plot(cumul_groomBout_categ)
    time_in_bout_categ = discretize(time_in_bout,5); %figure; plot(time_in_bout_categ)

    % All groom variables
    % % %     groom_variables=[cumul_groom, cumul_groomBout, ...
    % % %         total_groom,total_groomBout,...
    % % %         cumul_give, pGive, cumul_giveBouts,...
    % % %         cumul_get, pGet, cumul_getBouts,...
    % % %         recip_groom, recip_groomBout, cumul_total_time];

    % % %     groom_var_names={'Cumulative relative groom time', 'Cumulative relative groom bout',...
    % % %         'Total groom time', 'Total number of groom bouts',...
    % % %         'Cumulative groom give time', 'Proportion of groom give time', 'Cumulative groom give bouts',...
    % % %         'Cumulative groom receive time', 'Proportion of groom receive time', 'Cumulative groom receive bouts',...
    % % %         'Equality groom time', 'Equality groom bouts','Time'}

    % Only groom variagbles which do not necessarily increase monotonically
    % with time
    % % %     groom_variables=[cumul_groom, cumul_groomBout, ...
    % % %         pGive, pGet, ...
    % % %         recip_groom, recip_groomBout, cumul_total_time];
    % % %
    % % %     groom_var_names={'Cumulative relative groom time', 'Cumulative relative groom bout',...
    % % %         'Proportion of groom give time','Proportion of groom receive time', ...
    % % %         'Equality groom time', 'Equality groom bouts','Total groom time'};

    time_boutNum_combo = cumul_groomBout.*cumul_total_time;
    timInBout_num_combo = cumul_groomBout.*time_in_bout;

    %Final set of variables
%     groom_variables=[cumul_groomTime, cumul_groomBout, ...
%         cumul_total_time];
    groom_variables=[cumul_groomTime, cumul_groomBout, ...
        cumul_total_time, time_boutNum_combo, timInBout_numold_combo];
%     groom_variables=[cumul_groomTime_categ, cumul_groomBout_categ, ...
%         time_in_bout_categ];

%     groom_var_names={'Cumulative relative groom time', 'Cumulative relative groom bout',...
%         'Time in session'};
    groom_var_names={'Cumulative relative groom time', 'Cumulative relative groom bout',...
        'Time in session', 'Time X NetBout', 'TimeInBout x Net Bout'};


% %     figure; 
% %     %plot(cumul_groomBout)
% %     plot(cumul_groomBout.*time_in_bout)


    %% Predict observed grooming reciprocity from neural data

    figure(s); hold on
    for var=1:length(groom_var_names)


        neural_data = zscore(Spike_count_raster);
        groom_var = groom_variables(:,var);
        %idx_of_interest = ismember(behavior_labels,length(behav_categ)); %Rest
        idx_of_interest = ismember(behavior_labels,[7,8]); %Groom 
        %idx_of_interest = 1:length(behavior_labels); %All session


        %% Set up variables for regression
        X1= neural_data(idx_of_interest,:);
        %[~, score]= pca(neural_data(groom_idx,:)); X1=score(:,1:10);
        if var==6 %Don't z-score groom bout reciprocity because it yields all nans
            y_obs{s,var} = groom_var(idx_of_interest)';
        else
            y_obs{s,var} = zscore(groom_var(idx_of_interest)');
        end
        %figure; plot(y_obs{s,var})

        %Remove nans if any
        if any(isnan(y_obs{s,var}))
            non_nan_idx = ~isnan(y_obs{s,var});
            y_obs{s,var}=y_obs{s,var}(non_nan_idx); X1=X1(non_nan_idx,:);
            %             y_est=y_est(non_nan_idx);
        end


        %% fit a line to the observed data to measure how linear/monotonic
        %the variable is
        c = polyfit(1:length(y_obs{s,var}),y_obs{s,var},1); y_est = polyval(c,1:length(y_obs{s,var}));
        %figure; hold on; plot(y_obs,'LineWidth',2); plot(y_est, '--r')
        mae_linear_approx(s,var) = mean(abs(y_est-y_obs{s,var}));
        rsq_linear_approx(s,var) = corr(y_obs{s,var}', y_est').^2;


        %% Predict observed reciprocity from neural data
        %Fit ridge regression with 5-fold cross-validation (in blocks)
        [yfit]=log_crossValModel_groom(X1, y_obs{s,var}, [1:size(X1,2)], [1:size(X1,2)], [1:size(X1,2)], 5);
        Rsq_obs(s,var) = corr(y_obs{s,var}', yfit').^2; % Get rsquared
        residuals{s,var}= y_obs{s,var}' - yfit';
        yfit_plot{s,var} = yfit;

        % % %         %Fit ridge regression to linear appoximation
        % % %         [yfit_line]=log_crossValModel(X1, y_obs, [1:size(X1,2)], [1:size(X1,2)], [1:size(X1,2)], 5);
        % % %         Rsq_LinearApprox(s,var) = corr(y_est', yfit_line').^2; % Get rsquared

% % %                 %Plot fit
% % %                 if var==2
% % %                     gap= mean(yfit(1:500))- mean(y_obs{s,var}(1:500));
% % %                     yfit = yfit-gap;
% % %                 end
% % %                 yfit_smooth = movmean(yfit,50);
% % %                 subplot(2,3,var);hold on; plot(yfit); plot(y_obs{s,var},'LineWidth',3);
% % %                 if var==1
% % %                     legend('Prediction','Smoothed Prediction','Real','Location','best')
% % %                 end
% % %                 title([groom_var_names{var} ', cvRsq = ' num2str(Rsq_obs(s,var))])
% % %                 %disp([groom_var_names{var} ', cvRsq = ' num2str(Rsq_obs(s,var))])


        %% Predict observed reciprocity from each neuron individually
% % %         for n=1:size(X1,2)
% % % 
% % %             %Shuffle all except neuron of interest
% % %             X1_Full_Contrib = X1;
% % %             idx_to_shuffle = find(~ismember(1:size(X1,2), n));
% % %             for id=1:length(idx_to_shuffle)
% % %                 X1_Full_Contrib(:,id)=X1(randperm(size(X1,1)), id);
% % %             end
% % %             %Fit ridge regression with 5-fold cross-validation (in blocks)
% % %             %to single neuron model (full contribution)
% % %             [yfit_per_neuron]=log_crossValModel(X1_Full_Contrib, y_obs{s,var}, [1:size(X1,2)], [1:size(X1,2)], [1:size(X1,2)], 5);
% % %             FullRsq_per_neuron{s}(var,n) = corr(y_obs{s,var}', yfit_per_neuron').^2; % Get rsquared
% % % 
% % %             %Shuffle all except neuron of interest
% % %             X1_Unq_Contrib = X1;
% % %             idx_to_shuffle = find(ismember(1:size(X1,2), n));
% % %             for id=1:length(idx_to_shuffle)
% % %                 X1_Unq_Contrib(:,id)=X1(randperm(size(X1,1)), id);
% % %             end
% % %             %Fit ridge regression with 5-fold cross-validation (in blocks)
% % %             %to single neuron model (full contribution)
% % %             [yfit_per_neuron]=log_crossValModel(X1_Unq_Contrib, y_obs{s,var}, [1:size(X1,2)], [1:size(X1,2)], [1:size(X1,2)], 5);
% % %             UnqRsq_per_neuron{s}(var,n) = abs(Rsq_obs(s,var) - corr(y_obs{s,var}', yfit_per_neuron').^2); % Get rsquared
% % % 
% % %             if rem(n,10) == 0
% % %                 disp(['neuron ' num2str(n) '/' num2str(size(X1,2)) 'done.'])
% % %             end
% % % 
% % % 
% % %         end
        % %         %         figure; histogram(FullRsq_per_neuron{s}(var,:)); [maxFit idx_neuron]=max(FullRsq_per_neuron{s}(var,:)); figure; hold on; plot(y_obs); plot(movmean(X1(:,idx_neuron),50))
        % %         %         figure; histogram(UnqRsq_per_neuron{s}(var,:)); [maxFit idx_neuron]=max(UnqRsq_per_neuron{s}(var,:)); figure; hold on; plot(y_obs); plot(movmean(X1(:,idx_neuron),50))
        % %
        disp(['Variable:' num2str(var) '/' num2str(length(groom_var_names)) ' done.'])
    end

% % % % %     %Model explanation plots
% % % % %     figure;
% % % % %     plot(cumul_groomTime); ylim([-400 400])
% % % % %     figure
% % % % %     plot(cumul_groomBout); ylim([-3 3]); %xlim([0 6000])
% % % % %     figure; 
% % % % % %     cmap = [1 1 1; 0.8 0.8 0.8; 0 0 1];
% % % % % %     idx_plot = ones(size(idx_of_interest)); idx_plot(ismember(behavior_labels,length(behav_categ)))=2; idx_plot(idx_of_interest)=3;
% % % % %     cmap = [0.8 0.8 0.8; 0 0.5 1];
% % % % %     idx_plot = ones(size(idx_of_interest)); idx_plot(idx_of_interest)=2;
% % % % %     imagesc(idx_plot'); colormap(cmap)
% % % % %     
% % % % %     figure; 
% % % % %     plot(cumul_total_time)
% % % % %     figure; hold on; ni=1:10:10*12;c=1;
% % % % %     for n=1:9:100
% % % % %         plot(ni(c)+neural_data(:,n))
% % % % %         c=c+1;
% % % % %     end
% % % % %     ylim([0 120])


    chan=chan+1;
    disp(['Channel:' num2str(chan) ' done.'])
    
    end %end of channel loop
   
    disp(['Session ' num2str(s) ' done.'])

end

cd('~/Dropbox (Penn)/Datalogger/Results/All_sessions/Grooming_results/')
% save('GroomVar_Reg.mat','Rsq_obs','rsq_linear_approx','FullRsq_per_neuron','UnqRsq_per_neuron','session_range','a_sessions','h_sessions','groom_var_names','residuals')
% load('GroomVar_Reg.mat')

Rsq_obs(Rsq_obs==0)=NaN;
obs_recip(obs_recip ==0)=NaN;
obs_recip_bout(obs_recip_bout ==0)=NaN;

%Investigate the relationship between being monotonic and linear fit
var=2; groom_var_names{var}
figure; hold on
scatter(rsq_linear_approx(:,var), Rsq_obs(:,var),'filled')
xlabel('How monotonic (linear) the variable is')
ylabel('cvRsq')
[rho p]=corr(rsq_linear_approx(:,var), Rsq_obs(:,var),'rows','complete');[rho p]

soi=session_range; %find(Rsq_obs(:,var)<0.05);
for s=1:length(soi)%session_range
    figure; plot(y_obs{soi(s),var}); title(num2str(soi(s)))
end

% % %Investigate relationship between being equal in time and fit
% % var=3; groom_var_names{var}
% % figure; hold on
% % scatter(obs_recip, Rsq_obs(:,var))
% % xlabel('How monotonic (linear) the variable is')
% % ylabel('cvRsq')
% % [rho p]=corr(obs_recip', Rsq_obs(:,var),'rows','complete');[rho p]
% % 
% % %Investigate relationship between being equal in bouts and fit
% % var=2; groom_var_names{var}
% % figure; hold on
% % scatter(obs_recip_bout, Rsq_obs(:,var))
% % xlabel('How monotonic (linear) the variable is')
% % ylabel('cvRsq')
% % [rho p]=corr(obs_recip_bout', Rsq_obs(:,var),'rows','complete');[rho p]


%Plot results
chan=3;
results = squeeze(Rsq_obs(session_range,:));
[~,idx_sorted]=sort(nanmedian(results));
results_sorted = results(:,idx_sorted);

figure;
x=[ones(1,size(results_sorted,1)); 2*ones(1,size(results_sorted,1)); 3*ones(1,size(results_sorted,1))];
boxchart(results_sorted)
%swarmchart(x', results_sorted)
xticklabels(groom_var_names(idx_sorted))
ylabel('cvRsq')
ylim([0 1])

figure; 
data = cell2mat(FullRsq_per_neuron(session_range))';
x = [ones(1,size(data,1)); 2*ones(1,size(data,1)); 3*ones(1,size(data,1))];
swarmchart(x', data)
ylabel('cvRsq per neuron')
xticks([1,2,3])
xticklabels({'Relative groom duration', 'Relative number of bouts', 'Time in session'})
median(data)
mean(data)

figure; hold on
histogram(data(:,1))
histogram(data(:,2))
histogram(data(:,3))
n=3105

%Compute p-val comparing groom and non-groom epochs
[h, p]=ttest(groom_epochs_results(:,3), results_sorted(:,3))
%Relative # bouts p= 0.0058
%Time in session p=0.014
[h p]=ttest(Rsq_obs(:,1), Rsq_obs(:,2))
%Diff between two types of models p=0.0027

