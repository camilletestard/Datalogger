%% log_crossValReg_groom

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
channel_flag = "all"; %Channels considered
with_NC =1; %0: NC is excluded; 1:NC is included; 2:ONLY noise cluster
randomsample=0;
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well isolated units
num_iter = 100;%Number of SVM iterations
smooth= 0; % 1: smooth the data; 0: do not smooth
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


    %% Get data with specified temporal resolution and channels
    [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, ...
        unit_count, groom_labels_all, brain_label, behavior_log, behav_categ_original]= ...
        log_GenerateDataToRes_function(filePath, temp_resolution, channel_flag, ...
        is_mac, with_NC, isolatedOnly, smooth, sigma, threat_precedence, exclude_sq);

    disp('Data Loaded')


    Spike_count_raster = Spike_rasters';
    session_length(s) = size(Spike_count_raster,1);
    behavior_labels = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
    behavior_labels(behavior_labels==find(behav_categ=="Proximity"))=length(behav_categ);
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

    %% Extract grooming variables

     %Initialize variables which will be used to compute the different
    %grooming metrics
    behav_lbls=behavior_labels;
    groom_behav_perSession = behavior_labels;
    groomGive_perSession = behavior_labels;
    groomGet_perSession = behavior_labels;
    total_time=ones(size(behavior_labels));

    %Groom give only
    groomGive_perSession(groom_behav_perSession~=7)=0; groomGive_perSession(groom_behav_perSession==7)=1;
    groomGiveBout_perSession=zeros(size(groomGive_perSession)); groomGiveBout_perSession(find(diff(groomGive_perSession)<0)+1)=1;
 
    %Groom receive only
    groomGet_perSession(groom_behav_perSession~=8)=0; groomGet_perSession(groom_behav_perSession==8)=1;
    groomGetBout_perSession=zeros(size(groomGive_perSession)); groomGetBout_perSession(find(diff(groomGet_perSession)<0)+1)=1;

    %Groom give or receive
    groom_behav_perSession(groom_behav_perSession~=8 & groom_behav_perSession~=7)=0;
    groom_behav_perSession(groom_behav_perSession==8)=1; groom_behav_perSession(groom_behav_perSession==7)=-1;
    groomBout_perSession=zeros(size(groomGive_perSession)); groomBout_perSession(groomGetBout_perSession==1)=1; groomBout_perSession(groomGiveBout_perSession==1)=-1;

    %Cumulative grooming in the session (+1 if you get groomed, -1 if you groom)
    cumul_groom_perSession = cumsum(groom_behav_perSession); figure; plot(cumul_groom_perSession)
    cumul_groomBout_perSession= cumsum(groomBout_perSession);  figure; plot(cumul_groomBout_perSession)

    %Total amount of grooming thus far in the session (cumulative)
    total_groom_perSession = cumsum(abs(groom_behav_perSession));  figure; plot(total_groom_perSession)
    total_groomBout_perSession = cumsum(abs(groomBout_perSession));  figure; plot(total_groomBout_perSession)

    %Total amount of groom give in the session thus far
    cumul_give_perSession =cumsum(groomGive_perSession); figure; plot(cumul_give_perSession)
    pGive_perSession=cumul_give_perSession./cumsum(total_time);figure; plot(pGive_perSession)
    cumul_giveBouts_perSession = cumsum(groomGiveBout_perSession); figure; plot(cumul_giveBouts_perSession)

    %Total amount of groom receive in the session thus far
    cumul_get_perSession =cumsum(groomGet_perSession); figure; plot(cumul_get_perSession)
    pGet_perSession=cumul_get_perSession./cumsum(total_time);figure; plot(pGet_perSession)
    cumul_getBouts_perSession = cumsum(groomGetBout_perSession); figure; plot(cumul_getBouts_perSession)
    %pGive_perSession=cumul_give_perSession./total_groom_perSession;figure; plot(pGive_perSession)

    %Total reciprocity of groom duration
    recip_groom_perSession = 1-(cumul_groom_perSession./total_groom_perSession);figure; plot(recip_groom_perSession)
    recip_groomBout_perSession = 1-(cumul_groomBout_perSession./total_groomBout_perSession);figure; plot(recip_groomBout_perSession)


    %% Predict observed grooming reciprocity from neural data

    neural_data = zscore(Spike_count_raster);
    groom_idx = ismember(behavior_labels,[7,8]);

    %Predict observed reciprocity from neural data
    X1= neural_data(groom_idx,:);
    %[~, score]= pca(neural_data(groom_idx,:)); X1=score(:,1:10);
    %y_obs = cumul_groom_obs(groom_idx)';
    y_obs = zscore(cumul_groom_obs(groom_idx)'); 
    %Test autocorrelation of variable you're trying to predict
    %[acf, lags]=autocorr(y_obs,'NumLags',100); %HIGHLY autocorrelated
    %y_obs = zscore(cumul_groom_slidingWindow_obs(groom_idx)');
    figure; plot(y_obs)
    %y_obs =y_obs(1:1500); X1=X1(1:1500,:);

%     %Estimate autoregression?
%     Mdl = arima(2,0,1);
%     EstMdl = estimate(Mdl,y_obs)


    %Ridge regression?
    [yfit]=log_crossValModel(X1, y_obs, [1:size(X1,2)], [1:size(X1,2)], [1:size(X1,2)], 5);


%     mdl = fitlm(X1,y_obs);
%     yfit = predict(mdl,X1); Rsq_obs = corr(y_obs, yfit).^2; 
    
    Rsq_obs = corr(y_obs, yfit').^2; %
    yfit_smooth = movmean(yfit,50);
    figure;hold on; plot(yfit); plot(yfit_smooth,'LineWidth',3); plot(y_obs,'LineWidth',3); legend('Prediction','Smoothed Prediction','Real')

    %Correlate observed reciprocity to individual firing rate 
    groom_give = nan(size(groom_all)); groom_give(groom_all==1)=2;
    groom_receive = nan(size(groom_all)); groom_give(groom_all==-1)=2.1;
    groom_give=groom_give(idx_all); groom_receive=groom_receive(idx_all);
  
    [~, pca_score]= pca(neural_data(idx_all,:)); X1=score(:,1:10);
    for n= 1:size(X1,2)%1:10:200
        raw_activity(:,n) = pca_score(:,n);
        smooth_activity(:,n) = movmean(raw_activity(:,n),25);
        smooth_response = movmean(y_obs,50);
        [correl_per_neuron(n) p(n)] = corr(raw_activity(:,n),y_obs);

%         figure; hold on; 
%         scatter(1:length(groom_give),groom_give,10,'b','filled'); 
%         scatter(1:length(groom_receive),groom_receive,10,'c','filled'); 
%         plot(smooth_activity(:,n));
    end
    figure; histogram(correl_per_neuron)
    n_of_interest = find(abs(correl_per_neuron)>0.3);
    for n=1:length(n_of_interest)
    figure; hold on; plot(smooth_activity(:,n_of_interest(n))); plot(y_obs); title(num2str(correl_per_neuron(n_of_interest(n))))
    end


    %Randomized
    for rand = 1:1000
        X1= neural_data(idx_all,:);
        y = cumul_groom(idx_all(randsample(length(idx_all), length(idx_all))))'; %randomize response variable
        mdl = fitlm(X1,y);
        yfit = predict(mdl,X1);
        %values = crossval(@regf,X1,y);
        MAE_rand(rand) = mean(abs(yfit-y)); %mean(values(:,1));%
        adjMAE_rand(rand) = MAE_obs/range(y); %mean(values(:,2));%
        Rsq_rand(rand) = corr(y, yfit).^2; %mean(values(:,3)); %
        disp(rand)
    end
    figure; histogram(adjMAE_rand(1:100)); xline(adjMAE_obs,'--r')
    figure; histogram(Rsq_rand(1:100)); xline(Rsq_obs,'--r')


    %Randomized, keeping statistics of the behavior
    clear MAE_rand adjMAE_rand Rsq_rand correl_shuffle_real
    rand_iter = 100;
    Rsq_rand =nan(1,rand_iter);
    for rand = 1:rand_iter
        groom_bouts_perSession = allGroomBouts_sorted;
        n_bouts = size(groom_bouts_perSession ,1);
        duration_bouts = round(random('Exponential',120, 1,n_bouts));
        groom_bouts_perSession(:,3)=duration_bouts;
%         groom_bouts_perSession(:,1)=sort(randsample(1:size(Spike_count_raster,1),n_bouts));
%         groom_bouts_perSession(:,2)=groom_bouts_perSession(:,1) + groom_bouts_perSession(:,3);
        groom_bouts_perSession(:,6)=randsample(groom_bouts_perSession(:,6),length(groom_bouts_perSession(:,6)),'true');
        groom_all = zeros(1,size(neural_data,1));
        for bout = 1:size(groom_bouts_perSession,1)

            idx_bout{bout} = groom_bouts_perSession(bout,1):groom_bouts_perSession(bout,2);
            if groom_bouts_perSession(bout,6) ==8
                groom_all(idx_bout{bout})=ones(1,length(idx_bout{bout}));
            else
                groom_all(idx_bout{bout}) =ones(1,length(idx_bout{bout}))*(-1);
            end
        end
        idx_all = cell2mat(idx_bout);
        cumul_groom{rand} = cumsum(groom_all); %figure; plot(cumul_groom{rand} (idx_all))
        correl_shuffle_real(rand)=corr(cumul_groom{rand} (idx_all)', y_obs);

        % % %         %Simulate neural data
        % % %         mean_Hz = mean(Spike_count_raster);
        % % %         for n=1:size(Spike_count_raster,2)
        % % %             sim_neural_data(:,n) = poissrnd(mean_Hz(n),1,size(Spike_count_raster,1));
        % % % %             [acf, lags]=autocorr(sim_neural_data(:,n));
        % % % %             [acf, lags]=autocorr(Spike_count_raster(:,n));
        % % %             %figure; hold on; plot(movmean(sim_neural_data(:,n),100)); plot(movmean(Spike_count_raster(:,n),100))
        % % %         end
        % % %         sim_neural_data = zscore(sim_neural_data);

        %if correl_shuffle_real(rand)<0.1
            %return
            X1= neural_data(idx_all,:);
            %X1= sim_neural_data(idx_all,:); %use simulated neural data
            %IMPORTANT NOTE: this simulated neural data does NOT have temporal
            %auto-correlation

            %[~, score]= pca(neural_data(idx_all,:)); X1=score(:,1:10);
            y = zscore(cumul_groom{rand} (idx_all)'); %[acf, lags]=autocorr(y);
            [yfit]=log_crossValModel(X1, y, [1:size(X1,2)], [1:size(X1,2)], [1:size(X1,2)], 5);
            Rsq_rand(rand) = corr(y, yfit').^2; %mean(values(:,3)); %
            disp(rand)
        %end
        %clear cumul_groom yfit y
    end
    figure; histogram(Rsq_rand(~isnan(Rsq_rand)),10); xline(Rsq_obs,'--r','LineWidth',2)
    %figure; histogram(Rsq_rand(abs(correl_shuffle_real)<0.2),10); xline(Rsq_obs,'--r')
    figure; histogram(abs(correl_shuffle_real));
    figure; scatter(abs(correl_shuffle_real),Rsq_rand); [rho p]=corr(abs(correl_shuffle_real)',Rsq_rand')

    for t = find(Rsq_rand>0.4)
        figure; 
        plot(cumul_groom{t})
    end




    yfit_smooth = movmean(yfit,10);
    figure;hold on; plot(yfit); plot(yfit_smooth,'LineWidth',3); plot(y,'LineWidth',3); legend('Prediction','Real')

    figure; histogram(MAE); xline(50,'--r')
    sum(MAE<51)./length(MAE)
end


