%% Ron_DoSmth
% Have fun!

close all

%Set parameters
is_mac = 0; %For loading the data
is_ron = 1; %For setting path
with_partner =1;
temp_resolution = 1; %Temporal resolution of firing rate. 1sec
BRs = ["TEO", "vlPFC"]; %Channels considered (sets brain region TEO vlPFC all
with_NC =1; %0: Noise cluster (NC) is excluded; 1:NC is included; 2:ONLY noise cluster
isolatedOnly=0; %Only consider isolated units. 0=all units; 1=only well-isolated units


for br = 1:2 %Do separately by brain region to see if there is a difference
    
    channel_flag = BRs(br)
    
    subsample = 1; %Toggle whether to subsample
    
    
   

    for s = [1,2,3] %For now focus on first Amost sessions
        Results(br).sesh(s).brainregion = channel_flag;
        s
        %Set path
        if is_ron

            sessions = dir('C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_data'); sessions = sessions(3:end,:);
            filePath = ['C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_data\' sessions(s).name];
            savePath = ['C:\Users\ronwd\OneDrive\Documents\GitHub\Datalogger_results'];

        else


             %Set session list
        home = '~'; % set home directory
        cd([home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/'])
        sessions = dir('Ready to analyze output'); sessions = sessions(5:end,:);
        filePath = [home '/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/' sessions(s).name]; % Enter the path for the location of your Deuteron sorted neural .nex files (one per channel)
        savePath = [home '/Dropbox (Penn)/Datalogger/Results/' sessions(s).name '/SingleUnit_results/partner_vs_subject'];
        end
        %% Load data

        %Get data with specified temporal resolution and channels
        [Spike_rasters, labels, labels_partner, behav_categ, block_times, monkey, reciprocal_set, social_set, ME_final,unit_count, groom_labels_all]= log_GenerateDataToRes_Ron(filePath, temp_resolution, channel_flag, is_mac,is_ron, with_NC, isolatedOnly);
        % Output:
        %       1. Spike_rasters: n_neuron X time in seconds (if temp_resolution = 1)
        %       2. labels: time X 11 - check the function description for
        %       details. Essentially subject behavior is column 3.
        %       3. labels_partner: time X 11. Partner behavior is column 3.
        %       4. behav_categ: string vector; behavioral category do the numbers in
        %       labels correspond to.
        %       5. block_times: specifies the block times and ID
        %       6. monkey: subjecet monkey ID
        %       7. reciprocal_set: set of behaviors that are deeemed "reciprocal"
        %       8. social_set: set of behaviors that are deeemed "social"
        %       9. ME_final: Motion energy values
        %       10. unit_count: number of units in the session
        %       11. groom_labels_all: labels for grooming. Check function
        %       description for details

        session_length = size(Spike_rasters,2); % get session length
        
        if subsample
            
            n_2_use = 100; %for now just set using 100 neurons
            
            Spike_rasters = Spike_rasters(randperm(size(Spike_rasters,1),n_2_use),:);
            
            
        end
        
        Spike_count_raster = Spike_rasters';

        %Extract behavior labels for subject and partner
        behavior_labels_subject_init = cell2mat({labels{:,3}}'); %Extract unique behavior info for subject
        behavior_labels_partner_init = cell2mat({labels_partner{:,3}}'); %Extract unique behavior info for partner

        % Set proximity as rest
        behavior_labels_subject_init(behavior_labels_subject_init==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").
        behavior_labels_partner_init(behavior_labels_partner_init==find(behav_categ=="Proximity"))=length(behav_categ); %exclude proximity for now (i.e. mark as "undefined").

    

    %% Subdivide data for successive dimensionality tests
    %Going to just run successive tests based on what was outlined in the pitch
    %Can come back to clean this up to run more efficiently but for now just
    %want to confirm for this session that general idea is right.

    %Going to start with using PCA for speed and clarity (i.e. everyone gets
    %PCA)

    %Before doing anything else just get the inds for all the groups we want.

    groominds_give = groom_labels_all(:,1) == 7;
    groominds_receive = groom_labels_all(:,1) == 8;

    groominds_all = groominds_give + groominds_receive; %Since these are mutually exclusive and (as far as I can tell) define all groom times.
    groominds_all = groominds_all>0; %To change it back to logical indexing
    groominds_postthreat = groom_labels_all(:,3) == 1;
    groominds_reciprocal = groom_labels_all(:, 4) == 1;
    groominds_initiated = groom_labels_all(:,5) == 1;

    groominds_self = behavior_labels_subject_init == 24; %24 is label for self groom

    %Use Spike_count_raster as it is in obs x vars format these algorithms are
    %expecting

    %% Test 0.1: Basic checks - Dimensionality for grooming, the rest, and all
    %So far this does not match what I was predicting.  The deviation isn't
    %massive, but it is still surprising.  Doesn't entirely sink this train of
    %thought, so will contibue with other test, but need to reflect on what
    %this means.

    threshold = 50;

    [eigenvec, score, eigenval, ~, var_exp] = pca(Spike_count_raster);

    Results(br).sesh(s).all.eigenvec = eigenvec;
    Results(br).sesh(s).all.score = score;
    Results(br).sesh(s).all.eigenval = eigenval;
    Results(br).sesh(s).all.var_exp = var_exp;


    Data_grooming = Spike_count_raster(groominds_all,:);
    Data_self = Spike_count_raster(groominds_self,:);
    Data_groomgive = Spike_count_raster(groominds_give,:);
    Data_groomrec = Spike_count_raster(groominds_receive,:);


    [eigenvec, score, eigenval, ~, var_exp] = pca(Data_grooming);

    Results(br).sesh(s).groomall.eigenvec = eigenvec;
    Results(br).sesh(s).groomall.score = score;
    Results(br).sesh(s).groomall.eigenval = eigenval;
    Results(br).sesh(s).groomall.var_exp = var_exp;

    [eigenvec, score, eigenval, ~, var_exp] = pca(Data_groomgive);

    Results(br).sesh(s).groomgive.eigenvec = eigenvec;
    Results(br).sesh(s).groomgive.score = score;
    Results(br).sesh(s).groomgive.eigenval = eigenval;
    Results(br).sesh(s).groomgive.var_exp = var_exp;

    [eigenvec, score, eigenval, ~, var_exp] = pca(Data_groomrec);

    Results(br).sesh(s).groomrec.eigenvec = eigenvec;
    Results(br).sesh(s).groomrec.score = score;
    Results(br).sesh(s).groomrec.eigenval = eigenval;
    Results(br).sesh(s).groomrec.var_exp = var_exp;

    [eigenvec, score, eigenval, ~, var_exp] = pca(Data_self);

    Results(br).sesh(s).therest.eigenvec = eigenvec;
    Results(br).sesh(s).therest.score = score;
    Results(br).sesh(s).therest.eigenval = eigenval;
    Results(br).sesh(s).therest.var_exp = var_exp;

    %Get number of eigen values needed to cross threshold of variance explained

    Results(br).sesh(s).all.majdim = find(cumsum(Results(br).sesh(s).all.var_exp)>threshold,1,'first');
    Results(br).sesh(s).therest.majdim = find(cumsum(Results(br).sesh(s).therest.var_exp)>threshold,1,'first');
    Results(br).sesh(s).groomall.majdim = find(cumsum(Results(br).sesh(s).groomall.var_exp)>threshold,1,'first');
    Results(br).sesh(s).groomgive.majdim = find(cumsum(Results(br).sesh(s).groomgive.var_exp)>threshold,1,'first');
    Results(br).sesh(s).groomrec.majdim = find(cumsum(Results(br).sesh(s).groomrec.var_exp)>threshold,1,'first');

    %sorted_axis = 1:length(var_exp);%changed mind don't like this./length(var_exp); %quick trick for setting axes for eigen plots
    figure
    plot(1:length(Results(br).sesh(s).all.var_exp),cumsum(Results(br).sesh(s).all.var_exp),'k')
    title([channel_flag ' ' num2str(s)])
    ylabel('var explained')
    xlabel('Eigenvectors used')
    ylim([0,110])
    hold on
    plot(1:length(Results(br).sesh(s).groomall.var_exp),cumsum(Results(br).sesh(s).groomall.var_exp),'b')
    plot(1:length(Results(br).sesh(s).groomgive.var_exp),cumsum(Results(br).sesh(s).groomgive.var_exp),'b--')
    plot(1:length(Results(br).sesh(s).groomrec.var_exp),cumsum(Results(br).sesh(s).groomrec.var_exp),'b.-')
    plot(1:length(Results(br).sesh(s).therest.var_exp),cumsum(Results(br).sesh(s).therest.var_exp),'r')
    plot(1:length(var_exp), threshold * ones(size(var_exp)),'k--') %put dashed black line for majority of variance
    %number beside is number of PCs needed to cross threshold var explained
    legend(['All: ' num2str(Results(br).sesh(s).all.majdim)],...
        ['Grooming: ' num2str(Results(br).sesh(s).groomall.majdim)],...
        ['Groom Give: ' num2str(Results(br).sesh(s).groomgive.majdim)],...
        ['Groom Receive:' num2str(Results(br).sesh(s).groomrec.majdim)],...
        ['Groom Self: ' num2str(Results(br).sesh(s).therest.majdim)],...
        'Threshold')



    %% Test 1.1 Compare dimensionality: post threat, initiated, reciprocal

    %Dimensionality is effectively the same across these conditions as well as
    %compared to considering all grooming.  Kind of interesting to see; need to
    %think more on what it means

    Data_groompt = Spike_count_raster(groominds_postthreat,:);
    Data_groomrecip = Spike_count_raster(groominds_reciprocal,:);
    Data_groomini = Spike_count_raster(groominds_initiated,:);

    [eigenvec, score, eigenval, ~, var_exp] = pca(Data_groompt);

    Results(br).sesh(s).groompt.eigenvec = eigenvec;
    Results(br).sesh(s).groompt.score = score;
    Results(br).sesh(s).groompt.eigenval = eigenval;
    Results(br).sesh(s).groompt.var_exp = var_exp;

    [eigenvec, score, eigenval, ~, var_exp] = pca(Data_groomrecip);

    Results(br).sesh(s).groomrecip.eigenvec = eigenvec;
    Results(br).sesh(s).groomrecip.score = score;
    Results(br).sesh(s).groomrecip.eigenval = eigenval;
    Results(br).sesh(s).groomrecip.var_exp = var_exp;

    [eigenvec, score, eigenval, ~, var_exp] = pca(Data_groomini);

    Results(br).sesh(s).groomini.eigenvec = eigenvec;
    Results(br).sesh(s).groomini.score = score;
    Results(br).sesh(s).groomini.eigenval = eigenval;
    Results(br).sesh(s).groomini.var_exp = var_exp;


    Results(br).sesh(s).groompt.majdim = find(cumsum(Results(br).sesh(s).groompt.var_exp)>threshold,1,'first');
    Results(br).sesh(s).groomrecip.majdim = find(cumsum(Results(br).sesh(s).groomrecip.var_exp)>threshold,1,'first');
    Results(br).sesh(s).groomini.majdim = find(cumsum(Results(br).sesh(s).groomini.var_exp)>threshold,1,'first');


    figure
    plot(1:length(Results(br).sesh(s).groomall.var_exp),cumsum(Results(br).sesh(s).groomall.var_exp),'b')
    title([channel_flag ' ' num2str(s)])
    ylabel('var explained')
    xlabel('Eigenvectors used')
    hold on
    plot(1:length(Results(br).sesh(s).groompt.var_exp),cumsum(Results(br).sesh(s).groompt.var_exp),'bx')
    plot(1:length(Results(br).sesh(s).groomrecip.var_exp),cumsum(Results(br).sesh(s).groomrecip.var_exp),'b--')
    plot(1:length(Results(br).sesh(s).groomini.var_exp),cumsum(Results(br).sesh(s).groomini.var_exp),'b.-')

    plot(sorted_axis, threshold * ones(size(sorted_axis)),'k--') %put dashed black line for majority of variance
    %number beside is number of PCs needed to cross threshold var explained
    legend(['Grooming: ' num2str(Results(br).sesh(s).groomall.majdim)],...
        ['Groom Post Threat: ' num2str(Results(br).sesh(s).groompt.majdim)],...
        ['Groom Reciprocal: ' num2str(Results(br).sesh(s).groomrecip.majdim)],...
        ['Groom Initated: ' num2str(Results(br).sesh(s).groomini.majdim)],...
        'Threshold')

    %% Test 1.2 Grooming lies in different subspaces based on context: post threat, initiated, reciprocal
    %First do simple test with average angle between PCs
    %Next as there was more to that paper than I thought and my first tries are
    %not exactly working as expected; implement the analyses from the below
    %paper to see what that gets us and if some of their processing steps are
    %necessary to see this structure.

    %Not eigenvectors come out of PCA in matlab already normalized to have unit
    %vector norm (i.e. norm(eigenvec(:,1)) = 1

    %Just use maximum number need to cross threshold of variance explained from above.
    num_eigenvec = max([Results(br).sesh(s).groompt.majdim,Results(br).sesh(s).groomrecip.majdim, Results(br).sesh(s).groomini.majdim]); 


    collected_eigens = cell(3,1);

    collected_eigens{1} = Results(br).sesh(s).groompt.eigenvec(:,1:num_eigenvec);
    collected_eigens{2} = Results(br).sesh(s).groomrecip.eigenvec(:,1:num_eigenvec);
    collected_eigens{3} = Results(br).sesh(s).groomini.eigenvec(:,1:num_eigenvec);

    angles_between = cell(3,1); %n choose k with 3 is 3!/2! (3-2)! = 3
    avg_angle = nan(3,1);

    comparison = [1 2; 1 3; 2 3]; %better ways to do this but for now just manually put in all of the comparisons we want to do

    % for comp = 1:length(angles_between)
    %     
    %     angles_between{comp} = acos(dot(collected_eigens{comparison(comp,1)},collected_eigens{comparison(comp,2)}));
    %     
    %     avg_angle(comp) = rad2deg(mean(angles_between{comp}));
    %     
    %     
    % end
    % 
     some_labels = {'post-threat to reciprocal','post-threat to initiated','reciprocal to initated'}';
    % 
    % table(avg_angle,some_labels,'VariableNames',{'avg angle', 'comparison'})

    %Now trying the principal angles definition that is used in the below paper
    %and seems to be more standard in the field given some reading.  Take the
    %SVD of the product of the two subspaces

    angles_betweensvd = cell(3,1);
    principal_angle = nan(size(angles_betweensvd));

    for comp = 1:length(angles_between)

        [~,S,~]= svd(collected_eigens{comparison(comp,1)}' * collected_eigens{comparison(comp,2)}); %want to do SVD on PC x PC matrix, hence transpose

        angles_betweensvd{comp} = rad2deg(acos(diag(S))); %singular vlaues of this matrix relate to the angles

        principal_angle(comp) =     angles_betweensvd{comp}(end); %svd in matlab does descending order apparently.  Larget SV is 1st principal angle
    end

    t1 = table(principal_angle,some_labels,'VariableNames',{'principal angle', 'comparison'})

    Results(br).sesh(s).groomcontext.pa = t1;
    %% Test 1.3 sanity check on overlap between grooming and other behaviors
    %Replicating above analysis for grooming and groom self
    
    

    %% Test(s) 2-? Replicating Geometry of sequence working memory in macaque prefrontal cortex analyses
    %Our analyses don't exactly align due to the difference in our behavior vs
    %their task, but going to try to more or fully implement what they did and see
    %where it gets us.

    end
end