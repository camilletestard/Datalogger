%% Example code for how to fit MLE models to poission neuron data

%Running Notes: Will probably put this into SimNeuronChecks later, but
%separate for now for development.  Also this is just slapped together,
%will come back for more interesting simulations later on.

%% Set duration and other important parameters
%2022-03-02 NOTE: don't use all this yet
    %Timing parameters

    t_res = .001; % Start with milisecond resolution

    T = 300; %Start with just 5 minutes of recording

    t = 0:t_res:T;

    num_samples = length(t);


    %Behavior simulation parameters

    n_behaviors = 3;

    boi = 1:n_behaviors;

    unidis = 1; %For now simulate as if each behavior occurs with equal frequency to avoid having to subsample

    %Neuronal simulation parameters

    n_neurons = 5;

    as_process = 1; %Toggle whether to just simulate poisson neurons with a single lambda
    %Or to stimulate a poission process where lambdas are drawn for each frame
    %(think GLM).  Process will take longer to simulate I think

    fullrand = 0; %Toggle whether have neurons that are tuned to different catagories or not.

    %Binning parameters

    %100ms, 500ms, 1s, 2s, 5s

    binning = [.1 .5 1 2 5]*1/t_res; %convert time to indices
    
    %% A quite crap simulations
    %This doesn't work to illustrate the point I want yet...
    %Need to change how tuning is implemented.
    
    %Very basic situation, each behavior last 1/3 of session
    behav_vec = [ones(T/n_behaviors,1); ones(T/n_behaviors,1)*2; ones(T/n_behaviors,1)*3];
    
    spikes = nan(length(behav_vec),n_neurons);
    
    chunk = T/n_behaviors; %just for convenience for simulations
    
    gt_tuning = zeros(n_behaviors,n_neurons); %ground truth tuning
    
    fr_base = 10;
    tuning_bump = 5; %spike count bump to lambda (additive)
    
    for n = 1:n_neurons
        
        gt_tuning(randi(3),n) = 1; %Randomly tune neuron to one behavior
        
        for b = 1:n_behaviors
            
            lambda_temp = fr_base + tuning_bump*gt_tuning(b,n);
            
            gt_tuning(b,n) = lambda_temp; %save lambda value instead of 1 or 0
            
            %just being lazy otherwise this would be indexed by when
            %behavior occured
            spikes(1+chunk*(b-1):chunk*b,n) = poissrnd(gt_tuning(b,n),chunk,1);            
            
        end
        
        
    end
    
    expected_baseline = mean(gt_tuning); %expect baseline below to be average of underlying lambdas
    %lambda_hat_behavs below should be close to gt_tunning
    
    %% Actual method stuff
    %Fit Poisson
    %hat because we are getting it from the fit
    
   lambda_hat_baseline = nan(1,n_neurons); %Just doing whole session as baseline
   lambda_hb_ci = nan(2,n_neurons); %confidence interval for above
   lambda_hat_behav = nan(n_behaviors,n_neurons); %One for each behavior for each neuron
   lambda_hbeh_ci = nan(2,n_behaviors,n_neurons); %confidence interval for above
    for n = 1:n_neurons
        
        [lambda_hat_baseline(1,n), lambda_hb_ci(:,n)] = poissfit(spikes(:,n));
        
        for b = 1:n_behaviors
        [lambda_hat_behav(b,n), lambda_hbeh_ci(:,b,n)] = poissfit(spikes(1+chunk*(b-1):chunk*b,n));
        end
    
    end
    
    
  %Manual wald test (pulled from my previous work)
    
  %Review logic behind (i.e. my manuscript) this but go from 95% ci to se using below
  lambda_hb_se = diff(lambda_hb_ci)/2/1.96;
  lambda_hbeh_se = diff(lambda_hbeh_ci)/2/1.96; %Get width of confidence interval, scale back to se
  
  %Wiki https://en.wikipedia.org/wiki/Wald_test for formula until I explain
  %via manuscript.  sqW is wald statistic on asym z
  
  sqW = nan(size(lambda_hat_behav));
  
 %This almost certainly doesn't need to be for loops but I am a bit tired
 %to make something cute with the matrix math
 
 for n = 1:n_neurons
     for b = 1:n_behaviors
         
         sqW(b,n) = (lambda_hat_behav(b,n) - lambda_hat_baseline(1,n))./sqrt(lambda_hbeh_se(1,b,n).^2 + lambda_hb_se(1,n).^2);
         
     end
         
 end
  
  %Wald stat follows Z distribtution - need to check again the validity of
  %abs...remember that is in the text somewhere
  
  p_values = 1 - normcdf(abs(sqW),0,1);
  
  cutoff = 0.001; %set as desired
  

  gt_tuning
  lambda_hat_behav
  
  expected_baseline
  lambda_hat_baseline
 
  
  h = p_values < cutoff %return behaviors that are signficantly different from baseline
  
  figure
  imagesc(abs(sqW))
  %set(gca,'YDir','normal')
  xticks([1:n_neurons])
  xticklabels({'n1','n2','n3','n4','n5'})
  yticks([1:n_behaviors])
  yticklabels({'b1', 'b2', 'b3'})
  title('tuning per behavior')
    