mu1 = 15; % average firing rate of subpolutation 1
mu2 = 2; % average firing rate of subpopulation 2
frac1 = 0.5; % average fraction of neurons in subpopulation 1
gamma_shape = 5; % parameter controlling distribution of rates in each population (higher = more gaussian, lower (with minimum 1) = more exponential)

n = 100;  % num neurons
t = 2000; % num trials


figure;
subplot(1,2,1); hold on;
subplot(1,2,2); hold on;

for run=1:10
    rates = zeros(1, n);
    % Make some simulated neurons, each with a particular rate. Roughly frac1*n of them will have a
    % rate close to mu1, and (1-frac1)*n of them will have a rate close to mu2.
    rates1 = gamrnd(gamma_shape, mu1/gamma_shape, n, 1);
    rates2 = gamrnd(gamma_shape, mu2/gamma_shape, n, 1);
    is_population_1 = rand(n, 1) < frac1;
    rates(is_population_1) = rates1(is_population_1);
    rates(~is_population_1) = rates2(~is_population_1);
    
    subplot(1,2,1); hold on;
    histogram(rates, 0:2:50);
    
    % Generate poisson "spiking" data for n neurons on t trials
    spikes = poissrnd(repmat(rates, t, 1));
    
    c = cov(spikes);
    e = eig(c);
    
    subplot(1,2,2); hold on;
    plot(flipud(e));
end

subplot(1,2,1);
xlabel('rate');
title('(bimodal) histogram of rates');

subplot(1,2,2);
set(gca, 'yscale', 'log');
xlabel('rank');
ylabel('eigenvalue');
title('eigenspectrum');
