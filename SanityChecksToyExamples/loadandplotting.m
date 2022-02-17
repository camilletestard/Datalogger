%% Script to load in and plot results from various sanity checks

%% Tunesweep for varied number of neurons

%Will set this up for loading in all of the number of neurons sessions
%later  For now setting up plotting and relying on manual loading

labels = arrayfun(@num2str,binning*.001,'UniformOutput',false)

labels{length(labels)+1} = 'Chance';

%LinSep results

LinSep = horzcat(Sweep_Results.LinSep{:});
LinSep(1,:) = []; %Delete the first column as just full of nans

figure %make this better later, but fine to illustrate the points for now
plot(sweep,LinSep','.', 'MarkerSize',10)
hold on
plot(sweep,33*ones(length(sweep),1),'k--')
title(['Simulated Linear Decoder Performance for b = ' num2str(n_behaviors) ' and N =' num2str(n_neurons)])
xlabel('Increase in probability of firing (tunning)')
ylabel('% Correct')
legend(labels)

%Multinomial regression results
MNR = horzcat(Sweep_Results.MNR{:});
MNR(1,:) = []; %Delete the first column as just full of nans
figure %make this better later, but fine to illustrate the points for now
plot(sweep,MNR','.', 'MarkerSize',10)
hold on
plot(sweep,33*ones(length(sweep),1),'k--')
title(['Simulated MNR Performance for b = ' num2str(n_behaviors) ' and N =' num2str(n_neurons)])
xlabel('Increase in probability of firing (tunning)')
ylabel('% Correct')
legend(labels)
