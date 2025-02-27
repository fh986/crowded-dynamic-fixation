% plot RMSE with bootstrapping
% input: subj_rmse.mat
% 3 matrices with # of rows being: number of subjects
% and # of columns being: all time points collected
% assuming each participant is tested twice and the data are paired
% 1st row = 1st participant sess 1; 2nd row: 1st participant sess 2,etc.

clear all;
% close all;
clc;

%% set up
load("subj_rmse.mat");

assert(size(subj_stationary_rmse,1) == size(subj_dynamic_rmse,1));
assert(size(subj_stationary_rmse,2) == size(subj_dynamic_rmse,2));

assert(size(subj_dynamic_rmse,1) == size(subj_flies_rmse,1));
assert(size(subj_dynamic_rmse,2) == size(subj_flies_rmse,2));

numSess = size(subj_stationary_rmse,1);
numSubj = numSess/2; % assuming each participant is tested twice
numTP = size(subj_stationary_rmse,2);
%% average for each participant

avg_within_subj_stationary = NaN(numSubj,numTP);
avg_within_subj_dynamic = NaN(numSubj,numTP);
avg_within_subj_flies = NaN(numSubj,numTP);

for ii = 1:numSubj

    idx_sess1 = 2*ii-1;
    idx_sess2 = 2*ii;

    stationary_rmse1 = subj_stationary_rmse(idx_sess1,:);
    stationary_rmse2 = subj_stationary_rmse(idx_sess2,:);
    dynamic_rmse1 = subj_dynamic_rmse(idx_sess1,:);
    dynamic_rmse2 = subj_dynamic_rmse(idx_sess2,:);
    flies_rmse1 = subj_flies_rmse(idx_sess1,:);
    flies_rmse2 = subj_flies_rmse(idx_sess2,:);

    avg_within_subj_stationary(ii,:) = mean([stationary_rmse1;stationary_rmse2]);
    avg_within_subj_dynamic(ii,:) = mean([dynamic_rmse1;dynamic_rmse2]);
    avg_within_subj_flies(ii,:) = mean([flies_rmse1;flies_rmse2]);

end

    
%% sanity check. plot lines without bootstrapping
% 
dictionary_lineColors = dictionary({'Stationary','Dynamic','Flies'},{[0.4940 0.1840 0.5560],[0.8500 0.3250 0.0980],[0 0.4470 0.7410]});
dictionary_lineStyle = dictionary({'Stationary','Dynamic','Flies'},{'-','--','-.'});

avg_stationary_std = geomean(avg_within_subj_stationary);
avg_dynamic_std = geomean(avg_within_subj_dynamic);
avg_flies_std = geomean(avg_within_subj_flies);

figure;clf;

% set up the graph 

xlabel('Time (s)');
ylabel('Fixation error RMSE (deg)');
set(gca,'Fontsize',18)
xlim([-0.15 0.3])
ylim([0 2])

createPatch(0,0.15);

hold on;
l1 = plot(xValues,avg_stationary_std,cell2mat(dictionary_lineStyle({'Stationary'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Stationary'})),'DisplayName','Stationary Fixation');
l2 = plot(xValues,avg_dynamic_std,cell2mat(dictionary_lineStyle({'Dynamic'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Dynamic'})),'DisplayName','Dynamic Fixation');
l3 = plot(xValues,avg_flies_std,cell2mat(dictionary_lineStyle({'Flies'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Flies'})),'DisplayName','Crowded Dynamic Fixation');

legend('Location','northwest');
title('Average (N = 20)')
hold off;

% ratio
avg_dynamicVstationary = avg_dynamic_std ./ avg_stationary_std;
avg_fliesVstationary = avg_flies_std ./ avg_stationary_std;

figure;clf;

% set up the graph 

xlabel('Time (s)');
ylabel('Ratio of fixation error RMSE (relative to stationary)');
set(gca,'Fontsize',18)
xlim([-0.15 0.3])
ylim([0 1.2])

createPatch(0,0.15);

hold on;
l4 = plot(xValues,avg_dynamicVstationary,cell2mat(dictionary_lineStyle({'Dynamic'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Dynamic'})),'DisplayName','Dynamic Fixation');
l5 = plot(xValues,avg_fliesVstationary,cell2mat(dictionary_lineStyle({'Flies'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Flies'})),'DisplayName','Crowded Dynamic Fixation');

legend('Location','northwest');
title('Average (N = 20)')
hold off;


%% bootstrapping

[avg_stationary,lowerCI_stationary,upperCI_stationary] = bootstrapMean(avg_within_subj_stationary);
[avg_dynamic,lowerCI_dynamic,upperCI_dynamic] = bootstrapMean(avg_within_subj_dynamic);
[avg_flies,lowerCI_flies,upperCI_flies] = bootstrapMean(avg_within_subj_flies);

within_subj_dynamicVstationary = avg_within_subj_dynamic ./ avg_within_subj_stationary;
within_subj_fliesVstationary = avg_within_subj_flies ./ avg_within_subj_stationary;
[avg_dynamic_ratio,lowerCI_dynamic_ratio,upperCI_dynamic_ratio] = bootstrapMean(within_subj_dynamicVstationary);
[avg_flies_ratio,lowerCI_flies_ratio,upperCI_flies_ratio] = bootstrapMean(within_subj_fliesVstationary);


%%
figure;clf;

% set up the graph 

xlabel('Time (s)');
ylabel('Fixation Error RMSE (deg)');
set(gca,'Fontsize',18)
xlim([-0.15 0.3])
% xlim([min(xValues),max(xValues)])
ylim([0 2])
% yline(0.15,'k--','LineWidth',2)

createPatch(0,0.15);

hold on;
plot(xValues,avg_stationary,cell2mat(dictionary_lineStyle({'Stationary'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Stationary'})),'DisplayName','Stationary Fixation');
fill([xValues, fliplr(xValues)], [upperCI_stationary, fliplr(lowerCI_stationary)], ...
    cell2mat(dictionary_lineColors({'Stationary'})),'FaceAlpha', 0.2, ...
    'EdgeColor', 'none','HandleVisibility','off'); % Plot confidence intervals

plot(xValues,avg_dynamic,cell2mat(dictionary_lineStyle({'Dynamic'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Dynamic'})),'DisplayName','Dynamic Fixation');
fill([xValues, fliplr(xValues)], [upperCI_dynamic, fliplr(lowerCI_dynamic)], ...
    cell2mat(dictionary_lineColors({'Dynamic'})),'FaceAlpha', 0.2, ...
    'EdgeColor', 'none','HandleVisibility','off'); % Plot confidence intervals

plot(xValues,avg_flies,cell2mat(dictionary_lineStyle({'Flies'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Flies'})),'DisplayName','Crowded Dynamic Fixation');
fill([xValues, fliplr(xValues)], [upperCI_flies, fliplr(lowerCI_flies)],  ...
    cell2mat(dictionary_lineColors({'Flies'})),'FaceAlpha', 0.2, ...
    'EdgeColor', 'none','HandleVisibility','off'); % Plot confidence intervals


% legend('Location','northwest');
% title('Average (N = 20)')
hold off;
%%
%%
figure;clf;

% set up the graph 

xlabel('Time (s)');
ylabel('Ratio of fixation error RMSE (relative to Stationary)');
set(gca,'Fontsize',18)
xlim([-0.15 0.3])
% xlim([min(xValues),max(xValues)])
ylim([0 2])
% yline(0.15,'k--','LineWidth',2)

createPatch(0,0.15);

hold on;
plot(xValues,avg_dynamic_ratio,cell2mat(dictionary_lineStyle({'Dynamic'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Dynamic'})),'DisplayName','Dynamic Fixation');
fill([xValues, fliplr(xValues)], [upperCI_dynamic_ratio, fliplr(lowerCI_dynamic_ratio)], ...
    cell2mat(dictionary_lineColors({'Dynamic'})),'FaceAlpha', 0.2, ...
    'EdgeColor', 'none','HandleVisibility','off'); % Plot confidence intervals

plot(xValues,avg_flies_ratio,cell2mat(dictionary_lineStyle({'Flies'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Flies'})),'DisplayName','Crowded Dynamic Fixation');
fill([xValues, fliplr(xValues)], [upperCI_flies_ratio, fliplr(lowerCI_flies_ratio)],  ...
    cell2mat(dictionary_lineColors({'Flies'})),'FaceAlpha', 0.2, ...
    'EdgeColor', 'none','HandleVisibility','off'); % Plot confidence intervals


% legend('Location','northwest');
% title('Average (N = 20)')
hold off;

%% Below, we plot the average ratio over the View period

idx_view_start = find(xValues == 0);
[tp, idx_view_end] = min(abs(xValues - 0.15));

view_withinSubj_stationary = avg_within_subj_stationary(:, idx_view_start:idx_view_end);
view_withinSubj_dynamic = avg_within_subj_dynamic(:, idx_view_start:idx_view_end);
view_withinSubj_flies = avg_within_subj_flies(:, idx_view_start:idx_view_end);

ratio_subj_dynamic_stationary = mean(view_withinSubj_dynamic ./ view_withinSubj_stationary, 2);
ratio_subj_flies_stationary = mean(view_withinSubj_flies ./ view_withinSubj_stationary, 2);

[ratio_dynamicVstationary, lower_dynamicVstationary, higher_dynamicVstatinoary] = bootstrapMean(ratio_subj_dynamic_stationary);
[ratio_fliesVstationary, lower_fliesVstationary, higher_fliesVstatinoary] = bootstrapMean(ratio_subj_flies_stationary);

figure;
xlim([0 3])
ylim([0 0.8])
xticks([1 2])
xticklabels({'Dynamic', 'Crowded dynamic'})
set(gca, 'FontSize', 16)
hold on;
plot(1, ratio_dynamicVstationary, 'o','MarkerSize',10, ...
        'Color',cell2mat(dictionary_lineColors({'Dynamic'})), 'LineWidth',2)
plot(2, ratio_fliesVstationary, 'o','MarkerSize',10, ...
        'Color',cell2mat(dictionary_lineColors({'Flies'})), 'LineWidth',2)

errorbar(1, ratio_dynamicVstationary, ratio_dynamicVstationary - lower_dynamicVstationary, ...
        higher_dynamicVstatinoary - ratio_dynamicVstationary, ...
        'LineWidth',2,'CapSize',0,'Color',cell2mat(dictionary_lineColors({'Dynamic'})))
errorbar(2, ratio_fliesVstationary, ratio_fliesVstationary - lower_fliesVstationary, ...
        higher_fliesVstatinoary - ratio_fliesVstationary, ...
        'LineWidth',2,'CapSize',0,'Color',cell2mat(dictionary_lineColors({'Flies'})))


%%

%% functions


function [] = createPatch(onset,offset)

    % Create a patch to highlight the stimulus onset
    y_limits = ylim();
    stimulus_height = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)];
    stimulus_duration = [onset offset offset onset];
    patch(stimulus_duration, stimulus_height, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility','off')

end



function [meanError,lowerCI,upperCI] = bootstrapMean(data)

    % Assume your data matrix is stored in a variable named 'data'
    % data is a 20x43 matrix (20 participants, 43 time points)
    
    % Number of bootstrap samples
    numBootstraps = 1000;
    
    % Preallocate matrices to store bootstrap results
    bootMeans = zeros(numBootstraps, size(data, 2));
    
    % Perform bootstrapping for each time point
    for t = 1:size(data, 2)
        % Extract data for the current time point
        timePointData = data(:, t);
        
        % Perform bootstrapping to get means
        bootMeans(:, t) = bootstrp(numBootstraps, @geomean, timePointData);
    end
    
    % Calculate the 68% confidence intervals
    lowerCI = prctile(bootMeans, 16, 1); % Lower bound (16th percentile)
    upperCI = prctile(bootMeans, 84, 1); % Upper bound (84th percentile)

    % mean error
    meanError = geomean(data, 1);

end



    
 