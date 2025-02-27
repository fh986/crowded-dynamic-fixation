% takes in peeking summary and plot

clc;
clear all;
% close all;

%% specify variables

peek_criteria = [0.25 0.5 1 1.5 2];
numSessions = 40;
numTrials = 70;

%%
mydir = '/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Tracking Parameters Experiment/Final_Experiments/3_Stationary_Dynamic_Flies/aaa_Stationary_Dynamic_Flies_Codes';

d = dir(sprintf('%s/PeekingSummary_*.mat',mydir));
    
files = {d.name};

peek_summary_stationary = NaN(length(peek_criteria),numSessions);
peek_summary_dynamic = NaN(length(peek_criteria),numSessions);
peek_summary_flies = NaN(length(peek_criteria),numSessions);


for f = 1 :length(files)

    load(cell2mat(files(f)));

    peek_summary_stationary(f,:) = subj_stationary_count';
    peek_summary_dynamic(f,:) = subj_dynamic_count';
    peek_summary_flies(f,:) = subj_flies_count';

end

%%

% get percentages
peek_summary_stationary = peek_summary_stationary ./ numTrials .* 100;
peek_summary_dynamic = peek_summary_dynamic ./ numTrials .* 100;
peek_summary_flies = peek_summary_flies ./ numTrials .* 100;

% 
% avg_peek_stationary = mean(peek_summary_stationary,2);
% avg_peek_dynamic = mean(peek_summary_dynamic,2);
% avg_peek_flies = mean(peek_summary_flies,2);

med_peek_stationary = median(peek_summary_stationary,2);
med_peek_dynamic = median(peek_summary_dynamic,2);
med_peek_flies = median(peek_summary_flies,2);

%%
% bootstrap
numBootstraps = 10000;
confLevel = 0.68;

[means_stationary,lowerError_stationary,upperError_stationary] = bootstrapMed(peek_summary_stationary,numBootstraps,confLevel);
[means_dynamic,lowerError_dynamic,upperError_dynamic] = bootstrapMed(peek_summary_dynamic,numBootstraps,confLevel);
[means_flies,lowerError_flies,upperError_flies] = bootstrapMed(peek_summary_flies,numBootstraps,confLevel);
%% plot
CData = {[0.4940, 0.1840, 0.5560],[0.8500, 0.3250, 0.0980],[0, 0.4470, 0.7410]};

figure;

hold on;
errorbar(peek_criteria,means_stationary,lowerError_stationary,upperError_stationary,'-o','LineWidth',2,'Color',CData{1},'DisplayName','Stationary Fixation',...
    'CapSize',0,'MarkerSize',10)
errorbar(peek_criteria,means_dynamic,lowerError_dynamic,upperError_dynamic,'-square','LineWidth',2,'Color',CData{2},'DisplayName','Dynamic Fixation',...
    'CapSize',0,'MarkerSize',10)
errorbar(peek_criteria,means_flies,lowerError_flies,upperError_flies,'-d','LineWidth',2,'Color',CData{3},'DisplayName','Crowded Dynamic Fixation',...
    'CapSize',0,'MarkerSize',10)%,'CapSize',0
xlim([0 2]);ylim([0 50])
xticks([0 0.5 1 1.5 2])
% legend('Location','northeast')
set(gca,'FontSize',18)
xlabel('Fixation Tolerance (deg)')
ylabel('Percentage of Trials')
% title('Effect of Background on Breaking Fixation')

hold on;
scatter(1.55,9.3,'o','filled','SizeData',100, 'DisplayName', 'Stationary Fixation(Kurzawski et al., 2023)',...
    'MarkerFaceColor',CData{1},'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.6)
scatter(1.55,5.4,'square','filled','SizeData',100,'DisplayName', 'Dynamic Fixation(Kurzawski et al., 2023)', ...
    'MarkerFaceColor',CData{2},'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.6)

%% for plotting only fixation break towards correct directions
% %% plot
% CData = {[0.4940, 0.1840, 0.5560],[0.8500, 0.3250, 0.0980],[0, 0.4470, 0.7410]};
% 
% figure;
% 
% 
% hold on;
% errorbar(peek_criteria,means_stationary,lowerError_stationary,upperError_stationary,'-o','LineWidth',2,'Color',CData{1},'DisplayName','Stationary Fixation',...
%     'CapSize',0,'MarkerSize',8)
% errorbar(peek_criteria,means_dynamic,lowerError_dynamic,upperError_dynamic,'-square','LineWidth',2,'Color',CData{2},'DisplayName','Dynamic Fixation',...
%     'CapSize',0,'MarkerSize',8)
% errorbar(peek_criteria,means_flies,lowerError_flies,upperError_flies,'-d','LineWidth',2,'Color',CData{3},'DisplayName','Crowded Dynamic Fixation',...
%     'CapSize',0,'MarkerSize',8)%,'CapSize',0
% xlim([0 2]);ylim([0 50])
% xticks([0 0.5 1 1.5 2])
% % legend('Location','northeast')
% set(gca,'FontSize',18)
% xlabel('Fixation Tolerance (deg)')
% ylabel('Percentage of Trials')
% % title('Effect of Background on Breaking Fixation')
% 


%% plot
% CData = {[0.4940, 0.1840, 0.5560],[0.8500, 0.3250, 0.0980],[0, 0.4470, 0.7410]};
% 
% figure;
% 
% hold on;
% errorbar(peek_criteria(2:end),means_stationary(2:end),lowerError_stationary(2:end),upperError_stationary(2:end),'-o','LineWidth',2,'Color',CData{1},'DisplayName','Stationary',...
%     'CapSize',0,'MarkerSize',8)
% errorbar(peek_criteria(2:end),means_dynamic(2:end),lowerError_dynamic(2:end),upperError_dynamic(2:end),'-square','LineWidth',2,'Color',CData{2},'DisplayName','Dynamic',...
%     'CapSize',0,'MarkerSize',8)
% errorbar(peek_criteria(2:end),means_flies(2:end),lowerError_flies(2:end),upperError_flies(2:end),'-d','LineWidth',2,'Color',CData{3},'DisplayName','Flies',...
%     'CapSize',0,'MarkerSize',8)%,'CapSize',0
% xlim([0 2]);ylim([0 20])
% legend('Location','northeast')
% set(gca,'FontSize',18)
% xlabel('Fixation Tolerance (deg)')
% ylabel('Percentage of Trials')
% % title('Effect of Background on Breaking Fixation')
% 
% hold on;
% scatter(1.55,9.3,'o','filled','SizeData',80, 'DisplayName', 'Stationary (Kurzawski et al., 2023)',...
%     'MarkerFaceColor',CData{1},'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.6)
% scatter(1.55,5.4,'square','filled','SizeData',90,'DisplayName', 'Dynamic (Kurzawski et al., 2023)', ...
%     'MarkerFaceColor',CData{2},'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.6)

%% statistical test

% % put data into correct formats for ANOVA
% 
% % for each threshold, we want to compare across conditions
% % concatenate peek summaries matrices
% all_summaries = [peek_summary_stationary,peek_summary_dynamic,peek_summary_flies];
% % Create a 1x120 cell array
% conditions = cell(1, 120);
% conditions(1:40) = {'stationary'};
% conditions(41:80) = {'dynamic'};
% conditions(81:120) = {'flies'};
%%
% [p1, tbl1, stats1] = conductANOVA1(all_summaries,conditions,1);
% [p2, tbl2, stats2] = conductANOVA1(all_summaries,conditions,2);
% [p3, tbl3, stats3] = conductANOVA1(all_summaries,conditions,3);
% [p4, tbl4, stats4] = conductANOVA1(all_summaries,conditions,4);
% [p5, tbl5, stats5] = conductANOVA1(all_summaries,conditions,5);
% 
% %%
% [p1, tbl1, stats1] = conductKW(all_summaries,conditions,1);
% [p2, tbl2, stats2] = conductKW(all_summaries,conditions,2);
% [p3, tbl3, stats3] = conductKW(all_summaries,conditions,3);
% [p4, tbl4, stats4] = conductKW(all_summaries,conditions,4);
% [p5, tbl5, stats5] = conductKW(all_summaries,conditions,5);

%% use Wilcoxon test to compare dynamic and flies
% 
% for ii = 1:5
%     peek_summary_d = peek_summary_dynamic(ii,:);
%     peek_summary_f = peek_summary_flies(ii,:);
% 
%     % Perform the Wilcoxon rank-sum test
%     p = ranksum(peek_summary_d, peek_summary_f)
% 
% end

%% functions

function [means,lowerError,upperError] = bootstrapMeans(dataMatrix,numBootstraps,confLevel)
% Data should have the following format: columns are different samples
% (participants/sessions) and rows are different conditions

    % Number of conditions
    numConditions = size(dataMatrix, 1);

    % Initialize array to store bootstrap means
    bootstrapMeans = cell(1, numConditions);
    
    % Perform bootstrapping for each condition using the bootstrap function
    for cond = 1:numConditions
        data = dataMatrix(cond, :);
        
        % Define a function handle to calculate the mean
        meanFunc = @(data) mean(data);
        
        % Generate bootstrap samples and calculate means
        bootMeans = bootstrp(numBootstraps, meanFunc, data);
        
        bootstrapMeans{cond} = bootMeans;
    end

    % Calculate the 68% confidence intervals
    confIntervals = zeros(numConditions, 2);
    means = zeros(1, numConditions);
    for cond = 1:numConditions
        sortedMeans = sort(bootstrapMeans{cond});
        lowerBound = sortedMeans(round((1 - confLevel) / 2 * numBootstraps));
        upperBound = sortedMeans(round((1 + confLevel) / 2 * numBootstraps));
        confIntervals(cond, :) = [lowerBound, upperBound];
        means(cond) = mean(dataMatrix(cond, :)); % Mean of the data
    end
    
    % Calculate error bars
    lowerError = means - confIntervals(:, 1)';
    upperError = confIntervals(:, 2)' - means;
    % errorBars = [lowerError; upperError];


end

function [medians,lowerError,upperError] = bootstrapMed(dataMatrix,numBootstraps,confLevel)
% Data should have the following format: columns are different samples
% (participants/sessions) and rows are different conditions

    % Number of conditions
    numConditions = size(dataMatrix, 1);

    % Initialize array to store bootstrap means
    bootstrapMeans = cell(1, numConditions);
    
    % Perform bootstrapping for each condition using the bootstrap function
    for cond = 1:numConditions
        data = dataMatrix(cond, :);
        
        % Define a function handle to calculate the mean
        meanFunc = @(data) median(data);
        
        % Generate bootstrap samples and calculate means
        bootMeans = bootstrp(numBootstraps, meanFunc, data);
        
        bootstrapMeans{cond} = bootMeans;
    end

    % Calculate the 68% confidence intervals
    confIntervals = zeros(numConditions, 2);
    medians = zeros(1, numConditions);
    for cond = 1:numConditions
        sortedMeans = sort(bootstrapMeans{cond});
        lowerBound = sortedMeans(round((1 - confLevel) / 2 * numBootstraps));
        upperBound = sortedMeans(round((1 + confLevel) / 2 * numBootstraps));
        confIntervals(cond, :) = [lowerBound, upperBound];
        medians(cond) = median(dataMatrix(cond, :)); % Mean of the data
    end
    
    % Calculate error bars
    lowerError = medians - confIntervals(:, 1)';
    upperError = confIntervals(:, 2)' - medians;
    % errorBars = [lowerError; upperError];


end



function [means,lowerError,upperError] = bootstrapGeoMean(dataMatrix,numBootstraps,confLevel)
% Data should have the following format: columns are different samples
% (participants/sessions) and rows are different conditions

    % Number of conditions
    numConditions = size(dataMatrix, 1);

    % Initialize array to store bootstrap means
    bootstrapMeans = cell(1, numConditions);
    
    % Perform bootstrapping for each condition using the bootstrap function
    for cond = 1:numConditions
        data = dataMatrix(cond, :);
        if any(data == 0)
            data(data == 0) = 1e-10;
        end
        
        % Define a function handle to calculate the mean
        meanFunc = @(data) geomean(data);
        
        % Generate bootstrap samples and calculate means
        bootMeans = bootstrp(numBootstraps, meanFunc, data);
        
        bootstrapMeans{cond} = bootMeans;
    end

    % Calculate the 68% confidence intervals
    confIntervals = zeros(numConditions, 2);
    means = zeros(1, numConditions);
    for cond = 1:numConditions
        sortedMeans = sort(bootstrapMeans{cond});
        lowerBound = sortedMeans(round((1 - confLevel) / 2 * numBootstraps));
        upperBound = sortedMeans(round((1 + confLevel) / 2 * numBootstraps));
        confIntervals(cond, :) = [lowerBound, upperBound];
        means(cond) = geomean(dataMatrix(cond, :)); % Mean of the data
    end
    
    % Calculate error bars
    lowerError = means - confIntervals(:, 1)';
    upperError = confIntervals(:, 2)' - means;
    % errorBars = [lowerError; upperError];


end

function [p, tbl, stats] = conductANOVA1(all_summaries,conditions,tolerance)
% thresholds contain all data points
% conditions contain categorizing variables
% Example data (thresholds and condition labels):
% thresholds = [1.2, 1.5, 1.7, 2.3, 2.1, 2.0, 1.8, 1.9, 2.2];  % Dependent variable
% conditions = {'A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'C'};  % Independent variable
    thresholds = all_summaries(tolerance,:);
    [p, tbl, stats] = anova1(thresholds, conditions);
    multcompare(stats)

% multcompare(stats);

end

function [p, tbl, stats] = conductKW(all_summaries,conditions,tolerance)
% thresholds contain all data points
% conditions contain categorizing variables
% Example data (thresholds and condition labels):
% thresholds = [1.2, 1.5, 1.7, 2.3, 2.1, 2.0, 1.8, 1.9, 2.2];  % Dependent variable
% conditions = {'A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'C'};  % Independent variable
    thresholds = all_summaries(tolerance,:);
    % [p, tbl, stats] = anova1(thresholds, conditions);
    [p, tbl, stats] = kruskalwallis(thresholds, conditions);
    multcompare(stats, 'CType', 'dunn-sidak')

% multcompare(stats);

end