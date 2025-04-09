% plot_peeking_summaries.m
% Author: Helen Hu
% Last modified: Mar 7th, 2025

% This file takes in multiple summary files and plot a figure
% showing the number of trials where participants peeked 
% as a function of criteria (fixation tolerance)

clc;
clear all;
close all;

%% specify variables

peek_criteria = [0.25 0.5 1 1.5 2];
numSessions = 40;
numTrials = 70;

%%
% Get current directory and move one level up
mydir = pwd;

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

%% Calculate medium percentage

% get percentages
peek_summary_stationary = peek_summary_stationary ./ numTrials .* 100;
peek_summary_dynamic = peek_summary_dynamic ./ numTrials .* 100;
peek_summary_flies = peek_summary_flies ./ numTrials .* 100;

med_peek_stationary = median(peek_summary_stationary,2);
med_peek_dynamic = median(peek_summary_dynamic,2);
med_peek_flies = median(peek_summary_flies,2);

%% Bootstrap
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
xlabel('Fixation tolerance (deg)')
ylabel('Percentage of trials')
% title('Effect of Background on Breaking Fixation')

hold on;
scatter(1.55,9.3,'o','filled','SizeData',100, 'DisplayName', 'Stationary Fixation(Kurzawski et al., 2023)',...
    'MarkerFaceColor',CData{1},'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.6)
scatter(1.55,5.4,'square','filled','SizeData',100,'DisplayName', 'Dynamic Fixation(Kurzawski et al., 2023)', ...
    'MarkerFaceColor',CData{2},'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.6)


%% functions

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

