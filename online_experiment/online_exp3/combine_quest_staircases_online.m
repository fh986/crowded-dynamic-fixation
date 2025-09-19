% plot_thresholds_online.m
% Author: Helen Hu

% This script acquires threshold data, one sessions for each
% participant, and plots a stacked histogram of all thresholds.


addpath(genpath('/Applications/Psychtoolbox'))

clc;
clear all;
% close all;

eccentricity = 10;
%% set up files
scriptDir = pwd;
% repoDir = fileparts(scriptDir);
mydir = [scriptDir '/online3_data'];
d = dir(sprintf('%s/*.csv',mydir));
files = {d.name};

for f = 1 :length(files)
    mainFile = dir(sprintf('%s/%s*.csv',mydir,files{f}(1:15)));
    mainFiles{f} = mainFile.name;
end

mainFiles = unique(mainFiles);
numSubj = length(mainFiles);

conditions = {'Stationary', 'Dynamic', 'Flies'};

% %% filter out subjects excluded for final analysis
% load('subj_exclude.mat');
% 
% fprintf('Number of subjects excluded in the final analysis: %d\n', length(subj_exclude))
% 
% 

%%

%% extract thresholds
results = table();  
row = 1;  
for subj = 1:numSubj

    mainOutput = readtable([mydir filesep mainFiles{subj}], 'VariableNamingRule', 'preserve');

    if ismember('trials.response', mainOutput.Properties.VariableNames)
        mainOutput.Properties.VariableNames{'trials.response'} = 'trials_response';
    else
        disp('------')
        disp('Column not found.');
        disp(mainFiles{subj})

    end

    % extract threshold
    for cond = 1:length(conditions)

        condition = conditions{cond};
        conditionName_test = sprintf('%s_Right1', condition);
        conditionName_retest = sprintf('%s_Right2', condition);

        trials_test = mainOutput(strcmp(mainOutput.conditionName, conditionName_test), :);
        trials_retest = mainOutput(strcmp(mainOutput.conditionName, conditionName_retest), :);

        trials_combined = [trials_test; trials_retest];
        trials_sent = trials_combined(strcmp(trials_combined.trialGivenToQuest, 'TRUE'), :);

        levelProposed = trials_sent.levelProposedByQUEST;
        responses = trials_sent.trials_response;
        assert(size(levelProposed,1) == size(responses, 1))

        threshold_predicted_log = quest_estimate_from_trials(levelProposed, responses);
        threshold_predicted = 10 ^ threshold_predicted_log;


        % get thresholds predicted by quest
        threshold_test = 10.^trials_test.questMeanAtEndOfTrialsLoop(~isnan(trials_test.questMeanAtEndOfTrialsLoop));
        threshold_retest = 10.^trials_retest.questMeanAtEndOfTrialsLoop(~isnan(trials_retest.questMeanAtEndOfTrialsLoop));

        if numel(threshold_test) == 1 && numel(threshold_retest) == 1
    
            bouma_test = threshold_test ./ eccentricity;
            bouma_retest = threshold_retest ./ eccentricity;
            bouma_predicted = threshold_predicted ./ eccentricity;
    
            % Add a row to the results table
            results(row,:) = table(subj, {condition}, bouma_test, bouma_retest, bouma_predicted, ...
                                    threshold_test, threshold_retest, threshold_predicted, ...
                                   'VariableNames', {'Subject', 'Condition', 'BoumaTest', 'BoumaRetest', 'BoumaPredicted', ...
                                   'ThresholdTest', 'ThresholdRetest', 'ThresholdPredicted'});
    
            row = row + 1;

        else

            disp('Warning: no thresholds.')
            disp(mainFiles{subj})

        end
    end
end

%% determine minimum threshold from Kurzawski et al., 2023, JOV

repoDir = fileparts(fileparts(scriptDir));
jov_raw_data = readtable(sprintf('%s/analysis_codes/JoV23Data.csv',repoDir));

% filter
jov_filtered_data = jov_raw_data(strcmp(jov_raw_data.FlankinDirection, 'radial'),:);
jov_filtered_data = jov_filtered_data(strcmp(jov_filtered_data.Task, 'crowding'),:);
jov_filtered_data = jov_filtered_data(~strcmp(jov_filtered_data.Meridian, 'Upper'),:);
jov_filtered_data = jov_filtered_data(~strcmp(jov_filtered_data.Meridian, 'Lower'),:);
jov_filtered_data = jov_filtered_data(strcmp(jov_filtered_data.Font, 'Sloan'),:);
jov_filtered_data = jov_filtered_data((jov_filtered_data.RadialEccentricity == 10),:);

fprintf('Number of thresholds from Kurzawski, after filtering: %d\n', size(jov_filtered_data,1))

jov_filtered_crowding_distance = jov_filtered_data.CrowdingDistance;
jov_filtered_eccentricity = jov_filtered_data.RadialEccentricity;

jov_filtered_bouma = jov_filtered_crowding_distance./jov_filtered_eccentricity;

minJOVbouma = round(min(jov_filtered_bouma), 3);

fprintf('Minimum Bouma factor from Kurzawski et al., 2023, JoV: %f\n', minJOVbouma)



%% filter out participant with estimated Bouma of > 10
% crowding distance would be 100 deg which is impossible

ind_too_large = find(results.BoumaPredicted > 10);
subj_too_large = unique(results(ind_too_large,:).Subject);
rows_to_delete = ismember(results.Subject, subj_too_large);
results_clean = results(~rows_to_delete, :);


%% plot histograms for Bouma factors
 
CData = {[0.4940, 0.1840, 0.5560],[0.8500, 0.3250, 0.0980],[0, 0.4470, 0.7410]};
CData2 = {[0.5600, 0.1600, 0.6000],[0.8700, 0.3700, 0.1300],[0, 0.5000, 0.8000]};

results_for_plotting = results_clean;
table_stationary = results_for_plotting(strcmp(results_for_plotting.Condition, 'Stationary'),:);
table_dynamic = results_for_plotting(strcmp(results_for_plotting.Condition, 'Dynamic'),:);
table_flies = results_for_plotting(strcmp(results_for_plotting.Condition, 'Flies'),:);

ylimit = [0 8];
% bwidth = 0.1;
xlimit = [0.01 1.21];

minData = 0.01;
maxData = 1.2;
numBins = 35; 
binEdges = logspace(log10(minData), log10(maxData), numBins);
binWidths = diff(binEdges);

% make sure of that minJOVdata is included as a bin edge
[~, index] = min(abs(binEdges - minJOVbouma));
offset = binEdges(index) - minJOVbouma;
binEdges = binEdges - offset;

figure;
subplot(3,1,1)
plotStackedHist(table_stationary.BoumaTest,table_stationary.BoumaRetest,binEdges,xlimit,ylimit, ...
    'Stationary fixation',cell2mat(CData(1)),minJOVbouma)


subplot(3,1,2)
plotStackedHist(table_dynamic.BoumaTest,table_dynamic.BoumaRetest,binEdges,xlimit,ylimit, ...
    'Dynamic fixation',cell2mat(CData(2)),minJOVbouma)


subplot(3,1,3)
plotStackedHist(table_flies.BoumaTest,table_flies.BoumaRetest,binEdges,xlimit,ylimit, ...
    'Crowded dynamic fixation',cell2mat(CData(3)),minJOVbouma)
xlabel('Bouma factor b')

% sgtitle(sprintf('Number of subjects: %d; Number of thresholds per condition: %d', height(table_flies), height(table_flies)*2))


%% Plot sess 1 and 2 combined


plotCondScatter(table_stationary, 'Stationary', minJOVbouma)

plotCondScatter(table_dynamic, 'Dynamic', minJOVbouma)

plotCondScatter(table_flies, 'Flies', minJOVbouma)


%% Plot hist with 1 and 2 combined


ylimit = [0 6];
% bwidth = 0.1;
xlimit = [0.01 1.21];

minData = 0.01;
maxData = 1.2;
numBins = 35; 
binEdges = logspace(log10(minData), log10(maxData), numBins);
binWidths = diff(binEdges);

% make sure of that minJOVdata is included as a bin edge
[~, index] = min(abs(binEdges - minJOVbouma));
offset = binEdges(index) - minJOVbouma;
binEdges = binEdges - offset;

figure;
subplot(3,1,1)
plotStackedHist(table_stationary.BoumaPredicted,[],binEdges,xlimit,ylimit, ...
    'Stationary fixation',cell2mat(CData(1)),minJOVbouma)


subplot(3,1,2)
plotStackedHist(table_dynamic.BoumaPredicted,[],binEdges,xlimit,ylimit, ...
    'Dynamic fixation',cell2mat(CData(2)),minJOVbouma)


subplot(3,1,3)
plotStackedHist(table_flies.BoumaPredicted,[],binEdges,xlimit,ylimit, ...
    'Crowded dynamic fixation',cell2mat(CData(3)),minJOVbouma)
xlabel('Bouma factor b')


%% Print out standard deviation and mode for each condition
% SD of log thresholds
stationary_sd = std(log10(table_stationary.BoumaPredicted));
fprintf('SD of log of stationary thresholds: %f\n', stationary_sd)

dynamic_sd = std(log10(table_dynamic.BoumaPredicted));
fprintf('SD of log of dynamic thresholds: %f\n', dynamic_sd)

flies_sd = std(log10(table_flies.BoumaPredicted));
fprintf('SD of log of flies thresholds: %f\n', flies_sd);

%%
% modes of thresholds
stationary_mode = histMode(table_stationary.BoumaPredicted, binEdges, 'mean');
fprintf('Mode of stationary thresholds: %f\n', stationary_mode)

dynamic_mode = histMode(table_dynamic.BoumaPredicted, binEdges, 'mean');
fprintf('Mode of dynamic thresholds: %f\n', dynamic_mode)

flies_mode = histMode(table_flies.BoumaPredicted, binEdges, 'mean');
fprintf('Mode of flies thresholds: %f\n', flies_mode)


%%

function mode_est = histMode(data, binEdges, tieStrategy)
% histMode Estimate the mode using histogram binning
%
%   mode_est = histMode(data, binEdges)
%       returns the midpoints of the fullest bin(s).
%
%   mode_est = histMode(data, binEdges, tieStrategy)
%       tieStrategy options:
%           'all'   - return all tied modal bins (default)
%           'first' - return the first tied bin
%           'last'  - return the last tied bin
%           'mean'  - return the average of all tied bins
%
% Example:
%   data = randn(100,1);
%   edges = -3:0.5:3;
%   mode_est = histMode(data, edges, 'mean');

    if nargin < 3
        tieStrategy = 'all';
    end

    % Compute histogram counts
    [counts, edges] = histcounts(data, binEdges);

    % Find the fullest bin(s)
    maxCount = max(counts);
    modalBins = find(counts == maxCount);

    % Bin midpoints
    mids = mean([edges(modalBins); edges(modalBins+1)], 1);

    % Handle ties according to chosen strategy
    switch tieStrategy
        case 'all'
            mode_est = mids;          % return all tied bins
        case 'first'
            mode_est = mids(1);       % return first
        case 'last'
            mode_est = mids(end);     % return last
        case 'mean'
            mode_est = mean(mids);    % average across all ties
        otherwise
            error('Unknown tieStrategy: %s', tieStrategy);
    end
end

function [t] = quest_estimate_from_trials(levels, responses, ...
    tGuess, tGuessSd, pThreshold, beta, delta, gamma)
% levels:    N×1 vector of stimulus intensities (in the SAME units you use for Quest;
%            typically log10 contrast). If you stored linear contrast c, use log10(c).
% responses: N×1 vector of 0/1 (0 = wrong, 1 = right).


    arguments
        levels (:,1) double
        responses (:,1) double
        tGuess (1,1) double = log10(5) % from EE spreadsheet
        tGuessSd (1,1) double = 2 % easyeyes default
        pThreshold (1,1) double = 0.7 % EE spreadsheet
        beta (1,1) double = 2.3 % EE default
        delta (1,1) double = 0.01 % EE default
        gamma (1,1) double = 1/9 % 1/number of choices 
    end

    % Create Quest state with your prior & psychometric parameters
    q = QuestCreate(tGuess, tGuessSd, pThreshold, beta, delta, gamma);
    q.normalizePdf = 1;

    % drop NaN rows
    valid = ~(isnan(levels) | isnan(responses));
    levels = levels(valid);
    responses = responses(valid) ~= 0; % force logical 0/1

    % Feed each recorded trial to Quest
    for i = 1:numel(levels)
        q = QuestUpdate(q, levels(i), responses(i));
    end

    % Read out the estimate
    t  = QuestMean(q);     % posterior mean (recommended)
    % sd = QuestSd(q);       % posterior sd
    % Alternatives:
    % t  = QuestMode(q);    % MAP estimate
    % t  = QuestQuantile(q);% median
end



function [] = plotStackedHist(data1,data2,binEdges,xlimit,ylimit,titletxt,color1,minJOVbouma)
        
    hold on;
    xlim(xlimit)
    ylim(ylimit)
    title(titletxt)
    set(gca,'FontSize',15)

    set(gca, 'XScale', 'log');
    xticks([0.02, 0.04, 0.06, 0.08, ...
            0.1, 0.2, 0.4, 0.6, 0.8, 1]);
    set(gca, 'XMinorTick', 'off');  % since we're setting everything manually
    xtickangle(45)

    ylabel('Frequency')
    yticks(min(ylimit):2:max(ylimit))
    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off');
    set(gca, 'TickDir', 'out');

    counts1 = histcounts(data1, binEdges);
    counts2 = histcounts(data2, binEdges);

    for ii = 1:length(binEdges)-1
        % Dataset 1
        xPatch = [binEdges(ii), binEdges(ii+1), binEdges(ii+1), binEdges(ii)]; 
        yPatch1 = [0, 0, counts1(ii), counts1(ii)]; 
        patch(xPatch, yPatch1, color1, 'EdgeColor', 'k', ...
            'EdgeAlpha', 1, 'FaceAlpha', 0.5);
    
        % Dataset 2 (stacked on top of dataset 1)
        yPatch2 = [counts1(ii), counts1(ii), counts1(ii) + counts2(ii), counts1(ii) + counts2(ii)];
        patch(xPatch, yPatch2, color1, 'EdgeColor', 'k', ...
            'EdgeAlpha', 1, 'FaceAlpha', 0.8);

    end
    xline(minJOVbouma,'r--','LineWidth',2);

end

function [] = plotScatter(x_val, x_label, y_val, y_label, titletext)
    figure;
    scatter(x_val,y_val,'filled','SizeData',60, ...
        'MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
    hold on;
    xline(0.096, 'r--', 'LineWidth',1)
    yline(0.096, 'r--', 'LineWidth',1)
    plot([0.01,10], [0.01, 10], 'r--', 'LineWidth',1)
    ylabel(y_label)
    xlabel(x_label)
    title(titletext)
    set(gca, 'YScale', 'log') 
    set(gca, 'XScale', 'log') 
    set(gca, 'FontSize', 18)
    axis padded;
    axis equal;
    axis square;


end

function [] = setScatter(minJOVbouma)
    xlimit = [0.01, 100];
    ylimit = [0.01, 100];
    xlim(xlimit); ylim(ylimit);
    xline(minJOVbouma, 'r--')
    yline(minJOVbouma, 'r--')
    plot(xlimit, ylimit, 'r--')
    ylabel('Bouma factor, predicted by Quest (70 trials)')
    set(gca, 'FontSize', 15)
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    axis square 

end


function [] = plotCondScatter(table_for_plotting, title, minJOVbouma)
    figure('Position',[100 100 1000 500])
    subplot(1,2,1);
    scatter(table_for_plotting.BoumaTest, table_for_plotting.BoumaPredicted, 'filled', ...
        'SizeData',60, 'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7)
    hold on;
    xlabel('Bouma factor, test (35 trials)')
    setScatter(minJOVbouma);

    hold off;
    
    subplot(1,2,2);
    scatter(table_for_plotting.BoumaRetest, table_for_plotting.BoumaPredicted, 'filled', ...
        'SizeData',60, 'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7)
    hold on;
    xlabel('Bouma factor, retest (35 trials)')
    setScatter(minJOVbouma);

    hold off;

    sgtitle(title, 'FontSize', 18)
end