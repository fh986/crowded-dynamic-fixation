% plot thresholds


clc;
clear all;
%close all;
 

addpath('/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Gaze_Package/')

%% set up files

mydir = '/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Tracking_Parameters_Experiment/Final_Experiments/3_Stationary_Dynamic_Flies/data_include_threshold_18ppl_paired';
    
[cursorFiles,mainFiles, eyelinkFiles] = getFiles(mydir);

numSubj = length(cursorFiles);

%% calculate thresholds
subj_stationary_left = nan(1,numSubj);
subj_stationary_right = nan(1,numSubj);
subj_dynamic_left = nan(1,numSubj);
subj_dynamic_right = nan(1,numSubj);
subj_flies_left = nan(1,numSubj);
subj_flies_right = nan(1,numSubj);

subj_exclude = [];

for subj = 1:numSubj

    mainOutput = readtable([mydir filesep mainFiles{subj}]);

    block_sequence = unique(mainOutput.blockShuffleGroups1,'stable');
    rm = cellfun("isempty",block_sequence);
    block_sequence(rm) = [];

    thresholds_raw = mainOutput.questMeanAtEndOfTrialsLoop;
    rm = isnan(thresholds_raw);
    thresholds_raw(rm) = [];
    thresholds_raw = 10.^thresholds_raw;

    idx_block_stationary = find(strcmp(block_sequence,'Stationary'));
    subj_stationary_right(subj) = thresholds_raw(2*idx_block_stationary-1);
    subj_stationary_left(subj) = thresholds_raw(2*idx_block_stationary);

    idx_block_dynamic = find(strcmp(block_sequence,'Dynamic'));
    subj_dynamic_right(subj) = thresholds_raw(2*idx_block_dynamic-1);
    subj_dynamic_left(subj) = thresholds_raw(2*idx_block_dynamic);

    idx_block_flies = find(strcmp(block_sequence,'Flies'));
    subj_flies_right(subj) = thresholds_raw(2*idx_block_flies-1);
    subj_flies_left(subj) = thresholds_raw(2*idx_block_flies);


end

%% calculate averages
firstRun = 1:2:length(subj_stationary_left);
secondRun = 2:2:length(subj_stationary_left);
allRuns = 1:length(subj_stationary_left);
smallThresholds = [1 13 18 26 36];

whichRun = smallThresholds;

avg_stationary_right = geomean(subj_stationary_right(whichRun));
avg_stationary_left = geomean(subj_stationary_left(whichRun));
avg_dynamic_right = geomean(subj_dynamic_right(whichRun));
avg_dynamic_left = geomean(subj_dynamic_left(whichRun));
avg_flies_right = geomean(subj_flies_right(whichRun));
avg_flies_left = geomean(subj_flies_left(whichRun));

%% plot
% 
CData = {[0.4940, 0.1840, 0.5560],[0.8500, 0.3250, 0.0980],[0, 0.4470, 0.7410]};

conditions = [        {'Stationary Left'},...
    {'Stationary Right' },...
    {'Dynamic Left'   },...
    {'Dynamic Right'    },...
    {'Crowded Dynamic Left'     },...
    {'Crowded Dynamic Right'     }];



%% plot with bootstrapping, bootstrap after averaging across 2 sess for each subj, log thresholds
subj_stationary_left = subj_stationary_left(~isnan(subj_stationary_left));
subj_stationary_right = subj_stationary_right(~isnan(subj_stationary_right));subj_dynamic_left = subj_dynamic_left(~isnan(subj_dynamic_left));
subj_dynamic_right = subj_dynamic_right(~isnan(subj_dynamic_right));
subj_flies_left = subj_flies_left(~isnan(subj_flies_left));
subj_flies_right = subj_flies_right(~isnan(subj_flies_right));

subj_avg_stationary_left = avgOverSubj(subj_stationary_left);
subj_avg_stationary_right = avgOverSubj(subj_stationary_right);
subj_avg_dynamic_left = avgOverSubj(subj_dynamic_left);
subj_avg_dynamic_right = avgOverSubj(subj_dynamic_right);
subj_avg_flies_left = avgOverSubj(subj_flies_left);
subj_avg_flies_right = avgOverSubj(subj_flies_right);

   

thresholds = {
    subj_avg_stationary_left,...
    subj_avg_stationary_right,...
    subj_avg_dynamic_left, ...
    subj_avg_dynamic_right, ...
    subj_avg_flies_left, ...
    subj_avg_flies_right
};

% Number of bootstrap samples
numBootstraps = 10000;

% Confidence level (68%)
confLevel = 0.68;

% Initialize array to store bootstrap means
bootstrapMeans = cell(1, length(thresholds));


% Perform bootstrapping for each condition using the bootstrap function
for cond = 1:length(thresholds)
    data = thresholds{cond};

    % Take the logarithm of the thresholds
    logData = log(data);

    % Define a function handle to calculate the mean
    meanFunc = @(data) mean(data);

    % Generate bootstrap samples and calculate means
    bootMeans = bootstrp(numBootstraps, meanFunc, logData);

    bootstrapMeans{cond} = bootMeans;
end

% Calculate the 68% confidence intervals
confIntervals = zeros(length(thresholds), 2);
meanLogThresholds = zeros(1, length(thresholds));
for cond = 1:length(thresholds)
    sortedMeans = sort(bootstrapMeans{cond});
    lowerBound = sortedMeans(round((1 - confLevel) / 2 * numBootstraps));
    upperBound = sortedMeans(round((1 + confLevel) / 2 * numBootstraps));
    confIntervals(cond, :) = [lowerBound, upperBound];
    meanLogThresholds(cond) = mean(log(thresholds{cond})); % Mean of log-transformed data
end

% Convert mean thresholds and error bars back from log scale
meanThresholds = exp(meanLogThresholds);
lowerError = meanThresholds - exp(confIntervals(:, 1)');
upperError = exp(confIntervals(:, 2)') - meanThresholds;
errorBars = [lowerError; upperError];


%%
% Define jitter amount
jitterAmount = 0.05; 
jitter = (rand(1,3) - 0.5) * jitterAmount * 2; 

xPositions = [1, 2];

figure;
hold on;
numGroups = 3;

for ii = 1:numGroups
    idx1 = (ii - 1) * 2 + 1;
    idx2 = idx1 + 1;
    
    xJittered = xPositions + jitter(ii);

    plot(xJittered, [meanThresholds(idx1), meanThresholds(idx2)], 'o-', ...
        'Color', CData{ii}, 'LineWidth', 3, 'MarkerSize', 9);
    errorbar(xJittered(1), meanThresholds(idx1), lowerError(idx1), upperError(idx1), ...
        'LineStyle', 'none', 'Color', CData{ii}, 'LineWidth', 2, 'CapSize', 0);
    errorbar(xJittered(2), meanThresholds(idx2), lowerError(idx2), upperError(idx2), ...
        'LineStyle', 'none', 'Color', CData{ii}, 'LineWidth', 2, 'CapSize', 0);
end

set(gca, 'YScale', 'log');
set(gca, 'FontSize', 18);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Left', 'Right'});

xlim([0 3]); 
ylim([1, 6]);

% Labels
ylabel('Crowding thresholds (deg)');

hold off;



%% determine minimum threshold from Kurzawski et al., 2023, JOV

jov_raw_data = readtable('JoV23Data.csv');

% filter
jov_filtered_data = jov_raw_data(strcmp(jov_raw_data.FlankinDirection, 'radial'),:);
jov_filtered_data = jov_filtered_data(strcmp(jov_filtered_data.Task, 'crowding'),:);
jov_filtered_data = jov_filtered_data(~strcmp(jov_filtered_data.Meridian, 'Upper'),:);
jov_filtered_data = jov_filtered_data(~strcmp(jov_filtered_data.Meridian, 'Lower'),:);
jov_filtered_data = jov_filtered_data(strcmp(jov_filtered_data.Font, 'Sloan'),:);
jov_filtered_data = jov_filtered_data((jov_filtered_data.RadialEccentricity == 10),:);

fprintf('Number of thresholds: %d\n', size(jov_filtered_data,1))

jov_filtered_crowding_distance = jov_filtered_data.CrowdingDistance;
jov_filtered_eccentricity = jov_filtered_data.RadialEccentricity;

jov_filtered_bouma = jov_filtered_crowding_distance./jov_filtered_eccentricity;

minJOVbouma = round(min(jov_filtered_bouma), 3);

fprintf('Minimum Bouma factor from Kurzawski et al., 2023, JoV: %f\n', minJOVbouma)

%% plot histograms for Bouma factors

subj_stationary_leftB = subj_stationary_left./10;
subj_stationary_rightB = subj_stationary_right./10;
subj_dynamic_leftB = subj_dynamic_left./10;
subj_dynamic_rightB = subj_dynamic_right./10;
subj_flies_leftB = subj_flies_left./10;
subj_flies_rightB = subj_flies_right./10;

% subj_stationary_leftB = subj_stationary_left./10*0.7;
% subj_dynamic_leftB = subj_dynamic_left./10*0.7;
% subj_flies_leftB = subj_flies_left./10*0.7;
%%
CData2 = { ...
    [0.5600, 0.1600, 0.6000], ...  % Purple variation
    [0.8700, 0.3700, 0.1300], ...  % Orange variation
    [0, 0.5000, 0.8000] ...        % Blue variation
};

ylimit = [0 10];
% bwidth = 0.1;
xlimit = [0.03 1];

minData = 0.05;
maxData = 1;
numBins = 35; 
binEdges = logspace(log10(minData), log10(maxData), numBins);
binWidths = diff(binEdges);

% make sure of that minJOVdata is included as a bin edge
[~, index] = min(abs(binEdges - minJOVbouma));
offset = binEdges(index) - minJOVbouma;
binEdges = binEdges - offset;

figure;
subplot(3,1,1)
plotStackedHist(subj_stationary_rightB,subj_stationary_leftB,binEdges,xlimit,ylimit, ...
    'Stationary fixation',cell2mat(CData(1)),minJOVbouma)


subplot(3,1,2)
plotStackedHist(subj_dynamic_rightB,subj_dynamic_leftB,binEdges,xlimit,ylimit, ...
    'Dynamic fixation',cell2mat(CData(2)),minJOVbouma)


subplot(3,1,3)
plotStackedHist(subj_flies_rightB,subj_flies_leftB,binEdges,xlimit,ylimit, ...
    'Crowded dynamic fixation',cell2mat(CData(3)),minJOVbouma)
xlabel('Bouma factor b')


%%
function [] = plotStackedHist(data1,data2,binEdges,xlimit,ylimit,titletxt,color1,minJOVbouma)
        
    hold on;
    xlim(xlimit)
    ylim(ylimit)
    xscale('log')
    set(gca,'FontSize',15)
    title(titletxt)
    ylabel('Frequency')
    yticks(min(ylimit):2:max(ylimit))
    
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
    
        disp(xPatch)
        disp(yPatch1)
        disp(yPatch2)
    end
    xline(minJOVbouma,'r--','LineWidth',2);

end

function [] = plotHistLog(data,binEdges,xlimit,ylimit,titletxt,color)


    a = histogram(data,'BinEdges',binEdges);
    a.FaceColor = color; 
    a.EdgeColor = color;
    a.FaceAlpha = 0.6;
    a.EdgeAlpha = 0.6;
    ylim(ylimit)
    xlim(xlimit)
    set(gca,'XScale','log')
    set(gca,'FontSize',15)
    title(titletxt)
    ylabel('Frequency')


end

function [subj_avg_threshold] = avgOverSubj(session_threshold)

    numSubj = length(session_threshold)/2;
    subj_avg_threshold = NaN(1,numSubj);
    for ii = 1:numSubj
        sess1_threshold = session_threshold(ii*2-1);
        sess2_threshold = session_threshold(ii*2);
        subj_avg_threshold(ii) = mean([sess1_threshold,sess2_threshold]);
    end
end