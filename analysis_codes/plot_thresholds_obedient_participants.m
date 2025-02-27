% takes in peeking summary and plot

clc;
clear all;
close all;

%% specify variables

numSessions = 40;
numTrials = 70;

criterion = '1.510000';
addpath('/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Gaze_Package/')

%%
mydir = '/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Tracking Parameters Experiment/Final_Experiments/3_Stationary_Dynamic_Flies/aaa_Stationary_Dynamic_Flies_Codes';

d = dir(sprintf('%s/_PeekingSummary_%s.mat',mydir,criterion));
    
files = {d.name};

assert(length(files) == 1);

load(cell2mat(files));

assert(length(subj_stationary_count) == length(subj_dynamic_count));
assert(length(subj_flies_count) == length(subj_dynamic_count));


%% 

numSessions = length(subj_stationary_count);
sessions_include = [];

peek_summary = [subj_stationary_count'; subj_dynamic_count'; subj_flies_count'];
assert(size(peek_summary, 1) == 3);
assert(size(peek_summary, 2) == numSessions);

for session = 1:numSessions
    
    counts = peek_summary(:, session);
    
    if all(counts == 0)
        sessions_include = [sessions_include, session];
    end
end


disp(sessions_include)


%% set up files

mydir = '/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Tracking Parameters Experiment/Final_Experiments/3_Stationary_Dynamic_Flies/data_include_forGaze_20ppl_paired';
    
[cursorFiles,mainFiles, eyelinkFiles] = getFiles(mydir);

assert(length(cursorFiles) == numSessions);

%% calculate thresholds
subj_stationary_left = [];
subj_stationary_right = [];
subj_dynamic_left = [];
subj_dynamic_right = [];
subj_flies_left = [];
subj_flies_right = [];

subj_exclude = [];
for subj = sessions_include



    mainOutput = readtable([mydir filesep mainFiles{subj}]);

    block_sequence = unique(mainOutput.blockShuffleGroups1,'stable');
    rm = cellfun("isempty",block_sequence);
    block_sequence(rm) = [];

    thresholds_raw = mainOutput.questMeanAtEndOfTrialsLoop;
    rm = isnan(thresholds_raw);
    thresholds_raw(rm) = [];
    thresholds_raw = 10.^thresholds_raw;

    idx_block_stationary = find(strcmp(block_sequence,'Stationary'));
    subj_stationary_right = array_append(subj_stationary_right, thresholds_raw(2*idx_block_stationary-1));
    subj_stationary_left = array_append(subj_stationary_left, thresholds_raw(2*idx_block_stationary));

    idx_block_dynamic = find(strcmp(block_sequence,'Dynamic'));
    subj_dynamic_right = array_append(subj_dynamic_right, thresholds_raw(2*idx_block_dynamic-1));
    subj_dynamic_left = array_append(subj_dynamic_left, thresholds_raw(2*idx_block_dynamic));

    idx_block_flies = find(strcmp(block_sequence,'Flies'));
    subj_flies_right = array_append(subj_flies_right, thresholds_raw(2*idx_block_flies-1));
    subj_flies_left = array_append(subj_flies_left, thresholds_raw(2*idx_block_flies));


    % a couple of reasons to exclude a threshold:
    % 1. less than 30 trials were given to QUEST
    trialGivenBool = mainOutput.trialGivenToQuest;
    num_trialsNotGiven = sum(strcmp(trialGivenBool,'FALSE'));
    if num_trialsNotGiven > 5
        fprintf('Exclude (#1): %d \n',subj)
        subj_exclude = [subj_exclude,subj];
    end

    % 2. SD of threshold is larger than 0.1
    sd_raw = mainOutput.questSDAtEndOfTrialsLoop;
    rm = isnan(sd_raw);
    sd_raw(rm) = [];
    if any(sd_raw > 0.1)
        fprintf('Exclude (#2): %d \n',subj)
        subj_exclude = [subj_exclude,subj];

        stationary_blocks = [2*idx_block_stationary-1,2*idx_block_stationary];
        dynamic_blocks = [2*idx_block_dynamic-1,2*idx_block_dynamic];
        flies_blocks = [2*idx_block_flies-1,2*idx_block_flies];
        idx_largeSD = find(sd_raw > 0.1);
        if any(ismember(idx_largeSD,stationary_blocks))
            fprintf('Exclude (#2): %d -- Stationary',subj)
        end
        if any(ismember(idx_largeSD,dynamic_blocks))
            fprintf('Exclude (#2): %d -- Dynamic',subj)
        end
        if any(ismember(idx_largeSD,flies_blocks))
            fprintf('Exclude (#2): %d -- Flies',subj)
        end

    end
    

end

%%
subj_exclude_unique = unique(subj_exclude);

subj_exclude_unique = unique([subj_exclude_unique, 19, 20]); % excluding subject ML because of a temporary software bug

for ii = 1:length(subj_exclude_unique)

    rm = find(sessions_include == subj_exclude_unique(ii));

    subj_stationary_left(rm) = [];
    subj_stationary_right(rm) = [];
    subj_dynamic_left(rm) = [];
    subj_dynamic_right(rm) = [];
    subj_flies_left(rm) = [];
    subj_flies_right(rm) = [];
    sessions_include(rm) = [];
    
end

% %% calculate averages
% firstRun = 1:2:length(subj_stationary_left);
% secondRun = 2:2:length(subj_stationary_left);
% allRuns = 1:length(subj_stationary_left);
% smallThresholds = [1 13 18 26 36];
% 
% whichRun = smallThresholds;
% 
% avg_stationary_right = geomean(subj_stationary_right(whichRun));
% avg_stationary_left = geomean(subj_stationary_left(whichRun));
% avg_dynamic_right = geomean(subj_dynamic_right(whichRun));
% avg_dynamic_left = geomean(subj_dynamic_left(whichRun));
% avg_flies_right = geomean(subj_flies_right(whichRun));
% avg_flies_left = geomean(subj_flies_left(whichRun));
disp(sessions_include)

subj_stationary_left = avg_btw_subj(sessions_include,subj_stationary_left);
subj_stationary_right = avg_btw_subj(sessions_include,subj_stationary_right);
subj_dynamic_left = avg_btw_subj(sessions_include,subj_dynamic_left);
subj_dynamic_right = avg_btw_subj(sessions_include,subj_dynamic_right);
subj_flies_left = avg_btw_subj(sessions_include,subj_flies_left);
subj_flies_right = avg_btw_subj(sessions_include,subj_flies_right);

%% plot
% 
CData = {[0.4940, 0.1840, 0.5560],[0.8500, 0.3250, 0.0980],[0, 0.4470, 0.7410]};

conditions = [        {'Stationary Left'},...
    {'Stationary Right' },...
    {'Dynamic Left'   },...
    {'Dynamic Right'    },...
    {'Crowded Dynamic Left'     },...
    {'Crowded Dynamic Right'     }];



%% plot with bootstrapping, log thresholds
thresholds = {
    subj_stationary_left,...
    subj_stationary_right,...
    subj_dynamic_left, ...
    subj_dynamic_right, ...
    subj_flies_left, ...
    subj_flies_right
};

% Number of bootstrap samples
numBootstraps = 10000;

% Confidence level (68%)
confLevel = 0.68;

% Initialize array to store bootstrap means
bootstrap_Means = cell(1, length(thresholds));


% Perform bootstrapping for each condition using the bootstrap function
for cond = 1:length(thresholds)
    data = thresholds{cond};

    % Take the logarithm of the thresholds
    logData = log(data);

    % Define a function handle to calculate the mean
    meanFunc = @(data) mean(data);

    % Generate bootstrap samples and calculate means
    bootMeans = bootstrp(numBootstraps, meanFunc, logData);

    bootstrap_Means{cond} = bootMeans;
end

% Calculate the 68% confidence intervals
confIntervals = zeros(length(thresholds), 2);
meanLogThresholds = zeros(1, length(thresholds));
for cond = 1:length(thresholds)
    sortedMeans = sort(bootstrap_Means{cond});
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
% Plot the bar graph with error bars
figure;
hold on;
plot([1 2], [meanThresholds(1) meanThresholds(2)], 'o-', 'Color', CData{1}, 'LineWidth', 3,'MarkerSize',9, 'DisplayName','Stationary');
plot([1 2], [meanThresholds(3) meanThresholds(4)], 'o-', 'Color', CData{2}, 'LineWidth', 3,'MarkerSize',9, 'DisplayName','Dynamic');
plot([1 2], [meanThresholds(5) meanThresholds(6)], 'o-', 'Color', CData{3}, 'LineWidth', 3,'MarkerSize',9, 'DisplayName','Crowded dynamic');

set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
set(gca, 'XTick', [1 2], 'XTickLabel', [{'Left'} {'Right'}]);

xlim([0 3])
% xtickangle(45);
% xlabel('Condition');
ylabel('Crowding thresholds (deg)');
% title('Effect of Background on Crowding Thresholds');


errorbar(1, meanThresholds(1), lowerError(1),upperError(1), 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2, 'HandleVisibility','off');
errorbar(2, meanThresholds(2), lowerError(2),upperError(2), 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2, 'HandleVisibility','off');
errorbar(1, meanThresholds(3), lowerError(3),upperError(3), 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2, 'HandleVisibility','off');
errorbar(2, meanThresholds(4), lowerError(4),upperError(4), 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2, 'HandleVisibility','off');
errorbar(1, meanThresholds(5), lowerError(5),upperError(5), 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2, 'HandleVisibility','off');
errorbar(2, meanThresholds(6), lowerError(6),upperError(6), 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2, 'HandleVisibility','off');
hold off;
ylim([1,6]);
% legend();



%% plot histograms
% 
% subj_stationary_leftB = subj_stationary_left./10;
% subj_stationary_rightB = subj_stationary_right./10;
% subj_dynamic_leftB = subj_dynamic_left./10;
% subj_dynamic_rightB = subj_dynamic_right./10;
% subj_flies_leftB = subj_flies_left./10;
% subj_flies_rightB = subj_flies_right./10;
% 
% 


%% functions

function [output_data] = avg_btw_subj(indices, threshold_data)
    % assumes that indices and threshold_data are sorted

    output_data = [];
    
    % Identify unique participants based on index structure
    unique_participants = unique(floor((indices - 1) / 2) + 1);
    
    for i = 1:length(unique_participants)
        participant = unique_participants(i);
        
        % Identify all indices belonging to the current participant
        participant_indices = find(ismember(indices, [(participant-1)*2+1, (participant-1)*2+2]));
        % Extract the corresponding data
        session_data = threshold_data(participant_indices);
        
        % Average if there are two sessions, otherwise take the single session value
        output_data = [output_data, mean(session_data)];
    end
    
end



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

function [new_array] = array_append(array, number)
    
    new_array = [array, number];
end
