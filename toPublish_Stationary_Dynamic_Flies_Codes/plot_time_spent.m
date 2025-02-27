% plot time used for each block (condition)
% Purpose was to test if people spent extra time on flies (potentially more
% difficult)

% clc;
clear all;
%close all;

addpath('/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Gaze_Package/')

%% set up files

mydir = '/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Tracking Parameters Experiment/Final_Experiments/3_Stationary_Dynamic_Flies/data_include_forGaze_20ppl_paired';
% mydir = '/Users/helenhu/Documents/MATLAB/Stationary_Dynamic_Flies_Codes';
    
[cursorFiles,mainFiles, eyelinkFiles] = getFiles(mydir);

numSubj = length(cursorFiles);

%% time used for the block and response time for each trial
subj_stationary_timeUsed = nan(1,numSubj);
subj_dynamic_timeUsed = nan(1,numSubj);
subj_flies_timeUsed = nan(1,numSubj);

subj_stationary_respTime = nan(1,numSubj);
subj_dynamic_respTime = nan(1,numSubj);
subj_flies_respTime = nan(1,numSubj);

subj_stationary_trackTime = nan(1,numSubj);
subj_dynamic_trackTime = nan(1,numSubj);
subj_flies_trackTime = nan(1,numSubj);

for subj = 1:numSubj

    easyeyes = readtable([mydir filesep cursorFiles{subj}]);
    mainOutput = readtable([mydir filesep mainFiles{subj}]);

    % block names, in the sequence as they appeared
    block_sequence = unique(mainOutput.blockShuffleGroups1,'stable');
    rm = cellfun("isempty",block_sequence);
    block_sequence(rm) = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. time used for each block
    % find beginning of each block
    trialSteps = easyeyes.trialStep;
    bool_blockInstruction = strcmp(trialSteps,'_instructionRoutineEachFrame');
    idx_blockInstruction = find(diff(bool_blockInstruction) == -1)+1;
    idx_blockInstruction = [1;idx_blockInstruction]; % add the first block
    assert(length(idx_blockInstruction) == 3);

    % find end of the experiment
    bool_blockLoopEnd = strcmp(trialSteps,'blocksLoopEnd');
    idx_blockLoopEnd = find(diff(bool_blockLoopEnd) == 1)+1;

    % concatenating indices for calculating time used for each block
    idx_blocks = [idx_blockInstruction;idx_blockLoopEnd];

    % time used for each block
    timestamps_block = easyeyes.posixTimeSec(idx_blocks);
    time_used_block_Sec = diff(timestamps_block);
    assert(length(time_used_block_Sec) == 3);

    idx_block_stationary = find(strcmp(block_sequence,'Stationary'));
    subj_stationary_timeUsed(subj) = time_used_block_Sec(idx_block_stationary);

    idx_block_dynamic = find(strcmp(block_sequence,'Dynamic'));
    subj_dynamic_timeUsed(subj) = time_used_block_Sec(idx_block_dynamic);

    idx_block_flies = find(strcmp(block_sequence,'Flies'));
    subj_flies_timeUsed(subj) = time_used_block_Sec(idx_block_flies);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. response time for each block
    trialRoutineDuration = mainOutput.trialRoutineDurationFromBeginSec;
    rm = isnan(trialRoutineDuration);
    trialRoutineDuration(rm) = [];
    assert(length(trialRoutineDuration) == 210);

    avg_resp_time_Sec = [mean(trialRoutineDuration(1:70)),mean(trialRoutineDuration(71:140)),mean(trialRoutineDuration(141:210))];
    subj_stationary_respTime(subj) = avg_resp_time_Sec(idx_block_stationary);
    subj_dynamic_respTime(subj) = avg_resp_time_Sec(idx_block_dynamic);
    subj_flies_respTime(subj) = avg_resp_time_Sec(idx_block_flies);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. time used for tracking
    trialInstructionDuration = mainOutput.trialInstructionRoutineDurationFromBeginSec;
    rm = isnan(trialInstructionDuration);
    trialInstructionDuration(rm) = [];
    assert(length(trialInstructionDuration) == 210);

    avg_track_time_Sec = [mean(trialInstructionDuration(1:70)),mean(trialInstructionDuration(71:140)),mean(trialInstructionDuration(141:210))];
    subj_stationary_trackTime(subj) = avg_track_time_Sec(idx_block_stationary);
    subj_dynamic_trackTime(subj) = avg_track_time_Sec(idx_block_dynamic);
    subj_flies_trackTime(subj) = avg_track_time_Sec(idx_block_flies);
end


%% calculate averages

avg_stationary_timeUsed_min = mean(subj_stationary_timeUsed./60);
avg_dynamic_timeUsed_min = mean(subj_dynamic_timeUsed./60);
avg_flies_timeUsed_min = mean(subj_flies_timeUsed./60);
sem_stationary_timeUsed = sem(subj_stationary_timeUsed./60);
sem_dynamic_timeUsed = sem(subj_dynamic_timeUsed./60);
sem_flies_timeUsed = sem(subj_flies_timeUsed./60);

avg_stationary_respTime = mean(subj_stationary_respTime - 0.15);
avg_dynamic_respTime = mean(subj_dynamic_respTime - 0.15);
avg_flies_respTime = mean(subj_flies_respTime - 0.15);
sem_stationary_respTime = sem(subj_stationary_respTime - 0.15);
sem_dynamic_respTime = sem(subj_dynamic_respTime - 0.15);
sem_flies_respTime = sem(subj_flies_respTime - 0.15);

avg_stationary_trackTime = mean(subj_stationary_trackTime);
avg_dynamic_trackTime = mean(subj_dynamic_trackTime);
avg_flies_trackTime = mean(subj_flies_trackTime);
sem_stationary_trackTime = sem(subj_stationary_trackTime);
sem_dynamic_trackTime = sem(subj_dynamic_trackTime);
sem_flies_trackTime = sem(subj_flies_trackTime);


%% calculate percentage:

timeUsed_dynamic_vs_stationary = (avg_dynamic_timeUsed_min - avg_stationary_timeUsed_min)/avg_stationary_timeUsed_min*100
timeUsed_flies_vs_dynamic = (avg_flies_timeUsed_min - avg_dynamic_timeUsed_min)/avg_dynamic_timeUsed_min*100
timeUsed_flies_vs_stationary = (avg_flies_timeUsed_min - avg_stationary_timeUsed_min)/avg_stationary_timeUsed_min*100


respTime_dynamic_vs_stationary = (avg_dynamic_respTime - avg_stationary_respTime)/avg_stationary_respTime*100
respTime_flies_vs_dynamic = (avg_flies_respTime - avg_dynamic_respTime)/avg_dynamic_respTime*100
respTime_flies_vs_stationary = (avg_flies_respTime - avg_stationary_respTime)/avg_stationary_respTime*100


trackTime_dynamic_vs_stationary = (avg_dynamic_trackTime - avg_stationary_trackTime)/avg_stationary_trackTime*100
trackTime_flies_vs_dynamic = (avg_flies_trackTime - avg_dynamic_trackTime)/avg_dynamic_trackTime*100
trackTime_flies_vs_stationary = (avg_flies_trackTime - avg_stationary_trackTime)/avg_stationary_trackTime*100

%% plot
% 
CData = {[0.4940, 0.1840, 0.5560],[0.8500, 0.3250, 0.0980],[0, 0.4470, 0.7410]};

conditions = [        {'Stationary'},...
    {'Dynamic'   },...
    {'Crowded dynamic'     }];

markerType = 'o';%{'o','+','x','square','d','pentagram'};
% % plotPosition = linspace(-0.15,0.15,6);
markSize = 8;
lineWidth = 2;

%% plot time used for each block
plotPos = [1 2 3];

figure;
hold on;
plot(plotPos(1), avg_stationary_timeUsed_min, 'o', 'Color', CData{1}, 'LineWidth', 3,'MarkerSize',9);
plot(plotPos(2), avg_dynamic_timeUsed_min, 'o', 'Color', CData{2}, 'LineWidth', 3,'MarkerSize',9);
plot(plotPos(3), avg_flies_timeUsed_min, 'o', 'Color', CData{3}, 'LineWidth', 3,'MarkerSize',9);

errorbar(plotPos(1), avg_stationary_timeUsed_min, sem_stationary_timeUsed, 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
errorbar(plotPos(2), avg_dynamic_timeUsed_min, sem_dynamic_timeUsed, 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
errorbar(plotPos(3), avg_flies_timeUsed_min, sem_flies_timeUsed, 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
hold off;

set(gca, 'XTick', plotPos, 'XTickLabel', conditions);
set(gca,'FontSize',18)
xlim([0 4])
ylim([0 8])
title('Time used for each block')
ylabel('Time used (min)')



%% plot response time


figure;
hold on;
plot(plotPos(1), avg_stationary_respTime, 'o', 'Color', CData{1}, 'LineWidth', 3,'MarkerSize',9);
plot(plotPos(2), avg_dynamic_respTime, 'o', 'Color', CData{2}, 'LineWidth', 3,'MarkerSize',9);
plot(plotPos(3), avg_flies_respTime, 'o', 'Color', CData{3}, 'LineWidth', 3,'MarkerSize',9);

errorbar(plotPos(1), avg_stationary_respTime, sem_stationary_respTime, 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
errorbar(plotPos(2), avg_dynamic_respTime, sem_dynamic_respTime, 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
errorbar(plotPos(3), avg_flies_respTime, sem_flies_respTime, 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
hold off;

set(gca, 'XTick', plotPos, 'XTickLabel', conditions);
set(gca,'FontSize',18)
xlim([0 4])
ylim([0 2.5])
title('Response time per trial')
ylabel('Response time (sec)')



%% plot track time

figure;
hold on;
plot(plotPos(1), avg_stationary_trackTime, 'o', 'Color', CData{1}, 'LineWidth', 3,'MarkerSize',9);
plot(plotPos(2), avg_dynamic_trackTime, 'o', 'Color', CData{2}, 'LineWidth', 3,'MarkerSize',9);
plot(plotPos(3), avg_flies_trackTime, 'o', 'Color', CData{3}, 'LineWidth', 3,'MarkerSize',9);

errorbar(plotPos(1), avg_stationary_trackTime, sem_stationary_trackTime, 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
errorbar(plotPos(2), avg_dynamic_trackTime, sem_dynamic_trackTime, 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
errorbar(plotPos(3), avg_flies_trackTime, sem_flies_trackTime, 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
hold off;

set(gca, 'XTick', plotPos, 'XTickLabel', conditions);
set(gca,'FontSize',18)
xlim([0 4])
ylim([0 4])
title('Time spent tracking per trial')
ylabel('Time spent tracking (sec)')

%%
% for aa = whichRun%1:numSubj
%     plot(1,subj_stationary_left(aa),markerType,'MarkerSize',markSize,'MarkerEdgeColor','k','LineWidth',lineWidth)
%     plot(2,subj_stationary_right(aa),markerType,'MarkerSize',markSize,'MarkerEdgeColor','k','LineWidth',lineWidth)
%     plot(3,subj_dynamic_left(aa),markerType,'MarkerSize',markSize,'MarkerEdgeColor','k','LineWidth',lineWidth)
%     plot(4,subj_dynamic_right(aa),markerType,'MarkerSize',markSize,'MarkerEdgeColor','k','LineWidth',lineWidth)
%     plot(5,subj_flies_left(aa),markerType,'MarkerSize',markSize,'MarkerEdgeColor','k','LineWidth',lineWidth)
%     plot(6,subj_flies_right(aa),markerType,'MarkerSize',markSize,'MarkerEdgeColor','k','LineWidth',lineWidth)
% 
%     plot(1:6,[subj_stationary_left(aa),subj_stationary_right(aa),subj_dynamic_left(aa),subj_dynamic_right(aa),...
%         subj_flies_left(aa),subj_flies_right(aa)],'k-','LineWidth',2)
%     % pause(0.5)
% end
% 
% disp('')

% %% plot movie
% 
% markerType = 'o';%{'o','+','x','square','d','pentagram'};
% % % plotPosition = linspace(-0.15,0.15,6);
% markSize = 8;
% lineWidth = 2;
% % 
% fig = figure;
% % Create a VideoWriter object to write video to a file
% video = VideoWriter('plotting_process.avi');
% video.FrameRate = 1; % Set frame rate to 1 frame per second
% open(video);
% 
% 
% 
% b = bar(conditions,[avg_stationary_left,avg_stationary_right ,avg_dynamic_left,avg_dynamic_right,avg_flies_left,avg_flies_right]);
% 
% b.FaceColor = 'flat'; % This allows individual coloring of each bar
% b.EdgeColor = 'flat';
% b.CData(1,:) = cell2mat(CData(1)); 
% b.CData(2,:) = cell2mat(CData(1)); 
% b.CData(3,:) = cell2mat(CData(2)); 
% b.CData(4,:) = cell2mat(CData(2)); 
% b.CData(5,:) = cell2mat(CData(3)); 
% b.CData(6,:) = cell2mat(CData(3)); 
% b.FaceAlpha = 0.6;
% b.EdgeAlpha = 0.6;
% 
% set(gca,'FontSize',18)
% set(gca, 'xticklabel', conditions)
% ylim([0 8])
% title('Effect of Background on Crowding Thresholds')
% ylabel('Crowding Thresholds (deg)')
% frame = getframe(fig); % Capture the current figure as a frame
% writeVideo(video, frame); % Write the frame to the video
% 
% hold on;
% % errorbar(1, avg_noflies_peek, sem_noflies_peek, 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% % errorbar(2, avg_3flies_peek, sem_3flies_peek, 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% % errorbar(3, avg_5flies_peek, sem_5flies_peek, 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
% x = 1:6;
% y_data = {};
% for aa = whichRun%1:numSubj
%     y_data{end+1} = [subj_stationary_left(aa),subj_stationary_right(aa),subj_dynamic_left(aa),subj_dynamic_right(aa),...
%             subj_flies_left(aa),subj_flies_right(aa)];
% end
% 
% % Preallocate handles for the lines
% h = gobjects(length(whichRun), 1);
% 
% % Loop through each line, plotting and then erasing it
% for ii = 1:length(whichRun)
%     h(ii) = plot(x, y_data{ii}, 'o-','Color','#000000','LineWidth', 2);
%     frame = getframe(fig); % Capture the current figure as a frame
%     writeVideo(video, frame); % Write the frame to the video
%     pause(1); % Pause for 1 second to create a movie effect
%     if ii < length(whichRun)
%         delete(h(ii)); % Delete the current line if it's not the last one
%     end
% end
% 
% % Finally, plot all lines together
% hold on;
% for ii = 1:length(whichRun)
%     plot(x, y_data{ii}, 'o-','Color','#000000','LineWidth', 2);
% end
% hold off;
% frame = getframe(fig); % Capture the current figure as a frame
% writeVideo(video, frame); % Write the frame to the video
% 
% close(video);
% 
% disp('')
%% plot with error bar


% figure;
% b = bar(conditions,[avg_stationary_left,avg_stationary_right ,avg_dynamic_left,avg_dynamic_right,avg_flies_left,avg_flies_right]);
% 
% b.FaceColor = 'flat'; % This allows individual coloring of each bar
% b.EdgeColor = 'flat';
% b.CData(1,:) = cell2mat(CData(1)); 
% b.CData(2,:) = cell2mat(CData(1)); 
% b.CData(3,:) = cell2mat(CData(2)); 
% b.CData(4,:) = cell2mat(CData(2)); 
% b.CData(5,:) = cell2mat(CData(3)); 
% b.CData(6,:) = cell2mat(CData(3)); 
% b.FaceAlpha = 0.6;
% b.EdgeAlpha = 0.6;
% 
% set(gca,'FontSize',18)
% set(gca, 'xticklabel', conditions)
% ylim([0,5])
% title(sprintf('Effect of Background on Crowding Thresholds (N = %d)',numSubj))
% ylabel('Crowding Thresholds (Deg)')
% 
% hold on;
% errorbar(1, avg_stationary_left, sem_stationary_left, 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% errorbar(2, avg_stationary_right, sem_stationary_right, 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% errorbar(3, avg_dynamic_left, sem_dynamic_left, 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% errorbar(4, avg_dynamic_right, sem_dynamic_right, 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% errorbar(5, avg_flies_left, sem_flies_left, 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
% errorbar(6, avg_flies_right, sem_flies_right, 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);

% errorbar(1, avg_noflies_peek, sem_noflies_peek, 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% errorbar(2, avg_3flies_peek, sem_3flies_peek, 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% errorbar(3, avg_5flies_peek, sem_5flies_peek, 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);


%% plot median with bootstrapping, log thresholds
% subj_stationary_left = subj_stationary_left(~isnan(subj_stationary_left));
% subj_stationary_right = subj_stationary_right(~isnan(subj_stationary_right));
% subj_dynamic_left = subj_dynamic_left(~isnan(subj_dynamic_left));
% subj_dynamic_right = subj_dynamic_right(~isnan(subj_dynamic_right));
% subj_flies_left = subj_flies_left(~isnan(subj_flies_left));
% subj_flies_right = subj_flies_right(~isnan(subj_flies_right));
% 
% thresholds = {
%     subj_stationary_left,...
%     subj_stationary_right,...
%     subj_dynamic_left, ...
%     subj_dynamic_right, ...
%     subj_flies_left, ...
%     subj_flies_right
% };
% 
% % Number of bootstrap samples
% numBootstraps = 10000;
% 
% % Confidence level (68%)
% confLevel = 0.68;
% 
% % Initialize array to store bootstrap means
% bootstrapMeds = cell(1, length(thresholds));
% 
% 
% % Perform bootstrapping for each condition using the bootstrap function
% for cond = 1:length(thresholds)
%     data = thresholds{cond};
% 
%     % Take the logarithm of the thresholds
%     logData = log(data);
% 
%     % Define a function handle to calculate the median
%     medFunc = @(data) median(data);
% 
%     % Generate bootstrap samples and calculate medians
%     bootMeds = bootstrp(numBootstraps, medFunc, logData);
% 
%     bootstrapMeds{cond} = bootMeds;
% end
% 
% % Calculate the 68% confidence intervals
% confIntervals = zeros(length(thresholds), 2);
% medLogThresholds = zeros(1, length(thresholds));
% for cond = 1:length(thresholds)
%     sortedMeds = sort(bootstrapMeds{cond});
%     lowerBound = sortedMeds(round((1 - confLevel) / 2 * numBootstraps));
%     upperBound = sortedMeds(round((1 + confLevel) / 2 * numBootstraps));
%     confIntervals(cond, :) = [lowerBound, upperBound];
%     medLogThresholds(cond) = median(log(thresholds{cond})); % Mean of log-transformed data
% end
% 
% % Convert mean thresholds and error bars back from log scale
% medThresholds = exp(medLogThresholds);
% lowerError = medThresholds - exp(confIntervals(:, 1)');
% upperError = exp(confIntervals(:, 2)') - medThresholds;
% errorBars = [lowerError; upperError];
% 
% % Plot the bar graph with error bars
% figure;
% b = bar(medThresholds);
% set(gca, 'YScale', 'log');
% set(gca,'FontSize',18)
% set(gca, 'xticklabel', conditions)
% 
% b.FaceColor = 'flat'; 
% b.EdgeColor = 'flat';
% b.CData(1,:) = cell2mat(CData(1)); 
% b.CData(2,:) = cell2mat(CData(1)); 
% b.CData(3,:) = cell2mat(CData(2)); 
% b.CData(4,:) = cell2mat(CData(2)); 
% b.CData(5,:) = cell2mat(CData(3)); 
% b.CData(6,:) = cell2mat(CData(3)); 
% b.FaceAlpha = 0.6;
% b.EdgeAlpha = 0.6;
% 
% 
% 
% % xlabel('Condition');
% ylabel('Crowding Thresholds (deg)');
% title('Effect of Background on Crowding Thresholds');
% 
% hold on;
% errorbar(1, medThresholds(1), lowerError(1),upperError(1), 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% errorbar(2, medThresholds(2), lowerError(2),upperError(2), 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% errorbar(3, medThresholds(3), lowerError(3),upperError(3), 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% errorbar(4, medThresholds(4), lowerError(4),upperError(4), 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% errorbar(5, medThresholds(5), lowerError(5),upperError(5), 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
% errorbar(6, medThresholds(6), lowerError(6),upperError(6), 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
% hold off;
% ylim([1,5])
% 
% 
% 


% %% plot with bootstrapping, log thresholds
% subj_stationary_left = subj_stationary_left(~isnan(subj_stationary_left));
% subj_stationary_right = subj_stationary_right(~isnan(subj_stationary_right));
% subj_dynamic_left = subj_dynamic_left(~isnan(subj_dynamic_left));
% subj_dynamic_right = subj_dynamic_right(~isnan(subj_dynamic_right));
% subj_flies_left = subj_flies_left(~isnan(subj_flies_left));
% subj_flies_right = subj_flies_right(~isnan(subj_flies_right));
% 
% thresholds = {
%     subj_stationary_left,...
%     subj_stationary_right,...
%     subj_dynamic_left, ...
%     subj_dynamic_right, ...
%     subj_flies_left, ...
%     subj_flies_right
% };
% 
% % Number of bootstrap samples
% numBootstraps = 10000;
% 
% % Confidence level (68%)
% confLevel = 0.68;
% 
% % Initialize array to store bootstrap means
% bootstrapMeans = cell(1, length(thresholds));
% 
% 
% % Perform bootstrapping for each condition using the bootstrap function
% for cond = 1:length(thresholds)
%     data = thresholds{cond};
% 
%     % Take the logarithm of the thresholds
%     logData = log(data);
% 
%     % Define a function handle to calculate the mean
%     meanFunc = @(data) mean(data);
% 
%     % Generate bootstrap samples and calculate means
%     bootMeans = bootstrp(numBootstraps, meanFunc, logData);
% 
%     bootstrapMeans{cond} = bootMeans;
% end
% 
% % Calculate the 68% confidence intervals
% confIntervals = zeros(length(thresholds), 2);
% meanLogThresholds = zeros(1, length(thresholds));
% for cond = 1:length(thresholds)
%     sortedMeans = sort(bootstrapMeans{cond});
%     lowerBound = sortedMeans(round((1 - confLevel) / 2 * numBootstraps));
%     upperBound = sortedMeans(round((1 + confLevel) / 2 * numBootstraps));
%     confIntervals(cond, :) = [lowerBound, upperBound];
%     meanLogThresholds(cond) = mean(log(thresholds{cond})); % Mean of log-transformed data
% end
% 
% % Convert mean thresholds and error bars back from log scale
% meanThresholds = exp(meanLogThresholds);
% lowerError = meanThresholds - exp(confIntervals(:, 1)');
% upperError = exp(confIntervals(:, 2)') - meanThresholds;
% errorBars = [lowerError; upperError];
% %%
% % Plot the bar graph with error bars
% figure;
% hold on;
% plot([1 1.6], [meanThresholds(1) meanThresholds(2)], 'o-', 'Color', CData{1}, 'LineWidth', 3,'MarkerSize',9);
% % plot(1.6, meanThresholds(2), 'o', 'Color', CData{1}, 'LineWidth', 3,'MarkerSize',9);
% plot([3 3.6], [meanThresholds(3) meanThresholds(4)], 'o-', 'Color', CData{2}, 'LineWidth', 3,'MarkerSize',9);
% % plot(3.6, meanThresholds(4), 'o', 'Color', CData{2}, 'LineWidth', 3,'MarkerSize',9);
% plot([5 5.6], [meanThresholds(5) meanThresholds(6)], 'o-', 'Color', CData{3}, 'LineWidth', 3,'MarkerSize',9);
% % plot(5.6, meanThresholds(6), 'o', 'Color', CData{3}, 'LineWidth', 3,'MarkerSize',9);
% 
% set(gca, 'YScale', 'log');
% set(gca,'FontSize',18)
% set(gca, 'XTick', [1 1.6 3 3.6 5 5.6], 'XTickLabel', conditions);
% 
% xlim([0 7])
% xtickangle(45);
% % xlabel('Condition');
% ylabel('Crowding Thresholds (deg)');
% % title('Effect of Background on Crowding Thresholds');
% 
% 
% errorbar(1, meanThresholds(1), lowerError(1),upperError(1), 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% errorbar(1.6, meanThresholds(2), lowerError(2),upperError(2), 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% errorbar(3, meanThresholds(3), lowerError(3),upperError(3), 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% errorbar(3.6, meanThresholds(4), lowerError(4),upperError(4), 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% errorbar(5, meanThresholds(5), lowerError(5),upperError(5), 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
% errorbar(5.6, meanThresholds(6), lowerError(6),upperError(6), 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
% hold off;
% ylim([1,5])
% 
% %%
% % Plot the bar graph with error bars
% figure;
% hold on;
% plot([1 2], [meanThresholds(1) meanThresholds(2)], 'o-', 'Color', CData{1}, 'LineWidth', 3,'MarkerSize',9);
% plot([1 2], [meanThresholds(3) meanThresholds(4)], 'o-', 'Color', CData{2}, 'LineWidth', 3,'MarkerSize',9);
% plot([1 2], [meanThresholds(5) meanThresholds(6)], 'o-', 'Color', CData{3}, 'LineWidth', 3,'MarkerSize',9);
% 
% set(gca, 'YScale', 'log');
% set(gca,'FontSize',18)
% set(gca, 'XTick', [1 2], 'XTickLabel', [{'Left'} {'Right'}]);
% 
% xlim([0 3])
% % xtickangle(45);
% % xlabel('Condition');
% ylabel('Crowding Thresholds (deg)');
% % title('Effect of Background on Crowding Thresholds');
% 
% 
% errorbar(1, meanThresholds(1), lowerError(1),upperError(1), 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% errorbar(2, meanThresholds(2), lowerError(2),upperError(2), 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% errorbar(1, meanThresholds(3), lowerError(3),upperError(3), 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% errorbar(2, meanThresholds(4), lowerError(4),upperError(4), 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% errorbar(1, meanThresholds(5), lowerError(5),upperError(5), 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
% errorbar(2, meanThresholds(6), lowerError(6),upperError(6), 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
% hold off;
% ylim([1,5])
% 
% 
% 
% 
% %% plot minimum thresholds 
% % 
% % actualMins = [min(subj_stationary_left),...
% %     min(subj_stationary_right),...
% %     min(subj_dynamic_left), ...
% %     min(subj_dynamic_right), ...
% %     min(subj_flies_left), ...
% %     min(subj_flies_right)];
% % 
% % % Plot the bar graph with error bars
% % figure;
% % b = bar(actualMins);
% % set(gca, 'YScale', 'log');
% % set(gca,'FontSize',18)
% % set(gca, 'xticklabel', conditions)
% % 
% % b.FaceColor = 'flat'; 
% % b.EdgeColor = 'flat';
% % b.CData(1,:) = cell2mat(CData(1)); 
% % b.CData(2,:) = cell2mat(CData(1)); 
% % b.CData(3,:) = cell2mat(CData(2)); 
% % b.CData(4,:) = cell2mat(CData(2)); 
% % b.CData(5,:) = cell2mat(CData(3)); 
% % b.CData(6,:) = cell2mat(CData(3)); 
% % b.FaceAlpha = 0.6;
% % b.EdgeAlpha = 0.6;
% % 
% % 
% % ylim([0.3,3])
% % ylabel('Crowding Thresholds (deg)');
% % title('Effect of Background on Crowding Thresholds');
% 
% % hold on;
% % errorbar(1, actualMins(1), lowerError(1),upperError(1), 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% % errorbar(2, actualMins(2), lowerError(2),upperError(2), 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2);
% % errorbar(3, actualMins(3), lowerError(3),upperError(3), 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% % errorbar(4, actualMins(4), lowerError(4),upperError(4), 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2);
% % errorbar(5, actualMins(5), lowerError(5),upperError(5), 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
% % errorbar(6, actualMins(6), lowerError(6),upperError(6), 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2);
% % hold off;
% 
% %% plot histograms
% 
% subj_stationary_left = subj_stationary_left./10;
% subj_stationary_right = subj_stationary_right./10;
% subj_dynamic_left = subj_dynamic_left./10;
% subj_dynamic_right = subj_dynamic_right./10;
% subj_flies_left = subj_flies_left./10;
% subj_flies_right = subj_flies_right./10;
% 
% 
% %%
% ylimit = [0 11];
% % bwidth = 0.1;
% xlimit = [0.03 1];
% 
% minData = 0.05;
% maxData = 1;
% numBins = 35; % Define the number of bins you want
% binEdges = logspace(log10(minData), log10(maxData), numBins);
% 
% figure;
% subplot(3,1,1)
% plotHistLog([subj_stationary_left,subj_stationary_right],binEdges,xlimit,ylimit,'Stationary Fixation',cell2mat(CData(1)))
% hold on;xline(0.1,'r--','LineWidth',2);hold off;
% 
% subplot(3,1,2)
% plotHistLog([subj_dynamic_left,subj_dynamic_right],binEdges,xlimit,ylimit,'Dynamic Fixation',cell2mat(CData(2)))
% hold on;xline(0.1,'r--','LineWidth',2);hold off;
% 
% 
% subplot(3,1,3)
% plotHistLog([subj_flies_left,subj_flies_right],binEdges,xlimit,ylimit,'Crowded Dynamic Fixation',cell2mat(CData(3)))
% xlabel('Bouma Factor b')
% hold on;xline(0.1,'r--','LineWidth',2);hold off;
% 
% 
% 
% %% who have small thresholds?
% 
% 
% find(subj_stationary_left<1)
% find(subj_stationary_right<1)
% find(subj_dynamic_left<1)
% find(subj_dynamic_right<1)
% find(subj_flies_left<1)
% find(subj_flies_right<1)
% 
% smallThresholdSesh = [1 13 18 26 36];

%%
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
function [sem] = sem(data)
    n = length(data);
    std_dev = std(data);
    sem = std_dev / sqrt(n);
end

