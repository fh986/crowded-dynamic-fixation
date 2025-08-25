% plot_time_spent.m
% Author: Helen Hu
% Last edited: Apr 9th, 2025

% This script plots the time spent on each block (condition)

clc;
clear all;
close all;


%% set up files

% Get current directory and move one level up
scriptDir = pwd;
repoDir = fileparts(scriptDir);
% add some functions for gaze analysis
addpath(fullfile(repoDir,'Gaze_Package'));
% data files for gaze
mydir = fullfile(repoDir,'data_include_forGaze_20ppl_paired');
 
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

    easyeyes = readtable([mydir filesep cursorFiles{subj}],'VariableNamingRule','preserve');
    mainOutput = readtable([mydir filesep mainFiles{subj}],'VariableNamingRule','preserve');

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

timeUsed_dynamic_vs_stationary = (avg_dynamic_timeUsed_min - avg_stationary_timeUsed_min)/avg_stationary_timeUsed_min*100;
timeUsed_flies_vs_dynamic = (avg_flies_timeUsed_min - avg_dynamic_timeUsed_min)/avg_dynamic_timeUsed_min*100;
timeUsed_flies_vs_stationary = (avg_flies_timeUsed_min - avg_stationary_timeUsed_min)/avg_stationary_timeUsed_min*100;


respTime_dynamic_vs_stationary = (avg_dynamic_respTime - avg_stationary_respTime)/avg_stationary_respTime*100;
respTime_flies_vs_dynamic = (avg_flies_respTime - avg_dynamic_respTime)/avg_dynamic_respTime*100;
respTime_flies_vs_stationary = (avg_flies_respTime - avg_stationary_respTime)/avg_stationary_respTime*100;


trackTime_dynamic_vs_stationary = (avg_dynamic_trackTime - avg_stationary_trackTime)/avg_stationary_trackTime*100;
trackTime_flies_vs_dynamic = (avg_flies_trackTime - avg_dynamic_trackTime)/avg_dynamic_trackTime*100;
trackTime_flies_vs_stationary = (avg_flies_trackTime - avg_stationary_trackTime)/avg_stationary_trackTime*100;

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

errorbar(plotPos(1), avg_stationary_timeUsed_min, sem_stationary_timeUsed, 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2, 'CapSize', 0);
errorbar(plotPos(2), avg_dynamic_timeUsed_min, sem_dynamic_timeUsed, 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2, 'CapSize', 0);
errorbar(plotPos(3), avg_flies_timeUsed_min, sem_flies_timeUsed, 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2, 'CapSize', 0);
hold off;

set(gca, 'XTick', plotPos, 'XTickLabel', conditions);
set(gca,'FontSize',18)
xlim([0 4])
ylim([0 8])
title('Time used for each block')
ylabel('Duration (min)')



%% plot response time


figure;
hold on;
plot(plotPos(1), avg_stationary_respTime, 'o', 'Color', CData{1}, 'LineWidth', 3,'MarkerSize',9);
plot(plotPos(2), avg_dynamic_respTime, 'o', 'Color', CData{2}, 'LineWidth', 3,'MarkerSize',9);
plot(plotPos(3), avg_flies_respTime, 'o', 'Color', CData{3}, 'LineWidth', 3,'MarkerSize',9);

errorbar(plotPos(1), avg_stationary_respTime, sem_stationary_respTime, 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2, 'CapSize', 0);
errorbar(plotPos(2), avg_dynamic_respTime, sem_dynamic_respTime, 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2, 'CapSize', 0);
errorbar(plotPos(3), avg_flies_respTime, sem_flies_respTime, 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2, 'CapSize', 0);
hold off;

set(gca, 'XTick', plotPos, 'XTickLabel', conditions);
set(gca,'FontSize',18)
xlim([0 4])
ylim([0 2.5])
title('Response time per trial')
ylabel('Duration (sec)')



%% plot track time

figure;
hold on;
plot(plotPos(1), avg_stationary_trackTime, 'o', 'Color', CData{1}, 'LineWidth', 3,'MarkerSize',9);
plot(plotPos(2), avg_dynamic_trackTime, 'o', 'Color', CData{2}, 'LineWidth', 3,'MarkerSize',9);
plot(plotPos(3), avg_flies_trackTime, 'o', 'Color', CData{3}, 'LineWidth', 3,'MarkerSize',9);

errorbar(plotPos(1), avg_stationary_trackTime, sem_stationary_trackTime, 'LineStyle', 'none', 'Color', CData{1}, 'LineWidth', 2, 'CapSize', 0);
errorbar(plotPos(2), avg_dynamic_trackTime, sem_dynamic_trackTime, 'LineStyle', 'none', 'Color', CData{2}, 'LineWidth', 2, 'CapSize', 0);
errorbar(plotPos(3), avg_flies_trackTime, sem_flies_trackTime, 'LineStyle', 'none', 'Color', CData{3}, 'LineWidth', 2, 'CapSize', 0);
hold off;

set(gca, 'XTick', plotPos, 'XTickLabel', conditions);
set(gca,'FontSize',18)
xlim([0 4])
ylim([0 4])
title('Time spent tracking per trial')
ylabel('Duration (sec)')

%%
function [sem] = sem(data)
    n = length(data);
    std_dev = std(data);
    sem = std_dev / sqrt(n);
end

