% plot_all_trials_xy.m
% Author: Helen Hu
% Last modified: February 28th, 2025

% This file is used to generate Figure 3 of Hu et al.
% It plots a participant's gaze positions during three fixation conditions:
% Stationary, Dynamic, and Crowded dynamic ("flies").

%%
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


%% gaze analyses

for subj = 33 % plots the participant shown in Figure 3 of Hu et al.

    easyeyes = readtable([mydir filesep cursorFiles{subj}],'VariableNamingRule','preserve');
    mainOutput = readtable([mydir filesep mainFiles{subj}],'VariableNamingRule','preserve');
    eyelink = readtable([mydir filesep eyelinkFiles{subj}],'VariableNamingRule','preserve');
    
    % count stimuli presentations
    trialRoutine = strcmp(easyeyes.targetBool,'TRUE');
    diff_trialRout = diff(trialRoutine);
    stim_on = find(diff_trialRout == 1)+1;
    trialStimInfo = easyeyes(stim_on,:);

    % count tracking
    trackRoutine = strcmp(easyeyes.crosshairBool,'TRUE');
    diff_trackRout = diff(trackRoutine);
    track_on = find(diff_trackRout == 1)+1;
    trackingInfo = easyeyes(track_on,:);

    % special cases
    if contains(eyelinkFiles{subj}, 'MarcoLai061024')
        track_on(198) = [];
    end
    if contains(eyelinkFiles{subj}, 'KevinHong061624')
        track_on(15) = [];
    end
    if(length(stim_on) ~= length(track_on))
        track_on = [1;track_on];
    end

    % check if the row numbers for the stimulus presentations and 
    % trial-start timestamps are reasonable
    assert(~any((stim_on - track_on)<=0));
   
    % pixel to deg calculations
    screenWidthCm = unique(mainOutput.screenWidthByObjectCm(~isnan(mainOutput.screenWidthByObjectCm)));
    screenWidthPx = unique(mainOutput.screenWidthPx(~isnan(mainOutput.screenWidthPx)));
    distance = 40; %cm
    [PixelPerDeg,px_to_deg] = convertPxDeg(screenWidthPx,screenWidthCm,distance);

    % gaze data analysis
    timeframe_starttrack = 0.4; % get rid of the first 400 ms
    gaze_x_mtx = NaN(length(stim_on),95); % 100 frames per sec
    gaze_y_mtx = NaN(length(stim_on),95);

    % determine what trials to include
    trials_include = 1:length(stim_on);
    if contains(eyelinkFiles{subj}, 'KevinHong061624')
        trials_include = 2:length(stim_on);
    end    

    for s = trials_include

        %%%%%%%%%%%%%%%%%%%%
        % Drift correction %
        %%%%%%%%%%%%%%%%%%%%
        % focus on the full tracking period and correct for drift
        stim_timestamp = easyeyes(stim_on(s),:).posixTimeSec; 
        track_timestamp = easyeyes(track_on(s),:).posixTimeSec;
        trial_ee_track = easyeyes.posixTimeSec > (stim_timestamp - timeframe_starttrack) & easyeyes.posixTimeSec < (stim_timestamp);
        currenttrial_ee_track = easyeyes(trial_ee_track,:);
        % correct for a time difference between eyelink and easyeyes data
        trial_el_track = eyelink.t1+14400 > (stim_timestamp -  timeframe_starttrack) & eyelink.t1+14400 < (stim_timestamp);
        currenttrial_el_track = eyelink(trial_el_track,:);
        framesPerSec = 60;
        eyedrift = calcDrift(currenttrial_ee_track,currenttrial_el_track,PixelPerDeg,framesPerSec);

        %%%%%%%%%%%%%%%%%%%%
        % Drift correction %
        %%%%%%%%%%%%%%%%%%%%
        % record the x and y positions of gaze      
        trial_ee_stim = find(easyeyes.posixTimeSec == stim_timestamp);
        currenttrial_ee_stimon = easyeyes(trial_ee_stim,:);
        nearpointPx = str2num(cell2mat(currenttrial_ee_stimon.nearpointXYPx));
        last_crosshairPx = str2num(cell2mat(currenttrial_ee_stimon.crosshairPositionXYPx)) - nearpointPx;
        trial_eyelink = eyelink.t1+14400 > (stim_timestamp - 0.5) & eyelink.t1+14400 < (stim_timestamp + 0.5);
        currenttrial_el = eyelink(trial_eyelink,:);
        timestamps = currenttrial_el.t1;
  
        gazePx = [];
        gazePx(:,1) = currenttrial_el.gazeXYPix_1  - last_crosshairPx(1)- eyedrift(1);
        gazePx(:,2) = -currenttrial_el.gazeXYPix_2  - last_crosshairPx(2)- eyedrift(2); % note the y-flip
        [newGazePx,nantrialBool] = ignore_blink(gazePx,PixelPerDeg,framesPerSec);
        if nantrialBool
            disp('Warning: nantrial')
        end

        gazeXDeg = newGazePx(:,1) ./ PixelPerDeg;
        gazeYDeg = newGazePx(:,2) ./ PixelPerDeg;

        % record just the period we are plotting
        relative_time = timestamps - stim_timestamp + 14400;
        [val,idx_zero_time] = min(abs(relative_time));
        idx_record = (idx_zero_time-47):(idx_zero_time+47);

        gaze_x_mtx(s,:) = gazeXDeg(idx_record);
        gaze_y_mtx(s,:) = gazeYDeg(idx_record);

        % record the location of the stimulus: left or right meridian
        condNames = mainOutput.conditionName;
        rm = cellfun('isempty',condNames);
        condNames(rm) = [];
        rm = [71 72 143 144 215 216];
        condNames(rm) = [];
        rightBool = contains(condNames,'Right');

    end

    % plot x and y components of gaze
    block_sequence = unique(mainOutput.blockShuffleGroups1,'stable');
    rm = cellfun("isempty",block_sequence);
    block_sequence(rm) = [];
    frame_interval = mean(diff(eyelink.t1(1:10*60)));
    xValues = linspace(-frame_interval*47,frame_interval*47,95);

    % def blocks
    firstBlockEnd = 70;
    secondBlockStart = 71;
    secondBlockEnd = 140;
    thirdBlockStart = 141;
    thirdBlockEnd = 210;
    % manually correct for special cases
    if contains(eyelinkFiles{subj}, 'MarcoLai061024')
        thirdBlockEnd = 209;
    end
    if contains(eyelinkFiles{subj}, 'KevinHong061624')
        firstBlockEnd = 69;
        secondBlockStart = 70;
        secondBlockEnd = 139;
        thirdBlockStart = 140;
        thirdBlockEnd = 209;
    end

%%
    plotOnSameFig(xValues,gaze_x_mtx(1:firstBlockEnd,:),gaze_y_mtx(1:firstBlockEnd,:),rightBool(1:firstBlockEnd),block_sequence(1))
    plotOnSameFig(xValues,gaze_x_mtx(secondBlockStart:secondBlockEnd,:),gaze_y_mtx(secondBlockStart:secondBlockEnd,:),rightBool(secondBlockStart:secondBlockEnd),block_sequence(2))
    plotOnSameFig(xValues,gaze_x_mtx(thirdBlockStart:thirdBlockEnd,:),gaze_y_mtx(thirdBlockStart:thirdBlockEnd,:),rightBool(thirdBlockStart:thirdBlockEnd,:),block_sequence(3))

    
end

%% functions

function [] = plotOnSameFig(xValues,gaze_x_mtx,gaze_y_mtx,rightBool,title)
    figure;clf;

    subplot(2,1,1);
    hold on;
    % set up the graph
    xlabel('Time (s)');
    ylabel('X Position (deg)');
    set(gca,'Fontsize',18)
    xlim([-0.15 0.3])
    ylim([-10 10])
    yline(10,'k--','LineWidth',2)
    yline(-10,'k--','LineWidth',2)

    % Create a patch to highlight the stimulus onset
    y_limits = ylim();
    onset = 0; % Example: stimulus onset at 2 seconds
    offset = 0.15; % Example: stimulus offset at 4 seconds
    stimulus_height = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)];
    stimulus_duration = [onset offset offset onset];
    patch(stimulus_duration, stimulus_height, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

    % plot gaze positions
    for trial = 1:size(gaze_x_mtx, 1)
        if rightBool(trial) == 1
            lineColor = 'r';
        else
            lineColor = 'g';
        end
        xx = plot(xValues,gaze_x_mtx(trial,:),lineColor,'LineWidth',2);
        xx.Color = [xx.Color, 0.3];
        
    end
    hold off;



    subplot(2,1,2);
    hold on;
    % set up figure
    xlabel('Time (s)');
    ylabel('Y Position (deg)');
    set(gca,'Fontsize',18)
    xlim([-0.15 0.3])      
    ylim([-10 10])

    % patch
    patch(stimulus_duration, stimulus_height, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

    % plot gaze
    for trial = 1:size(gaze_y_mtx, 1)
        if rightBool(trial) == 1
            lineColor = 'r';
        else
            lineColor = 'g';
        end
        yy = plot(xValues,gaze_y_mtx(trial,:),lineColor,'LineWidth',2);
        yy.Color = [yy.Color, 0.3];
    end
    yline(0,'k--','LineWidth',2)
    hold off;    


    sgtitle(title)
end


