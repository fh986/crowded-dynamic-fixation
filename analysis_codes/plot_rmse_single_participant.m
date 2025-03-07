% plot_rmse_single_participant.m
% Author: Helen Hu
% Last modified: March 6th, 2025

% This file calculates and plots the RMSE for each session. 
% 1) Drift: drift of eye tracking was corrected for each frame
% by subtracting the average gaze error during the last 400 of 
% the "track" period.
% 2) Gaze error for each trial: for each frame in a trial, 
% the distance between gaze and crosshair center was calculated
% >>> One number per frame (time point) of each trial
% 3) Gaze error RMSE across trials: for each time point, RMSE was taken
% across 70 trials for each condition
% >>> One number per frame (time point) of each condition
% 4) Plotting: gaze error RMSE over time was plotted for each condition


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

for subj = 33

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

    % account for special cases
    if contains(eyelinkFiles{subj}, 'ML1_M')
        track_on(198) = [];
    end
    if contains(eyelinkFiles{subj}, 'KH2_K')
        track_on(15) = [];
    end
    if(length(stim_on) ~= length(track_on))
        track_on = [1;track_on];
    end

    % check that the indices are reasonable
    assert(~any((stim_on - track_on)<=0));
   
    % pixel to deg
    screenWidthCm = unique(mainOutput.screenWidthByObjectCm(~isnan(mainOutput.screenWidthByObjectCm)));
    screenWidthPx = unique(mainOutput.screenWidthPx(~isnan(mainOutput.screenWidthPx)));
    distance = 40; %cm
    [PixelPerDeg,px_to_deg] = convertPxDeg(screenWidthPx,screenWidthCm,distance);

    timeframe_starttrack = 0.4; % only include the last 400 ms of tracking
    % for calculating eye drift
    recFrames = 21;
    gaze_err_mtx = NaN(length(stim_on),2*recFrames+1);

    % determine what trials to include
    trials_include = 1:length(stim_on);
    if contains(eyelinkFiles{subj}, 'KH2_K')
        trials_include = 2:length(stim_on);
    end
    
    for s = trials_include

        %%%%%%%%%%%%%%%%%%%%
        % Drift correction %
        %%%%%%%%%%%%%%%%%%%%
        % focus on the tracking period and correct for drift
        stim_timestamp = easyeyes(stim_on(s),:).posixTimeSec; 
        track_timestamp = easyeyes(track_on(s),:).posixTimeSec;
        trial_ee_track = easyeyes.posixTimeSec > (stim_timestamp - timeframe_starttrack) & easyeyes.posixTimeSec < (stim_timestamp);
        currenttrial_ee_track = easyeyes(trial_ee_track,:);
        % correct for a time difference between eyelink and easyeyes data
        trial_el_track = eyelink.t1+14400 > (stim_timestamp -  timeframe_starttrack) & eyelink.t1+14400 < (stim_timestamp);
        currenttrial_el_track = eyelink(trial_el_track,:);
        framesPerSec = 60;
        eyedrift = calcDrift(currenttrial_ee_track,currenttrial_el_track,PixelPerDeg,framesPerSec);

        %%%%%%%%%%%%%%%%%%
        % Gaze positions %
        %%%%%%%%%%%%%%%%%%
        % record the x and y positions of gaze  
        trial_ee_stim = easyeyes.posixTimeSec > (stim_timestamp - 0.5) & easyeyes.posixTimeSec < (stim_timestamp + 0.5);
        currenttrial_ee_stimon = easyeyes(trial_ee_stim,:);
        nearpointPx = str2num(cell2mat(currenttrial_ee_stimon.nearpointXYPx));
        crosshairPx = str2num(cell2mat(currenttrial_ee_stimon.crosshairPositionXYPx)) - nearpointPx;
        timestamps = currenttrial_ee_stimon.posixTimeSec;

        trial_eyelink = eyelink.t1+14400 > (stim_timestamp - 0.5) & eyelink.t1+14400 < (stim_timestamp + 0.5);
        currenttrial_el = eyelink(trial_eyelink,:);
        

        gazePx = [];
        for tt = 1:length(crosshairPx)
            [v,ind] = min(abs(currenttrial_el.t1 + 14400 - currenttrial_ee_stimon.posixTimeSec(tt)));
            gazePx(tt,1) = currenttrial_el.gazeXYPix_1(ind) - eyedrift(1);
            gazePx(tt,2) = -currenttrial_el.gazeXYPix_2(ind) - eyedrift(2);
        end


        [newGazePx,nantrialBool] = ignore_blink(gazePx,PixelPerDeg,framesPerSec);
        if nantrialBool
            disp(currenttrial_ee_stimon.conditionName(1))
        end

        diff_gaze_cross_Px = [];
        for tt = 1:length(crosshairPx)
            diff_gaze_cross_Px(tt) = norm(crosshairPx(tt) - newGazePx(tt));
        end

        diff_gaze_cross_Deg = diff_gaze_cross_Px./PixelPerDeg;

        relative_time = timestamps - stim_timestamp;
        [val,idx_zero_time] = min(abs(relative_time));
        idx_record = (idx_zero_time-recFrames):(idx_zero_time+recFrames);

        gaze_err_mtx(s,:) = diff_gaze_cross_Deg(idx_record);

        % left and right
        condNames = mainOutput.conditionName;
        rm = cellfun('isempty',condNames);
        condNames(rm) = [];
        rm = [71 72 143 144 215 216];
        condNames(rm) = [];
        rightBool = contains(condNames,'Right');

    end

    %% plot

    block_sequence = unique(mainOutput.blockShuffleGroups1,'stable');
    rm = cellfun("isempty",block_sequence);
    block_sequence(rm) = [];
    plot_sequence = {'Stationary','Dynamic','Flies'};
   
    frame_interval = mean(diff(easyeyes.posixTimeSec(1:10*60)));
    xValues = linspace(-frame_interval*recFrames,frame_interval*recFrames,2*recFrames+1);

    dictionary_lineColors = dictionary({'Stationary','Dynamic','Flies'},{'#7E2F8E','#D95319','#0072BD'});
    dictionary_lineStyle = dictionary({'Stationary','Dynamic','Flies'},{'-','--','-.'});

    figure;clf;
    hold on;
    % set up the graph
    xlabel('Time (s)');
    ylabel('Fixation error RMSE (deg)');
    set(gca,'Fontsize',18)
    xlim([-0.15 0.3])
    ylim([0 2])

    createPatch(0,0.15);

    for line = 1:3
        whichBlock = find(strcmp(plot_sequence{line},block_sequence));
        trialStart = (whichBlock-1)*70+1;
        trialEnd = whichBlock*70;
        % if contains(eyelinkFiles{subj}, 'ML1_M') && trialEnd == 210
        %     trialEnd = 209;
        % end
        % if contains(eyelinkFiles{subj}, 'KH2_K') 
        %     trialEnd = trialEnd-1;
        %     trialStart = trialStart+1;
        % end


        blockCond = block_sequence(whichBlock);
        if strcmp(blockCond,'Stationary')
            blockName = 'Stationary Fixation';
        elseif strcmp(blockCond,'Dynamic')
            blockName = 'Dynamic Fixation';
        elseif strcmp(blockCond,'Flies')
            blockName = 'Crowded Dynamic Fixation';
        else
            disp('WARNING');
        end



        plotRMSE(xValues,gaze_err_mtx(trialStart:trialEnd,:),cell2mat(dictionary_lineStyle(blockCond)),blockName,cell2mat(dictionary_lineColors(blockCond))) 
    end

    legend('Location','northwest');
    title(sprintf('Gaze RMSE over time, session #%d',subj))
    hold off;

    
end

%% functions

function [] = plotRMSE(xValues,gaze_mtx,lineStyle,condition,lineColor)
    % figure;clf;

    % y values: std across all trials over time
    squared_distances = gaze_mtx .^ 2;
    mean_squared_distances = nanmean(squared_distances, 1);
    rmse = sqrt(mean_squared_distances);

    % plot gaze positions
    ll = plot(xValues,rmse,lineStyle,'LineWidth',2.5,'DisplayName',condition,'Color',lineColor);
    ll.Color = [ll.Color,0.9];

end

function [] = createPatch(onset,offset)

    % Create a patch to highlight the stimulus onset
    y_limits = ylim();
    stimulus_height = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)];
    stimulus_duration = [onset offset offset onset];
    patch(stimulus_duration, stimulus_height, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility','off')

end
