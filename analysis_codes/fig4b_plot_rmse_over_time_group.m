% plot_rmse_over_time_group.m
% Author: Helen Hu
% Last modified: Mar 7th, 2025

% This file calculates and saves fixation error RMSE over time.
% We start by acquiring gaze positions and crosshair positions from outputs
% of the experiment.
% For each trial, gaze positions were corrected with a drift calculated 
% from the track period. I.e., gaze correction (drift) is different from
% trial to trial.
% Then, fixation error was calculated as the distance between gaze and
% crosshair positions.
% For each participant, gaze error RMSE was calculated for each condition.
% The time points for each trial were aligned relative to the start of the
% stimulus position, and a gaze error RMSE was calculated over all trials
% in a condition, for each condition.
% In the end, the gaze error RMSE is saved locally, and the geomean of the
% RMSE across participants is plotted for each condition.
% (see bootstrap_plot_rmse_group.m for bootstrapping and plotting the RMSE
% curves with error bars)


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
recFrames = 21;
recFrames_total = 2*recFrames +1;
subj_stationary_rmse = NaN(numSubj,recFrames_total);
subj_dynamic_rmse = NaN(numSubj,recFrames_total);
subj_flies_rmse = NaN(numSubj,recFrames_total);
subj_n_blink_trials = NaN(numSubj, 1);

for subj = 1:numSubj

    %opt = detectImportOptions('*_cursor.csv');% if use opts, offset by 1
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

    % special cases: some participants had a time-out
    if contains(eyelinkFiles{subj}, 'ML1_M')
        track_on(198) = [];
    end
    if contains(eyelinkFiles{subj}, 'KH2_K')
        track_on(15) = [];
    end
    if(length(stim_on) ~= length(track_on))
        track_on = [1;track_on];
    end

    assert(~any((stim_on - track_on)<=0));

    % pixel to deg
    screenWidthCm = unique(mainOutput.screenWidthByObjectCm(~isnan(mainOutput.screenWidthByObjectCm)));
    screenWidthPx = unique(mainOutput.screenWidthPx(~isnan(mainOutput.screenWidthPx)));
    distance = 40; %cm
    [PixelPerDeg,px_to_deg] = convertPxDeg(screenWidthPx,screenWidthCm,distance);


    timeframe_starttrack = 0.4; % for calculating eye drift, only include the last 0.4 sec 
    % special case
    trials_include = 1:length(stim_on);
    if contains(eyelinkFiles{subj}, 'KH2_K')
        trials_include = 2:length(stim_on);
    end

    % gaze analysis
    gaze_err_mtx = NaN(length(stim_on),recFrames_total); 
    n_blink_trials = 0;

    for s = trials_include

        %%%%%%%%%%%%%%%%%%%%
        % Drift correction %
        %%%%%%%%%%%%%%%%%%%%
        % focus on the full tracking period and correct for drift
        stim_timestamp = easyeyes(stim_on(s),:).posixTimeSec; 
        track_timestamp = easyeyes(track_on(s),:).posixTimeSec;
        trial_ee_track = easyeyes.posixTimeSec > (stim_timestamp -  timeframe_starttrack) & easyeyes.posixTimeSec < (stim_timestamp);
        currenttrial_ee_track = easyeyes(trial_ee_track,:);
    
        trial_el_track = eyelink.t1+14400 > (stim_timestamp -  timeframe_starttrack) & eyelink.t1+14400 < (stim_timestamp);
        currenttrial_el_track = eyelink(trial_el_track,:);

        framesPerSec = 60;
        eyedrift = calcDrift(currenttrial_ee_track,currenttrial_el_track,PixelPerDeg,framesPerSec);
        
        %%%%%%%%%%%%%%%%%%
        % Gaze positions %
        %%%%%%%%%%%%%%%%%%
        % x and y positions of gaze     
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
            disp('Warning: empty trial')
            fprintf('Session %d', subj)
            fprintf('Condition: %s', currenttrial_ee_stimon.conditionName(1))
        end
        if newGazePx ~= gazePx
            n_blink_trials = n_blink_trials + 1;
        end

        % calculate tracking error: the distance between crosshair and gaze
        diff_gaze_cross_Px = [];
        for tt = 1:length(crosshairPx)
            diff_gaze_cross_Px(tt) = norm(crosshairPx(tt,:) - newGazePx(tt,:));
        end
        diff_gaze_cross_Deg = diff_gaze_cross_Px./PixelPerDeg;

        % record
        relative_time = timestamps - stim_timestamp;
        [val,idx_zero_time] = min(abs(relative_time));
        idx_record = (idx_zero_time-recFrames):(idx_zero_time+recFrames);
        assert(all(idx_record > 0));
        gaze_err_mtx(s,:) = diff_gaze_cross_Deg(idx_record);

    end

    %% save rmse sequences
    block_sequence = unique(mainOutput.blockShuffleGroups1,'stable');
    rm = cellfun("isempty",block_sequence);
    block_sequence(rm) = [];
    plot_sequence = {'Stationary','Dynamic','Flies'};

    %stationary
    [trialStart,trialEnd] = findCondTrials(plot_sequence{1},block_sequence,eyelinkFiles{subj});
    subj_stationary_rmse(subj,:) = rmseOverTrials(gaze_err_mtx(trialStart:trialEnd,:));
    
    %dynamic
    [trialStart,trialEnd] = findCondTrials(plot_sequence{2},block_sequence,eyelinkFiles{subj});
    subj_dynamic_rmse(subj,:) = rmseOverTrials(gaze_err_mtx(trialStart:trialEnd,:));

    %flies
    [trialStart,trialEnd] = findCondTrials(plot_sequence{3},block_sequence,eyelinkFiles{subj});
    subj_flies_rmse(subj,:) = rmseOverTrials(gaze_err_mtx(trialStart:trialEnd,:));  

    % number of trials with blinks
    subj_n_blink_trials(subj) = n_blink_trials;
end


%% calculate average percentage of trials with blinks
avg_n_blink_trials = mean(subj_n_blink_trials);
fprintf('Average percentage of trials with blinks: %f%. \n', avg_n_blink_trials./210*100);


%% save gaze positions as a mat file to local folder

frame_interval = mean(diff(easyeyes.posixTimeSec(1:10*60)));
xValues = linspace(-frame_interval*recFrames,frame_interval*recFrames,recFrames_total);
save('subj_rmse.mat','subj_stationary_rmse','subj_dynamic_rmse','subj_flies_rmse','xValues');

%% plot
dictionary_lineColors = dictionary({'Stationary','Dynamic','Flies'},{'#7E2F8E','#D95319','#0072BD'});
dictionary_lineStyle = dictionary({'Stationary','Dynamic','Flies'},{'-','--','-.'});

avg_stationary_std = geomean(subj_stationary_rmse,1);
avg_dynamic_std = geomean(subj_dynamic_rmse,1);
avg_flies_std = geomean(subj_flies_rmse,1);

figure;clf;

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
title('Average RMSE (N = 20)')
hold off;


%% functions

function [trialStart,trialEnd] = findCondTrials(plot_cond,block_sequence,subj_name)
    whichBlock = find(strcmp(plot_cond,block_sequence));
    trialStart = (whichBlock-1)*70+1;
    trialEnd = whichBlock*70;
    if contains(subj_name, 'ML1_M') && trialEnd == 210
        trialEnd = 209;
    end
    if contains(subj_name, 'KH2_K') 
        trialEnd = trialEnd-1;
        trialStart = trialStart+1;
    end
end

function [rmse] = rmseOverTrials(gaze_mtx)
% this function takes in the gaze position for each time point and
% calculate the RMSE over the 70 trials in a condition for that time
% point
% Input: gaze_mtx: a matric containing gaze positions over trials.
%                  (rows are trials, columns are time points)

    squared_distances = gaze_mtx .^ 2;
    mean_squared_distances = nanmean(squared_distances, 1);
    rmse = sqrt(mean_squared_distances);

end

function [] = createPatch(onset,offset)

    % Create a patch to highlight the stimulus onset
    y_limits = ylim();
    stimulus_height = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)];
    stimulus_duration = [onset offset offset onset];
    patch(stimulus_duration, stimulus_height, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility','off')

end
