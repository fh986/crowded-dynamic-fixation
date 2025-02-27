% plot all trials 
% plot x gaze position not downsampled
clc;
clear all;
%close all;
addpath('/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Gaze_Package/')

%% set up files

% mydir = '/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Tracking Parameters Experiment/Final_Experiments/3_Stationary_Dynamic_Flies/aaa_Stationary_Dynamic_Flies_Codes';
mydir = '/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Tracking Parameters Experiment/Final_Experiments/3_Stationary_Dynamic_Flies/data_include_forGaze_20ppl_paired';
 
[cursorFiles,mainFiles, eyelinkFiles] = getFiles(mydir);

numSubj = length(cursorFiles);

%% gaze analyses

for subj = 39

    %opt = detectImportOptions('*_cursor.csv');% if use opts, offset by 1
    easyeyes = readtable([mydir filesep cursorFiles{subj}]);
    mainOutput = readtable([mydir filesep mainFiles{subj}]);
    eyelink = readtable([mydir filesep eyelinkFiles{subj}]);
    
    % count stimuli presentations
    trialRoutine = strcmp(easyeyes.targetBool,'TRUE');
    diff_trialRout = diff(trialRoutine);
    stim_on = find(diff_trialRout == 1)+1;
    % rm = find(strcmp(easyeyes.conditionName(stim_on),'fixate'));
    % stim_on(rm) = [];
    trialStimInfo = easyeyes(stim_on,:);

    % count tracking
    trackRoutine = strcmp(easyeyes.crosshairBool,'TRUE');
    diff_trackRout = diff(trackRoutine);
    track_on = find(diff_trackRout == 1)+1;
    % 
    % rm = find(strcmp(easyeyes.conditionName(track_on),'fixate'));
    % track_on(rm) = [];
    trackingInfo = easyeyes(track_on,:);

    if contains(eyelinkFiles{subj}, 'MarcoLai061024')
        track_on(198) = [];
    end

    if contains(eyelinkFiles{subj}, 'KevinHong061624')
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




    % count peek
    timeframe_stimon = 0.15;
    timeframe_starttrack = 0.4; 
    half_frames_record = 27;      
    gaze_x_mtx = NaN(length(stim_on),half_frames_record*2+1); % 60 frames per sec
    gaze_y_mtx = NaN(length(stim_on),half_frames_record*2+1);

    trials_include = 1:length(stim_on);
    if contains(eyelinkFiles{subj}, 'KevinHong061624')
        trials_include = 2:length(stim_on);
    end
    
    % trials_include = 148;

    for s = trials_include

        % first, focus on the full tracking period and correct for drift
        % amount for correction
        stim_timestamp = easyeyes(stim_on(s),:).posixTimeSec; 
        track_timestamp = easyeyes(track_on(s),:).posixTimeSec;
        trial_ee_track = easyeyes.posixTimeSec > (stim_timestamp - timeframe_starttrack) & easyeyes.posixTimeSec < (stim_timestamp);
        currenttrial_ee_track = easyeyes(trial_ee_track,:);
    
        trial_el_track = eyelink.t1+14400 > (stim_timestamp -  timeframe_starttrack) & eyelink.t1+14400 < (stim_timestamp);
        currenttrial_el_track = eyelink(trial_el_track,:);

        framesPerSec = 60;
        eyedrift = calcDrift(currenttrial_ee_track,currenttrial_el_track,PixelPerDeg,framesPerSec);
        %disp(eyedrift)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % then, analyze gaze      
        trial_ee_stim = find(easyeyes.posixTimeSec == stim_timestamp);
        currenttrial_ee_stimon = easyeyes(trial_ee_stim,:);
        trial_ee = easyeyes.posixTimeSec > (stim_timestamp - 0.5) & easyeyes.posixTimeSec < (stim_timestamp + 0.5);
        currenttrial_ee = easyeyes(trial_ee,:);

        nearpointPx = str2num(cell2mat(currenttrial_ee_stimon.nearpointXYPx));
        % last_crosshairPx = str2num(cell2mat(currenttrial_ee_stimon.crosshairPositionXYPx)) - nearpointPx;

        trial_eyelink = eyelink.t1+14400 > (stim_timestamp - 0.5) & eyelink.t1+14400 < (stim_timestamp + 0.5);
        currenttrial_el = eyelink(trial_eyelink,:);
        timestamps_el = currenttrial_el.t1;

        [nearpointPx,crosshairPx,cursorPx] = getXYPositions_ee(currenttrial_ee);

        gazePx = [];
        for frame = 1:size(currenttrial_ee,1)
            tp_ee = currenttrial_ee.posixTimeSec(frame,:)-14400;
            [~,idx_tp_el] = min(abs(timestamps_el - tp_ee));
            gazePx(frame,1) = currenttrial_el.gazeXYPix_1(idx_tp_el) - eyedrift(1) - crosshairPx(frame,1);
            gazePx(frame,2) = -currenttrial_el.gazeXYPix_2(idx_tp_el) - eyedrift(2)  - crosshairPx(frame,2);
        end

        % gazePx = [];
        % gazePx(:,1) = currenttrial_el.gazeXYPix_1  - last_crosshairPx(1)- eyedrift(1);
        % gazePx(:,2) = -currenttrial_el.gazeXYPix_2  - last_crosshairPx(2)- eyedrift(2);
        [newGazePx,nantrialBool] = ignoreBlink(gazePx,PixelPerDeg,framesPerSec);
        if nantrialBool
            disp('nantrial')
        end
        

        gazeXDeg = newGazePx(:,1) ./ PixelPerDeg;
        gazeYDeg = newGazePx(:,2) ./ PixelPerDeg;

        relative_time = currenttrial_ee.posixTimeSec - stim_timestamp;
        [val,idx_zero_time] = min(abs(relative_time));

        idx_record = (idx_zero_time-half_frames_record):(idx_zero_time+half_frames_record);

        gaze_x_mtx(s,:) = gazeXDeg(idx_record);
        gaze_y_mtx(s,:) = gazeYDeg(idx_record);

        % left and right
        condNames = mainOutput.conditionName;
        rm = cellfun('isempty',condNames);
        condNames(rm) = [];
        rm = [71 72 143 144 215 216];
        condNames(rm) = [];
        rightBool = contains(condNames,'Right');

    end

    % plot
    block_sequence = unique(mainOutput.blockShuffleGroups1,'stable');
    rm = cellfun("isempty",block_sequence);
    block_sequence(rm) = [];
   
    frame_interval = mean(diff(eyelink.t1(1:10*60)));
    xValues = linspace(-frame_interval*half_frames_record,frame_interval*half_frames_record,half_frames_record*2+1);

    % def blocks
    firstBlockEnd = 70;
    secondBlockStart = 71;
    secondBlockEnd = 140;
    thirdBlockStart = 141;
    thirdBlockEnd = 210;
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
    for trial = 1:size(gaze_x_mtx)
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
    for trial = 1:size(gaze_y_mtx)
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

function [newGazePx,nantrialBool] = ignoreBlink(gazePx,pxPerDeg,framesPerSec)
% detects blinks based on the criterion: saccade velocity > 2000 deg/sec
% replaces the 100 msecs before and after the blink with the gaze position
% in the frame after the blink
% if the blink happens at the end of the sequence, replace with the frame
% before the blink

   
    gazeDeg = gazePx/pxPerDeg;
    
    velocity = eyeVelocity(gazeDeg,framesPerSec);
    delFrames = round(0.1*framesPerSec);
    
    replaceIdx = [];
    for ii = 1:length(velocity)
    
       if any([velocity(ii) > 1000;  abs(gazeDeg(ii,1)) >20;abs(gazeDeg(ii,2)) >20])
            replaceIdx = [replaceIdx,ii-delFrames:ii+delFrames];
       end

    end
    
    replaceIdx(replaceIdx<=0) = [];
    replaceIdx(replaceIdx>size(gazeDeg,1)) = [];
    replaceIdx = unique(replaceIdx);

    % are there multiple blinks?
    difference = diff(replaceIdx);
    multipleBool = any(difference~=1);
    blinks = {};
    if ~multipleBool
        blinks{1} = replaceIdx;
    else
        breaks = find(difference~=1);
        breaks = [0,breaks,length(replaceIdx)];
        for blinkNumber = 1:length(breaks)-1
            blinks{blinkNumber} = replaceIdx(breaks(blinkNumber)+1:breaks(blinkNumber+1));
        end
    end

    newGazePx = gazePx;
    nantrialBool = 0;
    if ~isempty(replaceIdx)
        if length(replaceIdx) == length(gazePx)
            newGazePx = NaN(size(gazePx));
            nantrialBool = 1;
        else
            for blinkNumber = 1:length(blinks)
                thisBlink = blinks{blinkNumber};
                replacePosition = [];
                
                if ismember(size(gazeDeg,1),thisBlink)
                    replacePosition = gazePx(thisBlink(1)-1,:);
                else
                    replacePosition = gazePx(thisBlink(end)+1,:);
                end 
    
                for frame = thisBlink
                    newGazePx(frame,:) = replacePosition;
                end
            end
        end
    end


end
