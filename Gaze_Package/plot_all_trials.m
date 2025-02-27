% plot all trials 
clc;
clear all;
%close all;


%% set up files

mydir = '/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Gaze_Analyses/draft_codes';
    
[cursorFiles,mainFiles, eyelinkFiles] = getFiles(mydir);

numSubj = length(cursorFiles);

%% count cheat

for subj = 1:numSubj

    %opt = detectImportOptions('*_cursor.csv');% if use opts, offset by 1
    easyeyes = readtable([mydir filesep cursorFiles{subj}]);
    mainOutput = readtable([mydir filesep mainFiles{subj}]);
    eyelink = readtable([mydir filesep eyelinkFiles{subj}]);
    
    % count stimuli presentations
    trialRoutine = strcmp(easyeyes.targetBool,'TRUE');
    diff_trialRout = diff(trialRoutine);
    stim_on = find(diff_trialRout == 1)+1; 
    rm = find(strcmp(easyeyes.conditionName(stim_on),'fixate'));
    stim_on(rm) = [];
    % rm = find(strcmp(easyeyes.trialStep(stim_on),'_instructionRoutineEachFrame'));
    % stim_on(rm) = [];
    trialStimInfo = easyeyes(stim_on,:);

    % count tracking
    trackRoutine = strcmp(easyeyes.trialStep,'trialInstructionRoutineEachFrame');
    diff_trackRout = diff(trackRoutine);
    track_on = find(diff_trackRout == 1)+1; 
    rm = find(strcmp(easyeyes.conditionName(track_on),'fixate'));
    track_on(rm) = [];
    trackingInfo = easyeyes(track_on,:);

    assert(length(stim_on) == length(track_on));
    assert(~any((stim_on - track_on)<=0));
   
    % pixel to deg
    screenWidthCm = unique(mainOutput.screenWidthByObjectCm(~isnan(mainOutput.screenWidthByObjectCm)));
    screenWidthPx = unique(mainOutput.screenWidthPx(~isnan(mainOutput.screenWidthPx)));
    distance = 40; %cm
     
    [PixelPerDeg,px_to_deg] = convertPxDeg(screenWidthPx,screenWidthCm,distance);

    % count peek
    
    timeframe_stimon = 0.15;
    timeframe_starttrack = 0.5; % disregard the first 500 ms of tracking 
    % as the eye might be still looking for the crosshair
    disp(length(stim_on))
    for s = 150 : 210 

        % first, focus on the full tracking period and correct for drift
        stim_timestamp = easyeyes(stim_on(s),:).posixTimeSec; 
        track_timestamp = easyeyes(track_on(s),:).posixTimeSec;
        trial_ee_track = easyeyes.posixTimeSec > (track_timestamp + timeframe_starttrack) & easyeyes.posixTimeSec < (stim_timestamp);
        currenttrial_ee_track = easyeyes(trial_ee_track,:);
    
        trial_el_track = eyelink.t1+14400 > (track_timestamp + timeframe_starttrack) & eyelink.t1+14400 < (stim_timestamp);
        currenttrial_el_track = eyelink(trial_el_track,:);

        framesPerSec = 60;
        drift = calcDrift(currenttrial_ee_track,currenttrial_el_track,PixelPerDeg,framesPerSec);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        plotGaze(easyeyes,stim_timestamp,eyelink,PixelPerDeg,drift,trialStimInfo.conditionName(s));
        
    
    end

    
end




