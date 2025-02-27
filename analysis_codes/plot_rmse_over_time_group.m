% plot all trials 
% plot x gaze position not downsampled
clc;
clear all;
%close all;
addpath('/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Gaze_Package/')

%% set up files

% mydir = '/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Tracking Parameters Experiment/Final_Experiments/3_Stationary_Dynamic_Flies/data_include_forGaze_secondRun_20ppl';
mydir = '/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Tracking Parameters Experiment/Final_Experiments/3_Stationary_Dynamic_Flies/data_include_forGaze_20ppl_paired';
% mydir = '/Users/helenhu/Documents/MATLAB/Stationary_Dynamic_Flies_Codes';

[cursorFiles,mainFiles, eyelinkFiles] = getFiles(mydir);

numSubj = length(cursorFiles);

%% gaze analyses
recFrames = 21;
recFrames_total = 2*recFrames +1;
subj_stationary_rmse = NaN(numSubj,recFrames_total);
subj_dynamic_rmse = NaN(numSubj,recFrames_total);
subj_flies_rmse = NaN(numSubj,recFrames_total);

for subj = 1:numSubj

    %opt = detectImportOptions('*_cursor.csv');% if use opts, offset by 1
    easyeyes = readtable([mydir filesep cursorFiles{subj}]);
    mainOutput = readtable([mydir filesep mainFiles{subj}]);
    eyelink = readtable([mydir filesep eyelinkFiles{subj}]);
    
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
    % timeframe_stimon = 0.15;
    timeframe_starttrack = 0.4; % for calculating eye drift, only include the last 0.4 sec 

    gaze_err_mtx = NaN(length(stim_on),recFrames_total); 

    % special case: KH2 - skip the first trial
    trials_include = 1:length(stim_on);
    if contains(eyelinkFiles{subj}, 'KevinHong061624')
        trials_include = 2:length(stim_on);
    end
    for s = trials_include

        % first, focus on the full tracking period and correct for drift
        % amount for correction
        stim_timestamp = easyeyes(stim_on(s),:).posixTimeSec; 
        track_timestamp = easyeyes(track_on(s),:).posixTimeSec;
        trial_ee_track = easyeyes.posixTimeSec > (stim_timestamp -  timeframe_starttrack) & easyeyes.posixTimeSec < (stim_timestamp);
        currenttrial_ee_track = easyeyes(trial_ee_track,:);
    
        trial_el_track = eyelink.t1+14400 > (stim_timestamp -  timeframe_starttrack) & eyelink.t1+14400 < (stim_timestamp);
        currenttrial_el_track = eyelink(trial_el_track,:);

        framesPerSec = 60;
        eyedrift = calcDrift(currenttrial_ee_track,currenttrial_el_track,PixelPerDeg,framesPerSec);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % then, analyze gaze      
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


        [newGazePx,nantrialBool] = ignoreBlink(gazePx,PixelPerDeg,framesPerSec);
        if nantrialBool
            disp(subj)
            disp(currenttrial_ee_stimon.conditionName(1))
        end

        diff_gaze_cross_Px = [];
        for tt = 1:length(crosshairPx)
            diff_gaze_cross_Px(tt) = norm(crosshairPx(tt,:) - newGazePx(tt,:));
        end

        diff_gaze_cross_Deg = diff_gaze_cross_Px./PixelPerDeg;

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
end

%% save

frame_interval = mean(diff(easyeyes.posixTimeSec(1:10*60)));
xValues = linspace(-frame_interval*recFrames,frame_interval*recFrames,recFrames_total);

save('subj_rmse.mat','subj_stationary_rmse','subj_dynamic_rmse','subj_flies_rmse','xValues');



%% plot


dictionary_lineColors = dictionary({'Stationary','Dynamic','Flies'},{'#7E2F8E','#D95319','#0072BD'});
dictionary_lineStyle = dictionary({'Stationary','Dynamic','Flies'},{'-','--','-.'});

avg_stationary_std = mean(subj_stationary_rmse,1);
avg_dynamic_std = mean(subj_dynamic_rmse,1);
avg_flies_std = mean(subj_flies_rmse,1);

figure;clf;

% set up the graph 

xlabel('Time (s)');
ylabel('Fixation Error RMSE (deg)');
set(gca,'Fontsize',18)
xlim([-0.15 0.3])
ylim([0 2])

createPatch(0,0.15);

hold on;
l1 = plot(xValues,avg_stationary_std,cell2mat(dictionary_lineStyle({'Stationary'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Stationary'})),'DisplayName','Stationary Fixation');
% l1.Color = [l1.Color,0.9]; 
l2 = plot(xValues,avg_dynamic_std,cell2mat(dictionary_lineStyle({'Dynamic'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Dynamic'})),'DisplayName','Dynamic Fixation');
% l2.Color = [l2.Color,0.9]; 
l3 = plot(xValues,avg_flies_std,cell2mat(dictionary_lineStyle({'Flies'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Flies'})),'DisplayName','Crowded Dynamic Fixation');
% l3.Color = [l3.Color,0.9]; 

% for line = 1:numSubj
%     ll = plot(xValues,subj_stationary_std(line,:),cell2mat(dictionary_lineStyle({'Stationary'})),'LineWidth',2.5,'Color',cell2mat(dictionary_lineColors({'Stationary'})));
%     ll.Color = [ll.Color,0.9];   
% end

legend('Location','northwest');
% title('Average (N = 20)')
hold off;


%% functions

function [trialStart,trialEnd] = findCondTrials(plot_cond,block_sequence,subj_name)
    whichBlock = find(strcmp(plot_cond,block_sequence));
    trialStart = (whichBlock-1)*70+1;
    trialEnd = whichBlock*70;
    if contains(subj_name, 'MarcoLai061024') && trialEnd == 210
        trialEnd = 209;
    end
    if contains(subj_name, 'KevinHong061624') 
        trialEnd = trialEnd-1;
        trialStart = trialStart+1;
    end
end


function [rmse] = rmseOverTrials(gaze_mtx)

    % y values: std across all trials over time
    squared_distances = gaze_mtx .^ 2;
    mean_squared_distances = nanmean(squared_distances, 1);
    rmse = sqrt(mean_squared_distances);

end


function [] = plotStd(xValues,gaze_mtx,lineStyle,condition,lineColor)
    % figure;clf;

    % y values: std across all trials over time
    std_over_trials = NaN(1,size(gaze_mtx,2));
    for ts = 1:size(gaze_mtx,2)
        position_trials = gaze_mtx(:,ts);
        rm = isnan(position_trials);
        position_trials(rm) = [];
        std_over_trials(ts) = std(position_trials);
    end


    % plot gaze positions
    ll = plot(xValues,std_over_trials,lineStyle,'LineWidth',2.5,'DisplayName',condition,'Color',lineColor);
    ll.Color = [ll.Color,0.9];

end

function [] = createPatch(onset,offset)

    % Create a patch to highlight the stimulus onset
    y_limits = ylim();
    stimulus_height = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)];
    stimulus_duration = [onset offset offset onset];
    patch(stimulus_duration, stimulus_height, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility','off')

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
