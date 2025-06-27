% flies_noFlies_dot.m
% Author: Helen Hu
% Last updated: Apr 9th, 2025

clc;
clear all;
%close all;

%% set up files

mydir = pwd;
repoDir = fileparts(mydir);

d = dir(sprintf('%s/*.csv',mydir));

files = {d.name};

for f = 1 :length(files)
    
    cursor = dir(sprintf('%s/*%s*_cursor.csv',mydir,files{f}(1:3)));
    cursorFiles{f} = cursor.name;

    mainFile = dir(sprintf('%s/*%s*_main.csv',mydir,files{f}(1:3)));
    mainFiles{f} = mainFile.name;
    
    matlabFile = dir(sprintf('%s/*%s*_matlabOutput.csv',mydir,files{f}(1:3)));
    matlabFiles{f} = matlabFile.name;
end

cursorFiles = unique(cursorFiles);
mainFiles = unique(mainFiles);
matlabFiles = unique(matlabFiles);

assert(length(cursorFiles) == length(mainFiles));
assert(length(cursorFiles) == length(matlabFiles));

numSubj = length(cursorFiles);

%% clean up data

% Define the conditions
all_conditions = [{'GravFlies_Left1Deg'}];

    % {'GravFlies_Left1Deg' }
    % {'GravFlies_Left4Deg' },...
    % {'GravFlies_Right2Deg'},...    
    % {'NoFlies_Left4Deg'   },...
    % {'NoFlies_Right2Deg'  },...
    % ,...
    % {'GravFlies_Right8Deg'},...
    % {'GravFlies_noDot'    },...
    % {'NoFlies_Left1Deg'   },...
    % {'NoFlies_Right8Deg'  },...
    % {'NoFlies_noDot'      }
%%

all_tracking_per_frame = [];
colors_h_codes = rgb2hex(orderedcolors('gem'));


for subj = 1:numSubj 

    easyeyes = readtable([mydir filesep cursorFiles{subj}], 'VariableNamingRule','preserve');
    mainOutput = readtable([mydir filesep mainFiles{subj}], 'VariableNamingRule','preserve');
    eyelink = readtable([mydir filesep matlabFiles{subj}], 'VariableNamingRule','preserve');

    % pixel to deg
    screenWidthCm = unique(mainOutput.screenWidthByObjectCm(~isnan(mainOutput.screenWidthByObjectCm)));
    screenWidthPx = unique(mainOutput.screenWidthPx(~isnan(mainOutput.screenWidthPx)));
    distance = 40; %cm

    PixelPerDeg = screenWidthPx / screenWidthCm * distance * pi / 180;
    px_to_deg = 1/PixelPerDeg;


    % count stimuli presentations
    begin_idx = strcmp(easyeyes.trialStep,'trialInstructionRoutineBegin');
    easyeyes(find(begin_idx),:) = []; % fix software bug
    
    trialRoutine = strcmp(easyeyes.targetBool,'TRUE');
    diff_trialRout = diff(trialRoutine);
    stim_on = find(diff_trialRout == 1)+2;

    rm = find(strcmp(easyeyes.conditionName(stim_on),'fixate'));
    stim_on(rm) = [];
    rm = find(strcmp(easyeyes.conditionName(stim_on),'Flies_RL_Left8Deg'));
    stim_on(rm) = [];
    trialStimInfo = easyeyes(stim_on,:);
    

    % mark beginnings of trials
    trialInstructionBegin = strcmp(easyeyes.trialStep,'trialInstructionRoutineEachFrame');
    diff_trialInst = diff(trialInstructionBegin);
    track_begin = find(diff_trialInst == 1)+2;
    rm = find(strcmp(easyeyes.conditionName(track_begin),'fixate'));
    track_begin(rm) = [];
    rm = find(strcmp(easyeyes.conditionName(track_begin),'Flies_RL_Left8Deg'));
    track_begin(rm) = [];
    track_begin = [1;track_begin];
    trialBeginInfo = easyeyes(track_begin,:);
    if subj == 4
        track_begin(134) = [];
    end

    assert(length(stim_on) == length(track_begin));

    timeframe = 1.5; 
    trial_tracking_error_perFrame_Px = [];

    for block = 1:length(all_conditions)
        this_block = cell2mat(all_conditions(block));
        this_trials_idx = find(strcmp(trialStimInfo.conditionName,this_block));
        trials_per_cond = length(this_trials_idx);
        figure('Position', [100, 100, 800, 1600]);
        for trial = this_trials_idx'
    
            track_timestamp = easyeyes(track_begin(trial),:).posixTimeSec; 
            stim_timestamp = easyeyes(stim_on(trial),:).posixTimeSec; 

            trial_easyeyes = easyeyes.posixTimeSec > (track_timestamp + timeframe) & easyeyes.posixTimeSec < stim_timestamp; 
            currenttrial_ee = easyeyes(trial_easyeyes,:);

            if size(currenttrial_ee,1) > 0
                [~,trial_tracking_error_perFrame_Px] = calcTrackErrorForTrial(currenttrial_ee);
            else
                disp('Warning: empty trial')
            end

            % convert to deg
            trial_tracking_error_perFrame_Deg = trial_tracking_error_perFrame_Px .* px_to_deg;
            
            % plot
            subplot(trials_per_cond, 1, find(this_trials_idx == trial));
            xValues = currenttrial_ee.posixTimeSec - currenttrial_ee.posixTimeSec(1) + timeframe;
            err_line = plot(xValues, trial_tracking_error_perFrame_Deg);
            err_line.Color = colors_h_codes(1);
            err_line.LineWidth = 2;
            err_line.DisplayName = 'Cursor error (deg)';
            yline(0.15,'--','Color','r','LineWidth',2,'Alpha',0.8,'DisplayName','Hotspot radius')

            legend('Location','northeastoutside')
            xLimits = xlim;
            xlim([timeframe, xLimits(2)]);
            xlabel('Time (sec)')
            ylim([0 0.4])
            ylabel(sprintf('Cursor error\n(RMSE deg)'))
            set(gca,'FontSize',12)
            
            % concatenate with other frames in the condition
            all_tracking_per_frame = [all_tracking_per_frame,trial_tracking_error_perFrame_Deg];

            
        end
        sgtitle(sprintf('Subject %d, Condition: %s', subj, ...
            cell2mat(all_conditions(block))),'FontSize',15)

    end

end



%% percentage of frames with large cursor error
large_error_count = sum(all_tracking_per_frame > 0.15);
percent_large_err = large_error_count / length(all_tracking_per_frame)

xlimit = [0 0.4];
% ylimit = [0 2000];
cond_name = cell2mat(all_conditions(1));

figure;
h = histogram(all_tracking_per_frame);
h.FaceColor = colors_h_codes(6);
% h.EdgeColor = [0 0 0];
h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 

xline(0.15,'--','Color','r','LineWidth',2,'Alpha',0.8,'DisplayName','Hotspot radius')
xlim(xlimit)
% ylim(ylimit)
% xscale('log')
set(gca,'FontSize',15)
title(sprintf('Cursor error at each frame, conition: %s', cond_name))
xlabel('Cursor error (RMSE deg)')
ylabel('Frequency')
legend('location','northeast')
% yticks(min(ylimit):2:max(ylimit))
%%
function [sem] = sem(data)
    n = length(data);
    std_dev = std(data);
    sem = std_dev / sqrt(n);
end


function [trialError, tracking_errors_Px] = calcTrackErrorForTrial(currenttrial_ee)
    
        nearpointPx = str2num(cell2mat(currenttrial_ee.nearpointXYPx(1,:)));

        crosshairPx = [];
        cursorPx = [];

        for ii = 1:length(currenttrial_ee.crosshairPositionXYPx)

            aa = currenttrial_ee.crosshairPositionXYPx(ii);
            aa = cell2mat(aa);
            comma_idx = find(aa == ',', 1);
            if comma_idx == 2
                aa(comma_idx) = '';
            end
            crosshairPx(ii,:) = str2num(aa) - nearpointPx;

            aa = currenttrial_ee.cursorPositionXYPx(ii);
            aa = cell2mat(aa);
            comma_idx = find(aa == ',', 1);
            if comma_idx == 2
                aa(comma_idx) = '';
            end
            cursorPx(ii,:) = str2num(aa) - nearpointPx;

        end
           
        tracking_errors_Px = [];

        for x = 1 : length(crosshairPx)
            tracking_errors_Px(x) = norm(crosshairPx(x,:)-cursorPx(x,:));
        end
        trialError = median(tracking_errors_Px);
end