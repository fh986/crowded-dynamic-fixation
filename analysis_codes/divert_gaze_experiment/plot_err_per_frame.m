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
all_conditions = [       {['Grav' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    '' ...
    'Flies_Right8Deg'] }];

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
            disp(size(trial_tracking_error_perFrame_Deg))
            all_tracking_per_frame = [all_tracking_per_frame,trial_tracking_error_perFrame_Deg];

            
        end

    end



    % 
    % % average trials of the same condition & create dictionary
    % cur_conditions = unique(trialStimInfo.conditionName,'stable'); 
    % cond_tracking_error_Deg = NaN(length(cur_conditions),1);
    % cond_trial_tracking_error_Deg = {};
    % 
    % for ii = 1:length(cur_conditions)
    %     idx_cond = find(strcmp(trialStimInfo.conditionName,cur_conditions(ii)));
    %     this_cond = trial_tracking_error_Deg(idx_cond);
    % 
    %     rm = find(this_cond == -1000*px_to_deg);
    %     this_cond(rm) = [];
    % 
    %     this_cond = removeOutliersIQR(this_cond);
    % 
    %     cond_tracking_error_Deg(ii) = mean(this_cond); %%median and mean
    %     cond_trial_tracking_error_Deg{ii} = this_cond;
    % end
    % 
    % cur_dictionary_tracking_Deg = dictionary(cur_conditions,cond_tracking_error_Deg);
    % cur_map_trial_tracking_Deg = containers.Map(cur_conditions,cond_trial_tracking_error_Deg);
    % 
    % all_subj_tracking_Deg{subj} = cur_dictionary_tracking_Deg;
    % all_subj_trial_tracking_Deg{subj} = cur_map_trial_tracking_Deg;

end

%%
large_error_count = sum(all_tracking_per_frame > 0.15);
percent_large_err = large_error_count / length(all_tracking_per_frame)

figure;
hist(all_tracking_per_frame)
xline(0.15)
xlim([0 3])
%%
% %% average across subjects
% 
% dictionary_avg_nTimeouts = dictionary({},[]);
% dictionary_avg_tracking_Deg = dictionary({},[]);
% 
% 
% for ii = 1:length(all_conditions)
% 
%     cur_condName = all_conditions(ii);
% 
%     total_timeouts = 0;
%     total_tracking_error = 0;
%     subjCount_track = 0;
% 
%     for jj = 1:numSubj
% 
%         subj_dic_tracking_error = all_subj_tracking_Deg{jj};
% 
% 
%         if isKey(subj_dic_tracking_error,cur_condName)
%             errorValue = subj_dic_tracking_error(cur_condName);
% 
%             if ~isnan(errorValue)
%                 subjCount_track = subjCount_track + 1;
%                 total_tracking_error = total_tracking_error + errorValue;
%             end
%         end
% 
%         dictionary_avg_nTimeouts = insert(dictionary_avg_nTimeouts,cur_condName,total_timeouts/numSubj);
%         dictionary_avg_tracking_Deg = insert(dictionary_avg_tracking_Deg,cur_condName,total_tracking_error/subjCount_track);
%     end
% end

% %% SEM for each mean
% 
% % concatenate individual data points
% map_allTrials_tracking_Deg = containers.Map('KeyType', 'char', 'ValueType', 'any');
% dictionary_SEM_tracking_Deg = dictionary({},[]);
% for ii = 1:length(all_conditions)
% 
%     cur_condName = all_conditions(ii);
% 
%     cur_trials = [];
% 
%     for jj = 1:numSubj      
%         cur_map_trial_tracking_Deg = all_subj_trial_tracking_Deg{jj};
%         cur_trials = [cur_trials;cur_map_trial_tracking_Deg(cell2mat(cur_condName))];
%     end
%     map_allTrials_tracking_Deg(cell2mat(cur_condName)) = cur_trials;
% 
% 
%     % calculate SEM
%     cur_sem = sem(cur_trials);
%     dictionary_SEM_tracking_Deg(cur_condName) = cur_sem;
% 
% end



% %% plot 
% 
% dot_location_LR = [0,1,2,4,8];
% NF_condNames = {'NoFlies_noDot','NoFlies_Left1Deg','NoFlies_Right2Deg','NoFlies_Left4Deg','NoFlies_Right8Deg'};
% F_3_Grav_condNames = {'GravFlies_noDot','GravFlies_Left1Deg','GravFlies_Right2Deg','GravFlies_Left4Deg','GravFlies_Right8Deg'};
% 
% for ii = 1:length(dot_location_LR)
% 
%     NF_tracking_error_Deg(ii) = dictionary_avg_tracking_Deg(NF_condNames(ii));
%     F_3_Grav_tracking_error_Deg(ii) = dictionary_avg_tracking_Deg(F_3_Grav_condNames(ii));
% 
%     %errorbars
%     NF_SEM(ii) = dictionary_SEM_tracking_Deg(NF_condNames(ii));
%     F_3_Grav_SEM(ii) = dictionary_SEM_tracking_Deg(F_3_Grav_condNames(ii));
% 
% 
% end
% 
% 
% 
% %%
% dotcolors = {'#7E2F8E','#0072BD','#EDB120','#D95319'};
% plotPositions = [0,0,0];
% 
% 
% % tracking errors
% 
% figure;
% xlabel('Fixation distance from center (deg)');
% ylabel('Cursor tracking error (deg)')
% xlim([0,8.1])
% ylim([0 0.3]);
% set(gca,'FontSize',18)
% hold on;
% errorbar(dot_location_LR,NF_tracking_error_Deg,NF_SEM,'o-', ...
%     'LineWidth',2,'Color',cell2mat(dotcolors(4)), ...
%     'DisplayName','No Flies', 'CapSize',0);
% errorbar(dot_location_LR + 0.1,F_3_Grav_tracking_error_Deg,F_3_Grav_SEM,'o-', ...
%     'LineWidth',2,'Color',cell2mat(dotcolors(2)), ...
%     'DisplayName','3 Flies','CapSize',0);
% yline(0.15,'k--','LineWidth',2,'DisplayName','Hotspot Radius = 0.15')
% % legend('Location','North');
% set(gcf,'Position',[100 100 400 400])
% hold off;
% 
% % title('Effect of Fixation Locations on Tracking Error')

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