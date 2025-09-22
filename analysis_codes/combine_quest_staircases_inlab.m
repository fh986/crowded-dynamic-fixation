% combine_quest_staircases_inlab.m
% Author: Helen Hu

% This script acquires threshold data, one sessions for each
% participant, and plots a stacked histogram of all thresholds.


addpath(genpath('/Applications/Psychtoolbox'))

clc;
clear all;
% close all;

plot_staircase_bool = false;
eccentricity = 10;
%% set up files
% Get current directory and move one level up
scriptDir = pwd;
repoDir = fileparts(scriptDir);
% add some functions for gaze analysis
addpath(fullfile(repoDir,'Gaze_Package'));
% data files for gaze
mydir = fullfile(repoDir,'data_include_threshold_18ppl_paired');
   
[cursorFiles,mainFiles, eyelinkFiles] = getFiles(mydir);

numSess = length(cursorFiles);

conditions = {'Stationary', 'Dynamic', 'Flies'};

%% plot staircases
% 
% eccentricity = 10;
% 
% for subj = 20
% 
%     mainOutput = readtable([mydir filesep mainFiles{subj}]);
%     ntrials = 35;
%     cond_seq = {'Stationary_Right1','Stationary_Right2', ...
%         'Dynamic_Right1','Dynamic_Right2','Flies_Right1','Flies_Right2'};
%     lineColors = {'#7E2F8E','#7E2F8E','#D95319','#D95319','#0072BD','#0072BD'};
%     disp(mainFiles{subj})
% 
%     figure;clf
%     hold on;
%     ct = 1;
%     for condition = 1:length(cond_seq)
% 
%         t_block_cond = mainOutput(strcmp(mainOutput.conditionName,cond_seq(condition)),:);
%         h(ct) = plot(10.^t_block_cond.questMeanBeforeThisTrialResponse,'Color',cell2mat(lineColors(ct)),'LineWidth',2);
%         ct = ct + 1;
% 
%         bouma = 10 ^ t_block_cond.questMeanAtEndOfTrialsLoop(end)/eccentricity;
%         disp(bouma)
% 
%     end
%     h(2).LineStyle = '-.';
%     h(4).LineStyle = '-.';
%     h(6).LineStyle = '-.';
%     xlabel('Trials');
%     ylabel('Quest Mean Before This Trial');
%     xlim([1,35])
%     % ylim([0,15])
%     legend(h,cond_seq, 'Location','bestoutside');
%     hold off;
%     set(gca,'FontSize',18)
%     title(sprintf('Participant #%d',subj))
% 
% end


%% extract thresholds
results = table();  
row = 1;  
assert(rem(numSess,2) == 0, 'Fatal: the number of sessions is not a multiple of 2 (test and retest)')
for subj = 1:(numSess/2)

    ind_test = 2 * subj - 1;
    ind_retest = 2 * subj;

    mainOutput_test = readtable([mydir filesep mainFiles{ind_test}], 'VariableNamingRule', 'preserve');
    mainOutput_retest = readtable([mydir filesep mainFiles{ind_retest}], 'VariableNamingRule', 'preserve');

    if ismember('trials.response', mainOutput_test.Properties.VariableNames) && ...
            ismember('trials.response', mainOutput_retest.Properties.VariableNames)
        mainOutput_test.Properties.VariableNames{'trials.response'} = 'trials_response';
        mainOutput_retest.Properties.VariableNames{'trials.response'} = 'trials_response';

    else
        disp('------')
        disp('Column not found.');
        disp(mainFiles{ind_test})
        disp(mainFiles{ind_retest})

    end

    % extract threshold
    for cond = 1:length(conditions)

        condition = conditions{cond};

        conditionName_left = sprintf('%s_Left', condition);
        conditionName_right = sprintf('%s_Right', condition);

        trials_left_test = mainOutput_test(strcmp(mainOutput_test.conditionName, conditionName_left), :);
        trials_left_retest = mainOutput_retest(strcmp(mainOutput_retest.conditionName, conditionName_left), :);
        trials_right_test = mainOutput_test(strcmp(mainOutput_test.conditionName, conditionName_right), :);
        trials_right_retest = mainOutput_retest(strcmp(mainOutput_retest.conditionName, conditionName_right), :);

        trials_left_combined = [trials_left_test; trials_left_retest];
        trials_right_combined = [trials_right_test; trials_right_retest];

        trials_left_sent = trials_left_combined(strcmp(trials_left_combined.trialGivenToQuest, 'TRUE'), :);
        trials_right_sent = trials_right_combined(strcmp(trials_right_combined.trialGivenToQuest, 'TRUE'), :);

        levelProposed_left = trials_left_sent.levelProposedByQUEST;
        levelProposed_right = trials_right_sent.levelProposedByQUEST;

        responses_left = trials_left_sent.trials_response;
        responses_right = trials_right_sent.trials_response;

        assert(size(levelProposed_left,1) == size(responses_left, 1))
        assert(size(levelProposed_right,1) == size(responses_right, 1))
        
        [threshold_predicted_left_log, threshold_predicted_left_sd] = quest_estimate_from_trials(levelProposed_left, responses_left);
        
        threshold_predicted_left = 10 ^ threshold_predicted_left_log;
        [threshold_predicted_right_log, threshold_predicted_right_sd] = quest_estimate_from_trials(levelProposed_right, responses_right);
        threshold_predicted_right = 10 ^ threshold_predicted_right_log;


        % get thresholds predicted by quest
        threshold_left_test = 10.^trials_left_test.questMeanAtEndOfTrialsLoop(~isnan(trials_left_test.questMeanAtEndOfTrialsLoop));
        threshold_left_retest = 10.^trials_left_retest.questMeanAtEndOfTrialsLoop(~isnan(trials_left_retest.questMeanAtEndOfTrialsLoop));
        threshold_right_test = 10.^trials_right_test.questMeanAtEndOfTrialsLoop(~isnan(trials_right_test.questMeanAtEndOfTrialsLoop));
        threshold_right_retest = 10.^trials_right_retest.questMeanAtEndOfTrialsLoop(~isnan(trials_right_retest.questMeanAtEndOfTrialsLoop));
        if numel(threshold_left_test) == 1 && numel(threshold_left_retest) == 1 ...
                && numel(threshold_right_test) == 1 && numel(threshold_right_retest) == 1
    
            bouma_left_test = threshold_left_test ./ eccentricity;
            bouma_left_retest = threshold_left_retest ./ eccentricity;
            bouma_left_predicted = threshold_predicted_left ./ eccentricity;
            bouma_right_test = threshold_right_test ./ eccentricity;
            bouma_right_retest = threshold_right_retest ./ eccentricity;
            bouma_right_predicted = threshold_predicted_right ./ eccentricity;    
            % Add a row to the results table
            results(row,:) = table(subj, {condition}, bouma_left_test, bouma_left_retest, bouma_left_predicted, ...
                                    bouma_right_test, bouma_right_retest, bouma_right_predicted, ...
                                    threshold_predicted_left_sd, threshold_predicted_right_sd, ...
                                   'VariableNames', {'Subject', 'Condition', 'BoumaLeftTest', 'BoumaLeftRetest', 'BoumaLeftPredicted', ...
                                   'BoumaRightTest', 'BoumaRightRetest', 'BoumaRightPredicted', ...
                                   'sdLeft', 'sdRight'});
    
            row = row + 1;

        else

            disp('Warning: no thresholds!!!')
            disp(mainFiles{subj})

        end
    end
end

%% determine minimum threshold from Kurzawski et al., 2023, JOV

jov_raw_data = readtable(sprintf('%s/JoV23Data.csv',scriptDir));

% filter
jov_filtered_data = jov_raw_data(strcmp(jov_raw_data.FlankinDirection, 'radial'),:);
jov_filtered_data = jov_filtered_data(strcmp(jov_filtered_data.Task, 'crowding'),:);
jov_filtered_data = jov_filtered_data(~strcmp(jov_filtered_data.Meridian, 'Upper'),:);
jov_filtered_data = jov_filtered_data(~strcmp(jov_filtered_data.Meridian, 'Lower'),:);
jov_filtered_data = jov_filtered_data(strcmp(jov_filtered_data.Font, 'Sloan'),:);
jov_filtered_data = jov_filtered_data((jov_filtered_data.RadialEccentricity == 10),:);

fprintf('Number of thresholds from Kurzawski, after filtering: %d\n', size(jov_filtered_data,1))

jov_filtered_crowding_distance = jov_filtered_data.CrowdingDistance;
jov_filtered_eccentricity = jov_filtered_data.RadialEccentricity;

jov_filtered_bouma = jov_filtered_crowding_distance./jov_filtered_eccentricity;

minJOVbouma = round(min(jov_filtered_bouma), 3);

fprintf('Minimum Bouma factor from Kurzawski et al., 2023, JoV: %f\n', minJOVbouma)





%% plot histograms for Bouma factors
 
CData = {[0.4940, 0.1840, 0.5560],[0.8500, 0.3250, 0.0980],[0, 0.4470, 0.7410]};
CData2 = {[0.5600, 0.1600, 0.6000],[0.8700, 0.3700, 0.1300],[0, 0.5000, 0.8000]};

results_for_plotting = results;
table_stationary = results_for_plotting(strcmp(results_for_plotting.Condition, 'Stationary'),:);
table_dynamic = results_for_plotting(strcmp(results_for_plotting.Condition, 'Dynamic'),:);
table_flies = results_for_plotting(strcmp(results_for_plotting.Condition, 'Flies'),:);

ylimit = [0 10];
% bwidth = 0.1;
xlimit = [0.01 1.21];

minData = 0.01;
maxData = 1.2;
numBins = 35; 
binEdges = logspace(log10(minData), log10(maxData), numBins);
binWidths = diff(binEdges);

% make sure of that minJOVdata is included as a bin edge
[~, index] = min(abs(binEdges - minJOVbouma));
offset = binEdges(index) - minJOVbouma;
binEdges = binEdges - offset;

figure;
subplot(3,1,1)
plotStackedHist([table_stationary.BoumaRightTest;table_stationary.BoumaRightRetest],...
    [table_stationary.BoumaLeftTest;table_stationary.BoumaLeftRetest],...
    binEdges,xlimit,ylimit, ...
    'Stationary fixation',cell2mat(CData(1)),minJOVbouma)


subplot(3,1,2)
plotStackedHist([table_dynamic.BoumaRightTest;table_dynamic.BoumaRightRetest],...
    [table_dynamic.BoumaLeftTest;table_dynamic.BoumaLeftRetest],...
    binEdges,xlimit,ylimit, ...
    'Dynamic fixation',cell2mat(CData(2)),minJOVbouma)


subplot(3,1,3)
plotStackedHist([table_flies.BoumaRightTest;table_flies.BoumaRightRetest],...
    [table_flies.BoumaLeftTest;table_flies.BoumaLeftRetest],...
    binEdges,xlimit,ylimit, ...
    'Crowded dynamic fixation',cell2mat(CData(3)),minJOVbouma)
xlabel('Bouma factor b')

% sgtitle(sprintf('Number of subjects: %d; Number of thresholds per condition: %d', height(table_flies), height(table_flies)*2))


%% Plot sess 1 and 2 combined


plotCondScatter(table_stationary, 'Stationary', minJOVbouma, ...
                'BoumaLeftTest','BoumaLeftRetest','BoumaLeftPredicted')
plotCondScatter(table_stationary, 'Stationary', minJOVbouma, ...
                'BoumaRightTest','BoumaRightRetest','BoumaRightPredicted')

plotCondScatter(table_dynamic, 'Dynamic', minJOVbouma, ...
                'BoumaLeftTest','BoumaLeftRetest','BoumaLeftPredicted')
plotCondScatter(table_dynamic, 'Dynamic', minJOVbouma, ...
                'BoumaRightTest','BoumaRightRetest','BoumaRightPredicted')

plotCondScatter(table_flies, 'Flies', minJOVbouma, ...
                'BoumaLeftTest','BoumaLeftRetest','BoumaLeftPredicted')
plotCondScatter(table_flies, 'Flies', minJOVbouma, ...
                'BoumaRightTest','BoumaRightRetest','BoumaRightPredicted')


%% Plot hist with 1 and 2 combined


ylimit = [0 8];
% bwidth = 0.1;
xlimit = [0.04 1.21];

minData = 0.01;
maxData = 1.2;
numBins = 40; 
binEdges = logspace(log10(minData), log10(maxData), numBins);
binWidths = diff(binEdges);

% make sure of that minJOVdata is included as a bin edge
[~, index] = min(abs(binEdges - minJOVbouma));
offset = binEdges(index) - minJOVbouma;
binEdges = binEdges - offset;

figure;
subplot(3,1,1)
plotStackedHist(table_stationary.BoumaRightPredicted,table_stationary.BoumaLeftPredicted,binEdges,xlimit,ylimit, ...
    'Stationary fixation',cell2mat(CData(1)),minJOVbouma)


subplot(3,1,2)
plotStackedHist(table_dynamic.BoumaRightPredicted,table_dynamic.BoumaLeftPredicted,binEdges,xlimit,ylimit, ...
    'Dynamic fixation',cell2mat(CData(2)),minJOVbouma)


subplot(3,1,3)
plotStackedHist(table_flies.BoumaRightPredicted,table_flies.BoumaLeftPredicted,binEdges,xlimit,ylimit, ...
    'Crowded dynamic fixation',cell2mat(CData(3)),minJOVbouma)
xlabel('Bouma factor b')


%%
function [t, sd] = quest_estimate_from_trials(levels, responses, ...
    tGuess, tGuessSd, pThreshold, beta, delta, gamma)
% levels:    N×1 vector of stimulus intensities (in the SAME units you use for Quest;
%            typically log10 contrast). If you stored linear contrast c, use log10(c).
% responses: N×1 vector of 0/1 (0 = wrong, 1 = right).


    arguments
        levels (:,1) double
        responses (:,1) double
        tGuess (1,1) double = log10(5) % from EE spreadsheet
        tGuessSd (1,1) double = 2 % easyeyes default
        pThreshold (1,1) double = 0.7 % EE spreadsheet
        beta (1,1) double = 2.3 % EE default
        delta (1,1) double = 0.01 % EE default
        gamma (1,1) double = 1/9 % 1/number of choices 
    end

    % Create Quest state with your prior & psychometric parameters
    q = QuestCreate(tGuess, tGuessSd, pThreshold, beta, delta, gamma);
    q.normalizePdf = 1;

    % drop NaN rows
    valid = ~(isnan(levels) | isnan(responses));
    levels = levels(valid);
    responses = responses(valid) ~= 0; % force logical 0/1

    % Feed each recorded trial to Quest
    for i = 1:numel(levels)
        q = QuestUpdate(q, levels(i), responses(i));
    end

    % Read out the estimate
    t  = QuestMean(q);     % posterior mean (recommended)
    % sd = QuestSd(q);       % posterior sd
    % Alternatives:
    % t  = QuestMode(q);    % MAP estimate
    % t  = QuestQuantile(q);% median
    sd = QuestSd(q);
end





function [] = plotStackedHist(data1,data2,binEdges,xlimit,ylimit,titletxt,color1,minJOVbouma)
        
    hold on;
    xlim(xlimit)
    ylim(ylimit)
    title(titletxt)
    set(gca,'FontSize',15)

    set(gca, 'XScale', 'log');
    xticks([0.02, 0.04, 0.06, 0.08, ...
            0.1, 0.2, 0.4, 0.6, 0.8, 1]);
    set(gca, 'XMinorTick', 'off');  % since we're setting everything manually
    xtickangle(45)

    ylabel('Frequency')
    yticks(min(ylimit):2:max(ylimit))
    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off');
    set(gca, 'TickDir', 'out');

    counts1 = histcounts(data1, binEdges);
    counts2 = histcounts(data2, binEdges);

    for ii = 1:length(binEdges)-1
        % Dataset 1
        xPatch = [binEdges(ii), binEdges(ii+1), binEdges(ii+1), binEdges(ii)]; 
        yPatch1 = [0, 0, counts1(ii), counts1(ii)]; 
        patch(xPatch, yPatch1, color1, 'EdgeColor', 'k', ...
            'EdgeAlpha', 1, 'FaceAlpha', 0.5);
    
        % Dataset 2 (stacked on top of dataset 1)
        yPatch2 = [counts1(ii), counts1(ii), counts1(ii) + counts2(ii), counts1(ii) + counts2(ii)];
        patch(xPatch, yPatch2, color1, 'EdgeColor', 'k', ...
            'EdgeAlpha', 1, 'FaceAlpha', 0.8);

    end
    xline(minJOVbouma,'r--','LineWidth',2);

end

function [] = plotScatter(x_val, x_label, y_val, y_label, titletext)
    figure;
    scatter(x_val,y_val,'filled','SizeData',60, ...
        'MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
    hold on;
    xline(0.096, 'r--', 'LineWidth',1)
    yline(0.096, 'r--', 'LineWidth',1)
    plot([0.01,10], [0.01, 10], 'r--', 'LineWidth',1)
    ylabel(y_label)
    xlabel(x_label)
    title(titletext)
    set(gca, 'YScale', 'log') 
    set(gca, 'XScale', 'log') 
    set(gca, 'FontSize', 18)
    axis padded;
    axis equal;
    axis square;


end

function [] = setScatter(minJOVbouma)
    xlimit = [0.01, 100];
    ylimit = [0.01, 100];
    xlim(xlimit); ylim(ylimit);
    xline(minJOVbouma, 'r--')
    yline(minJOVbouma, 'r--')
    plot(xlimit, ylimit, 'r--')
    ylabel('Bouma factor, predicted by Quest (70 trials)')
    set(gca, 'FontSize', 15)
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    axis square 

end


function [] = plotCondScatter(table_for_plotting, title, minJOVbouma, ...
                                testVarName, retestVarName, predVarName)
    BoumaTest = testVarName;
    BoumaRetest = retestVarName;
    BoumaPredicted = predVarName;


    figure('Position',[100 100 1000 500])
    subplot(1,2,1);
    scatter(table_for_plotting{:,BoumaTest}, table_for_plotting{:,BoumaPredicted}, 'filled', ...
        'SizeData',60, 'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7)
    hold on;
    xlabel('Bouma factor, test (35 trials)')
    setScatter(minJOVbouma);

    hold off;
    
    subplot(1,2,2);
    scatter(table_for_plotting{:,BoumaRetest}, table_for_plotting{:,BoumaPredicted}, 'filled', ...
        'SizeData',60, 'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7)
    hold on;
    xlabel('Bouma factor, retest (35 trials)')
    setScatter(minJOVbouma);

    hold off;

    sgtitle(title, 'FontSize', 18)
end