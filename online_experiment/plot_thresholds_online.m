% plot_thresholds_online.m
% Author: Helen Hu

% This script acquires threshold data, one sessions for each
% participant, and plots a stacked histogram of all thresholds.


clc;
clear all;
close all;

plot_staircase_bool = false;
eccentricity = 10;
%% set up files
scriptDir = pwd;
% repoDir = fileparts(scriptDir);
mydir = [scriptDir '/online_data'];
d = dir(sprintf('%s/*.csv',mydir));
files = {d.name};

for f = 1 :length(files)
    mainFile = dir(sprintf('%s/%s*.csv',mydir,files{f}(1:15)));
    mainFiles{f} = mainFile.name;
end

mainFiles = unique(mainFiles);
numSubj = length(mainFiles);


    
%% plot staircase

conditions = {'Stationary', 'Dynamic', 'Flies'};
colors_left = rgb2hex(orderedcolors("gem"));
colors_right = rgb2hex(orderedcolors("glow"));

if plot_staircase_bool

    for subj = 1:numSubj
    
        mainOutput = readtable([mydir filesep mainFiles{subj}],'VariableNamingRule', 'preserve');
    
        figure;
        hold on;
        title(sprintf('Subject %d', subj));
        set (gca,'FontSize',15);
        legend('Location','northeast');
        xlabel('Trials')
        ylabel('Crowding distance(deg)')
        ylim([0 15])
        for cond = 1:length(conditions)
    
            condition = conditions{cond};
            ee_condition = mainOutput(strcmp(mainOutput.blockShuffleGroups1,condition),:);
    
            trials_left = ee_condition(contains(ee_condition.conditionName, 'Left'), :);
            trials_right = ee_condition(contains(ee_condition.conditionName, 'Right'), :);
    
            staircase_left = trials_left.questMeanBeforeThisTrialResponse;
            staircase_right = trials_right.questMeanBeforeThisTrialResponse;
    
            plot(10.^staircase_left,'-','Color',colors_left(cond),'LineWidth',2, ...
                'DisplayName',sprintf('%s, Left', condition))
            plot(10.^staircase_right,'--','Color',colors_right(cond),'LineWidth',2, ...
                'DisplayName',sprintf('%s, Right', condition))
            
        end
    
    end

end


%% extract thresholds
results = table();  
row = 1;  
for subj = 1:numSubj

    mainOutput = readtable([mydir filesep mainFiles{subj}], 'VariableNamingRule', 'preserve');

    for cond = 1:length(conditions)

        condition = conditions{cond};
        conditionName_right = sprintf('%s_Right', condition);
        conditionName_left = sprintf('%s_Left', condition);

        trials_right = mainOutput(strcmp(mainOutput.conditionName, conditionName_right), :);
        trials_left = mainOutput(strcmp(mainOutput.conditionName, conditionName_left), :);

        % get thresholds
        threshold_left = 10.^trials_left.questMeanAtEndOfTrialsLoop(~isnan(trials_left.questMeanAtEndOfTrialsLoop));
        threshold_right = 10.^trials_right.questMeanAtEndOfTrialsLoop(~isnan(trials_right.questMeanAtEndOfTrialsLoop));

        assert(numel(threshold_left) == 1)
        assert(numel(threshold_right) == 1)

        bouma_left = threshold_left ./ eccentricity;
        bouma_right = threshold_right ./ eccentricity;

        % get measurements for filtering
        % 1. questSD
        questSD_right = trials_right.questSDAtEndOfTrialsLoop(~isnan(trials_right.questSDAtEndOfTrialsLoop));
        questSD_left = trials_left.questSDAtEndOfTrialsLoop(~isnan(trials_left.questSDAtEndOfTrialsLoop));
        assert(numel(questSD_right) == 1)
        assert(numel(questSD_left) == 1)
        % 2. trialGivenToQuest
        % Check if 'trialGivenToQuest' exists in the table
        if ismember('trialGivenToQuest', trials_right.Properties.VariableNames)
            trialGivenToQuest_right = trials_right.trialGivenToQuest;
            numTrials_right = sum(strcmp(trialGivenToQuest_right, 'TRUE'));
        else
            numTrials_right = 0;
        end
        
        if ismember('trialGivenToQuest', trials_left.Properties.VariableNames)
            trialGivenToQuest_left = trials_left.trialGivenToQuest;
            numTrials_left = sum(strcmp(trialGivenToQuest_left, 'TRUE'));
        else
            numTrials_left = 0;
        end

        % Add a row to the results table
        results(row,:) = table(subj, {condition}, bouma_left, bouma_right, threshold_left, threshold_right, ...
                               questSD_left, questSD_right, numTrials_left, numTrials_right, ...
                               'VariableNames', {'Subject', 'Condition', 'BoumaLeft', 'BoumaRight', 'ThresholdLeft', 'ThresholdRight', ...
                               'questSD_left', 'questSD_right', 'numTrials_left', 'numTrials_right'});

        row = row + 1;

    end
end

%% filter data based on quest SD and number of trials
subj_large_questSD = unique(results.Subject(results.questSD_left > 0.1 | results.questSD_right > 0.1));
subj_not_enough_trials = unique(results.Subject(results.numTrials_left < 35 | results.numTrials_left < 35)); 
subj_exclude = unique([subj_large_questSD; subj_not_enough_trials]);

results_filtered = results(~ismember(results.Subject, subj_exclude), :);

num_subj_remain = height(results_filtered)/3;
fprintf('Plotting results for %d out of %d subjects... \n\n', num_subj_remain, numSubj)

%% determine minimum threshold from Kurzawski et al., 2023, JOV

repoDir = fileparts(scriptDir);
jov_raw_data = readtable(sprintf('%s/analysis_codes/JoV23Data.csv',repoDir));

% filter
jov_filtered_data = jov_raw_data(strcmp(jov_raw_data.FlankinDirection, 'radial'),:);
jov_filtered_data = jov_filtered_data(strcmp(jov_filtered_data.Task, 'crowding'),:);
jov_filtered_data = jov_filtered_data(~strcmp(jov_filtered_data.Meridian, 'Upper'),:);
jov_filtered_data = jov_filtered_data(~strcmp(jov_filtered_data.Meridian, 'Lower'),:);
jov_filtered_data = jov_filtered_data(strcmp(jov_filtered_data.Font, 'Sloan'),:);
jov_filtered_data = jov_filtered_data((jov_filtered_data.RadialEccentricity == 10),:);

% fprintf('Number of thresholds: %d\n', size(jov_filtered_data,1))

jov_filtered_crowding_distance = jov_filtered_data.CrowdingDistance;
jov_filtered_eccentricity = jov_filtered_data.RadialEccentricity;

jov_filtered_bouma = jov_filtered_crowding_distance./jov_filtered_eccentricity;

minJOVbouma = round(min(jov_filtered_bouma), 3);

fprintf('Minimum Bouma factor from Kurzawski et al., 2023, JoV: %f\n', minJOVbouma)



%% plot histograms for Bouma factors
 
CData = {[0.4940, 0.1840, 0.5560],[0.8500, 0.3250, 0.0980],[0, 0.4470, 0.7410]};
CData2 = {[0.5600, 0.1600, 0.6000],[0.8700, 0.3700, 0.1300],[0, 0.5000, 0.8000]};

results_for_plotting = results_filtered;
table_stationary = results_for_plotting(strcmp(results_for_plotting.Condition, 'Stationary'),:);
table_dynamic = results_for_plotting(strcmp(results_for_plotting.Condition, 'Dynamic'),:);
table_flies = results_for_plotting(strcmp(results_for_plotting.Condition, 'Flies'),:);

ylimit = [0 12];
% bwidth = 0.1;
xlimit = [0.03 1];

minData = 0.05;
maxData = 1;
numBins = 35; 
binEdges = logspace(log10(minData), log10(maxData), numBins);
binWidths = diff(binEdges);

% make sure of that minJOVdata is included as a bin edge
[~, index] = min(abs(binEdges - minJOVbouma));
offset = binEdges(index) - minJOVbouma;
binEdges = binEdges - offset;

figure;
subplot(3,1,1)
plotStackedHist(table_stationary.BoumaRight,table_stationary.BoumaLeft,binEdges,xlimit,ylimit, ...
    'Stationary fixation',cell2mat(CData(1)),minJOVbouma)


subplot(3,1,2)
plotStackedHist(table_dynamic.BoumaRight,table_dynamic.BoumaLeft,binEdges,xlimit,ylimit, ...
    'Dynamic fixation',cell2mat(CData(2)),minJOVbouma)


subplot(3,1,3)
plotStackedHist(table_flies.BoumaRight,table_flies.BoumaLeft,binEdges,xlimit,ylimit, ...
    'Crowded dynamic fixation',cell2mat(CData(3)),minJOVbouma)
xlabel('Bouma factor b')

sgtitle(sprintf('Number of subjects: %d', height(table_flies)))


%%
function [] = plotStackedHist(data1,data2,binEdges,xlimit,ylimit,titletxt,color1,minJOVbouma)
        
    hold on;
    xlim(xlimit)
    ylim(ylimit)
    title(titletxt)
    set(gca,'FontSize',15)

    set(gca, 'XScale', 'log');
    xticks([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, ...
            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]);
    set(gca, 'XMinorTick', 'off');  % since we’re setting everything manually

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

