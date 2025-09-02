% plot_thresholds_online.m
% Author: Helen Hu

% This script acquires threshold data, one sessions for each
% participant, and plots a stacked histogram of all thresholds.


clc;
clear all;
% close all;

plot_staircase_bool = false;
eccentricity = 10;
%% set up files
scriptDir = pwd;
% repoDir = fileparts(scriptDir);
mydir = [scriptDir '/online3_data'];
d = dir(sprintf('%s/*.csv',mydir));
files = {d.name};

for f = 1 :length(files)
    mainFile = dir(sprintf('%s/%s*.csv',mydir,files{f}(1:15)));
    mainFiles{f} = mainFile.name;
end

mainFiles = unique(mainFiles);
numSubj = length(mainFiles);

conditions = {'Stationary', 'Dynamic', 'Flies'};

    
%% plot staircase

% colors_left = rgb2hex(orderedcolors("gem"));
% colors_right = rgb2hex(orderedcolors("glow"));

if plot_staircase_bool
    
    for subj = 1:numSubj
    
        mainOutput = readtable([mydir filesep mainFiles{subj}]);
        ntrials = 35;
        cond_seq = {'Stationary_Right1','Stationary_Right2',...
            'Dynamic_Right1','Dynamic_Right2',...
            'Flies_Right1','Flies_Right2'};
        lineColors = {'#7E2F8E','#7E2F8E','#D95319','#D95319','#0072BD','#0072BD'};
        disp(mainFiles{subj})
    
        figure;clf
        hold on;
        ct = 1;
        for condition = 1:length(cond_seq)
    
            t_block_cond = mainOutput(strcmp(mainOutput.conditionName,cond_seq(condition)),:);
            h(ct) = plot(10.^t_block_cond.questMeanBeforeThisTrialResponse,'Color',cell2mat(lineColors(ct)),'LineWidth',2);
            ct = ct + 1;
            
            bouma = 10 ^ t_block_cond.questMeanAtEndOfTrialsLoop(end)/eccentricity;
            disp(bouma)
        
        end
        h(2).LineStyle = '-.';
        h(4).LineStyle = '-.';
        h(6).LineStyle = '-.';
        xlabel('Trials');
        ylabel('Quest Mean Before This Trial');
        xlim([1,35])
        ylim([0,15])
        legend(h,cond_seq);
        hold off;
        set(gca,'FontSize',18)
        title(sprintf('Participant #%d',subj))
    
    end
end


%% extract thresholds
results = table();  
row = 1;  
for subj = 1:numSubj

    mainOutput = readtable([mydir filesep mainFiles{subj}], 'VariableNamingRule', 'preserve');
    
    screenWidthCm = unique(mainOutput.screenWidthByObjectCm(~isnan(mainOutput.screenWidthByObjectCm)));
    assert(length(screenWidthCm) == 1, 'Incorrect screen width input')


    % extract threshold
    for cond = 1:length(conditions)

        condition = conditions{cond};
        conditionName_test = sprintf('%s_Right1', condition);
        conditionName_retest = sprintf('%s_Right2', condition);

        trials_test = mainOutput(strcmp(mainOutput.conditionName, conditionName_test), :);
        trials_retest = mainOutput(strcmp(mainOutput.conditionName, conditionName_retest), :);

        % get thresholds
        threshold_test = 10.^trials_test.questMeanAtEndOfTrialsLoop(~isnan(trials_test.questMeanAtEndOfTrialsLoop));
        threshold_retest = 10.^trials_retest.questMeanAtEndOfTrialsLoop(~isnan(trials_retest.questMeanAtEndOfTrialsLoop));

        if numel(threshold_test) == 1 && numel(threshold_retest) == 1
    
            bouma_test = threshold_test ./ eccentricity;
            bouma_retest = threshold_retest ./ eccentricity;
    
            % get measurements for filtering
            % 1. questSD
            questSD_test = trials_test.questSDAtEndOfTrialsLoop(~isnan(trials_test.questSDAtEndOfTrialsLoop));
            questSD_retest = trials_retest.questSDAtEndOfTrialsLoop(~isnan(trials_retest.questSDAtEndOfTrialsLoop));
            assert(numel(questSD_test) == 1)
            assert(numel(questSD_retest) == 1)
            % 2. trialGivenToQuest
            % Check if 'trialGivenToQuest' exists in the table
            if ismember('trialGivenToQuest', trials_test.Properties.VariableNames)
                trialGivenToQuest_test = trials_test.trialGivenToQuest;
                numTrials_test = sum(strcmp(trialGivenToQuest_test, 'TRUE'));
            else
                numTrials_test = 0;
            end
            
            if ismember('trialGivenToQuest', trials_retest.Properties.VariableNames)
                trialGivenToQuest_retest = trials_retest.trialGivenToQuest;
                numTrials_retest = sum(strcmp(trialGivenToQuest_retest, 'TRUE'));
            else
                numTrials_retest = 0;
            end
            % 3. distanceCalibratedByBlindspot
            all_trials_block = [trials_test; trials_retest];
            distanceByBlindspotCm = unique(all_trials_block.viewingDistanceByBlindSpotCm(~isnan(all_trials_block.viewingDistanceByBlindSpotCm)),'stable');
            assert(numel(distanceByBlindspotCm) == 1)

            % Add a row to the results table
            results(row,:) = table(subj, {condition}, bouma_test, bouma_retest, threshold_test, threshold_retest, ...
                                   questSD_test, questSD_retest, numTrials_test, numTrials_retest, ...
                                   screenWidthCm, distanceByBlindspotCm, ...
                                   'VariableNames', {'Subject', 'Condition', 'BoumaTest', 'BoumaRetest', 'ThresholdTest', 'ThresholdRetest', ...
                                   'questSD_Test', 'questSD_Retest', 'numTrials_test', 'numTrials_retest', ...
                                   'screenWidthCm', 'distanceByBlindspotCm'});
    
            row = row + 1;

        else

            disp('Warning: no thresholds.')
            disp(mainFiles{subj})

        end
    end
end

%% determine minimum threshold from Kurzawski et al., 2023, JOV

repoDir = fileparts(fileparts(scriptDir));
jov_raw_data = readtable(sprintf('%s/analysis_codes/JoV23Data.csv',repoDir));

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

%% filter data based on quest SD and number of trials
subj_exclude = unique(results.Subject(results.questSD_Test > 0.1 | results.questSD_Retest > 0.1));

fprintf('Number of subjects with large questSD: %d\n', length(subj_exclude))

results_filtered = results(~ismember(results.Subject, subj_exclude), :);
results_out = results(ismember(results.Subject, subj_exclude), :);

num_subj_remain = height(results_filtered)/3;
fprintf('Plotting results for %d out of %d subjects... \n\n', num_subj_remain, numSubj)


%% plot histograms for Bouma factors
 
CData = {[0.4940, 0.1840, 0.5560],[0.8500, 0.3250, 0.0980],[0, 0.4470, 0.7410]};
CData2 = {[0.5600, 0.1600, 0.6000],[0.8700, 0.3700, 0.1300],[0, 0.5000, 0.8000]};

results_for_plotting = results_filtered;
table_stationary = results_for_plotting(strcmp(results_for_plotting.Condition, 'Stationary'),:);
table_dynamic = results_for_plotting(strcmp(results_for_plotting.Condition, 'Dynamic'),:);
table_flies = results_for_plotting(strcmp(results_for_plotting.Condition, 'Flies'),:);

ylimit = [0 8];
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
plotStackedHist(table_stationary.BoumaTest,table_stationary.BoumaRetest,binEdges,xlimit,ylimit, ...
    'Stationary fixation',cell2mat(CData(1)),minJOVbouma)


subplot(3,1,2)
plotStackedHist(table_dynamic.BoumaTest,table_dynamic.BoumaRetest,binEdges,xlimit,ylimit, ...
    'Dynamic fixation',cell2mat(CData(2)),minJOVbouma)


subplot(3,1,3)
plotStackedHist(table_flies.BoumaTest,table_flies.BoumaRetest,binEdges,xlimit,ylimit, ...
    'Crowded dynamic fixation',cell2mat(CData(3)),minJOVbouma)
xlabel('Bouma factor b')

% sgtitle(sprintf('Number of subjects: %d; Number of thresholds per condition: %d', height(table_flies), height(table_flies)*2))



 
%%
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