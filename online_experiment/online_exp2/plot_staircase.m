% set up files
% for each exp:

% clc;
clear all;
%close all;

% addpath('/Users/helenhu/Documents/MATLAB/Gaze_Package')

%% set up files

mydir = [pwd '/online2_data'];

d = dir(sprintf('%s/*.csv',mydir));

files = {d.name};

for f = 1 :length(files)

    mainFile = dir(sprintf('%s/*%s*.csv',mydir,files{f}(1:8)));
    mainFiles{f} = mainFile.name;

end

mainFiles = unique(mainFiles);

numSubj = length(mainFiles);

%% 

eccentricity = 10;

for subj = 1:numSubj

    mainOutput = readtable([mydir filesep mainFiles{subj}]);
    ntrials = 35;
    cond_seq = {'Stationary_Right','Stationary_Left', ...
        'Dynamic_Right','Dynamic_Left','Flies_Right','Flies_Left'};
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

%%
% 
% for subj = 1:numSubj
% 
%     mainOutput = readtable([mydir filesep mainFiles{subj}]);
%     ntrials = 35;
%     cond_seq = {'Stationary Right','Stationary Left','Dynamic Right','Dynamic Left','Flies Right','Flies Left'};
%     lineColors = {'#7E2F8E','#D95319','#0072BD'};
% 
%     figure;clf
%     hold on;
%     ct = 1;
%     for b = 1:3 % blocks 1 to 3
% 
%         for s = 1:2 % 1 = right threshold, 2 = left threshold
% 
%             blockOfInterest = b;
%             staircase = s;
% 
%             t_block_cond = mainOutput(contains(mainOutput.staircaseName,sprintf('%i_%i',blockOfInterest,staircase)),:);
%             h(ct) = plot(10.^t_block_cond.questMeanBeforeThisTrialResponse(1:ntrials),'Color',cell2mat(lineColors(b)),'LineWidth',2);
%             ct = ct + 1;
% 
%         end
%     end
%     h(2).LineStyle = '-.';
%     h(4).LineStyle = '-.';
%     h(6).LineStyle = '-.';
%     xlabel('Trials');
%     ylabel('Quest Mean Before This Trial');
%     xlim([1,35])
%     % ylim([0,25])
%     legend(h,cond_seq);
%     hold off;
%     set(gca,'FontSize',18)
%     title(sprintf('Participant #%d',subj))
% 
% end
% 
