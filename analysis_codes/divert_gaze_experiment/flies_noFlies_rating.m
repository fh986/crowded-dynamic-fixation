% flies_noFlies_rating.m
% Author: Helen Hu
% Last updated: Apr 9th, 2025

clc;
clear all;
close all;

%% set up files

mydir = pwd;
repoDir = fileparts(mydir);
d = dir(sprintf('%s/*.csv',mydir));

files = {d.name};

for f = 1 :length(files)
   
    mainFile = dir(sprintf('%s/*%s*_main.csv',mydir,files{f}(1:3)));
    mainFiles{f} = mainFile.name;
end

mainFiles = unique(mainFiles);
numSubj = length(mainFiles);

%% clean up


% Define the conditions
all_conditions = [       {'Flies_RL_Left1Deg'  },...
    {'Flies_RL_Left4Deg'  },...
    {'Flies_RL_Right2Deg' },...
    {'Flies_RL_Right8Deg' },...
    {'Flies_RL_noDot'     },...
    {'GravFlies_Left1Deg' },...
    {'GravFlies_Left4Deg' },...
    {'GravFlies_Right2Deg'},...
    {'GravFlies_Right8Deg'},...
    {'GravFlies_noDot'    },...
    {'NoFlies_Left1Deg'   },...
    {'NoFlies_Left4Deg'   },...
    {'NoFlies_Right2Deg'  },...
    {'NoFlies_Right8Deg'  },...
    {'NoFlies_noDot'      }];



for subj = 1:numSubj 

    mainOutput = readtable([mydir filesep mainFiles{subj}], 'VariableNamingRule','preserve');

    cur_all_conditions = unique(mainOutput.conditionName,'stable');
    rm = cellfun(@isempty, cur_all_conditions);
    cur_all_conditions(rm) = [];
    rm = find(strcmp(cur_all_conditions,'fixate'));
    cur_all_conditions(rm) = [];
    answers = mainOutput.questionAndAnswerResponse(~isnan(mainOutput.questionAndAnswerResponse));    
    cur_dictionary_qNa = dictionary(cur_all_conditions,answers);

    all_subj_qNa{subj} = cur_dictionary_qNa;
end

%% average across subjects

dictionary_avg_qNa = dictionary({},[]);
dictionary_SEM_qNa = dictionary({},[]);
dictionary_CI_qNa = dictionary({},{});

for ii = 1:length(all_conditions)

    cur_condName = all_conditions(ii);

    cur_score = [];

    for jj = 1:numSubj      
        cur_qNa_dict = all_subj_qNa{jj};
        cur_score = [cur_score;cur_qNa_dict(cur_condName)];
    end

    cur_geomean = geomean(cur_score);
    dictionary_avg_qNa(cur_condName) = cur_geomean;

    cur_sem = sem(cur_score);
    dictionary_SEM_qNa(cur_condName) = cur_sem;

    cur_CI = bootstrap_mean(cur_score);
    dictionary_CI_qNa{cur_condName} = cur_CI;

end



%% plot 


dot_location_LR = [0,1,2,4,8];
NF_condNames = {'NoFlies_noDot','NoFlies_Left1Deg','NoFlies_Right2Deg','NoFlies_Left4Deg','NoFlies_Right8Deg'};
F_5_condNames = {'Flies_RL_noDot','Flies_RL_Left1Deg','Flies_RL_Right2Deg','Flies_RL_Left4Deg','Flies_RL_Right8Deg'};
F_3_Grav_condNames = {'GravFlies_noDot','GravFlies_Left1Deg','GravFlies_Right2Deg','GravFlies_Left4Deg','GravFlies_Right8Deg'};


for ii = 1:length(dot_location_LR)

    % left and right together
    NF_qna(ii) = dictionary_avg_qNa(NF_condNames(ii));
    F_3_Grav_qna(ii) = dictionary_avg_qNa(F_3_Grav_condNames(ii));

    %errorbars
    NF_SEM(ii) = dictionary_SEM_qNa(NF_condNames(ii));
    F_3_Grav_SEM(ii) = dictionary_SEM_qNa(F_3_Grav_condNames(ii));

    % CI
    NF_CI = dictionary_CI_qNa{NF_condNames(ii)};
    NF_lowerBound(ii) = dictionary_avg_qNa(NF_condNames(ii)) - NF_CI(1);
    NF_upperBound(ii) = NF_CI(2) - dictionary_avg_qNa(NF_condNames(ii));

    F_3_Grav_CI = dictionary_CI_qNa{F_3_Grav_condNames(ii)};
    F_3_Grav_lowerBound(ii) = dictionary_avg_qNa(F_3_Grav_condNames(ii)) - F_3_Grav_CI(1);
    F_3_Grav_upperBound(ii) = F_3_Grav_CI(2) - dictionary_avg_qNa(F_3_Grav_condNames(ii));


end


%%

dotcolors = {'#7E2F8E','#0072BD','#EDB120','#D95319'};
plotPositions = [-0.1,0,0.1,0.2];


figure;
xlabel('Location of dot (deg from center)');
ylabel('Perceived difficulty rating')
xlim([0 8.1])
ylim([0 7]);

set(gca,'FontSize',18)
hold on;

errorbar(dot_location_LR,NF_qna,NF_lowerBound, NF_upperBound,'o-','LineWidth',2,'Color',cell2mat(dotcolors(4)),'DisplayName','No Flies','CapSize',0);
errorbar(dot_location_LR + 0.1,F_3_Grav_qna,F_3_Grav_lowerBound, F_3_Grav_upperBound,'o-','LineWidth',2,'Color',cell2mat(dotcolors(2)),'DisplayName','3 Flies, Gravity','CapSize',0);

% legend('Location','north');
% title('Effect of Fixation Location and Flies on Perceived Difficulty')
set(gcf,'Position',[100 100 400 400])

hold off;

%%
function [sem] = sem(data)
    n = length(data);
    std_dev = std(data);
    sem = std_dev / sqrt(n);
end

function [CI_68] = bootstrap_mean(data)

    % Number of bootstrap samples
    num_bootstraps = 10000; 
    
    % Preallocate array for bootstrap means
    bootstrap_means = zeros(num_bootstraps, 1);
    
    % Bootstrapping: Resample data with replacement and compute means
    for ii = 1:num_bootstraps
        resample = data(randi(length(data), 1, length(data))); % Resample with replacement
        bootstrap_means(ii) = geomean(resample); % Compute mean of resampled data
    end
    
    % Compute 68% confidence interval (16th and 84th percentiles)
    CI_68 = prctile(bootstrap_means, [16, 84]);
    
end
