% plot all trials 
% plot x gaze position not downsampled
clc;
clear all;
%close all;
addpath('/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Gaze_Package/')

%% set up files

mydir = '/Users/fh986/Documents/MATLAB/Tracking_Analyses_Codes/Tracking Parameters Experiment/Final_Experiments/3_Stationary_Dynamic_Flies/aaa_Stationary_Dynamic_Flies_Codes/tracking_demo_30s';
    

d = dir(sprintf('%s/*.csv',mydir));
    
files = {d.name};

for f = 1 :length(files)
    
    cursor = dir(sprintf('%s/*%s*_cursor.csv',mydir,files{f}(1:3)));
    cursorFiles{f} = cursor.name;

    mainFile = dir(sprintf('%s/*%s*_main.csv',mydir,files{f}(1:3)));
    mainFiles{f} = mainFile.name;

end

cursorFiles = unique(cursorFiles);
mainFiles = unique(mainFiles);

assert(length(cursorFiles) == length(mainFiles));


numSubj = length(cursorFiles);

%% cursor analyses
recFrames = 1300;

for subj = 1%:numSubj

    %opt = detectImportOptions('*_cursor.csv');% if use opts, offset by 1
    easyeyes = readtable([mydir filesep cursorFiles{subj}]);
    mainOutput = readtable([mydir filesep mainFiles{subj}]);
   

    % count tracking
    trackRoutine = strcmp(easyeyes.trialStep,'trialInstructionRoutineEachFrame');
    diff_trackRout = diff(trackRoutine);
    track_on = find(diff_trackRout == 1)+1;
    trackingInfo = easyeyes(track_on,:);

    % count the end of tracking
    trialRoutine = strcmp(easyeyes.trialStep,'trialsLoopEnd');
    diff_trialRout = diff(trialRoutine);
    stim_on = find(diff_trialRout == 1)+1;
    trialStimInfo = easyeyes(stim_on,:);

    if(length(stim_on) ~= length(track_on))
        track_on = [1;track_on];
    end

    assert(~any((stim_on - track_on)<=0));
   
    % pixel to deg
    screenWidthCm = unique(mainOutput.screenWidthByObjectCm(~isnan(mainOutput.screenWidthByObjectCm)));
    [nearpointPx,crosshairPx,cursorPx] = getXYPositions_ee(easyeyes(1,:));

    screenWidthPx = nearpointPx(1) * 2;
    distance = 40; %cm
    [PixelPerDeg,px_to_deg] = convertPxDeg(screenWidthPx,screenWidthCm,distance);




    % count peek
    timeframe_stimon = 0.15;
    timeframe_starttrack = 4; 

    cursor_x_mtx = NaN(length(stim_on),recFrames); 
    cursor_y_mtx = NaN(length(stim_on),recFrames);
    crosshair_x_mtx = NaN(length(stim_on),recFrames); 
    crosshair_y_mtx = NaN(length(stim_on),recFrames);
    tracking_error = NaN(length(stim_on),recFrames);


    % trials_include = 1:length(stim_on);
    trials_include = 3:4;

   
    for s = trials_include

        % cursor positions
        track_timestamp = easyeyes(track_on(s),:).posixTimeSec; 
        stim_timestamp = easyeyes(stim_on(s),:).posixTimeSec; 
        trial_ee_stim = easyeyes.posixTimeSec > (track_timestamp + timeframe_starttrack) & easyeyes.posixTimeSec < (stim_timestamp);
        currenttrial_ee_stimon = easyeyes(trial_ee_stim,:);
        [nearpointPx,crosshairPx,cursorPx] = getXYPositions_ee(currenttrial_ee_stimon);
        assert(size(crosshairPx,1) == size(cursorPx,1));
        assert(size(crosshairPx,2) == size(cursorPx,2));

        timestamps = currenttrial_ee_stimon.posixTimeSec;


        cursorXDeg = cursorPx(:,1) ./ PixelPerDeg;
        cursorYDeg = cursorPx(:,2) ./ PixelPerDeg;
        crosshairXDeg = crosshairPx(:,1) ./ PixelPerDeg;
        crosshairYDeg = crosshairPx(:,2) ./ PixelPerDeg;


        idx_record = 1:recFrames;

        cursor_x_mtx(s,:) = cursorXDeg(idx_record);
        cursor_y_mtx(s,:) = cursorYDeg(idx_record);
        crosshair_x_mtx(s,:) = crosshairXDeg(idx_record);
        crosshair_y_mtx(s,:) = crosshairYDeg(idx_record);

        
        diff_cursor_cross_Deg = [];
        for tt = 1:length(cursorXDeg)
            diff_cursor_cross_Deg(tt) = norm(cursorXDeg(tt) - crosshairXDeg(tt));
        end
        tracking_error(s,:) = diff_cursor_cross_Deg(idx_record);



    end

    % plot
    frame_interval = mean(diff(easyeyes.posixTimeSec(1:10*60)));
    xValues = linspace(0,frame_interval*recFrames,recFrames);

    % for ii = trials_include
    %     plot_tracking_xy(xValues,crosshair_x_mtx(ii,:),crosshair_y_mtx(ii,:),cursor_x_mtx(ii,:),cursor_y_mtx(ii,:))
    %     sgtitle(sprintf('Subject %d, trial %d', subj, ii))
    % end
    % 
    % for ii = trials_include
    %     figure;clf;
    %     hold on;
    %     xlabel('Time (s)');
    %     ylabel('Cursor tracking error (deg)');
    %     set(gca,'Fontsize',18)
    %     xlim([0 26])
    %     ylim([0 1])
    %     yline(0.15,'r--','LineWidth',2)
    %     plot(xValues,tracking_error(ii,:),'-','Color',[0 0.4470 0.7410],'LineWidth',3);
    %     title(sprintf('Subject %d, trial %d', subj, ii))
    %     set(gcf, 'Position', [100, 100, 800, 250])
    %     hold off;
    % end

    for ii = trials_include
        plot_trace_error(xValues,crosshair_x_mtx(ii,:), cursor_x_mtx(ii,:), tracking_error(ii,:))
        sgtitle(sprintf('Subject %d, trial %d', subj, ii))

    end


    
end



%% functions

function [] = plot_tracking_xy(xValues,crosshair_x,crosshair_y, ...
    cursor_x, cursor_y)


    figure;clf;

    subplot(2,1,1);
    hold on;
    % set up the graph
    xlabel('Time (s)');
    ylabel('X Position (deg)');
    set(gca,'Fontsize',18)
    xlim([0 26])
    ylim([-1 1])

    plot(xValues,crosshair_x,'-','Color','#808080','LineWidth',2);
    cursor_x = plot(xValues,cursor_x,'-','LineWidth',3);
    cursor_x.Color = [0 0.4470 0.7410 0.7];
    
    hold off;



    subplot(2,1,2);
    hold on;
    % set up figure
    xlabel('Time (s)');
    ylabel('Y Position (deg)');
    set(gca,'Fontsize',18)
    xlim([0 26])
    ylim([-1 1])

    plot(xValues,crosshair_y,'-','Color','#808080','LineWidth',2);
    cursor_y = plot(xValues,cursor_y,'-','LineWidth',3);
    cursor_y.Color = [0 0.4470 0.7410 0.7];


    hold off;    


end

function [] = plot_trace_error(xValues,crosshair_x, cursor_x, error_deg)


    figure;clf;

    subplot(2,1,1);
    hold on;
    % set up the graph
    xlabel('Time (s)');
    ylabel('X Position (deg)');
    set(gca,'Fontsize',18)
    xlim([0 26])
    ylim([-1 1])

    plot(xValues,crosshair_x,'-','Color','#808080','LineWidth',2);
    cursor_x = plot(xValues,cursor_x,'-','LineWidth',3);
    cursor_x.Color = [0 0.4470 0.7410 0.7];
    set(gcf, 'Position', [100, 100, 800, 600])
    title('Cursor and crosshair positions')

    hold off;



    subplot(2,1,2);
    hold on;
    % set up figure
    xlabel('Time (s)');
    ylabel('Cursor tracking error (deg)');
    set(gca,'Fontsize',18)
    xlim([0 26])
    ylim([0 1])
    yline(0.15,'r--','LineWidth',2)
    plot(xValues,error_deg,'-','Color',[0 0.4470 0.7410],'LineWidth',3);
    title('Cursor tracking error')



    hold off;    


end
