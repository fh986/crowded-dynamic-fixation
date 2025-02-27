function [] = plotGaze(easyeyes,stim_timestamp,eyelink,pxPerDeg,gaze_correction,trialName)

    trial_ee_stim = easyeyes.posixTimeSec > (stim_timestamp - 0.5) & easyeyes.posixTimeSec < (stim_timestamp + 0.5);
    currenttrial_ee_stimon = easyeyes(trial_ee_stim,:);

    nearpointPx = str2num(cell2mat(currenttrial_ee_stimon.nearpointXYPx(1,:)));
    % last_crosshairPx = str2num(cell2mat(currenttrial_ee_stimon.crosshairPositionXYPx(1,:))) - nearpointPx;
    crosshairPx = str2num(cell2mat(currenttrial_ee_stimon.crosshairPositionXYPx)) - nearpointPx;

    time_ee = currenttrial_ee_stimon.posixTimeSec; 

    trial_eyelink = eyelink.t1+14400 > (stim_timestamp - 0.5) & eyelink.t1+14400 < (stim_timestamp + 0.5);
    currenttrial_el = eyelink(trial_eyelink,:);
    % gazePx = [];
    % gazePx(:,1) = currenttrial_el.gazeXYPix_1 - gaze_correction(1);
    % gazePx(:,2) = - currenttrial_el.gazeXYPix_2 - gaze_correction(2);

    for xxx = 1 : size(currenttrial_ee_stimon,1)

        % find gaze position for each timepoint we looked at for cursor
        % error
        [val,ind] = min(abs(currenttrial_el.t1+14400 - currenttrial_ee_stimon.posixTimeSec(xxx)));
        gazePx(xxx,1) = currenttrial_el.gazeXYPix_1(ind) - gaze_correction(1);
        gazePx(xxx,2) = currenttrial_el.gazeXYPix_2(ind) - gaze_correction(2);% usually we do a y-flip
        % but here we are subtracting from the crosshair and we have to
        % "flip" again so that the Y graph we see will reflect the position
        % of the Y, i.e., positive number = up, negative number = down with
        % respect to the crosshair

    end
    % 
    % disp(crosshairPx)
    % disp(gazePx./60)
    gazeError = gazePx - crosshairPx;
    gazeError = gazeError ./ pxPerDeg;
    % disp(gazeError)



    % plot
    % Define the stimulus onset and offset times
    stimulus_onset = 0; % Example: stimulus onset at 2 seconds
    stimulus_offset = 0.15; % Example: stimulus offset at 4 seconds
    
    xValues = currenttrial_ee_stimon.posixTimeSec - stim_timestamp;
 
    figure;
    subplot(2,1,1);
    plot(xValues,gazeError(:,1),'k','LineWidth',2)
    ylabel('X Position (Deg)')
    xlabel('Time (s)')
    ylim([-10.1,10.1])
    xlim([-0.5,0.5])
    yticks(-10:1:10)
    % Create a patch to highlight the stimulus duration
    stimulus_duration = [stimulus_onset stimulus_offset stimulus_offset stimulus_onset];
    y_limits = ylim(); % Get the current y-axis limits
    stimulus_height = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)];
    patch(stimulus_duration, stimulus_height, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    title(trialName)
    hold on;
    yline(10,'r--','LineWidth',2);
    yline(-10,'r--','LineWidth',2);
    hold off;



    % 
    subplot(2,1,2);
    plot(xValues,gazeError(:,2),'k','LineWidth',2)
    ylabel('Y Position (Deg)')
    xlabel('Time (s)')    
    ylim([-10,10])
    xlim([-0.5,0.5])
    yticks(-10:1:10)
    % % Create a patch to highlight the stimulus duration
    stimulus_duration = [stimulus_onset stimulus_offset stimulus_offset stimulus_onset];
    y_limits = ylim(); % Get the current y-axis limits
    stimulus_height = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)];
    patch(stimulus_duration, stimulus_height, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on;
    yline(10,'r--','LineWidth',2);
    yline(-10,'r--','LineWidth',2);
    hold off;

    % remove blinks
    % [gazePx,blinkBool] = removeBlink(gazePx,PixelPerDeg);


end
    
  