function [drift] = calcDrift(currenttrial_ee,currenttrial_el,PixelPerDeg,framesPerSec)
% assumes that blinks and saccades are removed
% correct drift by taking the mean of the gaze error during fixations
% and subtracting it off of all the frames

    [nearpointPx,crosshairPx,cursorPx] = getXYPositions_ee(currenttrial_ee);
    [gazePx_downsampled] = getXYPositions_el(currenttrial_el,currenttrial_ee);

    % find indices for fixations
    % framesPerSec = 60;
    [fixIdx,gaze_fixatePx] = findFixations(gazePx_downsampled,PixelPerDeg,framesPerSec);
    crosshair_fixatePx = crosshairPx;
    crosshair_fixatePx(fixIdx,:) = [];


    assert(size(gaze_fixatePx,1) == size(crosshair_fixatePx,1));
    assert(size(gaze_fixatePx,2) == size(crosshair_fixatePx,2));

    % takes the average of fixations as error
    drift = mean((gaze_fixatePx-crosshair_fixatePx),1);

    % subtract the error for gazeErrorPx
    % gazeErrorPx = gazeErrorPx - repmat(drift,[size(gazeErrorPx,1),1]);

end