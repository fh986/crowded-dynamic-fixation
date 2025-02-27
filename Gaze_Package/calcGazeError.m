function [gaze_errors_Deg] = calcGazeError(gazePx,eyedrift,last_crosshairPx,px_to_deg)

    gaze_errors_Px = [];
    for x = 1 : size(gazePx,1)
        gaze_correct = gazePx(x,:) - eyedrift;
        gaze_errors_Px(x) = norm(gaze_correct-last_crosshairPx);
    end
    gaze_errors_Deg = gaze_errors_Px .* px_to_deg;


end