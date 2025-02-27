function [trialError] = calcTrackErrorForTrial(crosshairPx,cursorPx)

        assert(length(crosshairPx) == length(cursorPx));

        tracking_errors_Px = [];

        for x = 1 : length(crosshairPx)
            tracking_errors_Px(x) = norm(crosshairPx(x,:)-cursorPx(x,:));
        end
        trialError = median(tracking_errors_Px);
        
end