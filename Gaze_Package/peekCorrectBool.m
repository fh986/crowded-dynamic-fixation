function [correctBool] = peekCorrectBool(conditionName,gazePx,last_crosshairPx,px_to_deg,eyedrift)

    correctBool = NaN;

    [gaze_errors_Deg] = calcGazeError(gazePx,eyedrift,last_crosshairPx,px_to_deg);
    [a,Idx] = max(gaze_errors_Deg);
    peekPosition = gazePx(Idx,:);
    xPosition_relative = peekPosition(1) - last_crosshairPx(1);

    % disp('--------------newTrial--------------')
    % disp(conditionName)
    if xPosition_relative > 0 % participant looked right
        if contains(conditionName,'Left')
            correctBool = 0;%disp('Incorrect')
        else
            correctBool = 1;%disp('Correct')
        end
    elseif xPosition_relative < 0 % participant looked left
        if contains(conditionName,'Left')
            correctBool = 1;%disp('Correct')
        else
            correctBool = 0;%disp('Incorrect')
        end
    end
    
    assert(~isnan(correctBool))
    
end