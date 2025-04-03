function pressed = return_Pressed(keybs)

[ keyIsDown, timeSecs, keyCode ] = KbCheck(keybs);
if keyIsDown
    keyPressed = KbName(keyCode);
    if iscell(keyPressed)
        for i = 1:length(keyPressed)
            pressed = ~isempty(strfind(keyPressed{i},'Return')) || ~isempty(strfind(keyPressed{i},'RETURN'));
        end
    else
        pressed = ~isempty(strfind(keyPressed,'Return')) || ~isempty(strfind(keyPressed,'RETURN'));
    end
    %     pressed = strcmp(keyPressed(1),'E')  ||  strcmp(keyPressed(1),'e') ;
else
    pressed = 0;
end
end