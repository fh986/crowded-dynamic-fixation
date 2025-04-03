function pressed = s_Pressed(keybs)

[ keyIsDown, timeSecs, keyCode ] = KbCheck(keybs);
if keyIsDown
    keyPressed = KbName(keyCode);
    if iscell(keyPressed)
        for i = 1:length(keyPressed)
            pressed = ~isempty(strfind(keyPressed{i},'s')) || ~isempty(strfind(keyPressed{i},'S'));
        end
    else
        pressed = ~isempty(strfind(keyPressed,'s')) || ~isempty(strfind(keyPressed,'S'));
    end
    %     pressed = strcmp(keyPressed(1),'E')  ||  strcmp(keyPressed(1),'e') ;
else
    pressed = 0;
end
end