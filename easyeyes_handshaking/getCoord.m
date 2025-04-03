function [x,y,tPsych,tEyeLink]=getCoord(wPtr,const)
% ----------------------------------------------------------------------
% [x,y,t]=getCoord(wPtr,const)
% ----------------------------------------------------------------------
% Goal of the function :
% Get gaze coordinates or mouse's in dummy mode.
% ----------------------------------------------------------------------
% Input(s) :
% wPtr = window pointer
% const = struct containacing  const.TEST and const.DOMEYE;
% ----------------------------------------------------------------------
% Output(s):
% x : X eye/mouse coordinate (horizontal)
% y : Y eye/mouse coordinate (vertical)
% t : machine time
% ----------------------------------------------------------------------
% Function created by Martin SZINTE (martin.szinte@gmail.com)
% edited by Nina HANNING (hanning.nina@gmail.com)
% Last update : 09 / 07 / 2009
% Project :     Eyelink toolbox
% Version :     1.0
% ----------------------------------------------------------------------

if const.TEST
    [x,y]=GetMouse(wPtr.main); % gaze position simulate by mouse position
    t = GetSecs;
else
    evt = Eyelink('newestfloatsample');
    x = evt.gx(const.recEye);
    y = evt.gy(const.recEye);
    tPsych = GetSecs;
%     tPsych = clock;
    tEyeLink = evt.time;
    if evt.gx(const.recEye) == -32768 || evt.gy(const.recEye) == -32768
        x = 0;
        y = 0;
    end
end
end