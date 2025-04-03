
sampling = 0.01;
mycounter = 0:sampling:30*60;

collect_eyasyEyesinfo = 5;
firstime = 0;
firstimte_msg = 0;

Eyelink('startrecording');
Eyelink('message',  num2str(GetSecs));

clear dataEyeLink
dataEyeLink.t1 = NaN(length(mycounter),1);
dataEyeLink.t2 = NaN(length(mycounter),1);
allX = NaN(length(mycounter),1);
allY = NaN(length(mycounter),1);

ct = 1;
t0 = GetSecs;
t = 0;

% while true
% 
% %     escape = escPressed(keybs); if escape,break;end
% 
% 
% 
%     t = GetSecs - t0;
%     if ~mod(round(t,3),sampling) && ~firstime
%         
%         [x,y,timeEyes,timePsych] = getCoord(scr,const);
%         dataEyeLink.t1(ct) = timeEyes;
%         dataEyeLink.t2(ct) = timePsych;
%         dataEyeLink.x(ct) = x;
%         dataEyeLink.y(ct) = y;
%         ct = ct + 1;
%         
%         firstime = 1;
%         
%     elseif mod(round(t,3),sampling)
%         firstime = 0;
%     end
%     
% %     if ~mod(round(t,3),collect_eyasyEyesinfo)&& ~firstimte_msg
% %         response = webread(append(url, 'easyeyes'));
% %     elseif mod(round(t,3),collect_eyasyEyesinfo)
% %         firstimte_msg = 0;
% % 
% %     end
% 
% %     if strcmp(response.message, 'Save')
% %         break;
% %     end
% end
% 


elapsedTime = 0;

while true
    
      escape = escPressed(keybs); if escape,break;end
      %response = webread(append(url, 'easyeyes'));

    t = GetSecs - t0;

    if t >= elapsedTime
        elapsedTime = elapsedTime + sampling;

        [x, y, timeEyes, timePsych] = getCoord(scr, const);
        dataEyeLink.t1(ct) = timeEyes; % t1 is from GetSecs
        dataEyeLink.t2(ct) = timePsych; % t2 is from eyelink file, the "time" column
        allX(ct) = x;
        allY(ct) = y;
        ct = ct + 1;
    end
    
    % Rest of your code...
end


% allY = allY * -1;

xyPix= [allX,allY]-[scr.x_mid scr.y_mid];


dataEyeLink.gazeXYPix=xyPix;
dataEyeLink.gazeXYDeg=xyPix./oo.pixPerDeg;

GetSecs2posixOffset = convertTo(datetime(clock),'posix') - GetSecs;
dataEyeLink.t1 = dataEyeLink.t1 + GetSecs2posixOffset;



function pressed = escPressed(keybs)

[ keyIsDown, timeSecs, keyCode ] = KbCheck(keybs);
if keyIsDown
    keyPressed = KbName(keyCode);
    if iscell(keyPressed)
        for i = 1:length(keyPressed)
            pressed = ~isempty(strfind(keyPressed{i},'esc')) || ~isempty(strfind(keyPressed{i},'ESC'));
        end
    else
        pressed = ~isempty(strfind(keyPressed,'esc')) || ~isempty(strfind(keyPressed,'ESC'));
    end
    %     pressed = strcmp(keyPressed(1),'E')  ||  strcmp(keyPressed(1),'e') ;
else
    pressed = 0;
end
end
