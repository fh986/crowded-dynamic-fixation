% easyeyesHandshaking.m
% these codes were provided by Jan Kurzawski
% and edited by Helen Hu
% To run these codes, Psychtoolbox and CriticalSpacing must be installed.


clc
clear

% To run the codes, please replace these two lines with your paths for
% Psychtoolbox and CriticalSpacing:
% addpath(genpath('.../Psychtoolbox/'))
% addpath(genpath('.../CriticalSpacing'))

a = cd;
if a(1)=='/'
    a = PsychHID('Devices');
    for i = 1:length(a), d(i) = strcmp(a(i).usageName, 'Keyboard'); end
    keybs = find(d);
else
    keybs = [];
end


name = input('Participant ID? ','s');


%%


oo.beginningTime=now;
timeVector=datevec(oo.beginningTime);



oo.dataFolder = './../data/easyeyes';

url = 'https://easyeyes-server-netlify.netlify.app/.netlify/functions/api/';

% Wait for easyeyes sending the "Record" message

%%%%%%%%%%%%%%%%%%
calibrate_eyeLink
%%%%%%%%%%%%%%%%%%

while true
%     response = webread(append(url, 'easyeyes'));
%     if startsWith(response.message, 'Record')
%         break;
%     end
 pressed = return_Pressed(keybs); if pressed,break;end
 
end


% Get file name after receiving the "Record" message
% response = webread(append(url, 'easyeyes'));
% fileName = extractAfter(response.message, 'Record ');
% disp(fileName);

% oo.dataFilename=sprintf('%s.%d.%d.%d.%d.%d.%d',name,round(timeVector));

record_eyelink



% Start recording and send message "Recording" to easyeyes
% response = webwrite(append(url, 'matlab'), 'Recording');
% disp(response);



% %Recording data and wait for easyeyes sending the "Save" message
% while true
%     response = webread(append(url, 'easyeyes'));
%     if strcmp(response.message, 'Save')
%         break;
%     end
% end

% while true
% %     response = webread(append(url, 'easyeyes'));
% %     if startsWith(response.message, 'Record')
% %         break;
% %     end
%  pressed = return_Pressed(keybs); if pressed,break;end
%  
% end


oo.dataFilename = 'test';
stop_eyeLink

% Save data and send message "Saved" to easyeyes

response = webwrite(append(url, 'matlab'), 'Saved');
disp(response);


dataEyeLink.gazeXYDeg(:,2) = dataEyeLink.gazeXYDeg(:,2)* -1;
dataEyeLink.gazeXYPix(:,2) = dataEyeLink.gazeXYPix(:,2)* -1;



writetable(struct2table(dataEyeLink), sprintf("%s_matlabOutput.csv",name)) % CURRENTLY OVERWRITING. WE SHOULD AT LEAST ADD FILENAME + PARTICIPANT ID. 




