clc
clear
close all;

    
% load('./../data/easyeyes/JanK.2023.6.15.21.27.23.mat')
    
GetSecsRef = datetime('2001-01-01');
EasyEyesRef = datetime('1970-01-01');

newtime = GetSecs/60/60/24/365


GetSecsTime = (mytable.t1(1));
EasyEyesTime = (1686878929.195);

GetSecsTime2 = GetSecsTime/60/60/24/365
EasyEyesTime2 = EasyEyesTime/1000/60/60/24/365



EasyEyes_updated = EasyEyesTime - 978307200

X = convertTo(GetSecsTime,'epochtime','Epoch','2001-01-01')



current_datetime = datetime('now');

% milliseconds_later = 175348;
% miliseconds_GetSecs = 978307200;
seconds_later = milliseconds_later / 1000;

updated_datetime = current_datetime - seconds(seconds_later);
% 
% updated_date = datetime(updated_datetime);
% current_datetime = datetime(current_datetime);
% 
% disp(['Original Date: ' char(current_datetime)]);
% disp(['Updated Date: ' char(updated_date)]);


X = convertTo(D,'epochtime','Epoch','2001-01-01')