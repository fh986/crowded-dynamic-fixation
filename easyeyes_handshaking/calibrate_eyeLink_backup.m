
addpath(genpath('/Applications/Psychtoolbox/'))
% cavaaaa = cd;
% if a=='/'
%     a = PsychHID('Devices');
%     for i = 1:length(a), d(i) = strcmp(a(i).usageName, 'Keyboard'); end
%     keybs = find(d);
% else
%     keybs = [];
% end

oo.screen = 1;
oo.flipScreenHorizontally = false;
oo.useFractionOfScreenToDebug=0;
oo.resolution=Screen('Resolution',oo.screen);
[screenWidthMm,screenHeightMm]=Screen('DisplaySize',oo.screen);
oo.measuredScreenWidthCm = []
oo.measuredScreenHeightCm = []
oo.viewingDistanceCm = 40;

%%
oo.el = 1;
oo.trackRight = 1;
if isfinite(oo.measuredScreenWidthCm)
    screenWidthCm=oo.measuredScreenWidthCm;
else
    screenWidthCm=screenWidthMm/10;
end
if isfinite(oo.measuredScreenHeightCm)
    screenHeightCm=oo.measuredScreenHeightCm;
else
    screenHeightCm=screenHeightMm/10;
end

window=OpenWindow(oo);

screenRect=Screen('Rect',window);
screenWidthPix=RectWidth(screenRect);
screenHeightPix=RectHeight(screenRect);
pixPerDeg=0.05*screenHeightPix/atan2d(0.05*screenHeightCm,oo.viewingDistanceCm);
pixPerCm=screenWidthPix/screenWidthCm;
screenRect=Screen('Rect',window);
oo.screenRect=screenRect;
oo.viewingDistanceCm=oo.viewingDistanceCm;
oo.pixPerDeg=pixPerDeg;
oo.pixPerCm=pixPerCm;
cal.screen=oo.screen;
computer=Screen('Computer');
computer.system=strrep(computer.system,'Mac OS','macOS'); % Modernize the spelling.
[cal.screenWidthMm,cal.screenHeightMm]=Screen('DisplaySize',cal.screen);



oi = 1;

if oo.el %%% EYELINK (initialize once before all blocks)
    
    
    scr.main = window;
    scr.x_mid = screenWidthPix/2;
    scr.y_mid = screenHeightPix/2;
    scr.scr_sizeX = screenWidthMm;
    scr.scr_sizeY = screenHeightMm;
    const.TEST = 0;
    const.colBG = [220 220 220];
    const.calTargetRad_deg = 0.35;
    const.radCalib_val = 6.00;
    const.calTargetWidth_deg = 0.25;
    if oo.trackRight
        const.recEye = 2; % track right eye
    else
        const.recEye = 1;
    end
    [el,const]=initEyeLink(scr,const,oo,oi);
    %      key=0;
    %     while key ~=  0
    %         key = EyelinkGetKey(el);		% dump any pending local keys
    %     end
    
    calibresult = EyelinkDoTrackerSetup(el); % calibrate
    if calibresult==el.TERMINATE_KEY
        return
    end
end
sca

%%


function [el,const]=initEyeLink(scr,const,oo,oi)

% ----------------------------------------------------------------------
% [el,error]=initEyeLink(scr,const)
% ----------------------------------------------------------------------
% Goal of the function :
% Initializes eyeLink-connection, creates edf-file
% and writes experimental parameters to edf-file
% ----------------------------------------------------------------------
% Input(s) :
% scr : window pointer.
% const : struct containing many constant configuration.
% ----------------------------------------------------------------------
% Output(s):
% el : eye-link structure.c
% const : struct containing edfFileName
% ----------------------------------------------------------------------
% Function created by Martin SZINTE (martin.szinte@gmail.com)
% edited by Nina HANNING (hanning.nina@gmail.com)
% Last update : 12 / 03 / 2018
% Project : PirateAtt
% Version : 2.0
% ----------------------------------------------------------------------

%% Define EDF file name :
const.edffilename = 'XX.edf';

%% Modify different defaults settings :
el=EyelinkInitDefaults(scr.main);
el.backgroundcolour = const.colBG(1); %GrayIndex(el.window);
el.msgfontcolour    = WhiteIndex(el.window);
el.imgtitlecolour   = WhiteIndex(el.window);
el.targetbeep       = 0;
el.feedbackbeep     = 0;

el.calibrationtargetcolour= BlackIndex(el.window);
el.displayCalResults = 1;
el.eyeimgsize=50;
EyelinkUpdateDefaults(el);
el.backgroundcolour = const.colBG;
el.calibrationtargetsize=const.calTargetRad_deg;        % radius (deg) of calibration target
el.calibrationtargetwidth=const.calTargetWidth_deg; 	% radius (deg) of inside bull's eye of calibration target
el.txtCol = 15;
el.bgCol  = 0;

%% Initialization of the connection with the Eyelink Gazetracker.
if ~const.TEST
    dummymode = 0;
else
    dummymode = 1;
end

if ~EyelinkInit(dummymode)
    fprintf('Eyelink Init aborted.\n');
    Eyelink('Shutdown');
    Screen('CloseAll');
    return;
end

%% open file to record data to
res = Eyelink('Openfile', const.edffilename);
if res~=0
    fprintf('Cannot create EDF file ''%s'' ', const.edffilename);
    Eyelink('Shutdown');
    Screen('CloseAll');
    return;
end

% Describe general information on the experiment :
Eyelink('command', 'add_file_preamble_text ''Experiment by NMH''');

% make sure we're still connected.
if Eyelink('IsConnected')~=1 && ~dummymode
    fprintf('Not connected. exiting');
    Eyelink('Shutdown');
    Screen('CloseAll');
    return;
end

%% Set up tracker personal configuration :
rand('state',sum(100*clock));
%angle = 0:pi/3:5/3*pi;
const.numCali_dots = 4; % 4, 6, or 12
if const.numCali_dots == 4
    angle = 0:pi/2:3/2*pi;
elseif const.numCali_dots == 6 || const.numCali_dots == 12
    angle = 0:pi/3:5/3*pi;
end



if oo(oi).viewingDistanceCm <= 100
    refDist_val = const.radCalib_val/10/2;
else
    refDist_val = const.radCalib_val/10/1.25;
end

% compute calibration target locations
[cx1,cy1] = pol2cart(angle,refDist_val);
[cx2,cy2] = pol2cart(angle+(pi/6),refDist_val*0.65);
cx = round(scr.x_mid + scr.x_mid*[0 cx1 cx2]);
cy = round(scr.y_mid + scr.x_mid*[0 cy1 cy2]);

% start at center, select randomly, end at center
crp = randperm(const.numCali_dots);%+1;                   %crp = randperm(12)+1;
% crp = 1:4;

%crp = [2,4,3,1] % B T L R
crp = [4,2,3,1]; % T B L R

%crp = [3,4,2,1] % L T B R   BAD


% c(1:2:(numel(crp)+2)*2) = [cx(1) cx(crp) cx(1)];        %c(1:2:28) = [cx(1) cx(crp) cx(1)];
% c(2:2:(numel(crp)+2)*2) = [cy(1) cy(crp) cy(1)];        %c(2:2:28) = [cy(1) cy(crp) cy(1)];
c(1:2:(numel(crp))*2) = [cx(crp+1)];
c(2:2:(numel(crp))*2) = [cy(crp+1)];
c = [cx(1),cy(1),c,cx(1),cy(1)];

% compute validation target locations (ca libration targets smaller radius)
[vx1,vy1] = pol2cart(angle,refDist_val*0.85);
[vx2,vy2] = pol2cart(angle+pi/6,refDist_val*0.5);

vx = round(scr.x_mid + scr.scr_sizeY/2*[0 vx1 vx2]);
vy = round(scr.y_mid + scr.y_mid*[0 vy1 vy2]);

% start at center, select randomly, end at center
vrp = randperm(const.numCali_dots);%+1;                   %vrp = randperm(12)+1;
%v(1:2:(numel(crp)+2)*2) = [vx(1) vx(vrp) vx(1)];        %v(1:2:28) = [vx(1) vx(vrp) vx(1)];
%v(2:2:(numel(crp)+2)*2) = [vy(1) vy(vrp) vy(1)];        %v(2:2:28) = [vy(1) vy(vrp) vy(1)];
v(1:2:(numel(vrp))*2) = [vx(vrp+1)];
v(2:2:(numel(vrp))*2) = [vy(vrp+1)];
v = [vx(1),vy(1),v,vx(1),vy(1)];


Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, scr.scr_sizeX-1, scr.scr_sizeY-1);
Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, scr.scr_sizeX-1, scr.scr_sizeY-1);

Eyelink('command', 'calibration_type = HV5'); % Eyelink('command', 'calibration_type = HV13');
Eyelink('command', 'generate_default_targets = NO');

Eyelink('command', 'randomize_calibration_order 0');
Eyelink('command', 'randomize_validation_order 0');
Eyelink('command', 'cal_repeat_first_target 0');
Eyelink('command', 'val_repeat_first_target 0');

%Eyelink('command', 'calibration_samples=14');
Eyelink('command', 'calibration_samples=6');
%Eyelink('command', 'calibration_sequence=0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13');
Eyelink('command', 'calibration_sequence=0, 1, 2, 3, 4, 5');
%Eyelink('command', sprintf('calibration_targets = %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i',c));
Eyelink('command', sprintf('calibration_targets = %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i',c));

%Eyelink('command', 'validation_samples=14');
Eyelink('command', 'validation_samples=6');
%Eyelink('command', 'validation_sequence=0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13');
Eyelink('command', 'validation_sequence=0, 1, 2, 3, 4, 5');
%Eyelink('command', sprintf('validation_targets = %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i',v));
Eyelink('command', sprintf('validation_targets = %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i',c));




% Set parser
Eyelink('command', 'saccade_velocity_threshold = 35');
Eyelink('command', 'saccade_acceleration_threshold = 9500');

Eyelink('command', 'file_event_filter = RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
Eyelink('command', 'file_sample_data = RIGHT,GAZE,AREA');
Eyelink('command', 'link_event_filter = RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
Eyelink('command', 'link_sample_data  = RIGHT,GAZE,AREA');

Eyelink('command', 'heuristic_filter = 1 1');



%% Set pupil Tracking model in camera setup screen  (no = centroid. yes = ellipse)
Eyelink('command', 'use_ellipse_fitter =  NO');

%% set sample rate in camera setup screen
Eyelink('command', 'sample_rate = %d',1000);


% Test mode of eyelink connection :
status = Eyelink('IsConnected');
switch status
    case -1
        fprintf(1, '\tEyelink in dummymode.\n\n');
    case  0
        fprintf(1, '\tEyelink not connected.\n\n');
    case  1
        fprintf(1, '\tEyelink connected.\n\n');
end

% make sure we're still connected.
if Eyelink('IsConnected')~=1 && ~dummymode
    fprintf('Not connected. exiting');
    Eyelink('Shutdown');
    Screen('CloseAll');
    return;
end

end