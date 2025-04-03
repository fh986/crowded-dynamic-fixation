clc
clear
close all
easyeyes = readtable('HandyRedLemon754_labEasyEyes4_0001_2023-06-23_09h20.47.962_stimulus.csv');
el = readtable('S12JO_matlabOutput.csv');
stim_on = find(contains(easyeyes.easyEyesFunction,'trialRoutineBegin'));

%73.6 x 45.7 
%UP DOWN LEFT RIGHT
ct = 1
for s = 1 : length(stim_on)
%    
%     subplot(4,35,ct)
    time = easyeyes.posixTimeSec(stim_on(s)) - 4*60*60;
    
    trial = el.t1 > time & el.t1 < (time + 0.150);
    
    scatter(el.gazeXYDeg_1(trial),el.gazeXYDeg_2(trial),'filled');
    xlim([-73.6    73.6])
    ylim([-45.7  45.7])
    xlim([-15 15])
    ylim([-15 15])
    title(sprintf('timepoints = %i, trial = %i',sum(trial),s))
        ct = ct + 1;
axis square
hold on
    
end
    
    
    
    
    

