function [blinkIdx,gazeErrorPx,blinkBool] = removeBlink(gazeErrorPx,pxPerDeg,framesPerSec)
% detects blinks based on the criterion: saccade velocity > 2000 deg/sec
% deletes the 50 msecs before and after the blink


    blinkBool = 0;
    
    gazeErrorDeg = gazeErrorPx/pxPerDeg;
    
    velocity = eyeVelocity(gazeErrorDeg,framesPerSec);
    delFrames = round(0.1*framesPerSec);
    
    rm = [];
    for ii = 1:length(velocity)
    
       if velocity(ii) > 2000 || gazeErrorDeg(ii,1) < -40
            rm = [rm,ii-delFrames:ii+delFrames];
       end

    end
    
    rm(rm<=0) = [];
    rm(rm>size(gazeErrorDeg,1)) = [];
    rm = unique(rm);
    blinkIdx = rm;
    gazeErrorPx(rm,:) = [];

    if ~isempty(rm)
        blinkBool = 1;
    end

end
