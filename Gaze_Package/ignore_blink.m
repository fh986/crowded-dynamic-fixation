function [newGazePx,nantrialBool] = ignore_blink(gazePx,pxPerDeg,framesPerSec)
% detects blinks based on the criterion: saccade velocity > 2000 deg/sec
% replaces the 100 msecs before and after the blink with the gaze position
% in the frame after the blink
% if the blink happens at the end of the sequence, replace with the frame
% before the blink

   
    gazeDeg = gazePx/pxPerDeg;
    
    velocity = eyeVelocity(gazeDeg,framesPerSec);
    delFrames = round(0.1*framesPerSec);
    
    replaceIdx = [];
    for ii = 1:length(velocity)
    
       if any([velocity(ii) > 1000;  abs(gazeDeg(ii,1)) >20;abs(gazeDeg(ii,2)) >20])
            replaceIdx = [replaceIdx,ii-delFrames:ii+delFrames];
       end

    end
    
    replaceIdx(replaceIdx<=0) = [];
    replaceIdx(replaceIdx>size(gazeDeg,1)) = [];
    replaceIdx = unique(replaceIdx);

    % are there multiple blinks?
    difference = diff(replaceIdx);
    multipleBool = any(difference~=1);
    blinks = {};
    if ~multipleBool
        blinks{1} = replaceIdx;
    else
        breaks = find(difference~=1);
        breaks = [0,breaks,length(replaceIdx)];
        for blinkNumber = 1:length(breaks)-1
            blinks{blinkNumber} = replaceIdx(breaks(blinkNumber)+1:breaks(blinkNumber+1));
        end
    end

    newGazePx = gazePx;
    nantrialBool = 0;
    if ~isempty(replaceIdx)
        if length(replaceIdx) == length(gazePx)
            newGazePx = NaN(size(gazePx));
            nantrialBool = 1;
        else
            for blinkNumber = 1:length(blinks)
                thisBlink = blinks{blinkNumber};
                replacePosition = [];
                
                if ismember(size(gazeDeg,1),thisBlink)
                    replacePosition = gazePx(thisBlink(1)-1,:);
                else
                    replacePosition = gazePx(thisBlink(end)+1,:);
                end 
    
                for frame = thisBlink
                    newGazePx(frame,:) = replacePosition;
                end
            end
        end
    end


end