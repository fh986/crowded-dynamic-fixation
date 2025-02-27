function [saccIdx,saccBool] = findSaccade(gazeErrorPx,pxPerDeg,framesPerSec)
% take a blink-removed trial and identify indexs of time at which saccades
% happen 

    saccBool = 0;
    
    gazeErrorDeg = gazeErrorPx/pxPerDeg;
    
    [velocity,acceleration] = eyeVelocity(gazeErrorDeg,framesPerSec);


    saccIdx = [];
    before = round(0.02*framesPerSec);
    after = round(0.15*framesPerSec); % identify the 20ms before and 150ms after the saccade initiation
    for ii = 1:length(acceleration)
       if velocity(ii+1) > 30 && acceleration(ii) > 8000 % Eyelink Default saccade detector
            saccIdx = [saccIdx,ii-before:ii+after]; 
       end
    end
    
    saccIdx(saccIdx<=0) = [];
    saccIdx(saccIdx>size(gazeErrorDeg,1)) = [];
    saccIdx = unique(saccIdx);


    if ~isempty(saccIdx)
        saccBool = 1;
    end


    % debug
    % if saccBool == 1
    %     disp('------')
    %     disp(gazeErrorDeg)
    %     disp(velocity)
    %     afterDeletion = gazeErrorDeg;
    %     afterDeletion(saccIdx,:) = [];
    %     if isempty(afterDeletion)
    %         disp('TRIAL EMPTY')
    %     end
    %     disp(afterDeletion)
    % end


end