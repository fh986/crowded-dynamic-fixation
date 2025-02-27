function [fixIdx,gazePx] = findFixations(gazePx_downsampled,PixelPerDeg,framesPerSec)

        [blinkIdx] = removeBlink(gazePx_downsampled,PixelPerDeg,framesPerSec);
        [saccIdx] = findSaccade(gazePx_downsampled,PixelPerDeg,framesPerSec);
        fixIdx = unique([blinkIdx,saccIdx]);

        gazePx = gazePx_downsampled;
        gazePx(fixIdx,:) = [];

end