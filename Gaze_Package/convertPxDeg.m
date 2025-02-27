function [pxPerDeg,pxToDeg] = convertPxDeg(screenWidthPx,screenWidthCm,distance)
% calculate conversion between pixel and deg given pixel resolution, screen
% width, and distance from screen

    pxToCm = screenWidthCm/screenWidthPx; % how long is one px in cm?

    pxToDeg = 2 * atan(pxToCm/(2*distance)) * 180 / pi; % how many visual angles is one px?

    pxPerDeg = 1/pxToDeg;

end