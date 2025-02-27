function [newData] = rescalePxDeg(oldData,pxToDeg)
% converts data from pixel to deg

    newData = oldData.*pxToDeg;

end