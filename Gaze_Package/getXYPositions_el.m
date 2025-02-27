function [gazePx] = getXYPositions_el(currenttrial_el,currenttrial_ee)
% downsample gaze Px as the frame rate for eyelink is higher than that for
% EasyEyes

    gazePx = [];
    for xxx = 1 : size(currenttrial_ee,1)
        [val,ind] = min(abs(currenttrial_el.t1+14400 - currenttrial_ee.posixTimeSec(xxx)));
        % correct for calibration error
        gazePx(xxx,1) = currenttrial_el.gazeXYPix_1(ind);
        gazePx(xxx,2) = - currenttrial_el.gazeXYPix_2(ind); % correct for y-flip
    end

end