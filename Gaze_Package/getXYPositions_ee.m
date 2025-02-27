function [nearpointPx,crosshairPx,cursorPx] = getXYPositions_ee(currenttrial_ee)

    nearpointPx = str2num(cell2mat(currenttrial_ee.nearpointXYPx(1,:)));

    crosshairPx = [];
    cursorPx = [];

    for ii = 1:length(currenttrial_ee.crosshairPositionXYPx)

        aa = currenttrial_ee.crosshairPositionXYPx(ii);
        aa = cell2mat(aa);
        comma_idx = find(aa == ',', 1);
        if comma_idx == 2
            aa(comma_idx) = '';
        end
        crosshairPx(ii,:) = str2num(aa) - nearpointPx;

        aa = currenttrial_ee.cursorPositionXYPx(ii);
        aa = cell2mat(aa);
        comma_idx = find(aa == ',', 1);
        if comma_idx == 2
            aa(comma_idx) = '';
        end
        cursorPx(ii,:) = str2num(aa) - nearpointPx;

    end       

end
