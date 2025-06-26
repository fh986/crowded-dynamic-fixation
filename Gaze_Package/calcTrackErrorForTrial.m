function [trialError, tracking_errors_Px] = calcTrackErrorForTrial(currenttrial_ee)
    
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
           
        tracking_errors_Px = [];

        for x = 1 : length(crosshairPx)
            tracking_errors_Px(x) = norm(crosshairPx(x,:)-cursorPx(x,:));
        end
        trialError = median(tracking_errors_Px);
end


% function [trialError, tracking_errors_Px] = calcTrackErrorForTrial(crosshairPx,cursorPx)
% 
%         assert(length(crosshairPx) == length(cursorPx));
% 
%         tracking_errors_Px = [];
% 
%         for x = 1 : length(crosshairPx)
%             tracking_errors_Px(x) = norm(crosshairPx(x,:)-cursorPx(x,:));
%         end
% 
%         trialError = median(tracking_errors_Px);
% 
% end

