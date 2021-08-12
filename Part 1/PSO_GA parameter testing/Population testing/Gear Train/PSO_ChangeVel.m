function [TransPos, TransVel]= PSO_ChangeVel(PopNum, PosPop, VelPop, ProbDim, Pbest, Gbest, w, c1, c2, DimMinMax)
    TransPos = zeros(PopNum, ProbDim);
    TransVel = zeros(PopNum, ProbDim);

    for RowNum = 1:PopNum
        % Calculate new velocity
        TempVel = ceil(VelPop(RowNum, 1:ProbDim)*w + c1*rand()*(Pbest(RowNum, 1:ProbDim) - PosPop(RowNum, 1:ProbDim)) + c2*rand()*(Gbest(1, 1:ProbDim) - PosPop(RowNum, 1:ProbDim)));
        
        % Check if within Vmin, Vmax
        for i = 1:ProbDim
            VelBound = ceil((DimMinMax(i, 2) - DimMinMax(i, 1))/2);
            if TempVel(i)<-VelBound
                 TempVel(i) = -VelBound;
            elseif TempVel(i)>VelBound
                 TempVel(i) = VelBound;
            end
        end
        
        % Change position
%         x = bi2de(PosPop(RowNum, :));
%         y = floor(x + TempVel(1));
%         if y > 31
%             y = 31;
%         elseif y < 0
%             y = 0;
%         end
%         TempPos = fliplr(de2bi(y, ProbDim));

        TempPos = ceil(PosPop(RowNum, :) + TempVel);
        
        % Save in respective arrays
        TransPos(RowNum, :) = TempPos;
        TransVel(RowNum, :) = TempVel;
    end
end
