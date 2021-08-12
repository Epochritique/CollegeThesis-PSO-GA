% Function that recalculates the velocity and particle locations
function [TransPos, TransVel]= PSO_ChangeVel(PopNum, PosPop, VelPop, ProbDim, Pbest, Gbest, w, c1, c2, DimMinMax)
    TransPos = zeros(PopNum, ProbDim);
    TransVel = zeros(PopNum, ProbDim);

    for RowNum = 1:PopNum
        % Calculate new velocity
        TempVel = VelPop(RowNum, 1:ProbDim)*w + c1*rand()*(Pbest(RowNum, 1:ProbDim) - PosPop(RowNum, 1:ProbDim)) + c2*rand()*(Gbest(1, 1:ProbDim) - PosPop(RowNum, 1:ProbDim));
        
        % Check if within Vmin, Vmax
        for i = 1:ProbDim
            VelBound = (DimMinMax(i, 2) - DimMinMax(i, 1));
            % Adjust to bounds
            if TempVel(i)<-VelBound
                 TempVel(i) = -VelBound;
            elseif TempVel(i)>VelBound
                 TempVel(i) = VelBound;
            end
        end

        % Maintain feasible population using PSO-Flyback
        TempPos = PosPop(RowNum, :) + TempVel;
        [~, TotExc] = PSO_GA_Eval(TempPos(:), ProbDim, (ProbDim+1)/2);
        if TotExc > 1e-5
            TempPos = PosPop(RowNum, :);
        end
        
        % Save in respective arrays
        TransPos(RowNum, :) = TempPos;
        TransVel(RowNum, :) = TempVel;
    end
end
