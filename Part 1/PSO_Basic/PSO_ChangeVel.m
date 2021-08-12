function [TransPos, TransVel]= PSO_ChangeVel(PopNum, PosPop, VelPop, ProbDim, Pbest, Gbest, w, c1, c2)
    TransPos = zeros(PopNum, ProbDim);
    TransVel = zeros(PopNum, ProbDim);

    for RowNum = 1:PopNum
        % Calculate new velocity
        TempVel = VelPop(RowNum, 1:ProbDim)*w + c1*rand()*(Pbest(RowNum, 1:ProbDim) - PosPop(RowNum, 1:ProbDim)) + c2*rand()*(Gbest(1, 1:ProbDim) - PosPop(RowNum, 1:ProbDim));
        
        % Change position
        x = bi2de(PosPop(RowNum, :));
        y = floor(x + TempVel(1));
        if y > 31
            y = 31;
        elseif y < 0
            y = 0;
        end
        TempPos = fliplr(de2bi(y, ProbDim));
        % TempPos = PosPop + TempVel;
        
        % Save in respective arrays
        TransPos(RowNum, :) = TempPos;
        TransVel(RowNum, :) = TempVel;
    end
end
