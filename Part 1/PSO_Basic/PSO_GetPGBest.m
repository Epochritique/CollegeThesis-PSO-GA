function [Pbest, Gbest]= PSO_GetPGBest(PopNum, PosPop, ProbDim, FitVal, Pbest, Gbest)
    % Get Pbest and Gbest
    for RowNum = 1:PopNum
        if Pbest(RowNum, ProbDim+1) > FitVal(RowNum,1)
            Pbest(RowNum, :) = [PosPop(RowNum,:) FitVal(RowNum, 1)];
        end
    end
    
    if Gbest(ProbDim+1) > min(FitVal)
        Gbest(ProbDim+1) = min(FitVal);
        for RowNum = 1:PopNum
            if Pbest(RowNum, ProbDim + 1) == Gbest(ProbDim+1)
                Gbest = Pbest(RowNum, :);
                break;
            end
        end
    end
end
