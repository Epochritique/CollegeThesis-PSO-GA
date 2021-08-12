function [Pbest, Gbest]= PSO_GetInitPGBest(PopNum, PosPop, ProbDim, FitVal)
    Pbest = zeros(PopNum, ProbDim+1);
    Gbest = zeros(1, ProbDim+1);
    
    % Get Initial Pbest and Gbest
    for RowNum = 1:PopNum
        Pbest(RowNum, :) = [PosPop(RowNum,:) FitVal(RowNum, 1)];
    end

    Gbest(ProbDim+1) = min(FitVal);
    for RowNum = 1:PopNum
        if Pbest(RowNum, ProbDim + 1) == Gbest(ProbDim+1)
            Gbest = Pbest(RowNum, :);
            break;
        end
    end
end
