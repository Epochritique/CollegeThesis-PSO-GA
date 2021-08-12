% Funtion that records the inital personal and global bests
function [Pbest, Gbest]= PSO_GetInitPGBest(PopNum, PosPop, ProbDim, FitVal)
    Pbest = zeros(PopNum, ProbDim+1);
    Gbest = zeros(1, ProbDim+1);
    
    % Record inital as Pbest
    for RowNum = 1:PopNum
        Pbest(RowNum, :) = [PosPop(RowNum,:) FitVal(RowNum, 1)];
    end
    
    % Record the Gbest with the best value 
    Gbest(ProbDim+1) = min(FitVal);
    for RowNum = 1:PopNum
        if Pbest(RowNum, ProbDim + 1) == Gbest(ProbDim+1)
            Gbest = Pbest(RowNum, :);
            break;
        end
    end
end
