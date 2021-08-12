% Function that records personal and global bests
function [Pbest, Gbest]= PSO_GetPGBest(PopNum, PosPop, ProbDim, FitVal, Pbest, Gbest)
    % Get Pbest and Gbest
    for RowNum = 1:PopNum
        % if the balue is more 'fit', replace the current recorded personal
        % best
        if Pbest(RowNum, ProbDim+1) > FitVal(RowNum,1)
            Pbest(RowNum, :) = [PosPop(RowNum,:) FitVal(RowNum, 1)];
        end
    end
    
    % if the balue is more 'fit', replace the current recorded global best
    if Gbest(ProbDim+1) > min(FitVal)
        for RowNum = 1:PopNum
            if Pbest(RowNum, ProbDim + 1) == min(FitVal)
                Gbest = Pbest(RowNum, :);
                break;
            end
        end
    end
end
