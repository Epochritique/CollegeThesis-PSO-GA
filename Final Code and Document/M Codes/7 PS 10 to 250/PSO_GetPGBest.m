% Function that records personal and global bests
function [Pbest, Gbest]= PSO_GetPGBest(PopNum, PosPop, ProbDim, FitVal, Pbest, Gbest)
    % Get Pbest and Gbest
    for RowNum = 1:PopNum
        [~, TotExc1] = PSO_GA_Eval(PosPop(RowNum,:), ProbDim, (ProbDim+1)/2);
        [~, TotExc2] = PSO_GA_Eval(Pbest(RowNum,1:ProbDim), ProbDim, (ProbDim+1)/2);
        
        % if the value is more 'fit', replace the current recorded personal
        % best
        if Pbest(RowNum, ProbDim+1) > FitVal(RowNum,1) && TotExc1 <= 0
            Pbest(RowNum, :) = [PosPop(RowNum,:) FitVal(RowNum, 1)];
        elseif TotExc2 > TotExc1
            Pbest(RowNum, :) = [PosPop(RowNum,:) FitVal(RowNum, 1)];
        end
    end
    

    % if the value is more 'fit', replace the current recorded global best
    if Gbest(ProbDim+1) > min(FitVal)
        for RowNum = 1:PopNum
            [~, TotExc1] = PSO_GA_Eval(PosPop(RowNum,:), ProbDim, (ProbDim+1)/2);
            [~, TotExc2] = PSO_GA_Eval(Gbest(1:ProbDim), ProbDim, (ProbDim+1)/2);
            if Pbest(RowNum, ProbDim + 1) == min(FitVal) && TotExc1 <=0
                Gbest = Pbest(RowNum, :);
                break;
            elseif TotExc2 > TotExc1
                Gbest = Pbest(RowNum, :);
                break;
            end
        end
    end
end
