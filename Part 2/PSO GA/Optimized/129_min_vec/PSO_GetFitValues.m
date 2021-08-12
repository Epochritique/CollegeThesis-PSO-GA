% Get the fitness values of the whole population
function FitVal = PSO_GA_GetFitValues(PopNum, PosPop, ProbDim, CustCnt)
    FitVal = zeros(PopNum, 1);
    % Evaluate each particle's value
    for RowNum = 1:PopNum
        FitVal(RowNum, 1) = PSO_GA_Eval(PosPop(RowNum, :), ProbDim, CustCnt);
    end
end
