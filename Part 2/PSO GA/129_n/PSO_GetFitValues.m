% Get the fitness values of the whole population
function FitVal = PSO_GetFitValues(NodeCnt, VecCnt, alp, PopNum, PosPop, ProbDim, DimMinMax)
    FitVal = zeros(PopNum, 1);
    % Evaluate each particle's value
    for RowNum = 1:PopNum
        FitVal(RowNum, 1) = PSO_GA_Eval(PosPop(RowNum, :), ProbDim, NodeCnt, VecCnt, alp, DimMinMax);
    end
end
