function FitVal = PSO_GetFitValues(PopNum, PosPop, ProbDim, DimMinMax)
    FitVal = zeros(PopNum, 1);
    Errs = zeros(PopNum, 1);
    % Evaluate each particle's value
    for RowNum = 1:PopNum
        [FitVal(RowNum, 1), Errs(RowNum,1)] = PSO_GA_Eval(PosPop(RowNum, :), ProbDim);
    end
    FitWorst = max(FitVal);
    for RowNum = 1:PopNum
        probs = 0;
        for i = 1:ProbDim
           if ~ (PosPop(RowNum, i) >= DimMinMax(i, 1) && PosPop(RowNum, i) <= DimMinMax(i, 2))
               probs=probs+1;% compute the number of boundary errors
           end
        end
        if(probs~=0) || Errs(RowNum,1) > 0
            FitVal(RowNum, 1) = FitVal(RowNum, 1)+Errs(RowNum,1)*250;
        end
    end
end
