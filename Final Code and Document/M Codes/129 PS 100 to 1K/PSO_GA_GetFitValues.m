% Get the fitness values of the whole population
function FitVal = PSO_GA_GetFitValues(PopNum, PosPop, ProbDim, CustCnt)
    FitVal = zeros(PopNum, 1);
    TotExc = zeros(PopNum, 1);
    % Evaluate each particle's value
    for RowNum = 1:PopNum
        [FitVal(RowNum, 1), TotExc(RowNum,1)]= PSO_GA_Eval(PosPop(RowNum, :), ProbDim, CustCnt);
    end
    FitWorst = max(FitVal);
    % Evaluate each particle's value
    for RowNum = 1:PopNum
        if(TotExc(RowNum, 1) > 0)
            FitVal(RowNum, 1) = FitWorst + (TotExc(RowNum, 1));
            %FitWorst + (TotExc(RowNum, 1));
            %         else
            %disp('hey');
        end
    end
end
