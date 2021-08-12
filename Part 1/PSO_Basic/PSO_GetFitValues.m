function FitVal = PSO_GetFitValues(PopNum, PosPop, ProbDim)
    FitVal = zeros(PopNum, 1);
    % Evaluate each particle's value
    for RowNum = 1:PopNum
        x = 0;
        for i = 1:ProbDim
            x = x + PosPop(RowNum, i)*2^(ProbDim-i);
        end
        FitVal(RowNum, 1) = x^2;
    end
    
%     for RowNum = 1:PopNum
%         FitVal(RowNum, 1) = PosPop(RowNum)^2;
%     end
end
