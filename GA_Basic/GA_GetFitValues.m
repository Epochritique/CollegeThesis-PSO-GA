function GA_FitVal = GA_GetFitValues(GA_PS, GA_Chroms, ProbDim)
    GA_FitVal = zeros(GA_PS, 1);
    % Evaluate each particle's value
    for RowNum = 1:GA_PS
        GA_FitVal(RowNum, 1) = GA_Eval(GA_Chroms(RowNum, :), ProbDim); 
    end
end
