function GA_FitVal = GA_GetFitValues(GA_PS, GA_Chroms, ProbDim, DimMinMax)
    GA_FitVal = zeros(GA_PS, 1);
    
    % Evaluate each particle's value
    for RowNum = 1:GA_PS
        GA_FitVal(RowNum, 1) = GA_Eval(GA_Chroms(RowNum, :), ProbDim);
    end
    
    FitWorst = max(GA_FitVal);
    for RowNum = 1:GA_PS
        probs = 0;
        
        for i = 1:ProbDim
           if ~ (GA_Chroms(RowNum, i) >= DimMinMax(i, 1) && GA_Chroms(RowNum, i) <= DimMinMax(i, 2))
               probs=probs+1;% compute the number of boundary errors
           end
        end

        if(probs~=0)
            GA_FitVal(RowNum, 1) = abs(FitWorst)+probs;
        end
    end
end
