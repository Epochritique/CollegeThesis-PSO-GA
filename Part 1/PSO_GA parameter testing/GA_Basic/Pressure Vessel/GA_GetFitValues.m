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

        % Check if it satisfies constraint equations
        g1 = -GA_Chroms(RowNum, 1) + 0.0193*GA_Chroms(RowNum, 3);
        g2 = -GA_Chroms(RowNum, 2) + 0.0095*GA_Chroms(RowNum, 3);
        g3 = -(pi*(GA_Chroms(RowNum, 3)^2)*GA_Chroms(RowNum, 4)) - ((4/3)*pi*(GA_Chroms(RowNum, 3)^3)) + 1296000;
        g4 = GA_Chroms(RowNum, 4) - 240;
        if ~(g1<=0)
            probs=probs+1;% compute the number of errors
        end
        if ~(g2<=0)     
            probs=probs+1;% compute the number of errors
        end
        if ~(g3<=0)    
            probs=probs+1;% compute the number of errors
        end
        if ~(g4<=0)    
            probs=probs+1;% compute the number of errors
        end
        if(probs~=0)
            GA_FitVal(RowNum, 1) = abs(FitWorst) + abs(g1) + abs(g2) + abs(g3) + abs(g4);
        end
    end
end
