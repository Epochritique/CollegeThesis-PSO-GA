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
        g1 = 85.334407 + 0.0056858*GA_Chroms(RowNum, 2)*GA_Chroms(RowNum, 5) + 0.0006262*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 4) - 0.0022053*GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 5);
        g2 = 80.51249 + 0.0071317*GA_Chroms(RowNum, 2)*GA_Chroms(RowNum, 5) + 0.0029955*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 2) - 0.0021813*GA_Chroms(RowNum, 3)^2;
        g3 = 9.300961 + 0.0047026*GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 5) + 0.0012547*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 3) + 0.0019085*GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 4);
        if ~(g1>=0 && g1<=92)
            probs=probs+1;% compute the number of errors
        end
        if ~(g2>=90 && g2<=110)
            probs=probs+1;% compute the number of errors
        end
        if ~(g3>=20 && g3<=25)
            probs=probs+1;% compute the number of errors
        end
        if(probs~=0)
            GA_FitVal(RowNum, 1) = FitWorst + g1 + g2 + g3;
        end
    end
end
