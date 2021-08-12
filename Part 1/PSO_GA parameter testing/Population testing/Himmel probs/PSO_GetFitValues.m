function FitVal = PSO_GetFitValues(PopNum, PosPop, ProbDim, DimMinMax)
    FitVal = zeros(PopNum, 1);
    
    % Evaluate each particle's value
    for RowNum = 1:PopNum
        FitVal(RowNum, 1) = PSO_GA_Eval(PosPop(RowNum, :), ProbDim);
    end
    FitWorst = max(FitVal);
    for RowNum = 1:PopNum
        probs = 0;
        
        for i = 1:ProbDim
           if ~ (PosPop(RowNum, i) >= DimMinMax(i, 1) && PosPop(RowNum, i) <= DimMinMax(i, 2))
               probs=probs+1;% compute the number of boundary errors
           end
        end

        % Check if it satisfies constraint equations
        g1 = 85.334407 + 0.0056858*PosPop(RowNum, 2)*PosPop(RowNum, 5) + 0.0006262*PosPop(RowNum, 1)*PosPop(RowNum, 4) - 0.0022053*PosPop(RowNum, 3)*PosPop(RowNum, 5);
        g2 = 80.51249 + 0.0071317*PosPop(RowNum, 2)*PosPop(RowNum, 5) + 0.0029955*PosPop(RowNum, 1)*PosPop(RowNum, 2) - 0.0021813*PosPop(RowNum, 3)^2;
        g3 = 9.300961 + 0.0047026*PosPop(RowNum, 3)*PosPop(RowNum, 5) + 0.0012547*PosPop(RowNum, 1)*PosPop(RowNum, 3) + 0.0019085*PosPop(RowNum, 3)*PosPop(RowNum, 4);
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
            FitVal(RowNum, 1) = FitWorst + g1 + g2 + g3;
        end
    end
end
