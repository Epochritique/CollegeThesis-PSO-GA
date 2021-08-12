function FitVal = PSO_GetFitValues(PopNum, PosPop, ProbDim, DimMinMax)
    FitVal = zeros(PopNum, 1);
    
    % Evaluate each particle's value
    for RowNum = 1:PopNum
        FitVal(RowNum, 1) = 0.6224*PosPop(RowNum, 1)*PosPop(RowNum, 3)*PosPop(RowNum, 4) + 1.7781*PosPop(RowNum, 2)*(PosPop(RowNum, 3)^2) + 3.1661*(PosPop(RowNum, 1)^2)*PosPop(RowNum, 4) + 19.84*(PosPop(RowNum, 1)^2)*PosPop(RowNum, 3);
    end
    FitWorst = max(FitVal);
    for RowNum = 1:PopNum
        probs = 0;
        
        for i = 1:ProbDim
           if ~ (PosPop(RowNum, i) >= DimMinMax(i, 1) && PosPop(RowNum, i) <= DimMinMax(i, 2))
               probs=probs+1;% compute the number of boundary errors
           end
        end
        
        g1 = -PosPop(RowNum, 1) + 0.0193*PosPop(RowNum, 3);
        g2 = -PosPop(RowNum, 2) + 0.0095*PosPop(RowNum, 3);
        g3 = -(pi*(PosPop(RowNum, 3)^2)*PosPop(RowNum, 4)) - ((4/3)*pi*(PosPop(RowNum, 3)^3)) + 1296000;
        g4 = PosPop(RowNum, 4) - 240;
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
            FitVal(RowNum, 1) = abs(FitWorst) + abs(g1) + abs(g2) + abs(g3) + abs(g4);
        end
    end
end
