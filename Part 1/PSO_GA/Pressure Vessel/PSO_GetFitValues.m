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
        
%         g1c=0;
%         g2c=0;
%         g3c=0;

        % Check if it satisfies constraint equations
%         g1 = 85.334407 + 0.0056858*PosPop(RowNum, 2)*PosPop(RowNum, 5) + 0.00026*PosPop(RowNum, 1)*PosPop(RowNum, 4) - 0.0022053*PosPop(RowNum, 3)*PosPop(RowNum, 5);
%         g2 = 80.51249 + 0.0071317*PosPop(RowNum, 2)*PosPop(RowNum, 5) + 0.0029955*PosPop(RowNum, 1)*PosPop(RowNum, 2) - 0.0021813*PosPop(RowNum, 3)^2;
%         g3 = 9.300961 + 0.0047026*PosPop(RowNum, 3)*PosPop(RowNum, 5) + 0.0012547*PosPop(RowNum, 1)*PosPop(RowNum, 3) + 0.0019085*PosPop(RowNum, 3)*PosPop(RowNum, 4);
        g1 = -PosPop(RowNum, 1) + 0.0193*PosPop(RowNum, 3);
        g2 = -PosPop(RowNum, 2) + 0.0095*PosPop(RowNum, 3);
        g3 = -(pi*(PosPop(RowNum, 3)^2)*PosPop(RowNum, 4)) - ((4/3)*pi*(PosPop(RowNum, 3)^3)) + 1296000;
        g4 = PosPop(RowNum, 4) - 240;
        if ~(g1<=0)
%             g1c=1;    
            probs=probs+1;% compute the number of errors
        end
        if ~(g2<=0)
%             g2c=1;     
            probs=probs+1;% compute the number of errors
        end
        if ~(g3<=0)
%             g3c=1;     
            probs=probs+1;% compute the number of errors
        end
        if ~(g4<=0)
%             g3c=1;     
            probs=probs+1;% compute the number of errors
        end
        if(probs~=0)
            FitVal(RowNum, 1) = abs(FitWorst) + abs(g1) + abs(g2) + abs(g3) + abs(g4);
%             FitVal(RowNum, 1) = FitWorst;
%             if(probs~=0)
%                FitVal(RowNum, 1) = FitVal(RowNum, 1) + abs(probs);
%             end
%             if g1c~=0
%                FitVal(RowNum, 1) = FitVal(RowNum, 1) + abs(g1);
%             end
%             if g2c~=0
%                FitVal(RowNum, 1) = FitVal(RowNum, 1) + abs(g2);
%             end
%             if g3c~=0
%                FitVal(RowNum, 1) = FitVal(RowNum, 1) + abs(g3);
%             end
        end
    end
end
