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
        
        if(probs~=0)
            FitVal(RowNum, 1) = abs(FitWorst)+probs;
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
