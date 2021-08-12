% Generate initial GA population
function GA_Chroms = GA_InitPop(GA_PS, ProbDim, DimMinMax)
    GA_Chroms = zeros(GA_PS, ProbDim);
    
    % Randomly Generate GA_PS chromosomes
    for RowNum = 1:GA_PS
        TempPos = PSO_GA_Cons_Feasible(ProbDim, DimMinMax);
        
        GA_Chroms(RowNum, :) = TempPos; % Add to positions
    end
end
