function GA_Chroms = GA_InitPop(GA_PS, ProbDim, DimMinMax)
    GA_Chroms = zeros(GA_PS, ProbDim);
    % Randomly Generate GA_PS chromosomes
    for RowNum = 1:GA_PS
        TempPos = zeros(1, ProbDim);
        for i = 1:ProbDim
           % TempPos(1, i) = randi(2)-1;
           TempPos(1, i) = randi([DimMinMax(i, 1), DimMinMax(i, 2)], 1);
        end
        GA_Chroms(RowNum, :) = TempPos; % Add to positions
    end
end
