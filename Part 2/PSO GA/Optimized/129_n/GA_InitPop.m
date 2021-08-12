function GA_Chroms = GA_InitPop(GA_PS, ProbDim, DimMinMax)
    GA_Chroms = zeros(GA_PS, ProbDim);
    % Randomly Generate GA_PS chromosomes
    for RowNum = 1:GA_PS
        TempPos = zeros(1, ProbDim);
        for i = 1:ProbDim
           % TempPos(1, i) = randi(2)-1;
           TempPos(1, i) = (DimMinMax(i,2) - DimMinMax(i,1))*rand(1,'double') + DimMinMax(i,1);
        end
        
        uniq=0;
        while uniq==0
            % Randomize Position
            for i = 1:ProbDim
               TempPos(1, i) = (DimMinMax(i, 2) - DimMinMax(i, 1))*rand(1,'double') + DimMinMax(i, 1);
            end
            % Check for Uniqueness
            if length(TempPos(1,:)) == length(unique(TempPos(1,:)))
                uniq=1;
            end
        end
        GA_Chroms(RowNum, :) = TempPos; % Add to positions
    end
end
