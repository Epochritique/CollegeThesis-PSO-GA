function GA_RouWheel = GA_CreateWheel(GA_PS, GA_FitVal)
    % Build Roulette Wheel for GA Selection
    Probabilities  = zeros(GA_PS, 1);
    
    % get sum
    SumFit = (sum(GA_FitVal));
    
    % get fraction of fitness values
    for i = 1:GA_PS
        Probabilities(i, 1) = (GA_FitVal(i)/SumFit); 
    end
    
    GA_FitRatio = zeros(GA_PS, 1);
    
    % get fraction of proababilities
    SumProb = sum(Probabilities);
    
    % get probability ratio
    for i = 1:GA_PS % get proper probabilities
        GA_FitRatio(i, 1) = Probabilities(i)/SumProb; 
    end
    
    % Generate roulette wheel
    GA_RouWheel = zeros(GA_PS, 1);
    for RowNum = 1:GA_PS
        for i = 1:RowNum
            GA_RouWheel(RowNum) = GA_RouWheel(RowNum) + GA_FitRatio(i);
        end
    end
end
