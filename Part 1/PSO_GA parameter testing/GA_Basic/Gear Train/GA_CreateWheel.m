function GA_RouWheel = GA_CreateWheel(GA_PS, GA_FitVal)
    % Build Roulette Wheel for GA Selection
    Probabilities  = zeros(GA_PS, 1);
    Sum = (sum(GA_FitVal));
    for i = 1:GA_PS % perturb
        Probabilities(i, 1) = (1/(GA_FitVal(i)/Sum)); 
    end
    GA_FitRatio = zeros(GA_PS, 1);
    Sum = sum(Probabilities);
    for i = 1:GA_PS % get proper probabilities
        GA_FitRatio(i, 1) = Probabilities(i)/Sum; 
    end
    GA_RouWheel = zeros(GA_PS, 1);
    for RowNum = 1:GA_PS
        for i = 1:RowNum
            GA_RouWheel(RowNum) = GA_RouWheel(RowNum) + GA_FitRatio(i);
        end
    end
end
