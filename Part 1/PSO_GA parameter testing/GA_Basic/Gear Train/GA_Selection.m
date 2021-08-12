function [Parent1, Parent2] = GA_Selection(GA_PS, GA_Chroms, ProbDim, GA_FitVal, GA_RouWheel)
    Parent1 = zeros(1, ProbDim);
    Parent2 = zeros(1, ProbDim);
    % Select Parent 1
    RanSelVal = rand(1,'double');
    for RowNum = 1:GA_PS
        if GA_RouWheel(RowNum) > RanSelVal
            Parent1(1, :) = GA_Chroms(RowNum, :);
            break;
        end
    end
    
    % Select Parent 2
    RanSelVal = rand(1,'double');
    for RowNum = 1:GA_PS
        if GA_RouWheel(RowNum) > RanSelVal
            Parent2(1, :) = GA_Chroms(RowNum, :);
            break;
        end
    end
end
