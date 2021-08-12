% GA Selection Operator
function [Parent1, Parent2] = GA_Selection(GA_PS, GA_Chroms, ProbDim, GA_FitVal, GA_RouWheel)
    Parent1 = zeros(1, ProbDim);
    Parent2 = zeros(1, ProbDim);
    
    for RowNum = 1:GA_PS  
        % Select Parent 1
        if GA_RouWheel(RowNum) > rand(1,'double')
            Parent1(1, :) = GA_Chroms(RowNum, :);
            break;
        end
        % Select Parent 2
        if GA_RouWheel(RowNum) > rand(1,'double')
            Parent2(1, :) = GA_Chroms(RowNum, :);
            break;
        end
    end
end
