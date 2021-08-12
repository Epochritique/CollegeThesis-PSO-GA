% GA Selection Operator
function [Parent1, Parent2] = GA_Selection(GA_PS, GA_Chroms, ProbDim, GA_FitVal, GA_RouWheel)
    Parent1 = zeros(1, ProbDim);
    Parent2 = zeros(1, ProbDim);
    
    p1=rand(1,'double');
    p2=rand(1,'double');
    p1c=0; %Checker
    p2c=0; %Checker
    for RowNum = 1:GA_PS  
        % Select Parent 1
        if GA_RouWheel(RowNum) > p1 && p1c==0
            Parent1(1, :) = GA_Chroms(RowNum, :);
            p1c=1;
        end
        % Select Parent 2
        if GA_RouWheel(RowNum) > p2 && p2c==0
            Parent2(1, :) = GA_Chroms(RowNum, :);
            p2c=1;
        end
        if p1c == 1 && p2c == 1
            break;
        end
    end
    
end
