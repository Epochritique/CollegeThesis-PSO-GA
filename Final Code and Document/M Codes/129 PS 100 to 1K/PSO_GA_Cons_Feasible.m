% This function constructs a feasible solution using an insertion heuristic
% based on the capacity of a vehicle.
function Temp = PSO_GA_Cons_Feasible(ProbDim, DimMinMax)
    Temp = zeros(1, ProbDim);
    Temp1 = zeros(1, (ProbDim+1)/2);
    % Randomize Positions of Nodes
    for i = 1:(ProbDim+1)/2
       Temp1(1, i) = (DimMinMax(i, 2) - DimMinMax(i, 1))*rand(1,'double') + DimMinMax(i, 1);
    end
    
    % Keep randomized
    Temp(1:(ProbDim+1)/2) = Temp1(:);
    
    % Load data
    load('demand.mat');
    vec_cap = 12;
    
    Y = zeros(1, (ProbDim+1)/2);    % Checker
    sorted_temp = sort(Temp1);      % Sorted arrangement
    seq_route = zeros(1,(ProbDim+1)/2);
    AccWst = 0;% Accumulated Waste
    VecCnt = 1;% Vehicle Count
    % Sort the sequence
    for i = 1:(ProbDim+1)/2
        for j = 1:(ProbDim+1)/2
            if(sorted_temp(i) == Temp1(j) && Y(j)~=1)
                AccWst = AccWst + dmnd(j+1); % collect waste
                % If amount will exceed, place a depot node in the sequence
                % before this collection site
                if vec_cap < AccWst  
                    Temp(((ProbDim+1)/2) + VecCnt) = Temp1(seq_route(i-1)) + 0.001;
                    AccWst = dmnd(j+1);
                    VecCnt = VecCnt + 1;    
                end
                seq_route(i) = j;
                Y(j) = 1;
                break;
            end
        end
    end
    % The remaining vehicles are set to maximum value
    if ((ProbDim-1)/2) >= VecCnt
        for i=((ProbDim+1)/2) + VecCnt:ProbDim
            Temp(i) = DimMinMax(1, 2);
        end
    end
end
