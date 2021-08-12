% This function constructs a feasible solution using an insertion heuristic
% based on the capacity of a vehicle.
function Temp = PSO_GA_Cons_Feasible(ProbDim, DimMinMax)
    Temp = zeros(1, ProbDim);
    Temp1 = zeros(1, (ProbDim+1)/2);
    % Randomize Positions of Nodes
    for i = 1:(ProbDim+1)/2
       Temp1(1, i) = (DimMinMax(i, 2) - DimMinMax(i, 1))*rand(1,'double') + DimMinMax(i, 1);
    end
    
    Temp(1:(ProbDim+1)/2) = Temp1(:);
    
    load('demand.mat');
    vec_cap = 12;
    Y = zeros(1, (ProbDim+1)/2);
    sorted_temp = sort(Temp1);
    seq_route = zeros(1,(ProbDim+1)/2);
    AccWst = 0;
    VecCnt = 1;
    % Sort the sequence
    for i = 1:(ProbDim+1)/2
        for j = 1:(ProbDim+1)/2
            if(sorted_temp(i) == Temp1(j) && Y(j)~=1)
                AccWst = AccWst + dmnd(j+1);
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
    if ((ProbDim-1)/2) >= VecCnt
        for i=((ProbDim+1)/2) + VecCnt:ProbDim
            Temp(i) = 1;
        end
    end
end
