function [FitVal, ErrNum] = PSO_GA_Eval(Chrom, ProbDim)
    [dist_ij, node_stats, vehicle_cap]=Get_Graph_Data();

    % Obtain Sequence
    seq_route = ObtainSequence(Chrom, ProbDim);
    seq_inds = seq_route+1;
    % Calculate total distance
    total_distance = 0;
    for i=1:(ProbDim+2)-1
        total_distance = total_distance + dist_ij(seq_inds(i),seq_inds(i+1));
    end
    
    arrival_time = zeros(1,ProbDim+2);
    total_time_early = 0;
    total_time_late = 0;
    total_amt_collected = zeros(3,1);
    vehicle_num = 0;
    for i=1:(ProbDim+2)
        if seq_inds(i) == 1 % Depot
            arrival_time(i) = 0;
            vehicle_num = vehicle_num+1;
        else
            arrival_time(i) = max([node_stats(3,seq_inds(i-1))-arrival_time(i-1), 0]) + arrival_time(i-1)+dist_ij(seq_inds(i-1), seq_inds(i))/50 + node_stats(2,seq_inds(i-1));
            total_amt_collected(vehicle_num) = total_amt_collected(vehicle_num) + node_stats(1,seq_inds(i));
        end
    end
    
    % Calculate total time the vehicle was late and early
    for i=1:(ProbDim+2)
        if arrival_time(i) < node_stats(3, seq_inds(i))
            total_time_early = total_time_early + node_stats(3, seq_inds(i)) - arrival_time(i);
        elseif arrival_time(i) > node_stats(4, seq_inds(i))  
            total_time_late = total_time_late + arrival_time(i) - node_stats(4, seq_inds(i));
        end
    end
    
    ErrNum = 0;
    amt_exceed = 0;
    for i=1:3
        if total_amt_collected(i) > vehicle_cap
            amt_exceed = amt_exceed + total_amt_collected(i) - vehicle_cap;
            ErrNum = ErrNum+1;
        end
    end
    
    FitVal = total_distance + (total_time_early)*50 + (total_time_late)*50;
end
