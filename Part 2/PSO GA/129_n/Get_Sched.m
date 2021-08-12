function [FitVal, ErrNum] = Get_Sched(Chrom, ProbDim)
    load('vec_cap.mat');
    load('nodes.mat');
    load('distance_data.mat');
    
    % Obtain Sequence
    seq_route = ObtainSequence(Chrom, ProbDim);
    seq_inds = seq_route+1;
    
    time_early = 0;
    time_late = 0;
    
    total_distance = 0;
    seq_curr = 1;
    vec = 1;
    
    day = 1;
    Y = sprintf('\nDay: %d',day);
    disp(Y);
    
    while seq_curr < 129
        timer = 0;
        Y = sprintf('\nVehicle: %d',vec);
        disp(Y);
        while timer < 9
            total_amt_col = 0;
            prev_node = 1;    
            while total_amt_col < vehicle_cap(vec)
                Y = sprintf('\nPath: %d-%d',prev_node-1, seq_inds(seq_curr)-1);
                disp(Y);
                timer = timer + dist_ij(prev_node, seq_inds(seq_curr))*(1/30);
                total_distance = total_distance + dist_ij(prev_node, seq_inds(seq_curr));
                if total_amt_col + node_stats(1, seq_inds(seq_curr)) > vehicle_cap(vec)
                    avail_space = vehicle_cap(vec) - total_amt_col;
                    timer = timer + node_stats(2, seq_inds(seq_curr))*(avail_space/node_stats(1, seq_inds(seq_curr)));
                    node_stats(1, seq_inds(seq_curr)) = node_stats(1, seq_inds(seq_curr)) - avail_space;
                    total_amt_col = vehicle_cap(vec);
                    prev_node = seq_inds(seq_curr);
                else
                    total_amt_col = total_amt_col + node_stats(1, seq_inds(seq_curr));
                    timer = timer + node_stats(2, seq_inds(seq_curr));
                    node_stats(1, seq_inds(seq_curr)) = 0;
                    prev_node = seq_inds(seq_curr);
                    seq_curr = seq_curr+1;
                end
                if seq_curr > 129
                    break;
                end
                disp(timer);
            end
            timer = timer + dist_ij(prev_node, 131)*(1/30) + dist_ij(131,1)*(1/30);
            total_distance = total_distance + dist_ij(prev_node, 131) + dist_ij(131,1);

            if seq_curr > 129
                break;
            end
        end
        if vec < 19
            vec = vec+1;
        else
            day = day+1;
            Y = sprintf('\nDay: %d',day);
            disp(Y);
            vec = 1;
        end
        if timer > 9
            time_early = time_early + (timer - 9);
        elseif timer < 9
            time_late = time_late + (9 - timer);
        end
    end
    
    ErrNum = (time_early + time_late)*30;
    FitVal = total_distance;
end
