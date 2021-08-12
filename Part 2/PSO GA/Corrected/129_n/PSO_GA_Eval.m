% Get particular fitness
function FitVal = PSO_GA_Eval(Chrom, ProbDim, NodeCnt, VecCnt, alp, DimMinMax)
    % load data
    load('vec_cap.mat');
    load('nodes.mat');
    load('distance_data.mat');
    
    % Obtain Sequence
    seq_route = ObtainSequence(Chrom, ProbDim, NodeCnt);
    seq_inds = seq_route+1;
    
    % Calculate total distance
    total_distance = 0;
    for i=1:length(seq_inds)-1
        total_distance = total_distance + dist_ij(seq_inds(i),seq_inds(i+1));
    end
    
    % Calculate total excess
    vec_cnt = 1;
    accum_amt = 0;
    total_excess = 0;
    for i=2:length(seq_inds)
        if seq_inds(i) == 1
            total_excess = abs(vec_cap - accum_amt);
            accum_amt = 0;
        else
            accum_amt = accum_amt + dmnd(seq_inds(i));
        end
    end
end