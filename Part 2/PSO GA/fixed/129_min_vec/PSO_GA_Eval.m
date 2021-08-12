% Get particular fitness
function FitVal = PSO_GA_Eval(Chrom, ProbDim, DimMinMax)
    % load data
    % load('vec_cap.mat');
    load('nodes.mat');
    load('distance_data.mat');
    load('vcl_cp');
    
    % Obtain Sequence
    seq_route = ObtainSequence(Chrom, ProbDim);
    seq_inds = seq_route+1;
    
    SeqCur = 1; % Current Node
    VecCnt = 0; % Vehicle Count
    TotExSp = 0;% Total Extra Space
    TotDist = 0;% Total Distance
    while SeqCur < ProbDim
        PreNod = 1;      % Previous Node set to depot
        VecCnt = VecCnt + 1;% Increase vehicle count
        AccWst = 0;         % Accumulated Waste
        while AccWst < vec_cap
            TotDist = TotDist + dist_ij(PreNod, seq_inds(SeqCur));
            if AccWst + node_stats(1, seq_inds(SeqCur)) > vec_cap
                TotExSp = TotExSp+(vex_cap-AccWst);
                break;
            else
                AccWst = AccWst + node_stats(1, seq_inds(SeqCur));
                PreNod = seq_inds(SeqCur);
                SeqCur = SeqCur + 1;
            end
            if SeqCur > ProbDim
                break;
            end
        end
        TotDist = TotDist + dist_ij(PreNod, 1);
    end
    FitVal = TotDist + (1500)*nthroot(TotExSp,3);
end