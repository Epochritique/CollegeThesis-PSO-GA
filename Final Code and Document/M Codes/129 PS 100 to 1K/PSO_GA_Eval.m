% Function for getting particle fitness
function [FitVal, TotExc] = PSO_GA_Eval(Chrom, ProbDim, CustCnt)
    % Load data
    load('demand.mat');
    vec_cap = 12;
    load('distance_data..mat');
    
    % Obtain Sequence
    seq_route = ObtainSequence(Chrom, ProbDim, CustCnt);
    seq_inds = seq_route+1;
    
    % Calculate total distance
    TotDist = 0;
    for i=1:length(seq_inds)-1
        TotDist = TotDist + dist_ij(seq_inds(i),seq_inds(i+1));
    end
    
    % Calculate total excess
    VecCnt = 0;% vehicle count
    AccWst = 0;% accumulated by vehicle
    TotAcc = 0;% total accumulated
    TotUnc = 0;% uncollected
    for i=2:length(seq_inds)
        if seq_inds(i) == 1
            AccWst = 0;
            VecCnt = VecCnt+1;
        else
            AccWst = AccWst + dmnd(seq_inds(i));
            if(vec_cap - AccWst)>=0
                TotAcc = TotAcc + dmnd(seq_inds(i));
            else 
                TotUnc = TotUnc + dmnd(seq_inds(i));
            end
        end
    end
    
     FitVal = (TotDist*0.27*46.20)-(((TotAcc)/vec_cap)*(400)) + 840*VecCnt;
     TotExc = (TotUnc/vec_cap)*(400);
end