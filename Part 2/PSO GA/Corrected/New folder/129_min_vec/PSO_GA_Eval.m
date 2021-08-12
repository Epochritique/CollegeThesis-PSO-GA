% Get particular fitness
function [FitVal, TotExc] = PSO_GA_Eval(Chrom, ProbDim, CustCnt)
    % load data
    % load('vec_cap.mat');
    
%     load('demand.mat');
%     load('distance_data.mat');
%     load('vcl_cp');
    
    dmnd  = [0 7.47122 2.10123 9.27597 2.61348];
    vec_cap = 12;
    load('distance_data..mat');
    %load('vcl_cp');
    
    % Obtain Sequence
    seq_route = ObtainSequence(Chrom, ProbDim, CustCnt);
    seq_inds = seq_route+1;
    
    % Calculate total distance
    TotDist = 0;
    for i=1:length(seq_inds)-1
        TotDist = TotDist + dist_ij(seq_inds(i),seq_inds(i+1));
    end
    
    % Calculate total excess
    VecCnt = 0;
    AccWst = 0;
    TotExc = 0;
    TotUnu = 0;
    TotAcc = 0;
    TotUnc = 0;
    for i=2:length(seq_inds)
        if seq_inds(i) == 1
            if(vec_cap - AccWst)>=0
                TotUnu = TotUnu + abs(vec_cap - AccWst);
            else
                TotExc = TotExc + abs(vec_cap - AccWst);
            end
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
     TotExc = (TotExc/vec_cap)*(400);
% 
%     disp(seq_route);
%     disp(TotDist*0.27*46.20);
%     disp((TotAcc/VecCnt)*0.2*(400));
%     disp(840*(VecCnt));
%     disp(FitVal);
%     disp(TotExc);
end