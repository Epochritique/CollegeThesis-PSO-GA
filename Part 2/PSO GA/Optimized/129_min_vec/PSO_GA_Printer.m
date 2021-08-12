% % Get particular fitness
function [inf, FitVal] = PSO_GA_Printer(Chrom, ProbDim, CustCnt)
    % load data
    % load('vec_cap.mat');
    % load('demand.mat');
    dmnd = [0 7.47122 2.10123 9.27597 2.61348]; %[0 3 4 5];
    vec_cap = 12;
    load('distance_data_shorts.mat');
%     
% dist_ij = [ 0       0.75    0.5   1  ;
%             0.75    0       1.5   2  ;
%             0.5     1.5    0     2.5;
%             1       2       2.5   0  ;];

    %load('distance_data.mat');
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
    inf = 0;
x     disp(seq_route);
%     disp(TotAcc);
%     disp(seq_route);
%     disp(TotDist*0.27*46.20);
%     disp((TotAcc/vec_cap)*(400));
%     disp(500*(VecCnt));
    %disp(TotExc);
    FitVal = (TotDist*0.27*46.20)-(((TotAcc)/vec_cap)*(400)) + 840*VecCnt;
    if (TotExc > 0)
        inf =1;
    end
    disp(FitVal);
    
%     % Obtain Sequence
%     seq_route = ObtainSequence(Chrom, ProbDim);
%     seq_inds = seq_route+1;
%     
%     SeqCur = 1; % Current Node
%     VecCnt = 0; % Vehicle Count
%     TotExSp = 0;% Total Extra Space
%     TotDist = 0;% Total Distance
%     while SeqCur <= ProbDim
%         PreNod = 1;      % Previous Node set to depot
%         VecCnt = VecCnt + 1;% Increase vehicle count
%         AccWst = 0;         % Accumulated Waste
%         Y = sprintf('Vehicle: %d\n Route:',VecCnt);
%         disp(Y);
%         disp(PreNod);
%         while AccWst < vec_cap
%             if AccWst + dmnd(seq_inds(SeqCur)) > vec_cap
%                 TotExSp = TotExSp+abs(vec_cap-AccWst);
%                 break;
%             else
%                 TotDist = TotDist + dist_ij(PreNod, seq_inds(SeqCur));
%                 AccWst = AccWst + dmnd(seq_inds(SeqCur));
%                 PreNod = seq_inds(SeqCur);
%                 SeqCur = SeqCur + 1;
%                 disp(PreNod);
%             end
%             if SeqCur > ProbDim
%                 TotExSp = TotExSp+abs(vec_cap-AccWst);
%                 break;
%             end
%         end
%         TotDist = TotDist + dist_ij(PreNod, 1);
%     end
%     disp(seq_inds);
%     disp(TotExSp);
%     FitVal = TotDist*0.27*46.20;
%     disp(FitVal);
end
% %%
% % Get particular fitness
% function FitVal = PSO_GA_Printer(Chrom, ProbDim)
%     % load data
%     % load('vec_cap.mat');
%     
%     dmnd  = [0 7.9856 10.4665 6.6545 4.0003 9.6666];
%     load('distance_data.mat');
%     vec_cap = 12;
%     
%     % Obtain Sequence
%     seq_route = ObtainSequence(Chrom, ProbDim);
%     seq_inds = seq_route+1;
%     SeqCur = 1;
%     TotDist = 0;
%     PreNod = 1;
%     while SeqCur <= ProbDim
%         TotDist = TotDist + dist_ij(PreNod, seq_inds(SeqCur));
%         PreNod = seq_inds(SeqCur);
%         SeqCur = SeqCur+1;
%     end
%     TotDist = TotDist + dist_ij(PreNod, 1);
%     disp(seq_inds);
%     disp(TotDist);
%     FitVal = TotDist;
% end