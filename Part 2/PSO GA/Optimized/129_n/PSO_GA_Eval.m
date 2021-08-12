function FitVal = PSO_GA_Eval(Chrom, ProbDim, NodeCnt, VecCnt, alp, DimMinMax)
    load('vec_cap.mat');
    load('nodes.mat');
    load('distance_data.mat');
    
    % Obtain Sequence
    seq_route = ObtainSequence(Chrom, ProbDim, NodeCnt);
    seq_inds = seq_route+1;
    
    X = zeros(NodeCnt+1,NodeCnt+1,VecCnt);
    Y = zeros(NodeCnt+1,VecCnt);
    A = zeros(NodeCnt+1,NodeCnt+1,VecCnt);
    
    % Obtain decision variable values
    vec = 0;
    for i=1:(ProbDim+1)
        if(seq_inds(i) == 1)
            vec = vec+1;
        end
        X(seq_inds(i),seq_inds(i+1),vec)=X(seq_inds(i),seq_inds(i+1),vec)+1;
        X(seq_inds(i+1),seq_inds(i),vec)=X(seq_inds(i+1),seq_inds(i),vec)+1;
        Y(seq_inds(i),vec) = 1;
    end
    vec = 1;
    for i=2:(ProbDim+1)
        if(seq_inds(i) == 0)
            vec = vec+1;
        end
        A(seq_inds(i),seq_inds(i+1),vec)=A(seq_inds(i-1),seq_inds(i),vec)+node_stats(1, seq_inds(i));
    end
    
    % Obtain total distance
    total_cost = 0;
    for l=1:VecCnt
        for i=1:NodeCnt+1
            for j=1:NodeCnt+1
                total_cost = total_cost + X(i,j,l)*dist_ij(i,j);
            end
        end
    end
    
    % Obtain total demand
    W=0;
    for l=1:VecCnt
        for i=1:NodeCnt+1
            W = W + node_stats(1, seq_inds(i));
        end
    end
    
    % Obtain total collected
    total_coll=0;
    for l=1:VecCnt
        for i=1:NodeCnt+1
            total_coll = total_coll + Y(i,l)*node_stats(1, seq_inds(i));
        end
    end
    
    % Obtain objective value
    total_cost = total_cost + (10000/1000)*nthroot(W-total_coll,3);
    FitVal=total_cost;
    
    violated_num=0;
    
    % Constraint 1
    for j=2:NodeCnt+1
        c=0;
        for l=1:VecCnt
            for i=2:NodeCnt+1
                c=c+X(i,j,l);
            end
        end
        if c~=1
            violated_num = violated_num + 1;
        end
    end
    
    % Constraint 2
    for l=1:VecCnt
        for j=2:NodeCnt+1
            for i=2:NodeCnt+1
                if X(i,j,l) ~= Y(i,l)
                    violated_num = violated_num+1;
                end
            end
        end
    end
    
    % Constraint 3
    for l=1:VecCnt
        for j=2:NodeCnt+1
            if X(1,j,l) ~= 1
                violated_num = violated_num+1;
            end
        end
    end
    
    % Constraint 4
    for l=1:VecCnt
        for i=2:NodeCnt+1
            if X(i,1,l) ~= 1
                violated_num = violated_num+1;
            end
        end
    end
    
    % Constraint 5
    for l=1:VecCnt
        for j=2:NodeCnt+1
            if A(1,j,l) ~= 1
                violated_num = violated_num+1;
            end
        end
    end
    
    % Constraint 6
    for l=1:VecCnt
        for i=1:NodeCnt+1
            for j=1:NodeCnt+1
                if A(i,j,l) > vehicle_cap(l)
                    violated_num = violated_num+1;
                end
            end
        end
    end
    
    % Constraint 7
%     for l=1:VecCnt
%         for i=1:NodeCnt+1
%             for h=1:NodeCnt+1
%                 for j=1:NodeCnt+1
%                     if A(j,h,l)-A(i,j,l) ~= Y(j,l)*node_stats(1,j)
%                         disp('hi');
%                         violated_num = violated_num+1;
%                     end
%                 end
%             end
%         end
%     end
    ErrNum = violated_num;
    
    probs = 0;
    for i = 1:ProbDim
       if ~ (Chrom(i) >= DimMinMax(i, 1) && Chrom(i) <= DimMinMax(i, 2))
           probs=probs+1;% compute the number of boundary errors
       end
    end
    if(probs~=0)> 0 || ErrNum > 0
        FitVal = FitVal+(probs+ErrNum)*alp;
    end
end