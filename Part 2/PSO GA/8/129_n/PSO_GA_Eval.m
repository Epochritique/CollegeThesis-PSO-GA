% Get particular fitness
function FitVal = PSO_GA_Eval(Chrom, ProbDim, NodeCnt, VecCnt, alp, DimMinMax)
    % load data
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
        Y(seq_inds(i),vec) = 1;
        if(i>1)
            if(seq_inds(i)==1)
                A(seq_inds(i),seq_inds(i+1),vec)=0;
            else
                A(seq_inds(i),seq_inds(i+1),vec)=A(seq_inds(i-1),seq_inds(i),vec)+node_stats(1, seq_inds(i));
            end
        end
    end
    
    % Obtain total distance
    total_cost = 0;
    for l=1:VecCnt
        total_cost = total_cost + sum(sum(X(:,:,l).*dist_ij(:,:)));
    end
    
    % Obtain total demand
    % W=sum(node_stats(1,:));
    
    % Obtain total collected
    % total_coll = sum(tanspose(sum(Y)).*node_stats(1, seq_inds(:)));
    
    violated_num=0;
    
    % Constraint 1
    for j=2:NodeCnt+1
        if sum(sum(X(:,j,:))) ~= 1
            violated_num = violated_num + 1;
        end
    end
    
    for i=1:VecCnt
        % Constraint 2
        if sum(sum(X(:, :, i))) ~= sum(sum(Y(:,i)))
            violated_num = violated_num+1;
        end
        % Constraint 3
        if sum(sum(X(1,:,i))) ~= 1
            violated_num = violated_num+1;
        end
        % Constraint 4
        if sum(sum(X(:,1,i))) ~= 1
            violated_num = violated_num+1;
        end
        % Constraint 5
        if sum(sum(A(1,:,i))) ~= 0
            violated_num = violated_num+1;
        end
        % Constraint 6
        for k=1:NodeCnt+1
            if sum(sum(A(k,:,i))) > vehicle_cap(i)
                violated_num = violated_num + 1;
            end
        end
        % Constraint 7
        for k=2:NodeCnt
            if sum(sum(A(k,:,i))) - sum(sum(A(:,k,i))) > node_stats(1,k)*sum(sum(X(:,k,i)))
                violated_num = violated_num + 1;
            end
        end
    end
    
    ErrNum = violated_num;
    disp(ErrNum);
    probs = 0;
    for i = 1:ProbDim
       if ~ (Chrom(i) >= DimMinMax(i, 1) && Chrom(i) <= DimMinMax(i, 2))
           probs=probs+1;% compute the number of boundary errors
       end
    end
    if(probs~=0)> 0 || ErrNum > 0
        FitVal = total_cost + (probs+ErrNum)*alp;
    end
end