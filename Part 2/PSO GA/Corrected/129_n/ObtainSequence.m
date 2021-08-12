function seq_route = ObtainSequence(Chrom, ProbDim, NodeCnt)
    % Obtain Sequence
    Y = zeros(1, ProbDim);
    sorted_chrom = sort(Chrom);
    seq_route = zeros(1,ProbDim);
    for i=1:ProbDim
        for j=1:ProbDim
            if(sorted_chrom(i) == Chrom(j) && Y(j)~=1)
                if j>NodeCnt
                    seq_route(i) = 0;
                else
                    seq_route(i) = j;
                end
                Y(j) = 1;
                break;
            end
        end
    end
    seq_route = [0 seq_route 0 1];
    
    seq_route_nocon = [];
    % Consecutive Number Remover
    for i=1:length(seq_route)-1
        if ~(seq_route(i) == 0 && seq_route(i+1) == 0)
            seq_route_nocon = [seq_route_nocon seq_route(i)];
        end
    end
    seq_route = seq_route_nocon;
end