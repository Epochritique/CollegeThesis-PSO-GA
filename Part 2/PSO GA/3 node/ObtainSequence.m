function seq_route = ObtainSequence(Chrom, ProbDim)
    % Obtain Sequence
    sorted_chrom = sort(Chrom);
    seq_route = zeros(1,ProbDim);
    for i=1:ProbDim
        for j=1:ProbDim
            if(sorted_chrom(i) == Chrom(j))
                if j>3 % 9, 10 must be 0
                    seq_route(i) = 0;
                else
                    seq_route(i) = j;
                end
            end
        end
    end
    seq_route = [0 seq_route 0];
end