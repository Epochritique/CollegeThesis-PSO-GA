function seq_route = ObtainSequence(Chrom, ProbDim)
    % Obtain Sequence
    sorted_chrom = sort(Chrom);
    seq_route = zeros(1,ProbDim);
    for i=1:ProbDim
        for j=1:ProbDim
            if(sorted_chrom(i) == Chrom(j))
                seq_route(i) = j;
            end
        end
    end
end