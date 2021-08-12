function Chrom = GA_Mutation(Chrom, DimMinMax, ProbDim)
    % Mutation
    MutPos = randi(ProbDim);
    % Chrom(1, MutPos) = (DimMinMax(MutPos, 2)- DimMinMax(MutPos, 1))*rand(1,'double') + DimMinMax(MutPos, 1);
    Chrom(1, MutPos) = randi(2)-1;
end
