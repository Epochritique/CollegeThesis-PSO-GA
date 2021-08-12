% Function for GA Mutation mechanism
function Chrom = GA_Mutation(Chrom, DimMinMax, ProbDim)
    % Select Random Gene to Mutate
    MutPos = randi(ProbDim);    
    % Non-Binary Mutation
    Chrom(1, MutPos) = (DimMinMax(MutPos, 2)- DimMinMax(MutPos, 1))*rand(1,'double') + DimMinMax(MutPos, 1);
end
