function FitVal = GA_Eval(Chrom, ProbDim)
    x = 0;
    for i = 1:ProbDim
        x = x + Chrom(i)*2^(ProbDim-i);
    end
    FitVal = x^2;
end
