function FitVal = PSO_GA_Eval(Chrom, ProbDim)
%     x = 0;
%     for i = 1:ProbDim
%         x = x + Chrom(i)*2^(ProbDim-i);
%     end
%     FitVal = x^2;
    %FitVal = 5.3578547*Chrom(3)^2 + 0.8356891*Chrom(1)*Chrom(5) + 37.293239*Chrom(1) - 40792.141;
    %FitVal = 0.6224*Chrom(1)*Chrom(3)*Chrom(4) + 1.7781*Chrom(2)*(Chrom(3)^2) + 3.1661*(Chrom(1)^2)*Chrom(4) + 19.84*(Chrom(1)^2)*Chrom(3);
    FitVal = 1.10471*(Chrom(1)^2)*Chrom(2) + 0.04811*Chrom(3)*Chrom(4)*(14+Chrom(2));
end
