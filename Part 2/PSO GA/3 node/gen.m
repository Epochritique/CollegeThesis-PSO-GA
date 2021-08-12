ar = [];

for i = 1:4
   ar = [ar i]; 
end

permuts = perms(ar);

        DimMinMax = zeros(4, 2);
        DimMinMax(1, :) = [-1 10];
        DimMinMax(2, :) = [-1 10];
        DimMinMax(3, :) = [-1 10];
        DimMinMax(4, :) = [-1 10];

FitVal = PSO_GetFitValues(max(size(permuts)), permuts, 4, DimMinMax);
for i = 1:max(size(permuts))
    disp(permuts(i,:));
    disp(FitVal(i));
end
% 
% for RowNum = 1:max(size(permuts))
%     [FitVal(RowNum, 1), Errs(RowNum,1)] = PSO_GA_Eval(permuts(RowNum, :), 3);
%     disp(permuts(RowNum, :));
%     disp(FitVal(RowNum, 1));
%     disp(Errs(RowNum, 1));
% end