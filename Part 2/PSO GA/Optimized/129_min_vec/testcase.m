format long;
p = perms([0.1 0.2 0.3 0.4 0.5 0.6 0.7]);

FitVal = zeros(1,5040);
inf = zeros(1,5040);
for i = 1:5040
    [int(1,i), FitVal(1,i)] = PSO_GA_Printer(p(i, :), 7, 4);
end
Fit = [];
for i = 1:5040
    if int(1,i) ~= 1
        Fit=[Fit, FitVal(1,i)];
    end
end
disp(unique(FitVal));
% ni = max(Fit);
% disp(ni);
% ni = min(Fit);
% disp(ni);
% x=0;

%alp = ni/(12*0.2);
% for i = 1:5040
%     if abs(FitVal(i)-min(FitVal))<1e-6
%         %disp(i);
%         x = x+1;
%     end
% end