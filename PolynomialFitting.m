function [a,output,Syx,Sy] = PolynomialFitting(xi,yi,m)
c = length(xi); d = length(yi); if c~=d; disp('Vector dimensions must be the same');return;end
plot(xi,yi, 'ro' ,'MarkerSize',3, 'LineWidth', 3, 'MarkerFaceColor', 'r');grid minor;

totalXi = zeros(2*m+1,1);  

for i = 0:2*m
    totalXi(i+1) = sum(xi.^i); 
end

constantVectors = zeros(m+1,1);
mainMatrix = zeros(m+1,m+1);

for i = 0:m
    constantsVectors(i+1) = sum(yi.*(xi.^i)); 
    mainMatrix(i+1,:) = totalXi((i+1):(i+1+m)); 
end

RootofMatrix = mainMatrix\constantsVectors;
x = -max(xi)-1:0.1:max(xi)+1;
y = 0;
S = 0;

for i = 0:m
    a = RootofMatrix;
    y = y + a(i+1,:)*(x.^i);
    S = S + a(i+1)*(xi.^i);
    Sr = sum((yi - S).^2);
    Syx = sqrt(Sr/(c-(m+1)));
    arithmM = sum(yi)/c;
    St = sum((yi-arithmM).^2);
    Sy = sqrt(St/(c-1));
    rsquare = (St - Sr)/St;
    output = rsquare;
end

hold on; plot(x,y);
fprintf('%s \n%f \n%f \n%f \n%f \n','Root of Matrix: ',a)
fprintf('%s %f\n','Estimated Standard Error: ',Syx)
fprintf('%s %f\n','Estimated Standard Deviation: ',Sy)
fprintf('%s %f\n','Estimated Stability Coefficient: ',output)
end