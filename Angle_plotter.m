clear;
filename = 'object.csv';
A = importdata(filename);
rho = zeros(length(A), 1);
sigma_s = zeros(length(A), 1);
psi = zeros(length(A), 1);
tetta = zeros(length(A), 1);
t = zeros(length(A), 1);
for i = 1 : length(A)
    t(i) = A(i,1);
    psi(i) = A(i, 5);
    tetta(i) = A(i,4);
    rho(i) = A(i,2);
    sigma_s(i) = A(i,3);
end

subplot(4,1,1);
plot(t, rho(:));
title('Rho')

subplot(4,1,2); 
plot(t, tetta(:));
title('Tetta')

subplot(4,1,3); 
plot(t, psi(:));
title('Psi')

subplot(4,1,4); 
plot(t, sigma_s(:));
title('Sigma_s')