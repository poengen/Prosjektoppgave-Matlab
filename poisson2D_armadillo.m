% Plotting solution obtained from armadillo

%% Defining variables
n = 100; %internal grid points
h = 1/(n+1); %intervals

%% The Grid
x = zeros(n+2,1);
for i = 2:n+2
    x(i) = (i-1)*h;
end

%% Analytical Solution
uHandle = @(x,y) sin(pi*x)*sin(pi*y);
uAnal = zeros(n+2,n+2);
for i = 1:n+2
    for j = 1:n+2
        uAnal(i,j) = uHandle(x(i),x(j));
    end
end

%% read solution vector u from text file
Usol = importdata('u.dat');

%% Plotting
figure(1);
subplot(1,2,1);
mesh(x,x,uAnal);
title('Analytical solution');
subplot(1,2,2);
mesh(x,x,Usol)
title('Numerical solution');