



%% Two dimensional steady state poisson problem
% with homogeneous dirichlet boundary
n = 100; %internal grid points
uHandle = @(x,y) sin(pi*x)*sin(pi*y); %analytical solution
fHandle = @(x,y) 2*pi^2*sin(pi*x)*sin(pi*y); %load vector
g1=@(x) 0; g2=@(y) 0; g3=@(x) 0; g4=@(y) 0;

%% Two dimensional steady state poisson problem
% with non-homogeneous dirichlet boundary
n = 100;
uHandle = @(x,y) x;
fHandle = @(x,y) 0;
g1=@(x) x; g2=@(y) 1; g3=@(x) x; g4=@(y) 0;

%% Run calculations
[x, uNum, uAnal, error_L2g, error_L2] = poisson2D_steady(n, uHandle, fHandle, g1, g2, g3, g4);

%% TESTING

%% Error plot

% = 100:(doble hver gang):10000 %(minst fem punkter)

n = [50,100,200,400];
GERROR = zeros(1,length(n));
LERROR = zeros(1,length(n));
for i = 1:length(n)
    
    [x, uNum, uAnal, error_L2g, error_L2] = poisson2D_steady(n(i), uHandle, fHandle, g1, g2, g3, g4);
    GERROR(i) = error_L2g;
    LERROR(i) = error_L2;
    
end

%% Error plot
figure(1);
he = 1./(n+1);
plot(he,GERROR,'r',he,LERROR,'b');
legend('H1-error','L2-error')
xlabel('h: grid space');
ylabel('error');

%% Loglog-plot of the error
figure(2);
he = 1./(n+1);

G2 = log(GERROR(length(n))); G1 = log(GERROR(1));
L2 = log(LERROR(length(n))); L1 = log(LERROR(1));
x2 = log(he(length(n))); x1 = log(he(1));
kG = (G2-G1)/(x2-x1)
kL = (L2-L1)/(x2-x1)

loglog(he,GERROR,'-*r',he,LERROR,'-*b')
grid on
legend('H1-error, k= ','L2-error, k= ')
xlabel('log(h)');
ylabel('log(error)');



%% Plotting
figure(3);
subplot(1,2,1);
mesh(x,x,uAnal);
title('Analytical solution');
xlabel('x-axis');
ylabel('y-axis');
subplot(1,2,2);
mesh(x,x,uNum)
title('Numerical solution');
xlabel('x-axis');
ylabel('y-axis');