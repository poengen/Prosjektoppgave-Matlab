% Solving the time dependent 2D poisson problem with given boundary
% conditions

%% Image initial condition
PIC = imread('init_pic2.jpg');
PIC = rgb2gray(PIC);
PIC = double(PIC);
PIC = 1/256*PIC;
n = length(PIC);
PIC = ones(n,n)-PIC;

%% Defining variables
%n = 30; %internal grid points
hx = 1/(n+1); %intervals

ht = 0.001; %time step
k = 500; %number of time steps
T = ht*k; %total time

%% The Grid
x = zeros(n+2,1); %includes the boundary
for i = 2:n+2
    x(i) = (i-1)*hx;
end

%% Matrices
I = eye(n^2);
e1 = ones(n,1);
Adiag = spdiags([-e1 4*e1 -e1], -1:1, n, n);
Asub = diag(-e1);
Asup = diag(-e1);

A = -1/hx^2*blktridiag(Adiag,Asub,Asup,n);
G = I-ht*A;

%% Initial conditions
p = zeros(n*n,k); %solution vector
p(:,1) = reshape(PIC,n*n,1); %initial condition

%% Apply boundary conditions

pSol = zeros(n+2,n+2,k);

%boundary functions
f1 = @(x,t) 0.5*sin(t); %y=0
f2 = @(y,t) 0.5*(1-y)*sin(t); %x=1
f3 = @(x,t) (1-x)*sin(t); %y=1
f4 = @(y,t) 0.5*(1+y)*sin(t); %x=0

%% Solve the linear system Au=f

for i = 1:k-1
    f = zeros(n*n,1);
    for j=1:n %boundary independent of time
        f(j) = f(j) + f1(x(j+1),ht*i);
        f(n^2-n+j) = f(n^2-n+j) + f3(x(j+1),i*ht);
        f(n*j) = f(n*j) + f2(x(j+1),i*ht);
        f(n*j-n+1) = f(n*j-n+1) + f4(x(j+1),i*ht);
        
        pSol(j+1,1,i+1) = f1(x(j+1),(i+1)*ht);
        pSol(n+2,j+1,i+1) = f2(x(j+1),(i+1)*ht);
        pSol(j+1,n+2,i+1) = f3(x(j+1),(i+1)*ht);
        pSol(1,j+1,i+1) = f4(x(j+1),(i+1)*ht);
    end
    %corners
    pSol(1,1,i+1) = f1(0,(i+1)*ht);
    pSol(n+2,1,i+1) = f1(1,(i+1)*ht);
    pSol(n+2,n+2,i+1) = f3(1,(i+1)*ht);
    pSol(1,n+2,i+1) = f3(0,(i+1)*ht);
    
    p(:,i+1) = G\(p(:,i)+ht/hx^2*f);
     
end

%% Composing the solution matrix Usol
for i = 1:k
    pSol(2:n+1,2:n+1,i) = reshape(p(:,i),n,n);
end
%% Plotting

for i=1:k
    mesh(x,x,pSol(:,:,i));
    xlabel('x');
    ylabel('y');
    axis([0 1 0 1 -0.6 1]);
    pause(0.5);
end