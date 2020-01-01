function [] =  testForSolutionOfOriginalBurger(R,nu,N, b, h)
close all,
if nargin < 1
    R = 50;     % Reynold's number
    nu = 1;   % viscousity term
    N = 400;    % number of space discretization
    b = 1.0;    % width of the channel    
    h = 0.01;   % time step
end
P = R*(nu)^2/(b^2); % pressure
dx = b/N;           % space step 
x = (0:dx:b)';
U = 1.0;
v = sin((pi/b)*x);
u = [v;U];          % Solution of uncontrolled problem

% L is the linear part of the system
L = toeplitz([-2; 1; zeros(N-1,1)]);
L = nu*(1.0/dx)^2*L; 
L(N+2, N+2) = -nu/b;
c = 32; 
r = 15*exp(1i*pi*((1:c)-.5)/c); % roots of unity
A = h*L;
E = expm(A); E2 = expm(A/2);
I = eye(N+2); Z = zeros(N+2);
f1= Z; f2 = Z; f3 = Z; Q = Z;
for j = 1:c
    z = r(j);
    zIA= inv(z*I-A);
    Q = Q + h*zIA*(exp(z/2)-1);
    f1 = f1 + h*zIA*(-4-z+ exp(z)*(4-3*z + z^2))/z^2;
    f2 = f2 + h*zIA*(2 + z+ exp(z)*(z-2))/z^2;
    f3 = f3 + h*zIA*(-4-3*z -z^2+exp(z)*(4-z))/z^2;
end
f1 = real(f1/c); f2 = real(f2/c); 
f3 = real(f3/c); Q = real(Q/c);
% D: forward difference formula matrix for first derivative 
D = zeros(N+1, N+1);
for i = 1: N-2
    D(i,(i:i+3)) = [-7 6 3 -2];
end
for j = 1: 3
    D(N-2 + j, N-2-3 + j: N-2 + j)= [-2 9 -18 11];
end
D = (1.0/(6*dx))*D;


% % Composite Simpson's rule integration matrix 
S = ones(1, N+1);    
S(1,2:N-1) = repmat([4 2], 1, N/2 -1);
S(1, N) = 4;
S = ((2.0*dx)/6.0)*S;

% Time-stepping loop
uu = u; tt = 0;
uuexact = u; UUerror = 0;
%vverror = 0; 
tmax = 10.0; nmax = tmax/h;
for n = 1:nmax
    t = n*h;    
    % Solving the uncontroled problem
    
    Nu = [(1/b)*(u(N+2)*u(1:N+1)) - D*(u(1:N+1).^2) ...
        + exp(-t)*sin((pi/b)*x).*((2*pi/b)*cos((pi/b).*x) - exp(-t)/b - 1 + nu*pi^2/b^2) ; ...
        (P/b) - (1.0/b^2)*S*(u(1:N+1).^2) + (1/(2*b))*exp(-2*t) + (nu/(b^2)-1)*exp(-t) - P/b];    
    alpha = E2*u + Q*Nu;    
    Na = [(1/b)*alpha(N+2)*alpha(1:N+1) - D*(alpha(1:N+1).^2)...
        + exp(-t)*sin((pi/b)*x).*((2*pi/b)*cos((pi/b).*x) - exp(-t)/b - 1 + nu*pi^2/b^2); ...
        (P/b) - (1.0/b^2)*S*(alpha(1:N+1).^2) + (1/(2*b))*exp(-2*t) + (nu/(b^2)-1)*exp(-t) - P/b];
    beta = E2*u + Q*Na;
    Nb = [(1/b)*beta(N+2)*beta(1:N+1) - D*(beta(1:N+1).^2)...
        + exp(-t)*sin((pi/b)*x).*((2*pi/b)*cos((pi/b).*x) - exp(-t)/b - 1 + nu*pi^2/b^2); ...
        (P/b) - (1.0/b^2)*S*(beta(1:N+1).^2) + (1/(2*b))*exp(-2*t) + (nu/(b^2)-1)*exp(-t) - P/b];
    gamma = E2*alpha + Q*(2*Nb-Nu);
    Nc = [(1/b)*gamma(N+2)*gamma(1:N+1) - D*(gamma(1:N+1).^2)...
        + exp(-t)*sin((pi/b)*x).*((2*pi/b)*cos((pi/b).*x) - exp(-t)/b - 1 + nu*pi^2/b^2); ...
        (P/b) - (1.0/b^2)*S*(gamma(1:N+1).^2) + (1/(2*b))*exp(-2*t) + (nu/(b^2)-1)*exp(-t) - P/b];
    
    u = E*u + f1*Nu + 2*f2*(Na + Nb) + f3*Nc;
    u(1) = 0;
    u(N+1) = 0; % Homogeneuous Dirichlet BC's
    uu = [uu, u]; 
    tt = [tt, t];
    
    % Exact solution
    vexact = exp(-t).*sin(pi*x);    Uexact = exp(-t);
    uexact = [exp(-t).*sin(pi*x);exp(-t)];
    uuexact = [uuexact, uexact];
    
    % Error
    Uerror = abs(u(N+2)-uexact(N+2))/abs(u(N+2));
    UUerror = [UUerror; Uerror];
  
end

uu;
UUerror;
v = uu(1:N+1,:);
U = uu(N+2,:);
vexact = uuexact(1:N+1,:);


Uexact = uuexact(N+2,:);
% % Plot results
plot(tt,U);
hold on 
plot(tt, Uexact);
xlabel time, ylabel("$U(t)$",'interpreter','latex')
legend("Numerical Sol.", "Exact Sol.")
print -deps epsFig
figure
plot(tt,UUerror)
xlabel time, ylabel Uerror
title("Relative error between numerical and exact solution for U")
print -deps epsFig
figure
mesh(x,tt,v')
xlabel space, ylabel time, zlabel("v(x,t) - numerical")
title("v(x,t) - the numerical solution")
print -deps epsFig
figure
mesh(x,tt,vexact')
xlabel space, ylabel time, zlabel("v(x,t) - exact")
title("v(x,t) - the exact solution")
print -deps epsFig
end

