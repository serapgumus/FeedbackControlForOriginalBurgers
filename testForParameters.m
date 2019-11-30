% testForParameters.m
function [vErr] =  testForParameters(parameterType, input)
if parameterType ==0 && input == 0
    mu =1;
    M = 2;
    nu = 0.25;
    R = 260;
    N = 200;
    b = 1;
    h = 0.01;
    tmax=10;
elseif parameterType == 1
    mu = input;
    M = 15;
    nu = 0.25;        
    R = 50;     
    N = 200;    
    b = 1.0;       
    h = 0.01;  
    tmax = 10;  
elseif parameterType == 2
    mu = 3.5;
    M = input;
    nu = 0.25;        
    R = 50;     
    N = 200;    
    b = 1.0;       
    h = 0.01;  
    tmax = 10; 
elseif  parameterType == 3
    mu = 3.5;
    M = 15;
    nu = input;        
    R = 50;     
    N = 200;    
    b = 1.0;       
    h = 0.01;  
    tmax = 10; 
elseif parameterType == 4
    mu = 3.5;
    M = 15;
    nu = 0.25;        
    R = input;     
    N = 200;    
    b = 1.0;       
    h = 0.01;  
    tmax = 10;   
elseif parameterType == 5
    mu = 3.5;
    M = 15;
    nu = 0.25;        
    R = 50;     
    N = input;    
    b = 1.0;       
    h = 0.01;  
    tmax = 10;   
 elseif parameterType == 6
    mu = 3.5;
    M = 15;
    nu = 0.25;        
    R = 50;     
    N = 200;    
    b = input;       
    h = 0.01;  
    tmax = 10; 
elseif parameterType == 7
    mu = 3.5;
    M = 15;
    nu = 0.25;        
    R = 50;     
    N = 200;    
    b = 1.0;       
    h = input;  
    tmax = 10; 
elseif parameterType == 8
    mu = 3.5;
    M = 15;
    nu = 0.25;        
    R = 50;     
    N = 200;    
    b = 1.0;       
    h = 0.01;  
    tmax = input; 
end
    

P = R*(nu)^2/(b^2); % pressure
dx = b/N;           % space step 
x = (0:dx:b)';
U = 0.2;
Uti = 0.3;
v = 0.4*sin(pi/b*x);
vti = 0.5*sin(pi/b*x);
u = [v;U];          % Solution of uncontrolled problem
uti = [vti; Uti];   % Solution of controlled problem

% L is the linear part of the system
L = toeplitz([-2; 1; zeros(N-1,1)]); % Second derivative by central dif.
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

% Feedback control operator
w = zeros(N+1, M); 
k = 1:M;
w(:,k) = sqrt(2/b)*sin(((k*pi)/b).*x); % eigenfunctions of -d_xx ...
% (second derivative w.r.t. x) under Dirichlet BC

% Time-stepping loop
uu = u; tt = 0;
uuti = uti;
nmax = tmax/h;
for n = 1:nmax
    t = n*h;    
    % Solving the uncontroled problem
    Nu = [(1/b)*(u(N+2)*u(1:N+1)) - D*(u(1:N+1).^2); ...
        (P/b) - (1.0/b^2)*S*(u(1:N+1).^2)];    
    alpha = E2*u + Q*Nu;    
    Na = [(1/b)*alpha(N+2).*alpha(1:N+1) - D*(alpha(1:N+1).^2); ...
        (P/b) - (1.0/b^2)*S*(alpha(1:N+1).^2)];
    beta = E2*u + Q*Na;
    Nb = [(1/b)*beta(N+2)*beta(1:N+1) - D*(beta(1:N+1).^2); ...
        (P/b) - (1.0/b^2)*S*(beta(1:N+1).^2)];
    gamma = E2*alpha + Q*(2*Nb-Nu);
    Nc = [(1/b)*gamma(N+2)*gamma(1:N+1) - D*(gamma(1:N+1).^2); ...
        (P/b) - (1.0/b^2)*S*(gamma(1:N+1).^2)];
    u = E*u + f1*Nu + 2*f2*(Na + Nb) + f3*Nc;
    
    % Solving the feedback control problem 
    Nuti = [(1/b)*(uti(N+2)*uti(1:N+1)) - D*(uti(1:N+1).^2) ...
        - mu*w*(S*((uti(1:N+1)-u(1:N+1)).*w))';...
        (P/b) - (1.0/b^2)*S*(uti(1:N+1).^2)];
    alphati = E2*uti + Q*Nuti;
    Nati = [(1/b)*alphati(N+2)*alphati(1:N+1) - D*(alphati(1:N+1).^2)...
        - mu*w*(S*((alphati(1:N+1)-alpha(1:N+1)).*w))'; ...
        (P/b) - (1.0/b^2)*S*(alphati(1:N+1).^2)];
    betati = E2*uti + Q*Nati;
    Nbti = [(1/b)*betati(N+2)*betati(1:N+1) - D*(betati(1:N+1).^2)...
        - mu*w*(S*((betati(1:N+1)-beta(1:N+1)).*w))'; ...
        (P/b) - (1.0/b^2)*S*(betati(1:N+1).^2)];
    gammati = E2*alphati + Q*(2*Nbti-Nuti);
    Ncti = [(1/b)*gammati(N+2)*gammati(1:N+1) - D*(gammati(1:N+1).^2)...
        - mu*w*(S*((gammati(1:N+1)-gamma(1:N+1)).*w))'; ...
        (P/b) - (1.0/b^2)*S*(gammati(1:N+1).^2)];
    uti = E*uti + f1*Nuti + 2*f2*(Nati + Nbti) + f3*Ncti;
    
    u(1) = 0; u(N+1) = 0; % Homogeneuous Dirichlet BC's
    uti(1) = 0; uti(N+1) = 0; % Homogeneuous Dirichlet BC's
    uu = [uu, u]; 
    uuti = [uuti, uti];
    vDifference = sqrt(S*((uti(1:N+1)-u(1:N+1)).^2)); %L2-norm of (vti-v)
    UDifference = abs(uti(N+2)-u(N+2));  % |Uti - U|   
    utiNorm = sqrt(S*(uti(1:N+1)).^2);
    uNorm = sqrt(S*(u(1:N+1)).^2);
    tt = [tt, t];

end
vErr = vDifference/utiNorm; % Relative error.
end
