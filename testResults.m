% testResults.m
function [] = testResults()
close all; clc; clear;
type = input("Enter the corresponding number for \nthe parameter you want to test: \n1 = control parameter (mu) \n2 = the number of Fourier modes (M) \n3 = viscosity parameter (nu) \n4 = Reynolds number (R) \n5 = the number of space discritization points (N) \n6 = the width of the channel (b) \n7 = time step (h) \n8 = maximum time value (tmax) \n9 = both control parameter (mu) and Fourier modes (M) \n\nparameter number = ");

if type == 1
    mu = (0.1:0.1:50)';
    muLength = length(mu);
    vErr1= zeros(muLength,1);
    for i = 1: muLength
        vErr1(i) =  testForParameters(1,mu(i));
    end
    semilogy(mu,vErr1)
    xlabel("$\bf{\mu}$ - control parameter", "FontSize", 12,'interpreter', 'latex')
    ylabel("$\bf{\|\tilde v(:, t) - v(:,t)\|_{L^2(0,b)}}$","FontSize", 12,'Interpreter', 'latex')
    print -deps epsFig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 2
    M = (1:1:40)';
    MLength = length(M);
    vErr2= zeros(MLength,1);
    for i = 1: MLength
        vErr2(i) =  testForParameters(2,M(i));
    end
    semilogy(M,vErr2)
    xlabel("$M$ - the number of Fourier modes","FontSize", 12,'Interpreter', 'latex')
    ylabel("$\bf{\|\tilde v(:, t) - v(:,t)\|_{L^2(0,b)}}$",'Interpreter', 'latex')
    print -deps epsFig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 3
    nu = (0:0.1:5)';
    nuLength = length(nu);
    vErr3= zeros(nuLength,1);
    for i = 1: nuLength
        vErr3(i) =  testForParameters(3,nu(i));
    end

    semilogy(nu,vErr3)
    xlabel("$\nu$ - viscosity", "FontSize", 12,'Interpreter', 'latex')
    ylabel("$\bf{\|\tilde v(:, t) - v(:,t)\|_{L^2(0,b)}}$", "FontSize", 12,'Interpreter', 'latex')
    print -deps epsFig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 4
    R = (10:10:1000)';
    RLength = length(R);
    vErr4= zeros(RLength,1);
    for i = 1: RLength
        vErr4(i) =  testForParameters(4,R(i));
    end
            

    semilogy(R,vErr4)
    xlabel("$R$ - Reynolds number", "FontSize", 12,'Interpreter', 'latex')
    ylabel("$\bf{\|\tilde v(:, t) - v(:,t)\|_{L^2(0,b)}}$", "FontSize", 12,'Interpreter', 'latex')
    print -deps epsFig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 5
    N = [40,100, 200, 400, 800, 1600, 2000, 3000, 4000]';
    NLength = length(N);
    vErr5= zeros(NLength,1);

    for i = 1: NLength
        vErr5(i) =  testForParameters(5,N(i));
    end

    semilogy(N,vErr5)
    xlabel("$N$ - the number of space discretization points", "FontSize", 12,'Interpreter', 'latex')
    ylabel("$\bf{\|\tilde v(:, t) - v(:,t)\|_{L^2(0,b)}}$", "FontSize", 12,'Interpreter', 'latex')
    print -deps epsFig


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 6
    b = (0.5:0.1:20)';
    bLength = length(b);
    vErr6= zeros(bLength,1);

    for i = 1: bLength
        vErr6(i) =  testForParameters(6,b(i));
    end
    semilogy(b,vErr6)
    xlabel("$b$ - width of the channel", "FontSize", 12,'Interpreter', 'latex')
    ylabel("$\bf{\|\tilde v(:, t) - v(:,t)\|_{L^2(0,b)}}$", "FontSize", 12,'Interpreter', 'latex')
    print -deps epsFig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 7
    h = [10^(-5), 10^(-4), 10^(-3), 0.005 ,10^(-2), 0.05, 0.1, 0.11,0.12,0.13,0.14,0.15,0.2]';
    hLength = length(h);
    vErr7= zeros(hLength,1);
    for i = 1: hLength
        vErr7(i) =  testForParameters(7,h(i));
    end

    loglog(h,vErr7)
    xlabel("$h$ - time step", "FontSize", 12,'Interpreter', 'latex')
    ylabel("$\bf{\|\tilde v(:, t) - v(:,t)\|_{L^2(0,b)}}$", "FontSize", 12,'Interpreter', 'latex')
    print -deps epsFig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 8
    tmax = (1:1:100)';
    tmaxLength = length(tmax);
    vErr8= zeros(tmaxLength,1);
    for i = 1: tmaxLength
        vErr8(i) =  testForParameters(8,tmax(i));
    end
    semilogy(tmax,vErr8)
    xlabel("$t_{max}$ - maximum time value", "FontSize", 12,'Interpreter', 'latex')
    ylabel("$\bf{\|\tilde v(:, t) - v(:,t)\|_{L^2(0,b)}}$", "FontSize", 12,'Interpreter', 'latex')
    print -deps epsFig
    
elseif type == 9
     mu = (0.25:0.25:10)';
     M = (1:1:30)';
     muLength = length(mu);
     MLength = length(M);
     vErr9 = zeros(muLength, MLength);
     for j = 1:MLength
         for i = 1: muLength
            vErr9(i,j) =  testForControlParameters(mu(i),M(j));
         end
     end
    [minVerr9,index] = min(vErr9);
    [minimum,index2] = min(minVerr9);
    i = index(index2);
    j = index2;
    min_arg_mu = mu(i)
    min_arg_M = M(j)
    minimum
end
            
end


