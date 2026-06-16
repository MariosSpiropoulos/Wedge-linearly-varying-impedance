clear all 
close all 
clc 

addpath('aux/')


alpha = 10 * pi/180;   % Wedge angle 
gamma = pi-alpha ; 
Z1 = Inf ;             % Bottom-face impedance 
Z2= 1+1i ;                % Top-face impedance 


theta = linspace(-gamma,gamma,100); % Vector of polar anlge 
N = length(theta);

[x,A]=chebdif(N,2);      % Chebyshev differentiation  

D2 = squeeze(A(:,:,2));

D = squeeze(A(:,:,1));

fac = pi-alpha;

x = x*fac;
D = D/fac;
D2 = D2/fac^2;

L = D2;          % Differential operator 
F = -eye(size(L)); 


% Applying boundary condtitions 
if(isfinite(Z2))
    L(1,:) = 0; L(1,:) = 1i*Z2*D(1,:); L(1,1) = L(1,1) + 1; F(1,:) = 0;
else
    L(1,:) = 0; L(1,:) = D(1,:); F(1,:) = 0;
end
if(isfinite(Z1))
    L(N,:) = 0; L(N,:) = -1i*Z1*D(N,:); L(N,N) = L(N,N) + 1; F(N,:) = 0;
else
    L(N,:) = 0; L(N,:) = -D(N,:); F(N,:) = 0;
end
% solving the eigenvalue problem: D2 \Theta = \lambda F \Theta (Appendix A)
% Notice that for the present problem \lambda = k_{\theta,n}^2

[V,lambda]=eig(L,F);
lambda = diag(lambda);
pos = find(isfinite(lambda));
lambda = lambda(pos); V = V(:,pos);

[~,ord]=sort(real(lambda));
lambda = lambda(ord);
V = V(:,ord);

kt = sqrt(lambda); % take the square root, keeping only the positive branch (convention)


%% Choosing the number of modes: Notice that Nmodes <= length(theta)

Nmodes = 50 ;         
nvals = kt(1:Nmodes); 


%% Trace the number of modes in the complex plane 
figure 
plot(real(nvals),imag(nvals),'ko','LineWidth',2)
grid on ; grid minor 
hold on 
set(gca,'FontSize',12);
set(groot,'defaulttextinterpreter','latex');
xlabel("$\Re \{k_{\theta,n} \}$",'Interpreter','latex')
ylabel("$\Im \{k_{\theta,n} \}$",'Interpreter','latex')
xlim([0.01 5.2])


