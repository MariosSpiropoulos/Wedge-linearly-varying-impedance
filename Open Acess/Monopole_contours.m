clear all
clc 
close all 

addpath('aux/')



% Wedge-parameters
alpha = 10 *pi/180;
gamma= pi - alpha; 

% Discretisation in the polar coordinate system 
theta_start = - pi + alpha ;
theta_end =  pi - alpha ;
Nr = 100 ; 
Ntheta =250;
R = linspace(0,12,Nr);
Theta = linspace(theta_start,theta_end,Ntheta);
[r,theta] = meshgrid(R,Theta);
% Source's location 
r0 = 1 ;
theta0 = 0 ; 

% Wavenumber and frequency 
k0 = 1; 
c0=343;
w0 = c0.*k0;










% Quick functions to compute the summation coeffs. 
an=@(kn,z) -(1i * pi/2)* cbesselh(kn,z);


dn=@(kn,z) -(1i * pi/2)* cbesselj(kn,z);

A = @(kn,alpha,Z2) 1 ./ ( gamma +  gamma.*(  sinc(4*kn.*(gamma)./pi)));




Z1   = Inf;                              % Rigid bottom-face 
Z2 = input('Enter impedance Z2= ');      % Arbitrary top-face impedance




Nmodes= 100 ; %number of summation terms 








% Finding the eigenvalues
nvals=theta_equation(alpha,Z1,Z2,Nmodes,theta);

G=zeros(size(r));




for n=1:length(nvals)

    kn=nvals(n);

    
     An = cos((kn).*(theta0+gamma)).*A(kn,alpha,Z2);

    
    G_r_leq_r0= an(kn,k0.*r0).* cbesselj(kn,k0*r(1,:)).*An.*cos(kn.*(theta(:,1)+gamma));  % inner solution

    
     G_r_geq_r0= dn(kn,k0.*r0) .* cbesselh(kn,k0*r(1,:)) .*An.*cos(kn.*(theta(:,1)+gamma)); % outter solution

   
     G = G +(r<r0) .* G_r_leq_r0 + (r>r0) .* G_r_geq_r0 ;

    clear arg1 arg2 An kn  G_r_leq_r0 G_r_geq_r0


end




% Cartesian coordiante system 
X=r.*cos(theta);
Y=r.*sin(theta);
    

%% plotting 
%% Add Wedge (10-degree sector)
Length = 15;      % Length of the triangle's sides from the origin
half_theta_wedge = alpha; 
% Calculate the vertices of the triangle
origin = [0, 0]; % Vertex at the origin
v1 = [Length * cos(half_theta_wedge), Length * sin(half_theta_wedge)];  % Vertex at +x angle
v2 = [Length * cos(-half_theta_wedge), Length * sin(-half_theta_wedge)]; % Vertex at -x angle

% Combine vertices into arrays
x_vertices = [origin(1), -v1(1), -v2(1)];
y_vertices = [origin(2), -v1(2), -v2(2)];


figure('Units', 'inches', 'Position', [1, 1, 6,3]);  % [left, bottom, width, height]
pcolor(X,Y,real(G))
caxis([-1 1].*0.1)
shading interp 
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
axis equal tight;
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
colormap(bluewhitered)
hold on
fill(x_vertices, y_vertices,0.5.*[1 1 1]) % Green wedge outline
xlim([min(min(X)) max(max(X))])
set(gca,'FontSize',12)
hold on
xlim([-10 10])
ylim([-10 10])
    

     

