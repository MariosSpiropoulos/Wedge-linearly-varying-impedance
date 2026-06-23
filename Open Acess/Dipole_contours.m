clear 
clc 
close all 
addpath('aux/')

%% wedge angle 
alpha = 10*pi/180 ; 
gamma= pi - alpha; 


%% polar coordinates
theta_start = - gamma ;
theta_end =  gamma ;
R = linspace(0,10,200);
Theta = linspace(theta_start,theta_end,200); 
[r,theta] = meshgrid(R,Theta);




%% Plotting the wedge 
%% Add Wedge (5-degree sector)
 % Parameters
length = 15;      % Length of the triangle's sides from the origin
half_theta_wedge = alpha; 
% Calculate the vertices of the triangle
origin = [0, 0]; % Vertex at the origin
v1 = [length * cos(half_theta_wedge), length * sin(half_theta_wedge)];  % Vertex at +x angle
v2 = [length * cos(-half_theta_wedge), length * sin(-half_theta_wedge)]; % Vertex at -x angle

% Combine vertices into arrays
x_vertices = [origin(1), -v1(1), -v2(1)];
y_vertices = [origin(2), -v1(2), -v2(2)];

% Combine vertices into arrays
x_vertices = [origin(1), -v1(1), -v2(1)];
y_vertices = [origin(2), -v1(2), -v2(2)];

%% impedance faces  
Z2=0;
Z1 = Inf; 


%% source location 
r0 = 1 ;
theta0 = 0 ; 

%% Acoustic wavelength, frequency
k0 = 1; 
c0=343;
w0 = c0.*k0;



Nmodes= 90 ; %  number of eigenvalues used 




%% Computing the dipole acoustic field

%% --------------------------------------------------------------------------------------------------------------------





an=@(kn,z) -(1i * pi/2)* cbesselh(kn,z);


dn=@(kn,z) -(1i * pi/2)* cbesselj(kn,z);


%% Finding the eigenvalues
nvals = theta_equation(alpha,Z1,Z2,Nmodes,Theta);   % compute eigenvalues

%% Initialize Green function 

dGdr0_1=zeros(size(r));
dGdr0_2=dGdr0_1;
dGdtheta0_1= dGdr0_1;

dGdtheta0_2= dGdr0_1;


for n=1:Nmodes

    kn=(nvals(n));
    

    An = cos(kn.*(theta(:,1)+gamma)) ./ ( gamma  +gamma.*  sinc(4*kn.*(gamma )./(pi))).';

    dH_dr0 = (k0/2) .*  (cbesselh(kn-1,k0.*r0) -  cbesselh(kn+1,k0.*r0)); 
    dJ_dr0 = (k0/2) .*  (cbesselj(kn-1,k0.*r0) -  cbesselj(kn+1,k0.*r0));


    dGdr0_1 = dGdr0_1-(1i * pi /2) .* cbesselj(kn,k0.*r(1,:)) .* dH_dr0 .*  An .*  cos((kn).*(theta0+gamma));
    dGdr0_2 =  dGdr0_2-(1i * pi /2) .* cbesselh(kn,k0.*r(1,:)) .* dJ_dr0 .*  An .*  cos((kn).*(theta0+gamma));

    dGdtheta0_1 = dGdtheta0_1+ (kn).* (1i * pi /2) .* cbesselj(kn,k0.*r(1,:)) .* cbesselh(kn,k0.*r0).*  An .*  sin((kn).*(theta0+gamma));
    dGdtheta0_2 =  dGdtheta0_2+(kn).*(1i * pi /2) .* cbesselh(kn,k0.*r(1,:)) .* cbesselj(kn,k0.*k0.*r0) .*  An .*  sin((kn).*(theta0+gamma));

    clear arg1 arg2 An kn


end


    dGdx_in = (r<r0).* (dGdr0_1 .*cos(theta0) - dGdtheta0_1.*sin(theta0)./r0);
    dGdx_out = (r>r0).* (dGdr0_2.*cos(theta0) - dGdtheta0_2.*sin(theta0)./r0);
    
   
    dGdy_in = (r<r0).* (dGdr0_1.*sin(theta0)  +  dGdtheta0_1.*cos(theta0)./r0);
    dGdy_out = (r>r0).* (dGdr0_2.*sin(theta0) +  dGdtheta0_2.*cos(theta0)./r0);

    dGdx   = dGdx_in + dGdx_out ;

    dGdy   = dGdy_in + dGdy_out ;



%% ---------------------------------------------------------------------------------------------------------------------



X=r.*cos(theta);
Y=r.*sin(theta);

figure('Units', 'inches', 'Position', [1, 1, 4, 8]);  % [left, bottom, width, height]
subplot(2,1,1)
pcolor(X,Y,real(dGdx))
 caxis([-1 1].*0.1)
if isinf(Z2)
 title("$\Re \left[\frac{\partial \hat G}{\partial x_0} \right]$",'Interpreter','latex')
end
 shading interp 
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');

colormap(bluewhitered)
hold on
 
fill(x_vertices, y_vertices,0.5.*[1 1 1]) % Green wedge outline
xlim([min(min(X)) max(max(X))])
set(gca,'FontSize',10)
hold on
axis equal

subplot(2,1,2)
pcolor(X,Y,real(dGdy))
caxis([-1 1].*0.1)
if isinf(Z2)
 title("$\Re \left[ \frac{\partial \hat G}{\partial y_0} \right]$",'Interpreter','latex')
end
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');

set(gca,'FontSize',15)
shading interp 
%     colormap(bluewhitered)
axis equal
set(gca,'TickLabelInterpreter','latex')
if abs(imag(Z2) ) >0
    impedance='complex';
    % colorbar
else 
    impedance=num2str(Z2);
end
hold on
fill(x_vertices, y_vertices,0.5.*[1 1 1]) % Green wedge outline
xlim([min(min(X)) max(max(X))])
set(gca,'FontSize',10)
hold on
axis equal
