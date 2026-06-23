clear 
clc
close all 

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




%% compute Quadrupole radiation field 
%% Finding the eigenvalues
nvals = theta_equation(alpha,Z1,Z2,Nmodes,theta);

%% initialisation of summation-terms
Grr_in = zeros(size(r));


InTerm_xx_1  = Grr_in;
InTerm_xx_2  = Grr_in;
InTerm_xx_3  = Grr_in;

OutTerm_xx_1  = Grr_in;
OutTerm_xx_2  = Grr_in;
OutTerm_xx_3  = Grr_in;

InTerm_yy_1  = Grr_in;
InTerm_yy_2  = Grr_in;
InTerm_yy_3  = Grr_in;

OutTerm_yy_1  = Grr_in;
OutTerm_yy_2  = Grr_in;
OutTerm_yy_3  = Grr_in;


for n=1:Nmodes

    kn=(nvals(n));
    
    %% Terms used for the inner region
    Xi = cbesselj(kn,k0.*r(1,:)) .* cos(kn.*(theta(:,1)+gamma)) ./ ( gamma  +   gamma.* sinc(4*kn.*gamma./(pi)));
  
    
    H_cal = k0/2.* (   k0/2*( cbesselh(kn-2,k0.*r0)  -  cbesselh(kn,k0.*r0)   )-...
        k0/2*(  cbesselh(kn,k0.*r0)-cbesselh(kn+2,k0.*r0))    );
    
    dHdr0 = k0/2*( cbesselh(kn-1,k0.*r0)  -  cbesselh(kn+1,k0.*r0) )  ;
    
    %% Terms used for the inner region
    
    Phi = cbesselh(kn,k0.*r(1,:)) .* cos(kn.*(theta(:,1)+gamma)) ./ ( gamma  +  gamma.*  sinc(4*kn.*gamma./(pi)));
  
    
    J_cal = k0/2.* (   k0/2*( cbesselj(kn-2,k0.*r0)  -  cbesselj(kn,k0.*r0)   )-...
        k0/2*(  cbesselj(kn,k0.*r0)-cbesselj(kn+2,k0.*r0))    );
    
    dJdr0 = k0/2*( cbesselj(kn-1,k0.*r0)  -  cbesselj(kn+1,k0.*r0) )  ;
    %% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    %%                                 Gxx
    %% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    %% Gxx for r<r0
    
        InTerm_xx_1 = InTerm_xx_1+...
            1i.*pi.*sin(theta0)^2/(2*r0^2).*(kn)^2.* Xi .* cbesselh(kn,k0*r0) .* cos((gamma+theta0).*(kn));
        
        InTerm_xx_2= InTerm_xx_2-...
            1i.*pi.*cos(theta0)^2/(2).*  Xi .* H_cal .* cos((gamma+theta0).*(kn));
    
        InTerm_xx_3= InTerm_xx_3+...
            1i.*pi.*sin(theta0).*cos(theta0)/(r0).* (kn).* Xi .* dHdr0 .* sin((gamma+theta0).*(kn));
     
        %% Gxx for r>r0
     OutTerm_xx_1 = OutTerm_xx_1+...
            1i.*pi.*sin(theta0)^2/(2*r0^2).*(kn)^2.* Phi .* cbesselj(kn,k0*r0) .* cos((gamma+theta0).*(kn));
        
        OutTerm_xx_2= OutTerm_xx_2-...
            1i.*pi.*cos(theta0)^2/(2).*  Phi .* J_cal .* cos((gamma+theta0).*(kn));
    
        OutTerm_xx_3= OutTerm_xx_3+...
            1i.*pi.*sin(theta0).*cos(theta0)/(r0).* (kn).* Phi .* dJdr0 .* sin((gamma+theta0).*(kn));
     %% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    %%                                 Gyy
    %% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    %% Gyy for r<r0
    
        InTerm_yy_1 = InTerm_yy_1+...
            1i.*pi.*cos(theta0)^2/(2*r0^2).*(kn)^2.* Xi .* cbesselh(kn,k0*r0) .* cos((gamma+theta0).*(kn));
        
        InTerm_yy_2= InTerm_yy_2-...
            1i.*pi.*sin(theta0)^2/(2).*  Xi .* H_cal .* cos((gamma+theta0).*(kn));
    
        InTerm_yy_3= InTerm_yy_3+...
            1i.*pi.*sin(theta0).*cos(theta0)/(r0).* (kn).* Xi .* dHdr0 .* sin((gamma+theta0).*(kn));
     
        %% Gyy for r>r0
        OutTerm_yy_1 = OutTerm_yy_1+...
            1i.*pi.*cos(theta0)^2/(2*r0^2).*(kn)^2.* Phi .* cbesselj(kn,k0*r0) .* cos((gamma+theta0).*(kn));
        
        OutTerm_yy_2= OutTerm_yy_2-...
            1i.*pi.*sin(theta0)^2/(2).*  Phi .* J_cal .* cos((gamma+theta0).*(kn));
    
        OutTerm_yy_3= OutTerm_yy_3+...
            1i.*pi.*sin(theta0).*cos(theta0)/(r0).* (kn).* Phi .* dJdr0 .* sin((gamma+theta0).*(kn));
    


end

Gxx = (r<r0).*(InTerm_xx_1+InTerm_xx_2+InTerm_xx_3)+(r>r0).*(OutTerm_xx_1+OutTerm_xx_2+OutTerm_xx_3);

Gyy = (r<r0).*(InTerm_yy_1+InTerm_yy_2+InTerm_yy_3)+(r>r0).*(OutTerm_yy_1+OutTerm_yy_2+OutTerm_yy_3);





% Convert to Cartesian for contour plotting
X = r.*cos(theta);
Y = r.*sin(theta);



% Plot contour
figure('Units', 'inches', 'Position', [1, 1, 9, 3]);
subplot(1,2,1)
pcolor(X, Y, real(Gxx));
shading interp
caxis(0.1.*[-1 1])
% colormap(bluewhitered)
  colormap(bluewhitered)

hold on
fill(x_vertices, y_vertices,0.5.*[1 1 1]) % Green wedge outline
xlim([min(min(X)) max(max(X))])
set(gca,'FontSize',10)
hold on
axis equal tight;

if isinf(Z2)
title("$\Re  \left[  \frac{\partial^2 \hat G}{\partial x_0^2 } \right] $",'Interpreter','latex');
end
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
set(gca,'FontSize',13)
set(gca,'TickLabelInterpreter','latex')
axis equal
subplot(1,2,2)
pcolor(X, Y, real(Gyy));
shading interp
caxis(0.1.*[-1 1])
%colormap(bluewhitered)
   colormap(bluewhitered)
set(gca,'FontSize',10)

axis equal tight;
if isinf(Z2)
title("$\Re  \left[  \frac{\partial^2 \hat G}{\partial y_0^2 } \right] $",'Interpreter','latex');
end
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
set(gca,'FontSize',13)
set(gca,'TickLabelInterpreter','latex')
  hold on
fill(x_vertices, y_vertices,0.5.*[1 1 1]) % Green wedge outline
xlim([min(min(X)) max(max(X))])
hold on
axis equal tight;

axis equal
set(gca,'TickLabelInterpreter','latex')
