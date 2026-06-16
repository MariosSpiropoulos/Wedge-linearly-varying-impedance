clear all
clc 
close all 

% This code produces the directivity plots of figure 3 shown in the
% preprint : xxx


addpath('aux/')




% Wedge parameters
alpha = 10 *pi/180;
gamma= pi - alpha; 
theta_start = - pi + alpha ;
theta_end =  pi - alpha ;

% Polar grid 
Ntheta =150;
r=10;
theta = linspace(theta_start,theta_end,Ntheta);

an=@(kn,z) -(1i * pi/2)* cbesselh(kn,z);


dn=@(kn,z) -(1i * pi/2)* cbesselj(kn,z);


A = @(kn,alpha,Z2) 1 ./ ( gamma +  gamma.*(  sinc(4*kn.*(gamma)./pi)));



Z1   = Inf; 
i_Z  = 1;


r0 = 1 ;
theta0 = 0 ; 
k0 = 1; 
c0=343;
w0 = c0.*k0;


Nmodes= 80 ; %number of summation terms (eigenvalues) 

ik=1;
 figure('Units', 'inches', 'Position', [1, 1, 9, 2.7]);  % [left, bottom, width, height]

for k0     = [0.2*pi pi/2 pi]               %  wavenumber-loop

   
    i_Z=1;
    for Z2 = [Inf 0 exp(1i*pi/3).*k0]       % top-face impedance loop 
    
    

    %% Finding the eigenvalues
    nvals=theta_equation(alpha,Z1,Z2,Nmodes,theta);   % getting eigenvalues
    
    G=zeros(size(theta));
    
    
    
    
    for n=1:length(nvals)
    
        kn=nvals(n);
    
        
         An = cos((kn).*(theta0+gamma)).*A(kn,alpha,Z2);
    
        
         G_r_leq_r0= an(kn,k0.*r0).* cbesselj(kn,k0*r).*An.*cos(kn.*(theta+gamma));
    
        
         G_r_geq_r0= dn(kn,k0.*r0) .* cbesselh(kn,k0*r) .*An.*cos(kn.*(theta+gamma));
    
       
         G = G +(r<r0) .* G_r_leq_r0 + (r>r0) .* G_r_geq_r0 ;
    
         clear arg1 arg2 An kn  G_r_leq_r0 G_r_geq_r0
    
    
    end
    
    
    
    
        
        X=r.*cos(theta);
        Y=r.*sin(theta);
        
        if i_Z==1
            G1 = G ;
        elseif i_Z==2
            G2 = G ;
        elseif i_Z==3
            G3 = G;
        end
        

      
        i_Z=i_Z+1 ; 
        clear G 
    
         
    end
%% Plotting 

%% Add Wedge
   % Parameters
    Length = 15;      % Length of the triangle's sides from the origin
    half_theta_wedge = alpha; 
    % Calculate the vertices of the triangle
    origin = [0, 0]; % Vertex at the origin
    v1 = [Length * cos(half_theta_wedge), Length * sin(half_theta_wedge)];  % Vertex at +x angle
    v2 = [Length * cos(-half_theta_wedge), Length * sin(-half_theta_wedge)]; % Vertex at -x angle
    
    % Combine vertices into arrays
    x_vertices = [origin(1), -v1(1), -v2(1)];
    y_vertices = [origin(2), -v1(2), -v2(2)];
    wedge_center =pi/2;                                           % Center of wedge (90 degrees)
    wedge_angle = 2*alpha;                                        % Convert degrees to radians
    wedge_r = 0.7 ;                                               % Wedge radius
    
    % Define wedge boundary in polar coordinates
    wedge_theta = linspace(wedge_center - wedge_angle/2, wedge_center + wedge_angle/2, 20);
    wedge_r_values = ones(size(wedge_theta)) * wedge_r.*exp(1i.*pi); % Constant radius
    
    % Convert to Cartesian
    [wedge_x, wedge_y] = cart2pol(x_vertices, y_vertices);
    
    % Center point
    wedge_x = [0, wedge_x, 0]; 
    wedge_y = [0, wedge_y, 0];
    
    % Plot the wedge using a filled area approximation
    Line_x = [wedge_x, wedge_x(1)]; % Close the shape
    Line_y = [wedge_y, wedge_y(1)];


    subplot(1,3,ik)
    polarplot(theta,abs(G1),'k-',theta,abs(G2),'r-',theta,abs(G3),'b-','LineWidth',2)
    rlim([0 1.2.*max([max(abs(G1(:))),max(abs(G2(:))),max(abs(G3(:)))])])

    hold on
    plot(Line_x,Line_y,'Color',0.5.*[1 1 1], 'LineWidth', 3) % Green wedge outline
     set(gca,'FontSize',11)
    hold on
    set(gca,'TickLabelInterpreter','latex')


    grid on ; 

    ik=ik+1;
    clear G1 G2 G3
end

