clear;
close all;

% phi = linspace(0,2*pi);
% theta = linspace(0,pi);
phi = (0 : 2*pi/100 : 2*pi) + 2*pi/200;
theta = (0 : pi/50: pi) + pi/100;
KD_a = pi;  % azimuth
KD_e = pi;  % elevation
N_a = 3;
N_e = 3;


x = zeros(length(theta), length(phi), N_a*N_e);
y = zeros(length(theta), length(phi), N_a*N_e);
z = zeros(length(theta), length(phi), N_a*N_e);

for ii = 0:1:N_e-1
    for jj=0:1:N_a-1
        
        ite = ii*N_a + jj + 1;
        
        %Note that the DFT matrix with
        %dimension Na x Ne is not equal to the Kronecker product of two DFT matrices
        %with dimension Na and Ne.
        
        %Use the Kronecker product of two DFT matrices with dimension Na
        %and Ne as the basis
        beta_a = 2*pi/N_a * jj;
        beta_e = 2*pi/N_e * ii;
        
        %Use DFT matrix with dimension Na x Ne as the basis
        %The radiant pattern may not be symmetric
%         beta_a = 2*pi/N_a * ite;
%         beta_e = 2*pi*ite/N_e/N_a;

        Psi_a = KD_a * sin(theta)' * sin(phi) + beta_a;
        Psi_e = KD_e * cos(theta)' * ones(1, length(phi)) + beta_e;
        AF = ( 1/N_a * sin(N_a * Psi_a/2)./sin(Psi_a/2) ) .*...
            ( 1/N_e * sin(N_e * Psi_e/2)./sin(Psi_e/2));
        
        AF = abs(AF);
        % converting to spherical coordinates
        for m=1:length(theta)
            for n=1:length(phi)
                x_all(m,n, ite)=AF(m,n)*sin(theta(m))*cos(phi(n));
                y_all(m,n, ite)=AF(m,n)*sin(theta(m))*sin(phi(n));
                z_all(m,n, ite)=AF(m,n)*cos(theta(m));
            end
        end
        
        %plotting routine
        figure,
        %         subplot(2,2,1);
        
        x = x_all(:, :, ite);
        y = y_all(:, :, ite);
        z = z_all(:, :, ite);
        
        surf(x,y,z)

        title('Radiation Pattern for Planar Antenna')
        xlabel('x-axis')
        ylabel('y-axis')
        zlabel('z-axis')
        xlim([-1 1])
        ylim([-1 1])
        zlim([-1 1])
        xl = xlim;
        yl = ylim;
        zl = zlim;
        hold on
        
        %         subplot(2,2,2)
        %         plot(x, z);
        mesh(zeros(size(x)) - xl(1),y,z)
        
        %         subplot(2,2,3)
        %         plot(y,z)
        mesh(x,zeros(size(y)) - yl(1),z)
        
        mesh(x,y,zeros(size(z)) + zl(1))

        hold off
    end;
end;

figure,
for ii = 0:1:N_e-1
    for jj=0:1:N_a-1
        ite = ii*N_a + jj + 1;
        subplot(N_e, N_a, ite);
        % Note that in the subplot, the figures are placed in the order as
        % follows:
        %[1 2
        % 3 4]
        x = x_all(:, :, ite);
        y = y_all(:, :, ite);
        z = z_all(:, :, ite);
        
        surf(x,y,z)

%         title('Radiation Pattern for Planar Antenna')
        xlabel('x-axis')
        ylabel('y-axis')
        zlabel('z-axis')
        xlim([-1 1])
        ylim([-1 1])
        zlim([-1 1])
        xl = xlim;
        yl = ylim;
        zl = zlim;
%         hold on
%         
%         %         subplot(2,2,2)
%         %         plot(x, z);
%         mesh(zeros(size(x)) - xl(1),y,z)
%         
%         %         subplot(2,2,3)
%         %         plot(y,z)
%         mesh(x,zeros(size(y)) - yl(1),z)
%         
%         mesh(x,y,zeros(size(z)) + zl(1))
% 
%         
%         hold off
    end;
end;

% subplot(221)
% surf(x,y,z)
% title('Radition Pattern for Planar Antenna')
% xlabel('x-axis--->')
% ylabel('y-axis--->')
% zlabel('z-axis--->')
%
% subplot(222)
% polar(theta', AF(:,1))
% title('x-z plane')
%
% subplot(223)
% polar(phi, AF(1,:))
% title('x-y plane')