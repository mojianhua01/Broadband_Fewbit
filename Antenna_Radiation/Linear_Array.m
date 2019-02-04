clear;
close all;

% phi = linspace(0,2*pi);
% theta = linspace(0,pi);
phi = (0 : 2*pi/100 : 2*pi) + 2*pi/200;

KD_a = pi;  % azimuth

N_a = 6;


for ii=0:1:N_a-1
    
    beta_a = 2*pi/N_a * ii;
    
    Psi_a = KD_a * cos(phi) + beta_a;
    
    AF = ( 1/N_a * sin(N_a * Psi_a/2)./sin(Psi_a/2) );
    
    % converting to spherical coordinates
    %plotting routine
    subplot(1, N_a, ii+1);
    polar(phi, abs(AF))
    title('Radiation Pattern')
    
end;
