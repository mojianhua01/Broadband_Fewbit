clear;
close all;

a=linspace(0,2*pi);     %????????
b=linspace(0,pi);
f=sin((cos(a).*sin(b)-1)*(14/2)*pi)./(sin((cos(a).*sin(b)-1)*pi/2)*14);
subplot(221);
polar(a,f.*sin(b));    %???
title('14????H??d=??/2,??=??');
y1=(f.*sin(a))'*cos(b);
z1=(f.*sin(a))'*sin(b);
x1=(f.*cos(a))'*ones(size(b));
subplot(223);
surf(x1,y1,z1);
axis equal    %???????????
title('14???????');

a=linspace(0,2*pi);
b=linspace(0,pi);
f=sin((cos(a).*sin(b)+1)*(6/2)*pi)./(sin((cos(a).*sin(b)+1)*pi/2)*6);
subplot(222);
polar(a,f.*sin(b));
title('6????H??d=??/2,??=??');
y1=(f.*sin(a))'*cos(b);
z1=(f.*sin(a))'*sin(b);
x1=(f.*cos(a))'*ones(size(b));
subplot(224);
surf(x1,y1,z1);
axis equal
title('6???????');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Usage:
%This program plots 3-D radiation Pattern of a Dipole Antenna
%All the parameters are entered in the M-File

clear all

%Defining variables in spherical coordinates
theta=[0:0.12:2*pi];%theta vector 
phi=[0:0.12:2*pi];%phi vector

l_lamda1=1/100;% length of antenna in terms of wavelengths
I0=1;% max current in antenna structure
n=120*pi; %eta
 
% evaluating radiation intensity(U)
U1=( n*( I0^2 )*( ( cos(l_lamda1*cos(theta-(pi/2))/2) - cos(l_lamda1/2) )./ sin(theta-(pi/2)) ).^2 )/(8*(pi)^2);

%converting to dB scale
U1_1=10*log10(U1);

%normalizing in order to make U vector positive
min1=min(U1_1);
U=U1_1-min1;
 
% expanding theta to span entire space
U(1,1)=0;
for n=1:length(theta)
    theta(n,:)=theta(1,:);  % ??theta(h,:)????????theta(1,:)???????meshgrid(x,y)??x??????theta(1,:) = theta
end

% expanding phi to span entire space
phi=phi'; %phi??????
for m=1:length(phi)
    phi(:,m)=phi(:,1); % ??phi(:,m)????????phi(:,1)???????meshgrid(x,y)??y??????
end
 
% expanding U to span entire space
for k=1:length(U)
    U(k,:)=U(1,:); % U???theta???????theta???????!
end
 
% converting to spherical coordinates 
[x,y,z]=sph2cart(phi,theta,U);

%plotting routine
surf(x,y,z)
colormap(copper)
title('Radition Pattern for Dipole Antenna (length=1.5lamda)')
xlabel('x-axis--->')
ylabel('y-axis--->')
zlabel('z-axis--->')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;clc;
theta=(0:pi/100:pi)+pi/100;
phi=0:pi/100:2*pi;
for m=1:length(theta)
    E(m)=cos(pi*cos(theta(m))/2)/sin(theta(m));
    for n=1:length(phi)
        x(m,n)=E(m)*sin(theta(m))*cos(phi(n));
        y(m,n)=E(m)*sin(theta(m))*sin(phi(n));
        z(m,n)=E(m)*cos(theta(m));
    end
end
figure,
surf(x,y,z);