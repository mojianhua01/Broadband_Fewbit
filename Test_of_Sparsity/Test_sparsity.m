clear;
close all;
N = 200;
x = eps + (-3/2:0.0001:3/2);

figure, plot(x, abs( sin(N*pi*x)./N./sin(pi*x) ));
hold on;
plot(x, min(1, abs(1/2/N./(x- round(x)))),'r--');
ylim([0 1.2])

theta_vector = eps + (0: pi/N : 2*pi);
a = ones(1,N)';
xlabel('x', 'fontsize', 14)
h = legend('$\left|\frac{\sin(N \pi x)}{N \sin (\pi x)}\right|$','$\min\left(1, \frac{1}{2N \langle x \rangle} \right)$');
set(h,'Interpreter','latex', 'fontsize', 14)

for i=1:1:length(theta_vector)
    theta = theta_vector(i);
    b = exp(1j*theta*(1:1:N))';
    result(i) = abs(a'*b/N);
end;

figure,
plot(theta_vector, result);
hold on;
plot(theta_vector, min(1, 1/N./sin(theta_vector/2)),'r')
ylim([0 2])

Nt = 100;
Nr = 100;
for i = 1:1:Nt
    for j = 1:1:Nr
        X(i,j) = min(1, 1/2/Nr/abs((i/Nr - round(i/Nr)))) * ...
            min(1, 1/2/Nt/abs((j/Nt - round(i/Nt))));
    end;
end;
figure, 
bar3(abs(X))


y = sort(X(:), 'descend');
figure, loglog(y)
temp = sqrt(6)/pi * norm(X(:));
hold on;  plot( temp * (1:1:Nt*Nr).^(-1), 'r--')

N = 64;
delta = (1/N)* 1/8;
f1 = abs(sin(N*pi*(delta - (-N:1:N)/N))./N./sin(pi*(delta - (-N:1:N)/N)) );
temp = delta - (-N:1:N)/N;
temp = temp - round(temp);
f2 = abs(sin(pi*N*(delta))) /2/N./abs(temp);

figure,
plot(-N:1:N, f1, '-');
hold on;
plot(-N:1:N, f2, '--r');
xlim([-N N])
ylim([0 1])
xlabel('$i$', 'fontsize', 14, 'Interpreter','latex')
h = legend('$\left|\frac{\sin \left(\pi N \left( \delta - \frac{i}{N} \right) \right)}{N \sin \left(\pi \left(\delta - \frac{i}{N} \right) \right)} \right|$',...
    '$\frac{\left| \sin \left(\pi N \delta \right) \right|}{2 N \langle \delta - \frac{i}{N} \rangle}$');
set(h,'Interpreter','latex', 'fontsize', 14)

Nt = 100;
Nr = 100;
delta_r = 1/16/Nr + 0.01;
delta_t = 1/2/Nt + 0.01;
for i = 1:1:Nt
    for j = 1:1:Nr
        X(i,j) = sin(pi*Nr*(delta_r - i/Nr))/Nr/sin(pi*(delta_r - i/Nr)) * ...
            sin(pi*Nt*(delta_t - j/Nt))/Nt/sin(pi*(delta_t - j/Nt));
    end;
end;
figure, 
bar3(abs(X))

figure,
y = sort(abs(X(:)), 'descend');
loglog(y);
% xlim([0 300])
% ylim([0 10 * y(1)])
% temp1 = 6* y(1);
temp2 = sin(pi*Nr*delta_r) * sin(pi*Nt*delta_t) *( 1 + log(Nt * Nr/sin(pi*Nr*delta_r)/sin(pi*Nt*delta_t))); 
temp3 = sqrt(6)/pi * norm(X(:));
hold on;  plot( temp3 * ( (1:1:Nt*Nr)).^(-1), 'r--')
% plot(  y(1)* exp(- lambertw((1:1:BS_ant*MS_ant)/4)), 'k--')
ylabel('sorted |x|')
xlabel('index')
set(gca, 'fontsize', 12)
%title(['#clusters=', num2str( Num_cluster), ', #subpaths=', num2str(Num_ray) ] )