% function [H_scaled, H_v] = func_mmWave_Channel_UPA(BS_ant_W, BS_ant_H, MS_ant_W,MS_ant_H, ...
%     Angle_spread_W, Angle_spread_H, Num_cluster, Num_ray, Narrowband_flag, Nd)
%Generate

%% setup for test
clear;
%clc;
close all;

rng(1)

BS_ant = 64;
MS_ant = 64;
Num_cluster = 2;
Num_ray = 10;
Angle_spread = pi/24; % pi/24;

% Simulation parameters
BS_ant_index = 0:1:BS_ant - 1;
MS_ant_index = 0:1:MS_ant - 1;
H = zeros(MS_ant,BS_ant);
KD = pi;  % Assuming: K=2pi/lambda, D=lambda/2

Num_cluster = 2;
Num_ray = 10;
Narrowband_flag = 0;
Nd = 16;

srrc = false; % use SRRC instead of RC time-domain pulse?
alf_rc = 0; % raised-cosine parameter [dflt=0.5]
n = [0:Nd-1]+eps;

for channel_index=1:1:1
    if Narrowband_flag == 1
        H = zeros(MS_ant , BS_ant);
        % Channel generation
        for L = 1:1:Num_cluster
            alpha(1,L) = 1;%(sqrt(1/2)*(randn(1,1)+1j*randn(1,1))); % gain
            AoA_phi_mean = pi*rand(1,1) - pi/2;   % azimuthal angle
            AoD_phi_mean = pi*rand(1,1) - pi/2;
            
            for ii = 1:1:Num_ray
                %         alpha( ii + (L-1) * Num_ray) = randn(1,2)*[1;1i]/sqrt(2); % gain
                alpha( ii + (L-1) * Num_ray) = 1;
                AoA_offset = (-1 + 2* rand(1,1))*Angle_spread;
                AoD_offset = (-1 + 2* rand(1,1))*Angle_spread;
                
                AoA_phi = AoA_phi_mean + AoA_offset;
                AoD_phi = AoD_phi_mean + AoD_offset;
                
                BS_steering_temp = sqrt(1/BS_ant)*exp(1j * KD * BS_ant_index * sin(AoD_phi));
                
                BS_steering(:,ii+(L-1)*Num_ray) = BS_steering_temp(:);
                
                MS_steering_temp = sqrt(1/MS_ant)*exp(1j * KD * MS_ant_index * sin(AoA_phi));
                
                MS_steering(:,ii+(L-1)*Num_ray) = MS_steering_temp(:);
                
                H = H + MS_steering * BS_steering' * exp(-AoA_offset) * exp(-AoD_offset) * alpha(L);
            end;
        end
        
        H = MS_steering* diag(alpha)*BS_steering';
        H = H./norm(H, 'fro')*sqrt(BS_ant * MS_ant);
        
        Hv = 1/sqrt(BS_ant * MS_ant) * dftmtx(MS_ant)'* H *dftmtx(BS_ant);

        
    else  %% Broadband channel
        H = zeros(MS_ant, BS_ant, Nd);
        % Channel generation
        for L = 1:1:Num_cluster
            alpha(1,L) = 1;%(sqrt(1/2)*(randn(1,1)+1j*randn(1,1))); % gain
            AoA_phi_mean = pi*rand(1,1) - pi/2;   % azimuthal angle
            AoD_phi_mean = pi*rand(1,1) - pi/2;
            
            for ii = 1:1:Num_ray
                %         alpha( ii + (L-1) * Num_ray) = randn(1,2)*[1;1i]/sqrt(2); % gain
                alpha( ii + (L-1) * Num_ray) = 1;
                AoA_offset = (-1 + 2* rand(1,1))*Angle_spread;
                AoD_offset = (-1 + 2* rand(1,1))*Angle_spread;
                
                AoA_phi = AoA_phi_mean + AoA_offset;
                AoD_phi = AoD_phi_mean + AoD_offset;
                
                BS_steering_temp = sqrt(1/BS_ant)*exp(1j * KD * BS_ant_index * sin(AoD_phi));
                
                BS_steering(:,ii+(L-1)*Num_ray) = BS_steering_temp(:);
                
                MS_steering_temp = sqrt(1/MS_ant)*exp(1j * KD * MS_ant_index * sin(AoA_phi));
                
                MS_steering(:,ii+(L-1)*Num_ray) = MS_steering_temp(:);
                
            end;
            
            Npre = min(2, Nd/2);
            Npst = min(2, Nd/2);
            
            delay(L) = Npre+rand(1)*(Nd-Npre-Npst); % delay (between Npre and Nd-Npst)
            if srrc
                pulse = ( (1-alf_rc)*sinc((n-delay(L))*(1-alf_rc)) ...
                    + cos(pi*(n-delay(l))*(1+alf_rc))*4*alf_rc/pi ...
                    )./(1-(4*alf_rc*(n-delay(L))).^2);
            else % rc pulse
                pulse = cos(pi*alf_rc*(n-delay(L)))./(1-(2*alf_rc*(n-delay(L))).^2) ...
                    .*sinc(n-delay(L)); % raised-cosine impulse response
            end
            pulse = pulse/norm(pulse);
            %plot(pulse,'.-'); title('pulse'); pause;
            
            for d=1:Nd % for each delay ...
                H(:,:,d) =  H(:,:,d) + MS_steering(:, (L-1)*Num_ray+1 : L*Num_ray)* ...
                    diag(alpha((L-1)*Num_ray+1 : L*Num_ray))*...
                    BS_steering(:, (L-1)*Num_ray+1:L*Num_ray)' * pulse(d);
                Hv(:,:,d) = 1/sqrt(BS_ant * MS_ant) * dftmtx(MS_ant)'* H(:,:,d) *dftmtx(BS_ant);
            end
        end
        
        
    end;
    H = H./norm(H(:)) * sqrt(BS_ant * MS_ant);
    Hv = Hv./norm(Hv(:)) * sqrt(BS_ant * MS_ant);
    Hv_all(:,:,:,channel_index) = Hv;
end;

% adopt the normalization proposed in Omar's paper where norm(H_v,
% 'fro')^2 = BS_ant * MS_ant
Hv_all = Hv_all./norm(Hv_all(:)) * sqrt(BS_ant * MS_ant)*sqrt(channel_index);

% for channel_index = 1:1:size(Hv_all, 4)
%     Hv = Hv_all(:,:,:, channel_index);
%     Hv_power(channel_index) = norm(Hv(:)).^2
% end;

% filename_save = sprintf('H_UPA_%d_%d_%d_%dclusters', BS_ant, MS_ant, Nd, Num_cluster);
% save(filename_save,'Hv_all')

% save H_UPA_64_64_16_4clusters.mat Hv_all;
% %

Narrowband_flag = 1; % force narrowband plots!
%Hv = sum(abs(Hv),3); H = sum(abs(H),3);
Hv = sqrt(sum(abs(Hv).^2,3)); H = sqrt(sum(abs(H).^2,3)); 
%Hv = sum(abs(Hv).^2,3); H = sum(abs(H).^2,3); % energy

if Narrowband_flag == 1
    figure(1)
    bar3(abs(Hv))
    set(gca, 'FontSize',16)
    grid on;
    xlabel('Tx angle')
    ylabel('Rx angle')
    zlabel('$\sqrt{\sum_{\ell} |[\mathbf{X}[\ell]]_{i,j}|^2}$', 'Interpreter', 'latex');
    xlim([1 BS_ant])
    ylim([1 MS_ant])
    
    figure(2)
    contour(1:1:BS_ant ,1:1:MS_ant, (abs(Hv)))
    set(gca, 'FontSize',16)
    grid on;
    xlabel('TX angle i')
    ylabel('RX angle j')
    zlabel('$|\mathbf{X}_{\mathrm{v}}|$', 'Interpreter', 'latex');
    
    figure(3)
    bar3(abs(H))
    set(gca, 'FontSize',16)
    xlabel('Tx antenna')
    ylabel('Rx antenna')
    zlabel('$\sqrt{\sum_{\ell} |[\mathbf{H}[\ell]]_{i,j}|^2}$', 'Interpreter', 'latex');
    xlim([1 BS_ant])
    ylim([1 MS_ant])
    
    
else % Broadband channel
    figure,
    bar3((abs(reshape(permute(Hv,[1,3,2]), MS_ant*Nd, BS_ant))))
    set(gca, 'FontSize',12)
    grid on;
    xlabel('Tx angle')
    xlim([0 16])
    ylim([0 64])
%     xticks([0 5 16])
    ylabel('Rx angle & delay')
    zlabel('$|\mathbf{X}_{i,j}|$', 'Interpreter', 'latex'); 
    figure,
    bar3((abs(reshape(permute(H,[1,3,2]), MS_ant*Nd, BS_ant))))
    set(gca, 'FontSize',12)
    grid on;
    xlabel('Tx antenna')
    xlim([0 16])
    ylim([0 64])
%     xticks([0 5 16])
    ylabel('Rx antenna & delay')
    zlabel('$|\mathbf{H}_{i,j}|$', 'Interpreter', 'latex'); 
end;

Hv_sample = Hv_all(:,:,:,1);
y = sort(abs(Hv_sample(:)), 'descend');
figure(4); clf
loglog(y, 'linewidth', 1);
% xlim([0 300])
ylim([0 10 * y(1)])
% hold on;  plot( 15* y(1) * (1:1:BS_ant*MS_ant*Nd).^(-1), 'r--', 'linewidth', 1.5)
coeff = sqrt(6)/pi;
% coeff = 1;
hold on;  plot( coeff * norm(y) * (1:1:BS_ant*MS_ant*Nd).^(-1), 'r--', 'linewidth', 1)
set(gca, 'fontsize', 12)
xlabel('index $i$', 'fontsize', 14, 'Interpreter','latex')
handle = legend('$|x_{(i)}|$', '$R i^{-1}$');
set(handle, 'Interpreter','latex', 'fontsize', 16)


