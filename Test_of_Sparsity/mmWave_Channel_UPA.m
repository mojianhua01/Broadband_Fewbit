% function [H_scaled, H_v] = func_mmWave_Channel_UPA(BS_ant_W, BS_ant_H, MS_ant_W,MS_ant_H, ...
%     Angle_spread_W, Angle_spread_H, Num_cluster, Num_ray, Narrowband_flag, Nd)
%Generate

%% setup for test
clear;
clc;
close all;

BS_ant_W = 8;
BS_ant_H = 8;
MS_ant_W = 8;
MS_ant_H = 8;
Angle_spread_W = pi/24; % pi/24;
Angle_spread_H = pi/24; % pi/24;
Num_cluster = 2;
Num_ray = 10;
Narrowband_flag = 1;
Nd = 16;

srrc = false; % use SRRC instead of RC time-domain pulse?
alf_rc = 0; % raised-cosine parameter [dflt=0.5]
n = [0:Nd-1]+eps;

% Simulation parameters
BS_ant = BS_ant_W * BS_ant_H;
MS_ant = MS_ant_W * MS_ant_H;
BS_ant_W_index=0:1:BS_ant_W-1;
BS_ant_H_index=0:1:BS_ant_H-1;
MS_ant_W_index=0:1:MS_ant_W-1;
MS_ant_H_index=0:1:MS_ant_H-1;
%     H=zeros(MS_ant,BS_ant);
KD=pi;  % Assuming: K=2pi/lambda, D=lambda/2
for channel_index=1:1:1
    if Narrowband_flag == 1
        H = zeros(MS_ant , BS_ant);
        % Channel generation
        for L = 1:1:Num_cluster
            alpha_mean(L) = randn(1,2)*[1;1i]/sqrt(2); % gain
            AoA_phi_mean(L) = 2*pi*rand(1,1) ;   % azimuthal angle
            AoA_theta_mean(L) = pi*rand(1,1) ; % polar angle
            AoD_phi_mean(L) = 2*pi*rand(1,1) ;
            AoD_theta_mean(L) = pi*rand(1,1) ;
            
            for ii = 1:1:Num_ray
                AoA_phi( ii + (L-1) * Num_ray) = AoA_phi_mean(L) + (-1 + 2* rand(1,1))*Angle_spread_W;      % azimuthal angle spread
                AoA_theta( ii + (L-1) * Num_ray) = AoA_theta_mean(L) + (-1 + 2* rand(1,1))*Angle_spread_H;  % polar angle spread
                AoD_phi( ii + (L-1) * Num_ray) = AoD_phi_mean(L) + (-1 + 2* rand(1,1))*Angle_spread_W;
                AoD_theta( ii + (L-1) * Num_ray) = AoD_theta_mean(L) + (-1 + 2* rand(1,1))*Angle_spread_H;
                alpha( ii + (L-1) * Num_ray) = alpha_mean(L); % gain
                
                BS_steering_temp = sqrt(1/BS_ant)*...
                    exp(1j*KD*BS_ant_W_index*sin(AoD_phi( ii + (L-1) * Num_ray))*sin(AoD_theta( ii + (L-1) * Num_ray))).' ...
                    * exp(1j*KD*BS_ant_H_index*cos(AoD_theta( ii + (L-1) * Num_ray)));
                
                BS_steering(:,ii+(L-1)*Num_ray) = BS_steering_temp(:);
                
                MS_steering_temp = sqrt(1/MS_ant)*...
                    exp(1j*KD*MS_ant_W_index*sin(AoA_phi( ii + (L-1) * Num_ray))*sin(AoA_theta( ii + (L-1) * Num_ray))).' ...
                    * exp(1j*KD*MS_ant_H_index*cos(AoA_theta( ii + (L-1) * Num_ray)));
                
                MS_steering(:,ii+(L-1)*Num_ray) = MS_steering_temp(:);
                
            end;
        end
        
        H = MS_steering* diag(alpha)*BS_steering';
        H = H./norm(H, 'fro')*sqrt(BS_ant * MS_ant);
        
        %         Hv = 1/sqrt(BS_ant * MS_ant) * dftmtx(MS_ant)'* H *dftmtx(BS_ant);
        Hv = 1/sqrt(BS_ant * MS_ant) * kron( dftmtx(MS_ant_W), dftmtx(MS_ant_H) )' * H ...
            * kron( dftmtx(BS_ant_W), dftmtx(BS_ant_H) );
        
    else  %% Broadband channel
        H = zeros(MS_ant, BS_ant, Nd);
        % Channel generation
        for L = 1:1:Num_cluster
            
            AoA_phi_mean = 2*pi*rand(1,1) ;   % azimuthal angle
            AoA_theta_mean = pi*rand(1,1) ; % polar angle
            AoD_phi_mean = 2*pi*rand(1,1) ;
            AoD_theta_mean = pi*rand(1,1) ;
            
            for ii = 1:1:Num_ray
                alpha( ii + (L-1) * Num_ray) = randn(1,2)*[1;1i]/sqrt(2); % gain
                AoA_phi = AoA_phi_mean + (-1 + 2* rand(1,1))*Angle_spread_W;      % azimuthal angle spread
                AoA_theta = AoA_theta_mean + (-1 + 2* rand(1,1))*Angle_spread_H;  % polar angle spread
                AoD_phi = AoD_phi_mean + (-1 + 2* rand(1,1))*Angle_spread_W;
                AoD_theta = AoD_theta_mean + (-1 + 2* rand(1,1))*Angle_spread_H;
                
                BS_steering_temp = sqrt(1/BS_ant)*exp(1j*KD*BS_ant_W_index*sin(AoD_phi)*sin(AoD_theta))' * exp(1j*KD*BS_ant_H_index*cos(AoD_theta));
                
                BS_steering(:,ii+(L-1)*Num_ray) = BS_steering_temp(:);
                
                MS_steering_temp = sqrt(1/MS_ant)*exp(1j*KD*MS_ant_W_index*sin(AoA_phi)*sin(AoA_theta))' * exp(1j*KD*MS_ant_H_index*cos(AoA_theta));
                
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
                Hv(:,:,d) = 1/sqrt(BS_ant * MS_ant) * kron(dftmtx(MS_ant_W), dftmtx(MS_ant_H))'* H(:,:,d) * kron(dftmtx(BS_ant_W), dftmtx(BS_ant_H));
                %                 Hv(:,:,d) = 1/sqrt(BS_ant * MS_ant) * dftmtx(MS_ant)'* H(:,:,d) *dftmtx(BS_ant);
            end
        end
        
        
    end;
    Hv_all(:,:,:,channel_index) = Hv;
end;

% adopt the normalization proposed in Omar's paper where norm(H_v,
% 'fro')^2 = BS_ant * MS_ant
Hv_all = Hv_all./norm(Hv_all(:)) * sqrt(BS_ant * MS_ant)*sqrt(channel_index);

% for channel_index = 1:1:size(Hv_all, 4)
%     Hv = Hv_all(:,:,:, channel_index);
%     Hv_power(channel_index) = norm(Hv(:)).^2
% end;

filename_save = sprintf('H_UPA_%d_%d_%d_%dclusters', BS_ant, MS_ant, Nd, Num_cluster);
save(filename_save,'Hv_all')

% save H_UPA_64_64_16_4clusters.mat Hv_all;
% %

if Narrowband_flag == 1
    figure,
    bar3((abs(Hv)))
    set(gca, 'FontSize',16)
    grid on;
    xlabel('Tx angle')
    ylabel('Rx angle')
    zlabel('$|\mathbf{H}_{\mathrm{v}}|$', 'Interpreter', 'latex');
    
    figure,
    contour(1:1:BS_ant ,1:1:MS_ant, (abs(Hv)))
    set(gca, 'FontSize',16)
    grid on;
    xlabel('TX angle i')
    ylabel('RX angle j')
    zlabel('$|\mathbf{X}_{\mathrm{v}}|$', 'Interpreter', 'latex');
    
    figure(1)
    subplot(1,2,1)
    bar3(abs(H))
    xlabel('TX antenna i')
    ylabel('RX antenna j')
    set(gca, 'FontSize',14)
    zlabel('$|\mathbf{H}_{i,j}|$', 'Interpreter', 'latex');
    %title('Spatial Domain', 'FontSize', 20)
    
    subplot(1,2,2)
    bar3(abs(Hv))
    % stem3(1:1:BS_ant ,1:1:MS_ant, (abs(H_v)))
    set(gca, 'FontSize',14)
    grid on;
    xlabel('TX angle j')
    ylabel('RX agnle i')
    %title('Angular Domain', 'FontSize', 20)
    zlabel('$|\mathbf{X}_{i,j}|$', 'Interpreter', 'latex');
    
else % Broadband channel
    figure,
    bar3((abs(reshape(permute(Hv,[1,3,2]), MS_ant*Nd, BS_ant))))
    set(gca, 'FontSize',12)
    grid on;
    xlabel('Tx angle')
    ylabel('Rx angle & delay')
    zlabel('$|\mathbf{H}_{\mathrm{v}}|$', 'Interpreter', 'latex');
    
end;

Hv_sample = Hv_all(:,:,:,1);
y = sort(abs(Hv_sample(:)), 'descend');
figure,
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

