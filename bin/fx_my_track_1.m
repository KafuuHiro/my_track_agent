function fx_my_track_1()

%% Initialization
clc;clear;close all;
addpath './func/';

if_path_loss = 0;           % 0--No 1--Yes
if_imagesc_log = 0;         % 0--No 1--Yes
if_steer_vec_norm = 1;      % 0--No 1--Yes      % steering and beamforming vector normalization
select_field_type = 'near';  % 'near'--nearfield 'mid'--middlefield 'far'--farfield  % only valid when requirement=2
select_path_loss_type = 'free_space'; % 'free_space', 'urban', 'indoor', 'suburban', 'rural' % only valid when if_path_loss=1

% Constants
c_0 = physconst('LightSpeed');
jmath = sqrt(-1);
% pi = 3.1415926;

% System parameters
f_c = 5e9;                     % 無線周波数 (G)「Hz」
lambda = c_0/f_c;              % 波長 λ [m]
delta_f = 30e3;                % 子载波间隔 Δf [Hz]
N_sc = 1;                      % 子载波个数 % 32 % 8 % 2
B_rf = (N_sc-1)*delta_f;       % 射频带宽 B [Hz]
Tofdm   = 1 / delta_f;         % OFDM 符号时长 [s]（参考）
Tp      = 0.5e-3;              % 慢时间重复周期 [s] —— PRF = 1/Tp
N_tp    = 512;                % 慢时间点数 % 64 % 256 % 2048

% Adjustable
ang_res = 1;                   % Resolution/Grain of angle  5    1
L_win_AoA = 5;                 % Length of searching window for Rx-AoA estimation and tracking
L_win_AoD = 5;                 % Length of searching window for Tx-AoD estimation and tracking

% Structure
% Rx-1
Rx_1.pos_first = [1,0];
Rx_1.array_ang = 0; % deg
Rx_1.N = 16; % 8    16
Rx_1.d = lambda/2;
[Rx_1.pos_array,Rx_1.pos_center] = gen_gNB(Rx_1.pos_first,Rx_1.array_ang,Rx_1.N,Rx_1.d);

% Tx-1
Tx_1.pos_first = [0,0];
Tx_1.array_ang = 0; % deg
Tx_1.N = 8; % 4 8
Tx_1.d = lambda/2;
[Tx_1.pos_array,Tx_1.pos_center] = gen_gNB(Tx_1.pos_first,Tx_1.array_ang,Tx_1.N,Tx_1.d);

% Tar-1
if select_field_type == "near"
    Tar_1.pos_start = [0.5,1];            % near field
end
if select_field_type == "mid"
    Tar_1.pos_start = [0.5,5];            % middle field    % [0.5,6];
end
if select_field_type == "far"
    Tar_1.pos_start = [10,10*sqrt(3)];    % far field
end
Tar_1.sigma = 0.5; % 0.5
Tar_1.sigma_v = 0.3;      % m/s   2   5   5   1   0.2   0.3 0.4 0.5
Tar_1.sigma_theta = 30; % deg   1.5 4   30  30  30
% Tar_1.trajectory_vec = gen_trajectory_0(Tar_1.pos_start,N_tp,Tar_1.sigma);
Tar_1.trajectory_array = gen_trajectory(Tar_1.pos_start,N_tp,Tar_1.sigma_v,Tar_1.sigma_theta,Tp);
Tar_1.velocity_array = gen_velocity_array(Tar_1.trajectory_array,Tp);

% Alpha
alpha = ones(N_tp,1);

% Storage
ang_array_deg = (-90:ang_res:90).';
ang_array_rad = ang_array_deg./180.*pi;
Rx_1.DFT_codebook = exp(-1*jmath*2*pi*(0:Rx_1.N-1)'*sin(ang_array_rad.').*Rx_1.d./lambda); % Rx_1.DFT_codebook = exp(1*jmath*2*pi*(0:Rx_1.N-1)'*sin(ang_array_rad.')/2);
Tx_1.DFT_codebook = exp(-1*jmath*2*pi*(0:Tx_1.N-1)'*sin(ang_array_rad.').*Tx_1.d./lambda); % Tx_1.DFT_codebook = exp(1*jmath*2*pi*(0:Tx_1.N-1)'*sin(ang_array_rad.')/2);



%% Near field channel generation
% H_tilde_nr_nt__k_l = zeros(Rx_1.N,Tx_1.N,N_sc,N_tp);
% for k = 1:N_sc
%     for l = 1:N_tp
%         for nr = 1:Rx_1.N
%             for nt = 1:Tx_1.N
%             end
%         end
%     end
% end

% Tx_1 to Tar_1
[Downlink_1_1.Distance,Downlink_1_1.Velocity] = cal_dis_velocity(Tx_1.pos_array,Tar_1.trajectory_array,Tar_1.velocity_array);

% Tar_1 to Rx_1
[Uplink_1_1.Distance,Uplink_1_1.Velocity] = cal_dis_velocity(Rx_1.pos_array,Tar_1.trajectory_array,Tar_1.velocity_array);

% Link_1_1_1
D_nr_nt__l = zeros(Rx_1.N,Tx_1.N,N_tp);
V_nr_nt__l = zeros(Rx_1.N,Tx_1.N,N_tp);
for it = 1:N_tp
    for ntt = 1:Tx_1.N
        D_nr_nt__l(:,ntt,it) = Uplink_1_1.Distance(it,:).'+Downlink_1_1.Distance(it,ntt);
        V_nr_nt__l(:,ntt,it) = Uplink_1_1.Velocity(it,:).'+Downlink_1_1.Velocity(it,ntt);
    end
end

% Sensing channel
% Link_1_1_1
% Path loss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if if_path_loss == 1
    for l = 1:N_tp
        alpha(l,1) = cal_path_loss(D_nr_nt__l(1,1,l),f_c,select_path_loss_type);
    end
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_tilde_nr_nt__k_l = zeros(Rx_1.N,Tx_1.N,N_sc,N_tp);
if if_steer_vec_norm == 0
    for k = 1:N_sc
        for l = 1:N_tp
            H_tilde_nr_nt__k_l(:,:,k,l) = flip(alpha(l,1)*exp(-1*jmath*2*pi*(f_c+(k-1)*delta_f)/c_0*D_nr_nt__l(:,:,l)).*exp(jmath*2*pi*(l-1)*Tp/lambda*V_nr_nt__l(:,:,l)));
        end
    end
end
if if_steer_vec_norm == 1
    for k = 1:N_sc
        for l = 1:N_tp
            H_tilde_nr_nt__k_l(:,:,k,l) = flip(alpha(l,1)*exp(-1*jmath*2*pi*(f_c+(k-1)*delta_f)/c_0*D_nr_nt__l(:,:,l)).*exp(jmath*2*pi*(l-1)*Tp/lambda*V_nr_nt__l(:,:,l))./sqrt(Rx_1.N*Tx_1.N));
        end
    end
end

% A priori information calculation
% Link_1_1_1
Rx_1.ang_real_array_rad = cal_ang_center(Tar_1.trajectory_array,Rx_1.pos_center,'rad');
Tx_1.ang_real_array_rad = cal_ang_center(Tar_1.trajectory_array,Tx_1.pos_center,'rad');
Rx_1.ang_real_array_deg = cal_ang_center(Tar_1.trajectory_array,Rx_1.pos_center,'deg');
Tx_1.ang_real_array_deg = cal_ang_center(Tar_1.trajectory_array,Tx_1.pos_center,'deg');



%% Test
% ang_input_deg = -23
% beam_number = ang2beam_number(ang_array_deg,ang_input_deg,'deg')
% ang_output_deg = beam_number2ang(ang_array_deg,beam_number)

% ang_input_deg = -23
% beam_number = ang2beam_number(ang_array_deg,ang_input_deg)
% ang_output_deg = beam_number2ang(ang_array_deg,beam_number)

% ang2beam_number(ang_array_deg,-90)



%% ISAC System
h_eq_tilde_k_l = zeros(N_sc,N_tp);
h_eq_tilde_bir_bit__k_l = zeros(size(Rx_1.DFT_codebook,2),size(Tx_1.DFT_codebook,2),N_sc,N_tp);
% h_eq_tilde_bir_bit__k_l = zeros(size(Rx_1.DFT_codebook,2),size(Tx_1.DFT_codebook,2),N_sc,1);

% % Suppose the initial position of the target is known to the ISAC system by beamsweeping
% % A priori inf. calc.
% beam_number_AoA_initial = ang2beam_number(ang_array_deg,Rx_1.ang_real_array_deg(1,1))
% beam_number_AoD_initial = ang2beam_number(ang_array_deg,Tx_1.ang_real_array_deg(1,1))
% 
% % if if_steer_vec_norm == 1
% % Rx-BS
% w = Rx_1.DFT_codebook(:,beam_number_AoA_initial)./sqrt(Rx_1.N);
% % Tx-UE
% f = Tx_1.DFT_codebook(:,beam_number_AoD_initial)./sqrt(Tx_1.N);
% k = 1; l = 1;
% h_eq_tilde_k_l(k,l) = w'*H_tilde_nr_nt__k_l(:,:,k,l)*f;
% abs(h_eq_tilde_k_l(k,l))

% % Rx-BS
% w = Rx_1.DFT_codebook(:,beam_number_AoA_initial-1)./sqrt(Rx_1.N);
% % Tx-UE
% f = Tx_1.DFT_codebook(:,beam_number_initial_AoD)./sqrt(Tx_1.N);
% k = 1; l = 1;
% h_eq_tilde_k_l(k,l) = w'*H_tilde_nr_nt__k_l(:,:,k,l)*f;
% abs(h_eq_tilde_k_l(k,l))

% if if_steer_vec_norm == 0
% if if_steer_vec_norm == 1
%     l = 1;
%     for beam_index_tx = 1:size(Tx_1.DFT_codebook,2)
%         for beam_index_rx = 1:size(Rx_1.DFT_codebook,2)
%             % Rx-BS
%             w = Rx_1.DFT_codebook(:,beam_index_rx)./sqrt(Rx_1.N);
%             % Tx-UE
%             v = Tx_1.DFT_codebook(:,beam_index_tx)./sqrt(Tx_1.N);
%             % Equivalent channel
%             for k = 1:N_sc
%                 % for l = 1:size(Rx_1.DFT_codebook,2)*size(Tx_1.DFT_codebook,2) % 1:N_tp
%                     h_eq_tilde_k_l(k,l) = w'*H_tilde_nr_nt__k_l(:,:,k,l)*v;
%                     h_eq_tilde_bir_bit__k_l(beam_index_rx,beam_index_tx,k,l) = h_eq_tilde_k_l(k,l);
%                 % end
%             end
%             l = l+1;
%         end
%     end
% end
% Figure inf. calc.
% k_fix = 1;
% Distance_Rx1_Tx1 = norm(Rx_1.pos_first-Tx_1.pos_first);
% h_eq_tilde_bi = zeros(size(ang_array_deg,1),1);
% for beam_index = 1:size(ang_array_deg,1)
%     h_eq_tilde_bi(beam_index,1) = h_eq_tilde_bir_bit__k_l(beam_index,size(ang_array_deg,1)-beam_index+1,k_fix,1);
% end

% Initial beamsweeping
if if_steer_vec_norm == 1
    for l = 1:1
        for k = 1:N_sc
            for beam_number_tx = 1:size(Tx_1.DFT_codebook,2)
                for beam_number_rx = 1:size(Rx_1.DFT_codebook,2)
                    % Rx-BS
                    w = Rx_1.DFT_codebook(:,beam_number_rx)./sqrt(Rx_1.N);
                    % Tx-UE
                    f = Tx_1.DFT_codebook(:,beam_number_tx)./sqrt(Tx_1.N);
                    h_eq_tilde_bir_bit__k_l(beam_number_rx,beam_number_tx,k,l) = w'*H_tilde_nr_nt__k_l(:,:,k,l)*f;
                end
            end
        end
    end
    figure
    % subplot(1,2,1);
    imagesc(abs(h_eq_tilde_bir_bit__k_l(:,:,1,1)));colorbar;
    xlabel('beam number of AoD [Tx]');hold on;
    ylabel('beam number of AoA [Rx]');hold on;
    % subplot(1,2,2);
    % imagesc(ang_array_deg,ang_array_deg,abs(h_eq_tilde_bir_bit__k_l(:,:,1,1)));colorbar;
    % xlabel('AoD [deg]');hold on;
    % ylabel('AoA [deg]');hold on;
else
    % if_steer_vec_norm == 0 % or else
    for l = 1:1
        for k = 1:N_sc
            for beam_number_tx = 1:size(Tx_1.DFT_codebook,2)
                for beam_number_rx = 1:size(Rx_1.DFT_codebook,2)
                    % Rx-BS
                    w = Rx_1.DFT_codebook(:,beam_number_rx);
                    % Tx-UE
                    f = Tx_1.DFT_codebook(:,beam_number_tx);
                    h_eq_tilde_bir_bit__k_l(beam_number_rx,beam_number_tx,k,l) = w'*H_tilde_nr_nt__k_l(:,:,k,l)*f;
                end
            end
        end
    end
    figure
    % subplot(1,2,1);
    imagesc(abs(h_eq_tilde_bir_bit__k_l(:,:,1,1)));colorbar;
    xlabel('beam number of AoD [Tx]');hold on;
    ylabel('beam number of AoA [Rx]');hold on;
    % subplot(1,2,2);
    % imagesc(ang_array_deg,ang_array_deg,abs(h_eq_tilde_bir_bit__k_l(:,:,1,1)));colorbar;
    % xlabel('AoD [deg]');hold on;
    % ylabel('AoA [deg]');hold on;
    sgtitle('Initial beamsweeping result at k=1, l=1');
end

% Tracking process

% Suppose the initial position of the target is known to the ISAC system by beamsweeping
% A priori inf. calc.
beam_number_AoA_initial = ang2beam_number(ang_array_deg,Rx_1.ang_real_array_deg(1,1))
beam_number_AoD_initial = ang2beam_number(ang_array_deg,Tx_1.ang_real_array_deg(1,1))
if if_steer_vec_norm == 1
    % Rx-BS
    w = Rx_1.DFT_codebook(:,beam_number_AoA_initial)./sqrt(Rx_1.N);
    % Tx-UE
    f = Tx_1.DFT_codebook(:,beam_number_AoD_initial)./sqrt(Tx_1.N);
else
    % if_steer_vec_norm == 0 % or else
    % Rx-BS
    w = Rx_1.DFT_codebook(:,beam_number_AoA_initial);
    % Tx-UE
    f = Tx_1.DFT_codebook(:,beam_number_AoD_initial);
end
center_number_win_AoA = beam_number_AoA_initial;
center_number_win_AoD = beam_number_AoD_initial;
% Neibourhood searching
%  1    2    3    4    0
%  1   -1    2   -2    0
Rx_Tx_turn_indicator = 1; % (1)--Rx    (-1)--Tx
temp_maximum_h_eq_tilde = 0;
temp_number_AoA__maximum_h_eq_tilde = 0;
temp_number_AoD__maximum_h_eq_tilde = 0;
rec_center_number_win_AoA = zeros(ceil(N_tp/(L_win_AoA+L_win_AoD)),1);
rec_center_number_win_AoD = zeros(ceil(N_tp/(L_win_AoA+L_win_AoD)),1);
rec_channel_gain_max_period = zeros(ceil(N_tp/(L_win_AoA+L_win_AoD)),1);
if if_steer_vec_norm == 1
    for l = 1:N_tp
        for k = 1:N_sc
            if Rx_Tx_turn_indicator == 1 % --------Rx
                search_number_win_AoA = center_number_win_AoA + (-1)^(mod(l,L_win_AoA)+1)*ceil(mod(l,L_win_AoA)/2);
                if search_number_win_AoA < 1
                    search_number_win_AoA = 1;
                end
                if search_number_win_AoA > size(Rx_1.DFT_codebook,2)
                    search_number_win_AoA = size(Rx_1.DFT_codebook,2);
                end
                if center_number_win_AoD < 1
                    center_number_win_AoD = 1;
                end
                if center_number_win_AoD > size(Tx_1.DFT_codebook,2)
                    center_number_win_AoD = size(Tx_1.DFT_codebook,2);
                end
                % Rx-BS
                w = Rx_1.DFT_codebook(:,search_number_win_AoA);
                % Tx-UE
                f = Tx_1.DFT_codebook(:,center_number_win_AoD);
                h_eq_tilde_k_l(k,l) = w'*H_tilde_nr_nt__k_l(:,:,k,l)*f;
                h_eq_tilde_bir_bit__k_l(search_number_win_AoA,center_number_win_AoD,k,l) = h_eq_tilde_k_l(k,l);
                % maximum channel gain recording
                if abs(h_eq_tilde_k_l(k,l)) > temp_maximum_h_eq_tilde
                    temp_maximum_h_eq_tilde = abs(h_eq_tilde_k_l(k,l));
                    temp_number_AoA__maximum_h_eq_tilde = search_number_win_AoA;
                end
                if mod(l,L_win_AoA) == 0 % Potential error
                    center_number_win_AoA = temp_number_AoA__maximum_h_eq_tilde;
                    rec_center_number_win_AoA(ceil(l/(L_win_AoA+L_win_AoD)),1) = center_number_win_AoA;
                    temp_maximum_h_eq_tilde = 0;
                    Rx_Tx_turn_indicator = (-1)*Rx_Tx_turn_indicator;
                end
            else
                % Rx_Tx_turn_indicator == -1 % ----Tx
                search_number_win_AoD = center_number_win_AoD + (-1)^(mod(l,L_win_AoD)+1)*ceil(mod(l,L_win_AoD)/2);
                if search_number_win_AoD < 1
                    search_number_win_AoD = 1;
                end
                if search_number_win_AoD > size(Tx_1.DFT_codebook,2)
                    search_number_win_AoD = size(Tx_1.DFT_codebook,2);
                end
                if center_number_win_AoA < 1
                    center_number_win_AoA = 1;
                end
                if center_number_win_AoA > size(Rx_1.DFT_codebook,2)
                    center_number_win_AoA = size(Rx_1.DFT_codebook,2);
                end
                % Rx-BS
                w = Rx_1.DFT_codebook(:,center_number_win_AoA);
                % Tx-UE
                f = Tx_1.DFT_codebook(:,search_number_win_AoD);
                h_eq_tilde_k_l(k,l) = w'*H_tilde_nr_nt__k_l(:,:,k,l)*f;
                h_eq_tilde_bir_bit__k_l(center_number_win_AoA,search_number_win_AoD,k,l) = h_eq_tilde_k_l(k,l);
                % maximum channel gain recording
                if abs(h_eq_tilde_k_l(k,l)) > temp_maximum_h_eq_tilde
                    temp_maximum_h_eq_tilde = abs(h_eq_tilde_k_l(k,l));
                    temp_number_AoD__maximum_h_eq_tilde = search_number_win_AoD;
                end
                if mod(l,L_win_AoD) == 0 % Potential error
                    center_number_win_AoD = temp_number_AoD__maximum_h_eq_tilde;
                    rec_center_number_win_AoD(ceil(l/(L_win_AoA+L_win_AoD)),1) = center_number_win_AoD;
                    rec_channel_gain_max_period(ceil(l/(L_win_AoA+L_win_AoD)),1) = temp_maximum_h_eq_tilde; % special of AoD tracking process
                    temp_maximum_h_eq_tilde = 0;
                    Rx_Tx_turn_indicator = (-1)*Rx_Tx_turn_indicator;
                end
            end
        end
    end
else
    % if_steer_vec_norm == 0 % or else
end

% Trimming rec_center_number_win_AoA & rec_center_number_win_AoD & rec_channel_gain_max_period
if rec_center_number_win_AoA(size(rec_center_number_win_AoA,1),1) == 0
    rec_center_number_win_AoA(size(rec_center_number_win_AoA,1),1) = rec_center_number_win_AoA(size(rec_center_number_win_AoA,1)-1,1);
end
if rec_center_number_win_AoD(size(rec_center_number_win_AoD,1),1) == 0
    rec_center_number_win_AoD(size(rec_center_number_win_AoD,1),1) = rec_center_number_win_AoD(size(rec_center_number_win_AoD,1)-1,1);
end
if rec_channel_gain_max_period(size(rec_channel_gain_max_period,1),1) == 0
    rec_channel_gain_max_period(size(rec_channel_gain_max_period,1),1) = rec_channel_gain_max_period(size(rec_channel_gain_max_period,1)-1,1);
end

% Reconstruction of target trajectory in space domain
Tar_1.trajectory_array_est = zeros(size(rec_center_number_win_AoA,1),2);
for l_period = 1:size(rec_center_number_win_AoA,1)
    Tar_1.trajectory_array_est(l_period,:) = cal_ang_2_pos(Rx_1.pos_center,Tx_1.pos_center,beam_number2ang(ang_array_deg,rec_center_number_win_AoA(l_period,1)),beam_number2ang(ang_array_deg,rec_center_number_win_AoD(l_period,1)),'deg');
end
Tar_1.velocity_array_est = gen_velocity_array(Tar_1.trajectory_array_est,Tp);



%% Figure
figure
plot(Tx_1.pos_array(:,1),Tx_1.pos_array(:,2),'b:^');hold on;
plot(Rx_1.pos_array(:,1),Rx_1.pos_array(:,2),'r:sq');hold on;
plot(Tar_1.pos_start(:,1),Tar_1.pos_start(:,2),'g-o','MarkerSize',10);hold on;
plot(Tar_1.trajectory_array(:,1),Tar_1.trajectory_array(:,2),'k:o');hold on;
legend('Tx-1','Rx-1','Target starting position','Target trajectory');hold on;
xlabel('x [m]');hold on;
ylabel('y [m]');hold on;
title('Layout diagram')

figure
plot(vector2scalar(Tar_1.velocity_array),'k-o');hold on;
xlabel('slot index');hold on;
ylabel('velocity scalar [m/s]');hold on;
title('velocity trend');

% Tracking results of the designed ISAC system
figure
plot(Tar_1.trajectory_array(:,1),Tar_1.trajectory_array(:,2),'k-o');hold on;
plot(Tar_1.trajectory_array_est(:,1),Tar_1.trajectory_array_est(:,2),'b-sq');hold on;
plot(Tar_1.pos_start(:,1),Tar_1.pos_start(:,2),'g-o','MarkerSize',10);hold on;
legend('True trajectory','Estimated trajectory','Target starting position');hold on;
xlabel('x [m]');hold on;
ylabel('y [m]');hold on;
title('Trajectory tracking results');

figure
plot(vector2scalar(Tar_1.velocity_array(1:(L_win_AoA+L_win_AoD):N_tp)),'k-o');hold on;
plot(vector2scalar(Tar_1.velocity_array_est),'b-sq');hold on;
legend('True velocity','Estimated velocity');hold on;
xlabel('slot index');hold on;
ylabel('velocity scalar [m/s]');hold on;
title('Velocity tracking results');

figure
% plot(abs(h_eq_tilde_k_l));
plot(rec_channel_gain_max_period);
xlabel('slot index');hold on;
ylabel('channel gain');hold on;
title('Channel gain trend');

figure
plot(Rx_1.ang_real_array_deg(1:(L_win_AoA+L_win_AoD):N_tp),'r-o', 'LineWidth', 1, 'MarkerSize', 3, 'MarkerFaceColor', 'r');hold on;
plot(beam_number2ang(ang_array_deg,rec_center_number_win_AoA),'m-x', 'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', 'r');hold on;
plot(Tx_1.ang_real_array_deg(1:(L_win_AoA+L_win_AoD):N_tp),'b-o', 'LineWidth', 1, 'MarkerSize', 3, 'MarkerFaceColor', 'b');hold on;
plot(beam_number2ang(ang_array_deg,rec_center_number_win_AoD),'g-x', 'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', 'b');hold on;
xlabel('slot index');hold on;
ylabel('angle [deg]');hold on;
legend('Real AoA','Esti AoA','Real AoD','Esti AoD');hold on;
title('DoA estimation performance trend');



%% Functions
function [pos_array,pos_center] = gen_gNB(pos_start,array_ang,N,d)
    pos_array = zeros(N,2);
    for it = 1:N
        pos_array(it,:)=pos_start+(it-1)*d*[cosd(array_ang),sind(array_ang)];
    end
    pos_center = sum(pos_array)/N;
end

function trajectory_array = gen_trajectory_0(pos_start,points,sigma)
    % trajectory_vec = zeros(points,2);
    % trajectory_vec(1,:) = pos_start;
    % trajectory_vec(2,:) = 0.5*pos_start+0.5*sigma^2*randn(1,2);
    % for it = 3:points
    %     % trajectory_vec(it,:) = 0.5*trajectory_vec(it-1,:)+0.2*trajectory_vec(it-2,:)+0.3*sigma^2*randn(1,2);
    %     trajectory_vec(it,:) = 0.5*trajectory_vec(it-1,:)+0.5*sigma^2*randn(1,2);
    % end

    trajectory_array = zeros(points,2);
    trajectory_array(1,:) = pos_start;
    trajectory_array(2,:) = pos_start+sigma^2*randn(1,2);
    for it = 3:points
        % trajectory_vec(it,:) = 0.5*trajectory_vec(it-1,:)+0.2*trajectory_vec(it-2,:)+0.3*sigma^2*randn(1,2);
        trajectory_array(it,:) = trajectory_array(it-1,:)+sigma^2*randn(1,2);
    end
end

function trajectory_array = gen_trajectory(pos_start,points,sigma_v,sigma_theta,Tp)
    v_array=zeros(points,2);
    trajectory_array = zeros(points,2);
    trajectory_array(1,:) = pos_start;
    trajectory_array(2,:) = pos_start+(sigma_v*Tp)^2*randn(1,2);
    % v_vec(2,:) = trajectory_vec(2,:)-trajectory_vec(1,:);
    for it = 3:points
        v_array(it-1,:) = (trajectory_array(it-1,:)-trajectory_array(it-2,:))/Tp;
        omega = sigma_theta^2*randn(1,1);
        % v_vec(it,:) = 0.5*v_vec(it-1,:)+0.5*sigma_v^2*randn(1,1)*([cosd(omega),-sind(omega);sind(omega),cosd(omega)]*v_vec(it-1,:).').';
        v_array(it,:) = v_array(it-1,:)+sigma_v^2*randn(1,1)*([cosd(omega),-sind(omega);sind(omega),cosd(omega)]*[1;0]).';
        trajectory_array(it,:) = trajectory_array(it-1,:)+v_array(it,:)*Tp;
    end
end

function velocity_array = gen_velocity_array(trajectory_array,Tp)
    N_tp = size(trajectory_array,1);
    velocity_array = zeros(N_tp,2);
    for it = 2:N_tp
        velocity_array(it,:) = (trajectory_array(it,:)-trajectory_array(it-1,:))./Tp;
    end
    velocity_array(1,:) = velocity_array(2,:);
end

function scalar = vector2scalar(vector)
    scalar = zeros(size(vector,1),1);
    for it = 1:size(vector,1)
        scalar(it) = sqrt(vector(it,1)^2+vector(it,2)^2);
    end
end

function ang_real_array = cal_ang_center(tar_trajectory_array,gNB_pos_center,ang_string)
    % Use the centered vector of the RRU position, near field error as a result
    % Output: °[rad]
    fac = 1;
    ang_real_array = zeros(size(tar_trajectory_array,1),1);
    if ang_string == "rad"
        for it = 1:size(tar_trajectory_array,1)
            ang_real_array(it,1) = fac*atan((tar_trajectory_array(it,1)-gNB_pos_center(1,1))/(tar_trajectory_array(it,2)-gNB_pos_center(1,2)));
        end
    end
    if ang_string == "deg"
        for it = 1:size(tar_trajectory_array,1)
            ang_real_array(it,1) = fac*atand((tar_trajectory_array(it,1)-gNB_pos_center(1,1))/(tar_trajectory_array(it,2)-gNB_pos_center(1,2)));
        end
    end
end

% function beam_number = ang2beam_number(ang_array,ang,ang_string)
%     if ang_string == "deg"
%         parameter_rad_or_deg = 180;
%     end
%     if ang_string == "rad"
%         parameter_rad_or_deg = pi;
%     end
%     beam_number_rough = round((ang-ang_array(1,1))/parameter_rad_or_deg*size(ang_array,1))+1;
%     beam_number = beam_number_rough;
%     if ( beam_number_rough > 1 ) && ( beam_number_rough <= size(ang_array,1) )
%         if abs(ang_array(beam_number_rough-1,1)-ang) < abs(ang_array(beam_number_rough,1)-ang)
%             beam_number = beam_number_rough-1;
%         end
%     end
%     if beam_number_rough > size(ang_array,1)
%         beam_number = size(ang_array,1);
%     end
% end


function beam_number = ang2beam_number(ang_array,ang)
    parameter = ang_array(size(ang_array,1),1)-ang_array(1,1);
    beam_number_rough = round((ang-ang_array(1,1))/parameter*size(ang_array,1))+1;
    beam_number = beam_number_rough;
    if ( beam_number_rough > 1 ) && ( beam_number_rough <= size(ang_array,1) )
        if abs(ang_array(beam_number_rough-1,1)-ang) < abs(ang_array(beam_number_rough,1)-ang)
            beam_number = beam_number_rough-1;
        end
    end
    if beam_number_rough > size(ang_array,1)
        beam_number = size(ang_array,1);
    end
end

function ang = beam_number2ang(ang_array,beam_number)
    ang = ang_array(beam_number,1);
end



end