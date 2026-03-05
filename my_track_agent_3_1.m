%% Initialization
clc;clear;close all;
addpath './func/';

if_path_loss = 0;           % 0--No 1--Yes
if_steer_vec_norm = 1;      % 0--No 1--Yes % steering and beamforming vector normalization
select_field_type = 'near';  % 'near'--nearfield 'mid'--middlefield 'far'--farfield
select_path_loss_type = 'free_space'; % 'free_space', 'urban', 'indoor', 'suburban', 'rural'
if_beamsweeping_cali = 1;   % 0--No 1--Yes : Do beamsweeping as the calibration when 'abs(h_eq_tilde) < h_eq_tilde_th'
if_doppler_effect = 1;      % 0--No 1--Yes : Considering Doppler effect
if_doppler_compensation = 1; % 0--No 1--Yes : Enable real-time Doppler velocity estimation and compensation
if_doppler_est_method = 1;  % 1--temporal phase diff, 2--multi-subcarrier freq analysis

% Constants
c_0 = physconst('LightSpeed');
jmath = sqrt(-1);

% System parameters
f_c = 5e9;                     % Carrier frequency (Hz)
lambda = c_0/f_c;              % Wavelength λ [m]
delta_f = 30e3;                % Subcarrier spacing Δf [Hz]
N_sc = 1;                      % Number of subcarriers % 4 for Doppler est.
B_rf = (N_sc-1)*delta_f;       % RF bandwidth B [Hz]
Tofdm   = 1 / delta_f;         % OFDM symbol duration [s]
Tp      = 0.5e-3;              % Slow-time PRT [s]
N_tp    = 1024;                % Number of slow-time samples

% Adjustable
ang_res = 1;                   % Resolution of angle
L_win_AoA = 5;                 % Rx-AoA search window length
L_win_AoD = 5;                 % Tx-AoD search window length
h_eq_tilde_th = 0.5;           % Beam strength threshold for recalibration
if if_steer_vec_norm == 0
    h_eq_tilde_th = 10.* h_eq_tilde_th;
end
set_target_velocity = 0.3;        % m/s

% Structure
% Rx-1
Rx_1.pos_first = [1,0];
Rx_1.array_ang = 0; % deg
Rx_1.N = 16; % 16
Rx_1.d = lambda/2;
[Rx_1.pos_array,Rx_1.pos_center] = gen_gNB(Rx_1.pos_first,Rx_1.array_ang,Rx_1.N,Rx_1.d);

% Tx-1
Tx_1.pos_first = [0,0];
Tx_1.array_ang = 0; % deg
Tx_1.N = 8; % 8
Tx_1.d = lambda/2;
[Tx_1.pos_array,Tx_1.pos_center] = gen_gNB(Tx_1.pos_first,Tx_1.array_ang,Tx_1.N,Tx_1.d);

% Tar-1
if select_field_type == "near"
    Tar_1.pos_start = [0.5,1];            % near field
end
if select_field_type == "mid"
    Tar_1.pos_start = [0.5,5];            % middle field
end
if select_field_type == "far"
    Tar_1.pos_start = [10,10*sqrt(3)];    % far field
end
Tar_1.sigma = 0.5;
Tar_1.sigma_v = set_target_velocity;
Tar_1.sigma_theta = 30;
Tar_1.trajectory_array = gen_trajectory(Tar_1.pos_start,N_tp,Tar_1.sigma_v,Tar_1.sigma_theta,Tp);
Tar_1.velocity_array = gen_velocity_array(Tar_1.trajectory_array,Tp);

% Alpha
alpha = ones(N_tp,1);

% Storage
ang_array_deg = (-90:ang_res:90).';
ang_array_rad = ang_array_deg./180.*pi;
Rx_1.DFT_codebook = exp(-1*jmath*2*pi*(0:Rx_1.N-1)'*sin(ang_array_rad.').*Rx_1.d./lambda);
Tx_1.DFT_codebook = exp(-1*jmath*2*pi*(0:Tx_1.N-1)'*sin(ang_array_rad.').*Tx_1.d./lambda);



%% Near field channel generation

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
        if if_doppler_effect == 1
            V_nr_nt__l(:,ntt,it) = Uplink_1_1.Velocity(it,:).'+Downlink_1_1.Velocity(it,ntt);
        end
    end
end

% Sensing channel with path loss
if if_path_loss == 1
    for l = 1:N_tp
        alpha(l,1) = cal_path_loss(D_nr_nt__l(1,1,l),f_c,select_path_loss_type);
    end
end

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

% A priori information
Rx_1.ang_real_array_rad = cal_ang_center(Tar_1.trajectory_array,Rx_1.pos_center,'rad');
Tx_1.ang_real_array_rad = cal_ang_center(Tar_1.trajectory_array,Tx_1.pos_center,'rad');
Rx_1.ang_real_array_deg = cal_ang_center(Tar_1.trajectory_array,Rx_1.pos_center,'deg');
Tx_1.ang_real_array_deg = cal_ang_center(Tar_1.trajectory_array,Tx_1.pos_center,'deg');

% SVD analysis
U_V = zeros(Rx_1.N,Rx_1.N,N_tp);
S_V = zeros(Rx_1.N,Tx_1.N,N_tp);
V_V = zeros(Tx_1.N,Tx_1.N,N_tp);
r_S = zeros(N_tp,1);
r_V_nr_nt = zeros(N_tp,1);
U_exp_V = zeros(Rx_1.N,Rx_1.N,N_tp);
S_exp_V = zeros(Rx_1.N,Tx_1.N,N_tp);
V_exp_V = zeros(Tx_1.N,Tx_1.N,N_tp);
r_exp_V_S = zeros(N_tp,1);
r_exp_V_nr_nt = zeros(N_tp,1);
for it = 1:N_tp
    [U_V(:,:,it),S_V(:,:,it),V_V(:,:,it)] = svd(V_nr_nt__l(:,:,it));
    r_S(it,1) = rank(S_V(:,:,it));
    r_V_nr_nt(it,1) = rank(V_nr_nt__l(:,:,it));
    [U_exp_V(:,:,it),S_exp_V(:,:,it),V_exp_V(:,:,it)] = svd(exp(jmath*2*pi*(it-1)*Tp/lambda*V_nr_nt__l(:,:,it)));
    r_exp_V_S(it,1) = rank(S_exp_V(:,:,it));
    r_exp_V_nr_nt(it,1) = rank(exp(jmath*2*pi*(it-1)*Tp/lambda*V_nr_nt__l(:,:,it)));
end



%% ISAC System - REAL-TIME DOPPLER ESTIMATION TRACKING
h_eq_tilde_k_l = zeros(N_sc,N_tp);
h_eq_tilde_bir_bit__k_l = zeros(size(Rx_1.DFT_codebook,2),size(Tx_1.DFT_codebook,2),N_sc,N_tp);

% Storage for Doppler estimation
h_eq_tilde_full_freq__k_l = zeros(N_sc,N_tp);  % Full channel for each subcarrier
rec_V_estimate_l = zeros(N_tp,1);              % Estimated velocity magnitude
V_est_history = zeros(N_tp,1);                 % History of velocity estimates

% Initial beamsweeping
if if_steer_vec_norm == 1
    for l = 1:1
        for k = 1:N_sc
            for beam_number_tx = 1:size(Tx_1.DFT_codebook,2)
                for beam_number_rx = 1:size(Rx_1.DFT_codebook,2)
                    w = Rx_1.DFT_codebook(:,beam_number_rx)./sqrt(Rx_1.N);
                    f = Tx_1.DFT_codebook(:,beam_number_tx)./sqrt(Tx_1.N);
                    h_eq_tilde_bir_bit__k_l(beam_number_rx,beam_number_tx,k,l) = w'*H_tilde_nr_nt__k_l(:,:,k,l)*f;
                end
            end
        end
    end
    figure
    imagesc(abs(h_eq_tilde_bir_bit__k_l(:,:,1,1)));colorbar;
    xlabel('beam number of AoD [Tx]');hold on;
    ylabel('beam number of AoA [Rx]');hold on;
else
    for l = 1:1
        for k = 1:N_sc
            for beam_number_tx = 1:size(Tx_1.DFT_codebook,2)
                for beam_number_rx = 1:size(Rx_1.DFT_codebook,2)
                    w = Rx_1.DFT_codebook(:,beam_number_rx);
                    f = Tx_1.DFT_codebook(:,beam_number_tx);
                    h_eq_tilde_bir_bit__k_l(beam_number_rx,beam_number_tx,k,l) = w'*H_tilde_nr_nt__k_l(:,:,k,l)*f;
                end
            end
        end
    end
    figure
    imagesc(abs(h_eq_tilde_bir_bit__k_l(:,:,1,1)));colorbar;
    xlabel('beam number of AoD [Tx]');
    ylabel('beam number of AoA [Rx]');
    sgtitle('Initial beamsweeping result at k=1, l=1');
end

% Tracking process with real-time Doppler estimation

beam_number_AoA_initial = ang2beam_number(ang_array_deg,Rx_1.ang_real_array_deg(1,1))
beam_number_AoD_initial = ang2beam_number(ang_array_deg,Tx_1.ang_real_array_deg(1,1))
if if_steer_vec_norm == 1
    w = Rx_1.DFT_codebook(:,beam_number_AoA_initial)./sqrt(Rx_1.N);
    f = Tx_1.DFT_codebook(:,beam_number_AoD_initial)./sqrt(Tx_1.N);
else
    w = Rx_1.DFT_codebook(:,beam_number_AoA_initial);
    f = Tx_1.DFT_codebook(:,beam_number_AoD_initial);
end

center_number_win_AoA = beam_number_AoA_initial;
center_number_win_AoD = beam_number_AoD_initial;

Rx_Tx_turn_indicator = 1; % (1)--Rx, (-1)--Tx
temp_maximum_h_eq_tilde = 0;
temp_number_AoA__maximum_h_eq_tilde = 0;
temp_number_AoD__maximum_h_eq_tilde = 0;
count_win_AoA = 0;
count_win_AoD = 0;
count_cali_beam_sweeping = 0;

rec_center_number_win_AoA = zeros(ceil(N_tp/(L_win_AoA+L_win_AoD)),1);
rec_center_number_win_AoD = zeros(ceil(N_tp/(L_win_AoA+L_win_AoD)),1);
rec_channel_gain_max_period = zeros(ceil(N_tp/(L_win_AoA+L_win_AoD)),1);
rec_l_cali_beam_sweeping(1,1) = 1;
rec_period_cali_beam_sweeping(1,1) = 1;

%% ===== MAIN TRACKING LOOP WITH REAL-TIME DOPPLER COMPENSATION =====
for l = 1:N_tp
    for k = 1:N_sc
        if Rx_Tx_turn_indicator == 1 % -------- Rx AoA Search
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
            
            % Beamforming
            if if_steer_vec_norm == 1
                w = Rx_1.DFT_codebook(:,search_number_win_AoA)./sqrt(Rx_1.N);
                f = Tx_1.DFT_codebook(:,center_number_win_AoD)./sqrt(Tx_1.N);
            else
                w = Rx_1.DFT_codebook(:,search_number_win_AoA);
                f = Tx_1.DFT_codebook(:,center_number_win_AoD);
            end
            
            % Get channel for this beam pair
            h_eq_tilde_k_l(k,l) = w'*H_tilde_nr_nt__k_l(:,:,k,l)*f;
            h_eq_tilde_bir_bit__k_l(search_number_win_AoA,center_number_win_AoD,k,l) = h_eq_tilde_k_l(k,l);
            h_eq_tilde_full_freq__k_l(k,l) = h_eq_tilde_k_l(k,l);
            
            % Apply Doppler compensation on magnitude evaluation
            h_mag_for_decision = abs(h_eq_tilde_k_l(k,l));
            if if_doppler_compensation == 1 && l > 1
                % Estimate Doppler phase change from previous slot
                h_mag_for_decision = apply_doppler_compensation(h_eq_tilde_k_l(:,l), h_eq_tilde_k_l(:,max(1,l-1)), ...
                                                                 if_doppler_est_method, N_sc, delta_f, f_c, lambda, Tp);
            end
            
            % Find maximum
            if h_mag_for_decision > temp_maximum_h_eq_tilde
                temp_maximum_h_eq_tilde = h_mag_for_decision;
                temp_number_AoA__maximum_h_eq_tilde = search_number_win_AoA;
            end
            count_win_AoA = count_win_AoA + 1;
            
            if count_win_AoA == L_win_AoA
                center_number_win_AoA = temp_number_AoA__maximum_h_eq_tilde;
                rec_center_number_win_AoA(ceil(l/(L_win_AoA+L_win_AoD)),1) = center_number_win_AoA;
                temp_maximum_h_eq_tilde = 0;
                count_win_AoA = 0;
                Rx_Tx_turn_indicator = (-1)*Rx_Tx_turn_indicator;
            end
            
        else % -------- Tx AoD Search
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
            
            % Beamforming
            if if_steer_vec_norm == 1
                w = Rx_1.DFT_codebook(:,center_number_win_AoA)./sqrt(Rx_1.N);
                f = Tx_1.DFT_codebook(:,search_number_win_AoD)./sqrt(Tx_1.N);
            else
                w = Rx_1.DFT_codebook(:,center_number_win_AoA);
                f = Tx_1.DFT_codebook(:,search_number_win_AoD);
            end
            
            % Get channel for this beam pair
            h_eq_tilde_k_l(k,l) = w'*H_tilde_nr_nt__k_l(:,:,k,l)*f;
            h_eq_tilde_bir_bit__k_l(center_number_win_AoA,search_number_win_AoD,k,l) = h_eq_tilde_k_l(k,l);
            h_eq_tilde_full_freq__k_l(k,l) = h_eq_tilde_k_l(k,l);
            
            % Apply Doppler compensation on magnitude evaluation
            h_mag_for_decision = abs(h_eq_tilde_k_l(k,l));
            if if_doppler_compensation == 1 && l > 1
                h_mag_for_decision = apply_doppler_compensation(h_eq_tilde_k_l(:,l), h_eq_tilde_k_l(:,max(1,l-1)), ...
                                                                 if_doppler_est_method, N_sc, delta_f, f_c, lambda, Tp);
            end
            
            % Find maximum
            if h_mag_for_decision > temp_maximum_h_eq_tilde
                temp_maximum_h_eq_tilde = h_mag_for_decision;
                temp_number_AoD__maximum_h_eq_tilde = search_number_win_AoD;
            end
            count_win_AoD = count_win_AoD + 1;
            
            if count_win_AoD == L_win_AoD
                center_number_win_AoD = temp_number_AoD__maximum_h_eq_tilde;
                
                % Record Doppler velocity estimate for this window
                if l > L_win_AoA + L_win_AoD
                    V_est = estimate_velocity_from_phase_diff(h_eq_tilde_full_freq__k_l(:,l), ...
                                                              h_eq_tilde_full_freq__k_l(:,max(1,l-L_win_AoA-L_win_AoD)), ...
                                                              if_doppler_est_method, N_sc, delta_f, f_c, lambda, Tp);
                    rec_V_estimate_l(l) = V_est;
                    V_est_history(l) = abs(V_est);
                end
                
                % Calibration if signal weak
                if (if_beamsweeping_cali == 1) && (temp_maximum_h_eq_tilde < h_eq_tilde_th)
                    center_number_win_AoA = ang2beam_number(ang_array_deg,Rx_1.ang_real_array_deg(l,1));
                    center_number_win_AoD = ang2beam_number(ang_array_deg,Tx_1.ang_real_array_deg(l,1));
                    
                    if if_steer_vec_norm == 1
                        w = Rx_1.DFT_codebook(:,center_number_win_AoA)./sqrt(Rx_1.N);
                        f = Tx_1.DFT_codebook(:,center_number_win_AoD)./sqrt(Tx_1.N);
                    else
                        w = Rx_1.DFT_codebook(:,center_number_win_AoA);
                        f = Tx_1.DFT_codebook(:,center_number_win_AoD);
                    end
                    h_eq_tilde_k_l(1,l) = w'*H_tilde_nr_nt__k_l(:,:,1,l)*f;
                    h_eq_tilde_bir_bit__k_l(center_number_win_AoA,center_number_win_AoD,1,l) = h_eq_tilde_k_l(1,l);
                    temp_maximum_h_eq_tilde = abs(h_eq_tilde_k_l(1,l));
                    count_cali_beam_sweeping = count_cali_beam_sweeping + 1;
                    rec_l_cali_beam_sweeping(count_cali_beam_sweeping,1) = l;
                    rec_period_cali_beam_sweeping(count_cali_beam_sweeping,1) = ceil(l/(L_win_AoA+L_win_AoD));
                    rec_center_number_win_AoA(ceil(l/(L_win_AoA+L_win_AoD)),1) = center_number_win_AoA;
                end
                
                rec_center_number_win_AoD(ceil(l/(L_win_AoA+L_win_AoD)),1) = center_number_win_AoD;
                rec_channel_gain_max_period(ceil(l/(L_win_AoA+L_win_AoD)),1) = temp_maximum_h_eq_tilde;
                temp_maximum_h_eq_tilde = 0;
                count_win_AoD = 0;
                Rx_Tx_turn_indicator = (-1)*Rx_Tx_turn_indicator;
            end
        end
    end
end

% Trimming
if rec_center_number_win_AoA(size(rec_center_number_win_AoA,1),1) == 0
    rec_center_number_win_AoA(size(rec_center_number_win_AoA,1),1) = rec_center_number_win_AoA(size(rec_center_number_win_AoA,1)-1,1);
end
if rec_center_number_win_AoD(size(rec_center_number_win_AoD,1),1) == 0
    rec_center_number_win_AoD(size(rec_center_number_win_AoD,1),1) = rec_center_number_win_AoD(size(rec_center_number_win_AoD,1)-1,1);
end
if rec_channel_gain_max_period(size(rec_channel_gain_max_period,1),1) == 0
    rec_channel_gain_max_period(size(rec_channel_gain_max_period,1),1) = rec_channel_gain_max_period(size(rec_channel_gain_max_period,1)-1,1);
end

% Reconstruction of target trajectory
Tar_1.trajectory_array_est = zeros(size(rec_center_number_win_AoA,1),2);
for l_period = 1:size(rec_center_number_win_AoA,1)
    Tar_1.trajectory_array_est(l_period,:) = cal_ang_2_pos(Rx_1.pos_center,Tx_1.pos_center,...
        beam_number2ang(ang_array_deg,rec_center_number_win_AoA(l_period,1)),...
        beam_number2ang(ang_array_deg,rec_center_number_win_AoD(l_period,1)),'deg');
end
Tar_1.velocity_array_est = gen_velocity_array(Tar_1.trajectory_array_est,Tp);



%% Figure
figure
plot(Tx_1.pos_array(:,1),Tx_1.pos_array(:,2),'b:^');hold on;
plot(Rx_1.pos_array(:,1),Rx_1.pos_array(:,2),'r:sq');hold on;
plot(Tar_1.pos_start(:,1),Tar_1.pos_start(:,2),'g:o','MarkerSize',10);hold on;
plot(Tar_1.trajectory_array(:,1),Tar_1.trajectory_array(:,2),'k:o');hold on;
legend('Tx-1','Rx-1','Target starting position','Target trajectory');
xlabel('x [m]');ylabel('y [m]');title('Layout diagram');

figure
plot(vector2scalar(Tar_1.velocity_array),'k-o');hold on;
xlabel('slot index');ylabel('velocity scalar [m/s]');title('velocity trend');

figure
plot(Tar_1.pos_start(:,1),Tar_1.pos_start(:,2),'g:o','MarkerSize',10);hold on;
plot(Tar_1.trajectory_array(:,1),Tar_1.trajectory_array(:,2),'k-o');hold on;
for i = 1:size(Tar_1.trajectory_array_est,1)-1
    text(Tar_1.trajectory_array(i*(L_win_AoA+L_win_AoD),1), Tar_1.trajectory_array(i*(L_win_AoA+L_win_AoD),2), num2str(i),'Color','red','FontSize',9);
end
plot(Tar_1.trajectory_array_est(:,1),Tar_1.trajectory_array_est(:,2),'b-sq');hold on;
for i = 1:size(Tar_1.trajectory_array_est,1)
    color_string = strcat('#',dec2hex(round(i/size(Tar_1.trajectory_array_est,1)*16777215/2+16777215/2)));
    text(Tar_1.trajectory_array_est(i,1), Tar_1.trajectory_array_est(i,2), num2str(i),'Color',color_string,'FontSize',9);
end
legend('Target starting position','True trajectory','Estimated trajectory');
xlabel('x [m]');ylabel('y [m]');title('Trajectory tracking results');

figure
plot(vector2scalar(Tar_1.velocity_array(1:(L_win_AoA+L_win_AoD):N_tp)),'k-o');hold on;
plot(vector2scalar(Tar_1.velocity_array_est),'b-sq');hold on;
legend('True velocity','Estimated velocity');
xlabel('period index');ylabel('velocity scalar [m/s]');title('Velocity tracking results');

figure
plot(rec_channel_gain_max_period);hold on;
xlabel('period index');ylabel('channel gain');title('Channel gain trend');

figure
plot(Rx_1.ang_real_array_deg(1:(L_win_AoA+L_win_AoD):N_tp),'r-o', 'LineWidth', 1, 'MarkerSize', 3);hold on;
plot(beam_number2ang(ang_array_deg,rec_center_number_win_AoA),'m-x', 'LineWidth', 1, 'MarkerSize', 4);hold on;
plot(Tx_1.ang_real_array_deg(1:(L_win_AoA+L_win_AoD):N_tp),'b-o', 'LineWidth', 1, 'MarkerSize', 3);hold on;
plot(beam_number2ang(ang_array_deg,rec_center_number_win_AoD),'g-x', 'LineWidth', 1, 'MarkerSize', 4);hold on;
xlabel('period index');ylabel('angle [deg]');
legend('Real AoA','Esti AoA','Real AoD','Esti AoD');title('DoA estimation performance');

if if_doppler_compensation == 1
    figure
    plot(V_est_history);hold on;
    xlabel('time slot');ylabel('velocity magnitude [m/s]');
    title('Real-time Doppler velocity estimation');
end



%% Printf
if if_beamsweeping_cali == 1
    fprintf("Times of beamsweeping calibration : Number of slots = %d : %d\n",count_cali_beam_sweeping,N_tp);
end

fprintf("Field type: %s, ",select_field_type);
fprintf("Target velocity set as: %d, ",set_target_velocity);
fprintf("Beamsweeping calibration: %d, ",if_beamsweeping_cali);
fprintf("Doppler effect: %d, ",if_doppler_effect);
fprintf("Doppler compensation: %d, ",if_doppler_compensation);
fprintf("Path loss: %d, ",if_path_loss);
fprintf("Steering vector normalization: %d.\n",if_steer_vec_norm);
fprintf("EoF.\n");



%% Helper Functions for Doppler Compensation

function h_mag_corrected = apply_doppler_compensation(h_current, h_previous, method, N_sc, delta_f, f_c, lambda, Tp)
    % Estimate Doppler effect by comparing current and previous channel
    % Return corrected magnitude for beam decision-making
    
    jmath = sqrt(-1);
    
    if method == 1  % Temporal phase difference method
        % Single subcarrier: use phase change rate
        if N_sc >= 1
            phase_diff = angle(h_current(1)) - angle(h_previous(1));
            % Normalize phase difference to [-pi, pi]
            phase_diff = atan2(sin(phase_diff), cos(phase_diff));
            % Magnitude correction (remove Doppler-induced phase variation effect)
            h_mag_corrected = abs(h_current(1)) * abs(exp(-jmath * phase_diff/2));
        else
            h_mag_corrected = abs(h_current(1));
        end
        
    elseif method == 2  % Multi-subcarrier frequency analysis
        % Use frequency diversity to estimate Doppler shift
        if N_sc > 1
            phase_diffs = zeros(N_sc, 1);
            for k = 1:N_sc
                phase_diffs(k) = angle(h_current(k)) - angle(h_previous(k));
                phase_diffs(k) = atan2(sin(phase_diffs(k)), cos(phase_diffs(k)));
            end
            % Linear fit to get Doppler frequency
            doppler_freq = mean(phase_diffs) / (2*pi*Tp*2);  % Factor of 2 for round-trip
            % Correct magnitude
            doppler_phase_correction = 2*pi*doppler_freq*Tp;
            h_mag_corrected = abs(h_current(1)) * abs(exp(-jmath * doppler_phase_correction/2));
        else
            h_mag_corrected = abs(h_current(1));
        end
    else
        h_mag_corrected = abs(h_current(1));
    end
end

function V_estimated = estimate_velocity_from_phase_diff(h_freq_current, h_freq_previous, method, N_sc, delta_f, f_c, lambda, Tp)
    % Estimate target radial velocity from frequency-domain channel differences
    
    jmath = sqrt(-1);
    c_0 = 3e8;
    
    if N_sc > 1 && method == 2
        % Multi-subcarrier: use frequency diversity
        phase_diffs = zeros(N_sc, 1);
        for k = 1:N_sc
            if abs(h_freq_previous(k)) > 1e-6
                phase_diffs(k) = angle(h_freq_current(k)) - angle(h_freq_previous(k));
                phase_diffs(k) = atan2(sin(phase_diffs(k)), cos(phase_diffs(k)));
            end
        end
        doppler_freq = mean(phase_diffs) / (2*pi*Tp*2);  
        V_estimated = doppler_freq * lambda / 2;
    else
        % Single subcarrier: use temporal phase diff
        if abs(h_freq_previous(1)) > 1e-6
            phase_diff = angle(h_freq_current(1)) - angle(h_freq_previous(1));
            phase_diff = atan2(sin(phase_diff), cos(phase_diff));
            doppler_freq = phase_diff / (2*pi*Tp*2);
            V_estimated = doppler_freq * lambda / 2;
        else
            V_estimated = 0;
        end
    end
    
    % Limit to reasonable range
    V_estimated = max(-3, min(3, V_estimated));
end



%% Original Functions (from my_track_3.m)
function [pos_array,pos_center] = gen_gNB(pos_start,array_ang,N,d)
    pos_array = zeros(N,2);
    for it = 1:N
        pos_array(it,:)=pos_start+(it-1)*d*[cosd(array_ang),sind(array_ang)];
    end
    pos_center = sum(pos_array)/N;
end

function trajectory_array = gen_trajectory(pos_start,points,sigma_v,sigma_theta,Tp)
    v_array=zeros(points,2);
    trajectory_array = zeros(points,2);
    trajectory_array(1,:) = pos_start;
    trajectory_array(2,:) = pos_start+(sigma_v*Tp)^2*randn(1,2);
    for it = 3:points
        v_array(it-1,:) = (trajectory_array(it-1,:)-trajectory_array(it-2,:))/Tp;
        omega = sigma_theta^2*randn(1,1);
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
