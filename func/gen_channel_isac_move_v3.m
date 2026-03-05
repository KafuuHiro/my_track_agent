function [H_ul_total,H_dl_total,cal_perf,timing_error,G_hatAB,struct_Sense] = gen_channel_isac_move_v3(fc,nSC,M_ofdm,SCS,TA_interval,N_Ax,N_Bx,N_scatter,Pt_db)
% author : jqj
% date   : 4,Dec,204
% version1: 
% note   : 1. this versio is used to check process correct
%          2. add simulate of scatter location to generate true OTA channel
%          3. add comm. path loss and  radar path loss
%          4. add random of rf mismatch and Tx,Rx,Scatter location
%          5. To sim cell-free snr, must fix noise and change Power
%          6. To make Pt_db tansmit at Pt_db_eq, (to match jp_rzf_dl_mu_mimo_with_ac,noise = 1, 0dbm)
%          7. Beside 5&6, I do not normalize multi Tx-Scatter-Rx path (Pt<->All path),RZF do normalize
%          8. Vmax = 0 to simulate static channel
%          9. Restrict Vmax at [1~2]; 
%          10. Here we suppose scatter is  static but ue is moving!

%% system param
% fc = 26e9;          % carrier frequency
% nSC = 2048;         % number of subcarrier
% M_ofdm = 100;
% SCS = 120e3;    % subcarrier space
% TA_interval = 0.625e-3;      % signal repeat period


c0 = 3e8;
lambda = c0/fc;     % lamda
nComb = 1;          % gap of subcarrier, to do
nfft = 2048;        % nfft 
Tofdm = 1/SCS;  % ofdm

Radius = 100;
Vmax = 3;
Vmin = 4;
Rcs = 1;   % 1e-2 may cause scatter small


% N_Ax = 4;  % Can be seen as bs
% N_Bx = 2;  % Can be seen as ue
% N_scatter = 2;

bw = nSC*SCS;
n0 = -174;        % -174dbm/Hz
noise_db = 10*log10(bw) +  n0; 
% noise_pn = 10.^(noise_db/10);

% Pt_db = 0; % 0dbm
Pt_db_eq = Pt_db - noise_db;
Pt = 10.^(Pt_db_eq/10); % 
noise_pn = 10.^(0/10);  % note,To do so is to make normilize noise to match jp_rzf_dl_mu_mimo_with_ac    


%% scatter geometry, which used to generate necessary sensing params
N_path = N_scatter + 1; % if consider comm. + 1
[ap_pos, user_pos, scatter_pos]  = set_isac_bg(N_Ax,N_Bx,N_scatter,Radius);
TxRRU_pos = ap_pos;
RxRRU_pos = user_pos;
Struct_scatter = repmat(struct(), N_scatter, 1);
for is = 1:N_scatter
    Struct_scatter(is).pos = scatter_pos(is,:);
    % Struct_scatter(is).vec = rand(1, 2) * Vmax * 2 - Vmax;   
    % Struct_scatter(is).vec = [generateRandVec(Vmin,Vmax),generateRandVec(Vmin,Vmax)];
    Struct_scatter(is).vec = [0,0];

end

% 感知结构体-（1）
struct_Sense.ap_pos = ap_pos;
struct_Sense.user_pos = user_pos;
struct_Sense.Struct_scatter = Struct_scatter;


range_list = zeros(N_Bx,N_Ax,N_path);
vec_list = zeros(N_Bx,N_Ax,N_path);
alpha_list = zeros(N_Bx,N_Ax,N_path); 
for is = 1:N_scatter
    Obj_scatter = Struct_scatter(is);
    Obj_pos = Obj_scatter.pos;
    Obj_vec = Obj_scatter.vec;
    [range_u,velocity_u] = cal_dis_velocity(RxRRU_pos,Obj_pos,Obj_vec); % TxAP_obj
    [range_d,velocity_d] = cal_dis_velocity(TxRRU_pos,Obj_pos,Obj_vec); % RxAP_obj  
    for irx = 1:N_Bx
        for itx = 1:N_Ax 
            range = range_u(irx)+range_d(itx); % meter
            velocity = velocity_u(irx)+velocity_d(itx);  % m/s doppler = velocity/lambda;
            range_list(irx,itx,is) = range;
            vec_list(irx,itx,is) = velocity;
            alpha_list(irx,itx,is) = fbiradarloss(Pt,1,1,lambda,range_u(irx),range_d(itx),Rcs);
        end
    end
end
% if consider self-interfere , Comm
UE_vec = [generateRandVec(Vmin,Vmax),generateRandVec(Vmin,Vmax)];
for irx = 1:N_Bx
    for itx = 1:N_Ax 
        UE_pos = TxRRU_pos(itx,:);
        [range,velocity] = cal_dis_velocity(RxRRU_pos(irx,:),UE_pos,UE_vec); % RxAP_obj  
        % range = norm(RxRRU_pos(irx,:)-TxRRU_pos(itx,:),2); % meter
        % velocity = 0;  % m/s doppler = velocity/lambda;
        range_list(irx,itx,N_scatter+1) = range;
        vec_list(irx,itx,N_scatter+1) = velocity;
        alpha_list(irx,itx,N_scatter+1) = fcommloss(Pt,1,1,lambda,range);
    end
end

tau_list = range_list/c0;
fd_list = vec_list/c0*fc;


%% Generate H , G_B<-A , G_B<-A
%% Generate OTA h, it can be seen as sensing h
Hota_ce_slot_rx_tx_path = zeros(nSC,M_ofdm,N_Bx,N_Ax,N_path);
for iL = 1:N_path
    for isymbol = 1: M_ofdm
        for ibx = 1:N_Bx   
            for iax = 1:N_Ax 
                alpha = alpha_list(ibx,iax,iL);
                range = range_list(ibx,iax,iL);
                velocity = vec_list(ibx,iax,iL);
                Hota_ce_slot_rx_tx_path(:,isymbol,ibx,iax,iL) = alpha.*...
                                                                exp(-1j*2*pi*(0:nSC-1)*SCS*range/c0).*...
                                                                exp(1j*2*pi*(velocity/c0*fc)*(isymbol-1)*TA_interval);
            end
        end
    end
end
Hota_isac_ce_slot_rx_tx = sum(Hota_ce_slot_rx_tx_path,5);
HBA_dlota_isac_ce_slot_rx_tx = Hota_isac_ce_slot_rx_tx;
HAB_ulota_isac_ce_slot_rx_tx = permute(Hota_isac_ce_slot_rx_tx,[1,2,4,3]);

Hota_sense_ce_slot_rx_tx = sum(Hota_ce_slot_rx_tx_path(:,:,:,:,1:end-1),5);
HBA_dlota_sense_ce_slot_rx_tx = Hota_sense_ce_slot_rx_tx;
HAB_ulota_sense_ce_slot_rx_tx = permute(HBA_dlota_sense_ce_slot_rx_tx,[1,2,4,3]);


% 感知结构体-（2）!!
struct_Sense.HAB_ulota_isac_ce_slot_rx_tx = HAB_ulota_isac_ce_slot_rx_tx;   % 通感一体化理想空口信道
struct_Sense.HAB_ulota_sense_ce_slot_rx_tx = HAB_ulota_sense_ce_slot_rx_tx; % 单纯感知理想空口信道

% % % check if right
% % dt = TA_interval;
% % df = SCS;
% % N_tp_process = 100;
% % Nf = nSC;
% % dis_list = 1/df*c0*(0:Nf-1)/(Nf);
% % vec_list = 1/(dt*N_tp_process)/(fc)*c0*(-N_tp_process/2:N_tp_process/2-1);
% % % Div = H_ulota_isac_ce_slot_rx_tx_AssistedUe(:,1:N_tp_process,5,1);
% % Div = HAB_ulota_sense_ce_slot_rx_tx(:,1:N_tp_process,12,2);
% % Div_IFFT_dd = ifft(Div,[],1); 
% % Div_IFFT_rd = fftshift(fft(Div_IFFT_dd, [], 2),2); 
% % Div_IFFT_rd_db = 10*log10(abs(Div_IFFT_rd./max(abs(Div_IFFT_rd(:)))));
% % figure;
% % imagesc(vec_list,dis_list,Div_IFFT_rd_db,[-20,0]);colorbar;
% % xlabel('速度(m/s)');ylabel('距离/m');
% % ylim([0,300]);xlim([-10,10]);


%% Generate RF non-ideal factor 
maxTimingerr =  10e-9; % 10ns asy err
% maxEerr =  0.005e-6;     % 0.005 ppm err, incluence quiet big
maxEerr =  30e-9;

% RF mismatch model, which is from ~
% % Beta_At = zeros(N_Ax,1);
% % Beta_Ar = conj(Beta_At);
% % Beta_Bt = zeros(N_Bx,1);
% % Beta_Br = conj(Beta_Bt);

% contant for log-normal  log10 to ln
contant_lognormal = (log(10)/10)^2;
% variance of the natural logarithm of the amplitude of the RF gain of AP
var_t_BS = contant_lognormal .* 2;
var_r_BS = contant_lognormal .* 2;
var_t_UE = contant_lognormal .* 2;
var_r_UE = contant_lognormal .* 2;
% 收发端天线的初始相位和幅度，幅度遵循对数正态分布
tx_init_phase = rand(N_Ax,1)*2*pi;  % rf mismatch
rx_init_phase = rand(N_Ax,1)*2*pi;  % rf mismatch
tx_init_amp = exp( sqrt(var_t_BS) * randn(N_Ax,1) );  % rf mismatch
rx_init_amp = exp( sqrt(var_r_BS) * randn(N_Ax,1) );  % rf mismatch
Beta_At = tx_init_amp.*exp(1j*tx_init_phase);
Beta_Ar = rx_init_amp.*exp(1j*rx_init_phase);
% 收发端天线的初始相位和幅度，幅度遵循对数正态分布
tx_init_phase = rand(N_Bx,1)*2*pi;  % rf mismatch
rx_init_phase = rand(N_Bx,1)*2*pi;  % rf mismatch
tx_init_amp = exp( sqrt(var_t_UE) * randn(N_Bx,1) );  % rf mismatch
rx_init_amp = exp( sqrt(var_r_UE) * randn(N_Bx,1) );  % rf mismatch
Beta_Bt = tx_init_amp.*exp(1j*tx_init_phase);
Beta_Br = rx_init_amp.*exp(1j*rx_init_phase);


% asy model
tau_At = (2*rand(N_Ax,1)-1)*maxTimingerr;        
e_At =  (2*rand(N_Ax,1)-1)*maxEerr;      
phi_At  =  (2*rand(N_Ax,1)-1)*pi;
tau_Ar = tau_At;
e_Ar = e_At;
phi_Ar = phi_At;
tau_Bt = (2*rand(N_Bx,1)-1)*maxTimingerr;    
e_Bt = (2*rand(N_Bx,1)-1)*maxEerr;      
phi_Bt  = (2*rand(N_Bx,1)-1)*pi;
tau_Br = tau_Bt;
e_Br = e_Bt;
phi_Br = phi_Bt;


%% thus,which is important, To generae relative err~
tau_BA = zeros(N_Bx,N_Ax);
e_BA = zeros(N_Bx,N_Ax);
phi_BA = zeros(N_Bx,N_Ax);
for irx = 1:N_Bx   
    for itx = 1:N_Ax 
        tau_BA(irx,itx) = tau_At(itx) - tau_Br(irx);
        e_BA(irx,itx) = e_At(itx) - e_Br(irx);
        phi_BA(irx,itx) = phi_At(itx) - phi_Br(irx);
    end
end

tau_AB = zeros(N_Ax,N_Bx);
e_AB = zeros(N_Ax,N_Bx);
phi_AB = zeros(N_Ax,N_Bx);
for irx = 1:N_Ax   
    for itx = 1:N_Bx 
        tau_AB(irx,itx) = tau_Bt(itx) - tau_Ar(irx);
        e_AB(irx,itx) = e_Bt(itx) - e_Ar(irx);
        phi_AB(irx,itx) = phi_Bt(itx) - phi_Ar(irx);
    end
end


%% generate ideal rf mismatch and tf asy
H_At =  zeros(nSC,M_ofdm,N_Ax);
H_Ar =  zeros(nSC,M_ofdm,N_Ax);
H_Bt =  zeros(nSC,M_ofdm,N_Bx);
H_Br =  zeros(nSC,M_ofdm,N_Bx);
C_At = H_At;
C_Ar = H_Ar;
C_Bt = H_Bt;
C_Br = H_Br;
for isymbol = 1: M_ofdm
    for iax = 1:N_Ax 
        tau_t = tau_At(iax);tau_r = tau_Ar(iax);
        e_t = e_At(iax);e_r = e_Ar(iax);
        phi_t = phi_At(iax);phi_r = phi_Ar(iax);
        beta_t = Beta_At(iax);beta_r = Beta_Ar(iax);
        H_At(:,isymbol,iax) =  exp(1j*(phi_t)).*exp(-1j*2*pi*(0:nSC-1)*SCS*tau_t).*exp(1j*2*pi*(e_t*fc)*(isymbol-1)*TA_interval);
        H_Ar(:,isymbol,iax) =  exp(-1j*(phi_r)).*exp(1j*2*pi*(0:nSC-1)*SCS*tau_r) .*exp(-1j*2*pi*(e_r*fc)*(isymbol-1)*TA_interval);
        C_At(:,isymbol,iax) = beta_t.*H_At(:,isymbol,iax);
        C_Ar(:,isymbol,iax) = beta_r.*H_Ar(:,isymbol,iax);
    end
end


for isymbol = 1: M_ofdm
    for ibx = 1:N_Bx  
        tau_t = tau_Bt(ibx);tau_r = tau_Br(ibx);
        e_t = e_Bt(ibx);e_r = e_Br(ibx);     
        phi_t = phi_Bt(ibx);phi_r = phi_Br(ibx);
        beta_t = Beta_Bt(ibx);beta_r = Beta_Br(ibx);
        H_Bt(:,isymbol,ibx) =  exp(1j*(phi_t)).*exp(-1j*2*pi*(0:nSC-1)*SCS*tau_t).*exp(1j*2*pi*(e_t*fc)*(isymbol-1)*TA_interval);
        H_Br(:,isymbol,ibx) =  exp(-1j*(phi_r)).*exp(1j*2*pi*(0:nSC-1)*SCS*tau_r) .*exp(-1j*2*pi*(e_r*fc)*(isymbol-1)*TA_interval);
        C_Bt(:,isymbol,ibx) = beta_t.*H_Bt(:,isymbol,ibx);
        C_Br(:,isymbol,ibx) = beta_r.*H_Br(:,isymbol,ibx);
    end
end

%% generate ul channel and dl channel
% ISAC 通信感知一体信道
HBA_dlisac_ce_slot_rx_tx = zeros(nSC,M_ofdm,N_Bx,N_Ax); % A发B收
HAB_ulisac_ce_slot_rx_tx = zeros(nSC,M_ofdm,N_Ax,N_Bx); % B发A收
HBA_dlisac_ce_slot_rx_tx_noise = zeros(nSC,M_ofdm,N_Bx,N_Ax); % A发B收
HAB_ulisac_ce_slot_rx_tx_noise = zeros(nSC,M_ofdm,N_Ax,N_Bx); % B发A收

% 单纯感知信道
HBA_dlsnese_ce_slot_rx_tx = zeros(nSC,M_ofdm,N_Bx,N_Ax); % A发B收
HAB_ulsense_ce_slot_rx_tx = zeros(nSC,M_ofdm,N_Ax,N_Bx); % B发A收
HBA_dlsense_ce_slot_rx_tx_noise = zeros(nSC,M_ofdm,N_Bx,N_Ax); % A发B收
HAB_ulsense_ce_slot_rx_tx_noise = zeros(nSC,M_ofdm,N_Ax,N_Bx); % B发A收

for irx = 1:N_Bx   
    for itx = 1:N_Ax 
        HBA_dlisac_ce_slot_rx_tx(:,:,irx,itx) = C_Br(:,:,irx).*HBA_dlota_isac_ce_slot_rx_tx(:,:,irx,itx).*C_At(:,:,itx);
        HBA_dlsnese_ce_slot_rx_tx(:,:,irx,itx) = C_Br(:,:,irx).*HBA_dlota_sense_ce_slot_rx_tx(:,:,irx,itx).*C_At(:,:,itx);
        % another way ,add noise
        noise = sqrt(1/2) * sqrt(noise_pn) * (randn(nSC,M_ofdm) + 1j*randn(nSC,M_ofdm));
        HBA_dlisac_ce_slot_rx_tx_noise(:,:,irx,itx) = C_Br(:,:,irx).*HBA_dlota_isac_ce_slot_rx_tx(:,:,irx,itx).*C_At(:,:,itx)...
            + noise;
        HBA_dlsense_ce_slot_rx_tx_noise(:,:,irx,itx) = C_Br(:,:,irx).*HBA_dlota_sense_ce_slot_rx_tx(:,:,irx,itx).*C_At(:,:,itx)...
            + noise;
    end
end

for irx = 1:N_Ax   
    for itx = 1:N_Bx 
        HAB_ulisac_ce_slot_rx_tx(:,:,irx,itx) = C_Ar(:,:,irx).*HAB_ulota_isac_ce_slot_rx_tx(:,:,irx,itx).*C_Bt(:,:,itx);
        HAB_ulsense_ce_slot_rx_tx(:,:,irx,itx) = C_Ar(:,:,irx).*HAB_ulota_sense_ce_slot_rx_tx(:,:,irx,itx).*C_Bt(:,:,itx);
        % another way ,add noise
        noise = sqrt(1/2) * sqrt(noise_pn) * (randn(nSC,M_ofdm) + 1j*randn(nSC,M_ofdm));
        HAB_ulisac_ce_slot_rx_tx_noise(:,:,irx,itx) = C_Ar(:,:,irx).*HAB_ulota_isac_ce_slot_rx_tx(:,:,irx,itx).*C_Bt(:,:,itx)...
            + noise;
        HAB_ulsense_ce_slot_rx_tx_noise(:,:,irx,itx) = C_Ar(:,:,irx).*HAB_ulota_sense_ce_slot_rx_tx(:,:,irx,itx).*C_Bt(:,:,itx)...
            + noise;
    end
end

% 感知结构体-（3）!!
% struct_Sense.HAB_ulsense_ce_slot_rx_tx_noise = HAB_ulsense_ce_slot_rx_tx;
struct_Sense.HAB_ulsense_ce_slot_rx_tx_noise = HAB_ulsense_ce_slot_rx_tx_noise;

% % % check if right
% % figure;hold on
% % bs_select = 4; ue_select = 1;
% % plot(unwrap(angle(squeeze(HAB_ulisac_ce_slot_rx_tx(1,:,bs_select,ue_select)))),'k-v')
% % plot(unwrap(angle(squeeze(HAB_ulota_isac_ce_slot_rx_tx(1,:,bs_select,ue_select)))),'r-^')
% % legend('srs track','ota track') % attention,,here we use HAB_ulota_isac_ce_slot_rx_tx or HAB_ulota_sense_ce_slot_rx_tx


%% simulate channel precode to get non-ideal factor 
% step1:get HBA
% step2:mulit ./H to get argos
G_hatAB = zeros(nSC,M_ofdm,N_Ax,N_Bx);
for isc = 1: nSC
    for isymbol = 1: M_ofdm
        hB = reshape(HBA_dlisac_ce_slot_rx_tx(isc,isymbol,:,:),N_Bx,N_Ax);
        hA = reshape(HAB_ulisac_ce_slot_rx_tx(isc,isymbol,:,:),N_Ax,N_Bx);
        G_hatAB(isc,isymbol,:,:) = hA ./ hB.';
    end
end


G_hatAB_noise = zeros(nSC,M_ofdm,N_Ax,N_Bx);
for isc = 1: nSC
    for isymbol = 1: M_ofdm
        hB = reshape(HBA_dlisac_ce_slot_rx_tx_noise(isc,isymbol,:,:),N_Bx,N_Ax);
        hA = reshape(HAB_ulisac_ce_slot_rx_tx_noise(isc,isymbol,:,:),N_Ax,N_Bx);
        G_hatAB_noise(isc,isymbol,:,:) = hA ./ hB.';
    end
end

 
% To be consisdent as before  
H_ul_total = HAB_ulisac_ce_slot_rx_tx_noise;
% H_ul_total = HAB_ulisac_ce_slot_rx_tx; %  just for debug using 
H_dl_total = HBA_dlisac_ce_slot_rx_tx_noise;
% cal_perf = C_At./C_Ar;
cal_perf = C_Ar./C_At;
cal_perf = cal_perf./cal_perf(:,:,end); % 用最后一根天线的校准系数归一化（Attention,最后一根天线每个子载波每个时隙都变成了1!）
cal_perf = permute(cal_perf,[3 2 1]); % cal_perf维度为L_BS*numBs,Nt,Nf
H_ul_total = permute(H_ul_total,[3,4,1,2]);
H_dl_total = permute(H_dl_total,[3,4,1,2]);
timing_error = tau_AB;                
timing_error = timing_error(:,1)-timing_error(end,1);


% % % check if right
% % figure;hold on
% % bs_select = 4; ue_select = 1;
% % plot(unwrap(angle(squeeze(HAB_ulisac_ce_slot_rx_tx(1,:,bs_select,ue_select)))),'k-v')
% % plot(unwrap(angle(squeeze(HAB_ulota_isac_ce_slot_rx_tx(1,:,bs_select,ue_select)))),'r-^')
% % plot(unwrap(angle(squeeze(H_ul_total(bs_select,ue_select,1,:)))),'b-*')
% % legend('srs track','ota track','ul') % attention,,here we use HAB_ulota_isac_ce_slot_rx_tx or HAB_ulota_sense_ce_slot_rx_tx


end


function x = generateRandVec(vmin,vmax)
    if rand() < 0.5
    % 50%概率选择 [-2, -1]
        x = -vmax + (rand() * (vmax-vmin)); % 映射到 [-2, -1]
    else
    % 50%概率选择 [1, 2]
        x = vmin + (rand() * (vmax-vmin));  % 映射到 [1, 2]
    end
end
