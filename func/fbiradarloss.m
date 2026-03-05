function [Amp_r,P_r] = fbiradarloss(P_t,G_t,G_r,lamda,R1,R2,rcs)
% 双基雷达自由空间损耗公式

% P_t = 1;  % 发射总功率
% G_t = 1;  % 发射天线增益
% G_r = 1;  % 接收天线增益

P_r = P_t * G_t * G_r * lamda.^2 * rcs ./ ((4*pi).^3 .* R1.^2.* R2.^2);
Amp_r = sqrt(P_r);

end