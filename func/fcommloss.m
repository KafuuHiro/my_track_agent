function [Amp_r,P_r] = fcommloss(P_t,G_t,G_r,lamda,R)
% 通信自由空间损耗公式

% P_t = 1;  % 发射总功率
% G_t = 1;  % 发射天线增益
% G_r = 1;  % 接收天线增益

P_r = P_t * G_t * G_r * lamda.^2  / ((4*pi).^2 * R^2);
Amp_r = sqrt(P_r);

end