function [new_AoA, new_AoD, success] = coordinate_descent_tracking(...
    current_AoA, current_AoD, current_gain, threshold,...
    H_matrix, Rx, Tx, ang_array)
% 坐标下降法跟踪函数
% 输入:
%   current_AoA - 当前AoA波束编号
%   current_AoD - 当前AoD波束编号  
%   current_gain - 当前信道增益
%   threshold - 增益门限值
%   H_matrix - 当前信道矩阵
%   Rx - 接收端参数结构体
%   Tx - 发送端参数结构体
%   ang_array - 角度数组(度)
% 输出:
%   new_AoA - 新的AoA波束编号
%   new_AoD - 新的AoD波束编号
%   success - 是否成功标志

% 初始化参数
max_iterations = 5;  % 最大迭代次数
step_size = 1;       % 初始步长
min_step = 0.5;      % 最小步长
new_AoA = current_AoA;
new_AoD = current_AoD; 
success = false;

% 参数边界检查
if current_AoA < 1 || current_AoA > size(Rx.DFT_codebook,2) || ...
   current_AoD < 1 || current_AoD > size(Tx.DFT_codebook,2)
    return;
end

% 迭代优化
for iter = 1:max_iterations
    improved = false;
    
    % AoA方向优化
    for delta = [-step_size, 0, step_size]
        test_AoA = round(new_AoA + delta);
        test_AoA = max(1, min(test_AoA, size(Rx.DFT_codebook,2)));
        
        w = Rx.DFT_codebook(:,test_AoA)./sqrt(Rx.N);
        f = Tx.DFT_codebook(:,new_AoD)./sqrt(Tx.N);
        test_gain = abs(w'*H_matrix*f);
        
        if test_gain > current_gain
            new_AoA = test_AoA;
            current_gain = test_gain;
            improved = true;
        end
    end
    
    % AoD方向优化 
    for delta = [-step_size, 0, step_size]
        test_AoD = round(new_AoD + delta);
        test_AoD = max(1, min(test_AoD, size(Tx.DFT_codebook,2)));
        
        w = Rx.DFT_codebook(:,new_AoA)./sqrt(Rx.N);
        f = Tx.DFT_codebook(:,test_AoD)./sqrt(Tx.N);
        test_gain = abs(w'*H_matrix*f);
        
        if test_gain > current_gain
            new_AoD = test_AoD;
            current_gain = test_gain;
            improved = true;
        end
    end
    
    % 检查是否满足门限要求
    if current_gain > threshold
        success = true;
        return;
    end
    
    % 调整步长
    if ~improved
        step_size = max(step_size/2, min_step);
    end
end

% 最终检查是否满足门限
success = current_gain > threshold;
end
