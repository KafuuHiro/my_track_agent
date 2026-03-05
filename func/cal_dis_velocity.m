function [range_u,velocity_u] = cal_dis_velocity(RRU_pos,Obj_pos,Obj_velocity)
% 输入：
%   RRU_pos      : [N_node×2] 节点坐标（若为 1×2 也可）
%   Obj_pos      : [Q×2] 目标位置（可支持多个目标 Q>1；本例 Q=1）
%   Obj_velocity : [Q×2] 目标速度
% 输出：
%   range_u      : [Q×N_node] 目标到各节点的欧氏距离
%   velocity_u   : [Q×N_node] 目标速度在“目标→节点”视线方向上的径向分量（正负号由几何决定）

    assert(size(Obj_pos,1)==size(Obj_velocity,1),'输入维度错误\n')
    Q   = size(Obj_pos,1);
    N_U = size(RRU_pos,1);
    range_u   = zeros(Q,N_U);
    velocity_u= zeros(Q,N_U);

    % 计算距离
    for q=1:Q
        for i = 1:N_U
            diff = RRU_pos(i,:) - Obj_pos(q,:);
            range_u(q,i) = sqrt(sum(diff.^2));
        end
    end
    % 计算径向速度：v_rad = ((node - obj) · v_obj) / |node - obj|
    for q=1:Q
        for i = 1:N_U
            velocity_u(q,i) = ((RRU_pos(i,:)-Obj_pos(q,:))*Obj_velocity(q,:).')./range_u(q,i);
        end
    end
end