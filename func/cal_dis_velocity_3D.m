%% cal dis between AP and Target
% author:
% date：
% input: AP_Pos ,  Obj_pos
% output: range_u 
function [range_u,velocity_u] = cal_dis_velocity_3D(RRU_pos,Obj_pos,Obj_velocity) 
    assert(size(Obj_pos,1)==size(Obj_velocity,1),'输入维度错误\n')
    Q = size(Obj_pos,1); %number of reflections
    N_U = size(RRU_pos,1); %number of RXs
    range_u = zeros(Q,N_U);
    velocity_u = zeros(Q,N_U);
    for q=1:Q
        for i = 1:N_U
            range_u(q,i) = sqrt(sum((RRU_pos(i,:)-Obj_pos(q,:)).^2));
        end
    end
    for q=1:Q
        for i = 1:N_U
            velocity_u(q,i) = ((RRU_pos(i,:)-Obj_pos(q,:))*Obj_velocity(q,:).')./range_u(q,i);
            % velocity_u(q,i) = ((RRU_pos(i,:)-Obj_pos(q,:)).*Obj_velocity(q,:))./range_u(q,i);
        end
    end

end