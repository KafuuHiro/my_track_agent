function ang_doa = cal_doa(RRU_pos_center,Tar_pos,ang_string)
    % Author: Hua Geyang @ SEU
    % 
    % RRU_pos_center: 1x2
    % Tar_pos: 1x2
    % ang_doa: 1x1 (scalar)
    % ang_string: 'deg' or 'rad'
    % angle > 0: counterclockwise from the positive y-axis of the RRU
    % 
    % Test example: >> cal_doa([2,0],[1,1],'deg')----------------(45)
    % Test example: >> cal_doa([0,0],[1,1],'deg')---------------(-45)
    % Test example: >> cal_doa([4,0],[1,sqrt(3)],'deg')----------(60)
    % Test example: >> cal_doa([2,0],[1,1],'rad')--------------(pi/4)
    % Test example: >> cal_doa([4,0],[1,sqrt(3)],'rad')--------(pi/3)

    delta_x = Tar_pos(1,1) - RRU_pos_center(1,1);
    delta_y = Tar_pos(1,2) - RRU_pos_center(1,2);
    if ang_string == "deg"
        ang_doa = atand(-delta_x/delta_y); % angle from the positive y-axis
    end
    if ang_string == "rad"
        ang_doa = atan(-delta_x/delta_y); % angle from the positive y-axis
    end
end