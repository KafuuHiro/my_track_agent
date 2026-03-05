function Tar_pos = cal_ang_2_pos(Rx_pos_center,Tx_pos_center,ang_W,ang_F,ang_string)
    % Author: Hua Geyang @ SEU
    % 
    % R/Tx_pos_center: 1x2
    % ang_W/F : 1x1 (scalar)
    % Tar_loc: 1x2
    % 
    % Test example: >> cal_ang_2_pos([2,0],[0,0],45,-45,'deg')--------------[1,1]
    % Test example: >> cal_ang_2_pos([4,0],[0,0],60,-30,'deg')--------------[1,1]
    % Test example: >> cal_ang_2_pos([2,0],[0,0],pi/4,-pi/4,'rad')----[1,sqrt(3)]
    % Test example: >> cal_ang_2_pos([4,0],[0,0],pi/3,-pi/6,'rad')----[1,sqrt(3)]
    ang_W=-ang_W;ang_F=-ang_F; %???
    Tar_pos = zeros(1,2);
    if ang_string == "deg"
        ninety_degree = 90;
    end
    if ang_string == "rad"
        ninety_degree = pi/2;
    end

    if ang_string == "deg"
        a_1 = tand(ninety_degree + ang_F);
        b_1 = -1;
        c_1 = tand(ninety_degree + ang_F) * Tx_pos_center(1,1) - Tx_pos_center(1,2);

        a_2 = tand(ninety_degree + ang_W);
        b_2 = -1;
        c_2 = tand(ninety_degree + ang_W) * Rx_pos_center(1,1) - Rx_pos_center(1,2);
    end

    if ang_string == "rad"
        a_1 = tan(ninety_degree + ang_F);
        b_1 = -1;
        c_1 = tan(ninety_degree + ang_F) * Tx_pos_center(1,1) - Tx_pos_center(1,2);

        a_2 = tan(ninety_degree + ang_W);
        b_2 = -1;
        c_2 = tan(ninety_degree + ang_W) * Rx_pos_center(1,1) - Rx_pos_center(1,2);
    end

    A = [a_1, b_1; a_2, b_2];
    D = det(A);

    D_x = b_2 * c_1 - b_1 * c_2;
    D_y = a_1 * c_2 - a_2 * c_1;

    Tar_pos(1,1) = D_x / D;
    Tar_pos(1,2) = D_y / D;
end