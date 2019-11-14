function angle_end = changeWindAngle(angle_ini)
    % The initial angle is in degree
    % The output angle is in degree
    
    angle_ini = deg2rad(angle_ini);
    
    if (angle_ini<=0) && (angle_ini<=pi)
        angle_ini = angle_ini+pi;
    else
        angle_ini = angle_ini-pi;
    end

    angle_end = rad2deg(angle_ini);
    
end