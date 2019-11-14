function box_out = evaluateBoxNotVortex(s0, t, z, angles, v, radii, vent_height, generic_times)

% This function evaluates the volume of the box at different heights and
% provides all the info required for the plotting part. This is the "old"
% method: it means that you start with an initial height and then it
% propagates according to the rise of the plume.

% Note that for plotting the two circles you need 6 information:
% 1° circle: center1(i,:),normal1(i,:),radius1(i)
% 2° circle: center2(i,:),normal2(i,:),radius2(i)

    s0 = s0+vent_height;
    s1(1) = vent_height;
    s2(1) = s0;
    z1_mem(1) = vent_height;
    z2_mem(1) = s0*sin(angles(1)); 
    v1_mem(1) = interp1(z,v,z1_mem(1),'linear','extrap');
    v2_mem(1) = interp1(z,v,z2_mem(1),'linear','extrap');
    
    % Here I pass from the position on the z axis to the position along the s line
    s = z./sin(angles);
    
    % Correggiamo il tempo iniziale per il punto h2
    offset = interp1(s,t,s0,'linear','extrap');
    t1 = 0;
    t2 = offset;
    
    for i = 1:length(generic_times)-1
        dt = generic_times(i+1)-generic_times(i);
        
        s1(i+1) = interp1(t,s,t1+dt,'linear','extrap');
        z1_mem(i+1) = interp1(t,z,t1+dt,'linear','extrap');
        v1_mem(i+1) = interp1(z,v,z1_mem(i+1),'linear','extrap');
        
        s2(i+1) = interp1(t,s,t2+dt,'linear','extrap');
        z2_mem(i+1) = interp1(t,z,t2+dt,'linear','extrap');
        v2_mem(i+1) = interp1(z,v,z2_mem(i+1),'linear','extrap');
        
        t1 = t1+dt;
        t2 = t2+dt;
    end
   
    table = [];
    
    min_x = 0;    max_x = 0;    max_z = 0;
    max_x_plot = 0;    min_x_plot = 0;
    max_z_plot = 0;    min_z_plot = 0;
    
    for i = 1:length(generic_times)       
        angle_iter1(i) = interp1(z,angles,z1_mem(i),'linear','extrap');
        center_z1(i) = z1_mem(i); 
        center_x1(i) = z1_mem(i)*cos(angle_iter1(i))/sin(angle_iter1(i));
        center_y1(i) = 0;
        center1(i,:) = [center_x1(i) center_y1(i) center_z1(i)];
        normal1(i,:) = [cos(angle_iter1(i)) 0 sin(angle_iter1(i))];
        radius1(i) = interp1(z,radii,z1_mem(i),'linear','extrap');
        
        angle_iter2(i) = interp1(z,angles,z2_mem(i),'linear','extrap');
        center_z2(i) = z2_mem(i); 
        center_x2(i) = z2_mem(i)*cos(angle_iter2(i))/sin(angle_iter2(i));
        center_y2(i) = 0;
        center2(i,:) = [center_x2(i) center_y2(i) center_z2(i)];
        normal2(i,:) = [cos(angle_iter2(i)) 0 sin(angle_iter2(i))];
        radius2(i) = interp1(z,radii,z2_mem(i),'linear','extrap');      
    end
    
    box_out.center1 = center1;
    box_out.normal1 = normal1;
    box_out.radius1 = radius1;
    box_out.center2 = center2;
    box_out.normal2 = normal2;
    box_out.radius2 = radius2;
    
end