function box_out = evaluateBox(z, angles, radii, generic_times)

% This function evaluates the volume of the box at different heights and
% provides all the info required for the plotting part.
% Note that for plotting the two circles you need 6 information:
% 1° circle: center1(i,:),normal1(i,:),radius1(i)
% 2° circle: center2(i,:),normal2(i,:),radius2(i)

     % Here I pass from the position on the z axis to the position along the s line
    s = z./sin(angles);
    
    for i=1:length(s)
        s1(i) = s(i);
        r1(i) = r(i);
        angle_iter1(i) = angles(i);
        center_z1(i) = interp1(s, z, s1(i), 'linear','extrap');
        center_x1(i) = center_z1(i)*cos(angle_iter1(i))/sin(angle_iter1(i));
        center_y1(i) = 0;
        center1(i,:) = [center_x1(i) center_y1(i) center_z1(i)];
        normal1(i,:) = [cos(angle_iter1(i)) 0 sin(angle_iter1(i))];
               
        box_out.center1(i,:) = center1(i,:);
        box_out.normal1(i,:) = normal1(i,:);
        box_out.radius1(i) = r1(i);
        box_out.angle1(i) = angle_iter1(i);
    end       
    
end