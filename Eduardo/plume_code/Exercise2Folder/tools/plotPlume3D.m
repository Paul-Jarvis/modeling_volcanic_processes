function plotPlume3D(t, z, angles, radii)
% The entry values are in meters!    
   
    z = z.*1e-3; % From m to km
    radii = radii.*1e-3; % From m to km
    for i = 1:length(t)-20       
            center_z = z(i);
            center_x = z(i)*cos(angles(i))/sin(angles(i));
            center_y = 0;
            center = [center_x center_y center_z];
            normal = [cos(angles(i)) 0 sin(angles(i))];
            radius = radii(i);
        
            figure(234121)
            plotCircle3D(center,normal,radius)

            axis equal;
            hold on;
    end

end

function plotCircle3D(center,normal,radius)
% This function plots a circle

    theta=0:0.01:2*pi;
    v = null(normal);
    points = repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    plot3(points(1,:),points(2,:),points(3,:),'r-');
    xlabel('x');
    ylabel('y');
    zlabel('z');

end