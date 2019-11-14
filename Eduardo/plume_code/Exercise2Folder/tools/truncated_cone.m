function [x,y,z]= truncated_cone(volcano_height, radius, alpha)
% 
%     R_top = radius;
%     R_bottom = radius;
    
    R_top = radius;
    R_bottom = R_top+volcano_height*tan(alpha);
   
    r = linspace(R_top,R_bottom,5); % Anello circolare
    zboh = linspace(0,2500,5);
    theta = linspace(0,2*pi,10); % Angolo
    [r,theta] = meshgrid(r,theta)
    x = r.*cos(theta);
    y = r.*sin(theta);
%     z = linspace(0,volcano_height,10);
    z = r;
%     mesh(x,y,z);
  

%  [x,y,z] = cylinder(0:volcano_height);
% figure
% surfnorm(x,y,z)
% axis([-12 12 -12 12 -0.1 1])
end