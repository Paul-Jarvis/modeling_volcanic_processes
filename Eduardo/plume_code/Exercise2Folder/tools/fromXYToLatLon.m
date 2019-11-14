function [LAT, LON, Z] = fromXYToLatLon(lat0, lon0, x_in, y_in, z_in, theta_wind)
    % This function provides the latitude and longitude of a point
    % expressed in a system of reference of plume model. 
    % It needs:
    
    % lat0: vent latitude
    % lon0: vent longitude
    % x_in: x particle coordinate in the plume SOR
    % y_in: y particle coordinate in the plume SOR
    % z_in: z particle coordinate in the plume SOR
    % theta_wind: angle FROM which the wind blows. It is zero from north
    % and then clockwise (IN RADIANS)
    
    % output: LAT, LON. Real coordinates on the Earth

%         alpha = theta_wind-pi;
    alpha = 3/2*pi-theta_wind;
        
    x_m = x_in;
    y_m = y_in;
    z_m = z_in;

    x_s = x_m.*cos(alpha) - y_m.*sin(alpha);
    y_s = x_m.*sin(alpha) + y_m.*cos(alpha);
    z_s = z_m;
        
    dx = x_s;
    dy = y_s;
    dz = z_s;
        
    R = 6371*10^3;

    d = sqrt(dx.^2+dy.^2);
    b = mod(atan2d(dx,dy),360);

    LAT = asind(sind(lat0) .* cosd(d./R) + cosd(lat0) .* sin(d./R) .* cosd(b));
    LON = lon0 + atan2d( sind(b) .* sin(d./R) .* cosd(lat0), cos(d./R) - sind(lat0) .* sind(LAT));
    Z = z_in;
    
end