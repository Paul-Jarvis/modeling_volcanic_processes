function plot3D2DDem(plume_out, dem)

% This function plots the volcano and the aggregation on the same figure
% (using subplot)
    
    z_max = max(plume_out.z);
    z_all = plume_out.z;
    x_all = plume_out.x;
    y_all = zeros(length(x_all),1);
    
    idxs = find(plume_out.z<=plume_out.nbl_Height);
    
    z = plume_out.z(idxs);
    x = plume_out.x(idxs);
    y = zeros(length(x),1);
    angles = plume_out.angle(idxs);
    r = plume_out.r(idxs);   
    s = plume_out.s(idxs);
    
    for i=1:length(idxs)
        s1(i) = s(i);
        r1(i) = r(i);
        angle_iter1(i) = angles(i);
        center_z1(i) = z(i);
%         center_x1(i) = center_z1(i)*cos(angle_iter1(i))/sin(angle_iter1(i));
        center_x1(i) = x(i);
        center_y1(i) = 0;
        center1(i,:) = [center_x1(i) center_y1(i) center_z1(i)];
        normal1(i,:) = [cos(angle_iter1(i)) 0 sin(angle_iter1(i))];            
    end
             
    N_font = 20;
    
    
    figure(1)
    plot(plume_out.vertical_velox, plume_out.z, '-', 'LineWidth', 5)
    hold on
    plot(plume_out.ambientwind, plume_out.z, '-r', 'LineWidth', 5)
    xlabel('Velocity (m/s)')
    ylabel('Altitude (a.s.l.) (m)')
    legend('Plume', 'Wind')
    set(gca,'FontSize', N_font)
    
    figure(2)
    plot(plume_out.density, plume_out.z, '-', 'LineWidth', 5)
    hold on
    plot(plume_out.ambientdensity, plume_out.z, '-r', 'LineWidth', 5)
    xlabel('Density (kg/m^3)')
    ylabel('Altitude (a.s.l.) (m)')
    legend('Plume', 'Atmosphere')
    set(gca,'FontSize', N_font)
    
    figure(3)
    plot(plume_out.theta, plume_out.z, '-', 'LineWidth', 5)
    xlabel('Plume temperature (K)')
    ylabel('Altitude (a.s.l.) (m)')
    set(gca,'FontSize', N_font)
    
    figure(4)
    plot(plume_out.n_s, plume_out.z, '-', 'LineWidth', 5)
    hold on
    plot(plume_out.n_gas, plume_out.z, '-r', 'LineWidth', 5)
    xlabel('Mass fractions (-)')
    ylabel('Altitude (a.s.l.) (m)')
    legend('Solid', 'Gas (dry+vapor)')
    set(gca,'FontSize', N_font)
   
    figure(5)
    plot(plume_out.ambientwind, plume_out.z, '-', 'LineWidth', 5)
    xlabel('Wind velocity (m/s)')
    ylabel('Altitude (a.s.l.) (m)')
    set(gca,'FontSize', N_font)
    
    figure(6)
    plot(plume_out.ambientT, plume_out.z, '-', 'LineWidth', 5)
    xlabel('Atmosphere temperature (K)')
    ylabel('Altitude (a.s.l.) (m)')
    set(gca,'FontSize', N_font)
    
    % Color shading of the plume
    color_vector = linspace(0.15,0.90,length(z));
    color_agg = [1 0 0];
  
    % Dem
    dem.X(dem.X<0) = 360+(dem.X(dem.X<0));   
    dem_lon = dem.X(100:250,200:350);
    dem_lat = dem.Y(100:250,200:350);
    dem_z = dem.Z(100:250,200:350);
       
    hFig = figure(7);
%         set(hFig, 'Position', [55, 55, 1600, 900]);
%         set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    hold on
    
    lat0 = 37.75;                                                     % Vent latitude
    lon0 = 15;                                                    % Vent longitude
    theta_wind = 0;
    
    % Dem plot
    surf(dem_lon, dem_lat, dem_z)
    shading interp
    colormap copper
    camlight('left');
    lighting phong;
    hold on;
    
    for i = 1:length(idxs)
            
            plotCircle3D(center1(i,:), normal1(i,:), r1(i), 'r', lat0, lon0, theta_wind);
            zlabel('Height (m)');
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            view(48,28)
            set(gca,'xaxisLocation','top')
            zlim([0 z_max]);
       
% %        To create gif video
%        cd ../video
%             rootOUT = 'video';
%             fNamePPM = getFileName(rootOUT,i,'ppm');
%             print(hFig,'-dppm',fNamePPM)
%        cd ../tools
            
%        pause(0.5)
       
    end
   
    hold on
    [LAT_all, LON_all, Z_all] = fromXYToLatLon(lat0, lon0, x_all, y_all, z_all, theta_wind);
    plot3(LON_all, LAT_all, Z_all, '.')
   
end

function plotCircle3D(center, normal, radius, color, lat0, lon0, theta_wind)
    
% This function plots a circle    
    
    theta = 0:0.1:2*pi;
    v = null(normal);
    points = repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    grey = [0.4,0.4,0.4];
%     plot3(points(1,:),points(2,:),points(3,:),'Color',color);
    
    all_x = points(1,:);
    all_y = points(2,:);
    all_z = points(3,:);
    
    [LAT, LON, Z] = fromXYToLatLon(lat0, lon0, all_x, all_y, all_z, theta_wind);    
    
    % Conversion of the points from [x,y,z] to [lon,lat,z]    
    plot3(LON, LAT, Z, 'Color',color);
    
end

function plotCircle3DAgg(center,normal,radius,color)
% This function plots a circle

    theta=0:0.01:2*pi;
    v = null(normal);
    points = repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    fill3(points(1,:), points(2,:), points(3,:), color)
%     plot3(points(1,:),points(2,:),points(3,:),'r');
end

function fileName = getFileName(root,iter,ext)

            % ... create iteration number
            iterStr = num2str(iter);
            zeri = '_00000000';
            iterStr = strcat( zeri(1:(9-length(iterStr))),iterStr);
            fileName = strcat( root,iterStr,'.',ext);
            
end

function plotDataSurf(i, data)
  
    IO = data.DT_IO;
    
    
    empty = [1 1 1];
    full = 'c';
    % Faces definition on the diagonal 
    fac1 = [1 5 8 4; 5 6 2 1; 6 7 3 2; 8 7 3 4; 5 6 7 8];
     % Faces definition for the triangular elements
    fac2 = [1 4 6 3; 5 2 1 4; 6 5 2 3];
    fac3 = [1 2 3];
    
    % Elements on the diagonal 
    for counter = 1:data.n_bins_diagonal
            
%         z_single_element = log10(N_p_time(i,counter));
        z_single_element = data.wf(i,counter);

        vertixes_ground = data.vertixes(counter).table;
        vertixes_top = data.vertixes(counter).table;
        vertixes_ground(:,3) = 0;
        vertixes_top(:,3) = z_single_element;
        
        all_x = vertixes_ground(:,1);
        all_x(end+1:end+length(vertixes_ground)) = vertixes_top(:,1);
        
        all_y = vertixes_ground(:,2);
        all_y(end+1:end+length(vertixes_ground)) = vertixes_top(:,2);
        
        all_z = vertixes_ground(:,3);
        all_z(end+1:end+length(vertixes_ground)) = vertixes_top(:,3);
        
        all_vertixes_1 = [all_x all_y all_z];
        
        if (z_single_element)> 0.01;
            h1(counter) = patch('Vertices', all_vertixes_1 ,'Faces', fac1, 'FaceColor',full,'FaceAlpha',1, 'LineWidth', 2);
        else
            h1(counter) = patch('Vertices', all_vertixes_1 ,'Faces', fac1, 'FaceColor',empty,'FaceAlpha',1, 'LineWidth', 2);
        end
        hold on;
    
    end
              
    % Triangular elements
    for counter = data.n_bins_diagonal+1:data.n_bins
%         z_single_element = log10(N_p_time(i,counter));
        z_single_element = data.wf(i,counter);
        vertixes_ground = data.vertixes(counter).table;
        vertixes_top = data.vertixes(counter).table;
        vertixes_ground(:,3) = 0;
        vertixes_top(:,3) = z_single_element;
        
        all_x = vertixes_ground(:,1);
        all_x(end+1:end+length(vertixes_ground)) = vertixes_top(:,1);
        
        all_y = vertixes_ground(:,2);
        all_y(end+1:end+length(vertixes_ground)) = vertixes_top(:,2);
        
        all_z = vertixes_ground(:,3);
        all_z(end+1:end+length(vertixes_ground)) = vertixes_top(:,3);
        
        all_vertixes_2 = [all_x all_y all_z];
        
        if z_single_element>0.01
        
            h2(counter) = patch('Vertices', all_vertixes_2 ,'Faces', fac2, 'FaceColor',full,'FaceAlpha',1, 'LineWidth', 2);
            hold on;
            h3(counter) = patch('Vertices', vertixes_top ,'Faces', fac3, 'FaceColor',full,'FaceAlpha',1, 'LineWidth', 2);
            hold on
            
        else
            
            h2(counter) = patch('Vertices', all_vertixes_2 ,'Faces', fac2, 'FaceColor',empty,'FaceAlpha',1, 'LineWidth', 2);
            hold on;
            h3(counter) = patch('Vertices', vertixes_top ,'Faces', fac3, 'FaceColor',empty,'FaceAlpha',1, 'LineWidth', 2);
            hold on
            
        end
       
    end
               
    zlim([0 1])
    view([-82 68])
    pause(1)
    delete(h1)
    delete(h2)
    delete(h3)
            
end





function tick = plotBarColor(data, X, Y, phi_c, volumes_c_phi)
    
    x_axis = X(:,1);
    y_axis = Y(1,:);
    out_bar = bar3(data');
    zlim([0 1]);
    set(gcf,'units','normalized','outerposition',[0 0 0.95 0.95])    
    tick = get(gca,{'XTick','YTick'});
    
    x_tick_old = cell2mat(tick(1));
    
    y_tick_old = cell2mat(tick(2));
    
    for i=1:length(x_tick_old)
        id = 1;
        if (x_tick_old(i)>=1) && (x_tick_old(i)<=length(x_axis)) 
%             x_tick_new(i) = x_axis(x_tick_old(i));
            x_tick_new(i) = phi_c(i);
        else
            x_tick_new(i) = -999; 
        end
    end

    for i=1:length(y_tick_old)
        id = 1;
        if (y_tick_old(i)>= 1) && (y_tick_old(i)<=length(y_axis)) 
%             y_tick_new(i) = y_axis(y_tick_old(i));
            y_tick_new(i) = volumes_c_phi(i);
        else
            y_tick_new(i) = NaN; 
        end
    end
 
    set(gca,'XTickLabel',x_tick_new);
    set(gca,'YTickLabel',y_tick_new); 

    b = out_bar;
    for k = 1:length(b)
        zdata = b(k).ZData;
        b(k).CData = zdata;
        b(k).FaceColor = 'interp';
        
    end


end