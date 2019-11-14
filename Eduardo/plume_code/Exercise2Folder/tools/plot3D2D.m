function plot3D2D(plume_out, radius_vent, vent_height)

% This function plots the volcano and the aggregation on the same figure
% (using subplot)
    
    radius_vent = plume_out.radius_vent;
    vent_height = plume_out.vent_height;
    idx = length(plume_out.z)-100;
    z = plume_out.z(1:idx);
    angles = plume_out.angle(1:idx);
    r = plume_out.r(1:idx);   
    s = z./sin(angles);
    wf_plot = plume_out.wf_plot;
    X = plume_out.X;
    Y = plume_out.Y;
    phi_c = plume_out.phi_c;
    volumes_c_phi = plume_out.volumes_c_phi;
    
    min_x = 0;    max_x = 0;    max_z = 0;
    max_x_plot = 0;    min_x_plot = 0;
    max_z_plot = 0;    min_z_plot = 0;
    
    x_min_bar = X(1,1);  
    x_max_bar = X(end,1);
    y_min_bar = Y(1,1);
    y_max_bar = Y(1,end);
    z_min_bar = min(min(min(wf_plot(:,:,:))));
    z_max_bar = max(max(max(wf_plot(:,:,:))));
    
    for i=1:idx
        s1(i) = s(i);
        r1(i) = r(i);
        angle_iter1(i) = angles(i);
        center_z1(i) = z(i);
        center_x1(i) = center_z1(i)*cos(angle_iter1(i))/sin(angle_iter1(i));
        center_y1(i) = 0;
        center1(i,:) = [center_x1(i) center_y1(i) center_z1(i)];
        normal1(i,:) = [cos(angle_iter1(i)) 0 sin(angle_iter1(i))];
        
        % Automatically checking the boundaries of the plot
        angle_iter2(i) = angle_iter1(i);
        center_x2(i) = center1(i,1);    center_z2(i) = center1(i,3);
        r2(i) = r1(i);
        
        if (round(center_x2(i))<=min_x)
            min_x = center_x2(i);
            min_x_plot = min_x-abs(r2(i)*sin(angle_iter2(i)));
        end
        
        if (center_x2(i)>max_x)
            max_x = center_x2(i);
            max_x_plot = max_x+abs(r2(i)*sin(angle_iter2(i)));
        end
        
        if (center_z2(i)>max_z)
            max_z = center_z2(i);
            max_z_plot = max_z+abs(r2(i)*cos(angle_iter2(i)));
        end
        
    end
    
    max_y_plot = max(r2);    min_y_plot = -max(r2);    
        
    % Color shading of the plume
    color_vector = linspace(0.15,0.90,length(z));
    color_agg = [1 0 0];
  
    % Get points
   [Points, G, min_x, min_y] = drawVolcano(radius_vent, vent_height);    
   
   if (min_x<min_x_plot)
      min_x_plot = min_x-1e3; 
   end
   
   if (min_y<min_y_plot)
      min_y_plot = min_y-1e3; 
   end
      
   % Plotting part: first we plot the entire column on the left plot
   
   hFig = figure(3);
        set(hFig, 'Position', [55, 55, 1600, 900]);

   for i = 1:idx
        hold off;
       
             subplot(5,6,[1:3 7:9 13:15 19:21]);

            % Plot ground
%             axis equal
            fill3([min_x_plot max_x_plot max_x_plot min_x_plot],...
                  [min_y_plot min_y_plot max_y_plot max_y_plot],[0 0 0 0],[0.55 0.7 0.2]);
            hold on;
            % Plot volcano      
            fill3([Points.P8(1) Points.P1(1) Points.P9(1) Points.P16(1)],...
                  [Points.P8(2) Points.P1(2) Points.P9(2) Points.P16(2)],...
                  [Points.P8(3) Points.P1(3) Points.P9(3) Points.P16(3)],G);
            hold on;
            fill3([Points.P1(1) Points.P2(1) Points.P10(1) Points.P9(1)],...
                  [Points.P1(2) Points.P2(2) Points.P10(2) Points.P9(2)],...
                  [Points.P1(3) Points.P2(3) Points.P10(3) Points.P9(3)],G);
            hold on;
            fill3([Points.P2(1) Points.P3(1) Points.P11(1) Points.P10(1)],...
                  [Points.P2(2) Points.P3(2) Points.P11(2) Points.P10(2)],...
                  [Points.P2(3) Points.P3(3) Points.P11(3) Points.P10(3)],G);
            hold on;
            fill3([Points.P3(1) Points.P4(1) Points.P12(1) Points.P11(1)],...
                  [Points.P3(2) Points.P4(2) Points.P12(2) Points.P11(2)],...
                  [Points.P3(3) Points.P4(3) Points.P12(3) Points.P11(3)],G);
            hold on;  
            fill3([Points.P4(1) Points.P5(1) Points.P13(1) Points.P12(1)],...
                  [Points.P4(2) Points.P5(2) Points.P13(2) Points.P12(2)],...
                  [Points.P4(3) Points.P5(3) Points.P13(3) Points.P12(3)],G);
            hold on;
            fill3([Points.P5(1) Points.P6(1) Points.P14(1) Points.P13(1)],...
                  [Points.P5(2) Points.P6(2) Points.P14(2) Points.P13(2)],...
                  [Points.P5(3) Points.P6(3) Points.P14(3) Points.P13(3)],G);
            hold on;
            fill3([Points.P6(1) Points.P7(1) Points.P15(1) Points.P14(1)],...
                  [Points.P6(2) Points.P7(2) Points.P15(2) Points.P14(2)],...
                  [Points.P6(3) Points.P7(3) Points.P15(3) Points.P14(3)],G);
            hold on;
            fill3([Points.P7(1) Points.P8(1) Points.P16(1) Points.P15(1)],...
                  [Points.P7(2) Points.P8(2) Points.P16(2) Points.P15(2)],...
                  [Points.P7(3) Points.P8(3) Points.P16(3) Points.P15(3)],G);
            hold on;
             
            plotCircle3D(center1(i,:),normal1(i,:),r1(i),color_vector(i));
            
            xlim([min_x_plot max_x_plot]);
            ylim([min_y_plot max_y_plot]);
            zlim([min_z_plot max_z_plot]);
            xlabel('Horizontal distance (m)');
            zlabel('Height (m)');
            set(gca,'ytick',[])
            view([-12 6]);
            set(gca,'xaxisLocation','top')
                        
            subplot(5,6,[16:18 22:24 28:30]);
%             subplot(1,2,2)
            
            mat(:,:) = wf_plot(i,:,:);
            plotBarColor(mat, X, Y, phi_c, volumes_c_phi);
            view(29,42)
                        
            xlabel('Masses');
            ylabel('Volumes');
            stringa = 'Height (above vent) = ';
            tempo = num2str(sprintf('%1.0f',z(i)-vent_height));
            unit = ' m';
            str_time = strcat(stringa, tempo, unit);
            title(str_time);                      
        
           
       
%        To create gif video
       cd ../video
            rootOUT = 'video';
            fNamePPM = getFileName(rootOUT,i,'ppm');
            print(hFig,'-dppm',fNamePPM)
       cd ../tools
            
       drawnow
%        pause(0.5)
%        input('dddd')
       hold off
   end
   
end


function plotCircle3D(center,normal,radius,color)
% This function plots a circle

    theta=0:0.01:2*pi;
    v = null(normal);
    points = repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    fill3(points(1,:), points(2,:), points(3,:), [color, color, color])
    grey = [0.4,0.4,0.4];
%     plot3(points(1,:),points(2,:),points(3,:),'Color',grey);
    
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

function plotTotalBar(data3, X, Y, t)

    x_min = X(1,1);  
    x_max = X(end,1);
    y_min = Y(1,1);
    y_max = Y(1,end);
    z_min = min(min(min(data3(:,:,:))));
    z_max = max(max(max(data3(:,:,:))));
    
    for i=1:length(t)
        mat(:,:) = data3(i,:,:);
%         mat = data3;
%         plotBarColor(mat, X, Y, x_min, y_min, z_min, x_max, y_max, z_max);
        plotBarColor(mat, X, Y);
        view(-25,48)
        drawnow;
        pause(0.2)
%         input('sdsdsdds')
    end

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