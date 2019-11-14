function plot3DFinal(plume_out, radius_vent, vent_height)

% This function plots the volcano and the aggregation on the same figure
% (using subplot)
    
    idx = length(plume_out.z)-3;
    z = plume_out.z(1:idx);
    angles = plume_out.angle(1:idx);
    r = plume_out.r(1:idx);
    sizes = plume_out.aggregation.sizes; 
    dwdphi = plume_out.aggregation.dwdphi(1:idx,:);
    
    mass_fractions = plume_out.mass_fractions;
    
    min_x = 0;    max_x = 0;    max_z = 0;
    max_x_plot = 0;    min_x_plot = 0;
    max_z_plot = 0;    min_z_plot = 0;

    max_agg_x = max(max(sizes)); min_agg_x = min(min(sizes));
    max_agg_y = max(max(dwdphi)); min_agg_y = min(min(dwdphi));
    
% Here I pass from the position on the z axis to the position along the s line
    s = z./sin(angles);
    
    for i=1:length(s)
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
   
   hFig = figure(234123);
        set(hFig, 'Position', [55, 55, 1600, 900]);

%    for i = 1:length(z)   
%         hold off;
%        for j = 1:length(z)
%              subplot(5,6,[1:6 7:12 13:18]);
% %             subplot(4,5,[1:3 6:8 11:13 16:18])
%             % Plot ground
%             axis equal
%             fill3([min_x_plot max_x_plot max_x_plot min_x_plot],...
%                   [min_y_plot min_y_plot max_y_plot max_y_plot],[0 0 0 0],[0.55 0.7 0.2]);
%             hold on;
%             % Plot volcano      
%             fill3([Points.P8(1) Points.P1(1) Points.P9(1) Points.P16(1)],...
%                   [Points.P8(2) Points.P1(2) Points.P9(2) Points.P16(2)],...
%                   [Points.P8(3) Points.P1(3) Points.P9(3) Points.P16(3)],G);
%             hold on;
%             fill3([Points.P1(1) Points.P2(1) Points.P10(1) Points.P9(1)],...
%                   [Points.P1(2) Points.P2(2) Points.P10(2) Points.P9(2)],...
%                   [Points.P1(3) Points.P2(3) Points.P10(3) Points.P9(3)],G);
%             hold on;
%             fill3([Points.P2(1) Points.P3(1) Points.P11(1) Points.P10(1)],...
%                   [Points.P2(2) Points.P3(2) Points.P11(2) Points.P10(2)],...
%                   [Points.P2(3) Points.P3(3) Points.P11(3) Points.P10(3)],G);
%             hold on;
%             fill3([Points.P3(1) Points.P4(1) Points.P12(1) Points.P11(1)],...
%                   [Points.P3(2) Points.P4(2) Points.P12(2) Points.P11(2)],...
%                   [Points.P3(3) Points.P4(3) Points.P12(3) Points.P11(3)],G);
%             hold on;  
%             fill3([Points.P4(1) Points.P5(1) Points.P13(1) Points.P12(1)],...
%                   [Points.P4(2) Points.P5(2) Points.P13(2) Points.P12(2)],...
%                   [Points.P4(3) Points.P5(3) Points.P13(3) Points.P12(3)],G);
%             hold on;
%             fill3([Points.P5(1) Points.P6(1) Points.P14(1) Points.P13(1)],...
%                   [Points.P5(2) Points.P6(2) Points.P14(2) Points.P13(2)],...
%                   [Points.P5(3) Points.P6(3) Points.P14(3) Points.P13(3)],G);
%             hold on;
%             fill3([Points.P6(1) Points.P7(1) Points.P15(1) Points.P14(1)],...
%                   [Points.P6(2) Points.P7(2) Points.P15(2) Points.P14(2)],...
%                   [Points.P6(3) Points.P7(3) Points.P15(3) Points.P14(3)],G);
%             hold on;
%             fill3([Points.P7(1) Points.P8(1) Points.P16(1) Points.P15(1)],...
%                   [Points.P7(2) Points.P8(2) Points.P16(2) Points.P15(2)],...
%                   [Points.P7(3) Points.P8(3) Points.P16(3) Points.P15(3)],G);
%             hold on;
%             
%             if (i==j)
%                 plotCircle3DAgg(center1(j,:),normal1(j,:),r1(j),color_agg);
% %                 plotCircle3D(center1(j,:),normal1(j,:),r1(j),color_vector(j));
%                 hold on;
%             else
%                 plotCircle3D(center1(j,:),normal1(j,:),r1(j),color_vector(j));
%                 hold on;
%             end
%                        
%             xlim([min_x_plot max_x_plot]);
%             ylim([min_y_plot max_y_plot]);
%             zlim([min_z_plot max_z_plot]);
%             xlabel('Horizontal distance (m)');
% %             ylabel('y (m)');
%             zlabel('Height (m)');
%             title('Day: 4th of May ; Hour = [18 24]')
%             set(gca,'ytick',[])
%             view([-15 16]);
%             set(gca,'xaxisLocation','top')
%                         
%             subplot(5,6,[22:23 28:29]);
% %             subplot(4,5,[9:10 14:15]);
% %             plot(sizes,mass_fractions(i,:),'-');
%             bar(sizes,mass_fractions(i,:));
%             axis([min_agg_x max_agg_x min_agg_y max_agg_y]);
%             xlabel('phi');
%             ylabel('mass fraction');
%             stringa = 'Height (above vent) = ';
%             tempo = num2str(sprintf('%1.0f',z(i)-vent_height));
%             unit = ' m';
%             str_time = strcat(stringa, tempo, unit);
%             title(str_time);                      
%         
% %        drawnow
%        end
%        
% % %        To create gif video
% %        cd ../video
% %             rootOUT = 'video';
% %             fNamePPM = getFileName(rootOUT,i,'ppm');
% %             print(hFig,'-dppm',fNamePPM)
% %        cd ../tools
%             
%        drawnow
%    end
   
% Istogram and differences

   % Create axes
   axes1 = axes('Parent',hFig,'FontSize',20);
   % Uncomment the following line to preserve the X-limits of the axes
   xlim(axes1,[-6.5 11.5]);
   box(axes1,'on');
   hold(axes1,'on');

   bar(sizes,mass_fractions(end,:),'FaceColor', 'none', 'EdgeColor', 'red', 'LineWidth', 5);
   hold on
   bar(sizes,mass_fractions(1,:),'FaceColor', 'none', 'EdgeColor', 'blue', 'LineWidth', 5);
   xlabel('Phi', 'FontSize', 20);
   ylabel('mass fraction', 'FontSize', 20);
   hLeg = legend('Final', 'Initial');
   set(hLeg,'Position',[0.15 0.85 0.08 0.07],'FontSize',18);
   
%    a1 = gca;
%    pos = get(a1, 'Position')

    % Plot differences
    diff = (mass_fractions(end,:)-mass_fractions(1,:))./mass_fractions(1,:)*100;
    axes('Parent', hFig,'FontSize',20, 'Position', [.6 .6 .3 .3]);
    plot(sizes, diff, 'o-', 'LineWidth', 5)
    xlabel('phi', 'FontSize', 20);
    ylabel('% rel.', 'FontSize', 20);
    title('Mass fraction difference')
    box(axes1,'on');
    hold(axes1,'on');
%     axis(h, 'off', 'tight')
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