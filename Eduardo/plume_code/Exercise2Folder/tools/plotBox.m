function table = plotBox(box_out, plot_choice, radius_vent, vent_height)

% h0 is the initial length of the box
% velox is the velocity along the axis (not to be confused with velocity
% projection on z-axis).
% If plot_choice = 0: plot box
% If plot_choice = 1: plot box with aggregation
% If plot_choice = 2: plot plume
% Else: no plot

% Here we evaluate the output from aggregation and we define "generic_times"
%     dataAggregation = getDataAggregation;
%        
%     for i=1:length(dataAggregation)
%         sizes(i,:) = dataAggregation(i).sizes;
%         dwdphi(i,:) = dataAggregation(i).dwdphi;
%     end       
%     
%     max_agg_x = max(max(sizes)); min_agg_x = min(min(sizes));
%     max_agg_y = max(max(dwdphi)); min_agg_y = min(min(dwdphi));
  
    table = [];
    
    min_x = 0;    max_x = 0;    max_z = 0;
    max_x_plot = 0;    min_x_plot = 0;
    max_z_plot = 0;    min_z_plot = 0;
    
    center1 = box_out.center1;
    normal1 = box_out.normal1;
    radius1 = box_out.radius1;
    center2 = box_out.center2;
    normal2 = box_out.normal2;
    radius2 = box_out.radius2;
   
  % For to automatically design the boundaries of the plot
    for i = 1:length(box_out.times)          
        angle_iter2(i) = box_out.angle2(i);
        center_x2(i) = center2(i,1);    center_z2(i) = center2(i,3);
       
        if (round(center_x2)<=min_x)
            min_x = center_x2(i);
            min_x_plot = min_x-abs(radius2(i)*sin(angle_iter2(i)));
        end
        
        if (center_x2(i)>max_x)
            max_x = center_x2(i);
            max_x_plot = max_x+abs(radius2(i)*sin(angle_iter2(i)));
        end
        
        if (center_z2(i)>max_z)
            max_z = center_z2(i);
            max_z_plot = max_z+abs(radius2(i)*cos(angle_iter2(i)));
        end
    end
    
    max_y_plot = max(radius2);
    min_y_plot = -max(radius2);    
    
    % This useless for is just to have the borders of the plot in
    % authomatic
    color_vector = linspace(0.15,0.90,length(box_out.times));
 
    % Get points
   [Points, G, min_x, min_y] = drawVolcano(radius_vent, vent_height);    
   
   if (min_x<min_x_plot)
      min_x_plot = min_x; 
   end
   
   if (min_y<min_y_plot)
      min_y_plot = min_y; 
   end
   
     
   % ______________________________________________________________________
   % ______________________________________________________________________
   
   if (plot_choice == 0)
       % Plot box volume
       hFig = figure(234122)
       set(hFig, 'Position', [30, 30, 1500, 700]);
       
       for i = 1:length(box_out.times) 
           %input('Plot')
            hold off;
%             % Plot ground
%             fill3([min_x_plot max_x_plot max_x_plot min_x_plot],...
%                 [min_y_plot min_y_plot max_y_plot max_y_plot],[0 0 0 0],[0.55 0.7 0.2]);
%             hold on;
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
            plotCircle3D(center1(i,:),normal1(i,:),radius1(i),color_vector(i));
            hold on;
            plotCircle3D(center2(i,:),normal2(i,:),radius2(i),color_vector(i));
                
            xlim([min_x_plot max_x_plot]);
            ylim([min_y_plot max_y_plot]);
            zlim([min_z_plot max_z_plot]);
            
            view([-15 16]);
       
            drawnow
            hold off
        end
   % ______________________________________________________________________
   elseif (plot_choice == 1)
       % Plot box volume + aggregation
       
       hFig = figure(234123);
       set(hFig, 'Position', [30, 30, 1500, 700]);
       
       for i = 1:length(box_out.times)   
            hold off;
            subplot(4,5,[1:3 6:8 11:13 16:18])
%             % Plot ground
%             fill3([min_x_plot max_x_plot max_x_plot min_x_plot],...
%                 [min_y_plot min_y_plot max_y_plot max_y_plot],[0 0 0 0],[0.55 0.7 0.2]);
%             hold on;
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
            plotCircle3D(center1(i,:),normal1(i,:),radius1(i),color_vector(i));
            hold on;
            plotCircle3D(center2(i,:),normal2(i,:),radius2(i),color_vector(i));
                
            xlim([min_x_plot max_x_plot]);
            ylim([min_y_plot max_y_plot]);
            zlim([min_z_plot max_z_plot]);

            view([-15 16]);
            drawnow
            hold off
            subplot(4,5,[9:10 14:15])
%             data_plot_aggr = getDataAggregation;
%             plot(data_plot_aggr(i).sizes,data_plot_aggr(i).dwdphi,'-')
%             %semilogx(data_plot_aggr(i).sizes,data_plot_aggr(i).dwdphi,'-')
%             axis([min_agg_x max_agg_x min_agg_y max_agg_y]);
            hold on
          
       end
   % ______________________________________________________________________
   elseif (plot_choice == 2)
        % Plot plume
        figure(234121)
        for i = 1:length(box_out.times)
            center_z = z1_mem(i);
            center_x = z1_mem(i)*cos(angle_iter1(i))/sin(angle_iter1(i));
            center_y = 0;
            center = [center_x center_y center_z];
            normal = [cos(angle_iter1(i)) 0 sin(angle_iter1(i))];
            radius = radius1(i);
            
%             hold off;
            % Plot ground
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
            
            plotCircle3D(center,normal,radius,color_vector(i));
            xlim([min_x_plot max_x_plot]);
            ylim([min_y_plot max_y_plot]);
            zlim([min_z_plot max_z_plot]);
            view([-15 9]);
            drawnow
            hold on;
        end
   else
       % No plot
       disp('No plot for plume');
   end
      
end



function plotCircle3D(center,normal,radius,color)
% This function plots a circle

    theta=0:0.01:2*pi;
    v = null(normal);
    points = repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    fill3(points(1,:), points(2,:), points(3,:), [color, color, color])
%     plot3(points(1,:),points(2,:),points(3,:),'r-');
    xlabel('x');
    ylabel('y');
    zlabel('z');

end