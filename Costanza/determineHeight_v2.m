function [z_out_2] = determineHeight_v2(downwind_distance, crosswind_distance, density, pSizecm, scenario, data)

    % The following Matlab function returns the estimated plume height
    % starting from known maximum downwind and crosswind distances. 
    % This version is 2.0 updated to the 31/10/2018
    
    % INPUT VALUES:
    % - downwind_distance: maximum downwind distance [m]
    % - crosswind_distance: maximum crosswind distance [m]
    % - density: clast density [kg/m^3]
    % - pSizecm: clast size [cm]
    % - scenario: 1, small eruption; 2, intermediate; 3, large
    % - data: all data (to import as it is)
    
    % OUTPUT VALUE
    % - z_out_2: plume height in meters
    
    % EXAMPLE of USE
    % 1) Import data: press "Import Data" on the "Home" window and import the
    % file "data.mat"
    % 2) Call the function "determineHeight" in the Command Window, inserting 
    % the known downwind, crosswind ranges; particle density, particle size
    % and eruptive scenario. Take care in setting the right units.
    % 3) Press: determineHeight_v2(10000, 10000, 2500, 3.2, 2, data). 
    
    % Program written in Matlab by Eduardo Rossi
    % For any problems or questions please do not hesitate to contact me:
    % Eduardo.Rossi@unige.ch
    
    close all;
    
    x_down = downwind_distance;
    y_down = crosswind_distance;
    
    table_08 = data.out08;
    table_16 = data.out16;
    table_32 = data.out32;
    table_64 = data.out64;
    
    vector_density_08 = [2500 2000 1500 1000 500 250];                     % [kg/m^3]
    vector_density_16 = vector_density_08;
    vector_density_32 = vector_density_08;
    vector_density_64 = vector_density_08;
       
    vector_size_16 = [1.6 2 2.7 4 8 16];                                   % [cm]
    vector_size_32 = [3.2 4 5.3 7.9 15.8 31];
    vector_size_08 = [0.8 1.1 1.3 2.0 4.1 8.1];
    vector_size_64 = [6.4 8.0 10.8 16.1 32.1 64];
    

    if scenario == 1
        
        matrix_X_08 = table_08.X_scen_1;
        matrix_Y_08 = table_08.Y_scen_1;
        matrix_Z_08 = table_08.Z_scen_1_avg;
        
        matrix_X_16 = table_16.X_scen_1;
        matrix_Y_16 = table_16.Y_scen_1;
        matrix_Z_16 = table_16.Z_scen_1_avg;
        
        matrix_X_32 = table_32.X_scen_1;
        matrix_Y_32 = table_32.Y_scen_1;
        matrix_Z_32 = table_32.Z_scen_1_avg;
        
        matrix_X_64 = table_64.X_scen_1;
        matrix_Y_64 = table_64.Y_scen_1;
        matrix_Z_64 = table_64.Z_scen_1_avg;
        
    elseif scenario == 2
        
        matrix_X_08 = table_08.X_scen_2;
        matrix_Y_08 = table_08.Y_scen_2;
        matrix_Z_08 = table_08.Z_scen_2_avg;
        
        matrix_X_16 = table_16.X_scen_2;
        matrix_Y_16 = table_16.Y_scen_2;
        matrix_Z_16 = table_16.Z_scen_2_avg;
        
        matrix_X_32 = table_32.X_scen_2;
        matrix_Y_32 = table_32.Y_scen_2;
        matrix_Z_32 = table_32.Z_scen_2_avg;
        
        matrix_X_64 = table_64.X_scen_2;
        matrix_Y_64 = table_64.Y_scen_2;
        matrix_Z_64 = table_64.Z_scen_2_avg;
        
    elseif scenario == 3
        
        matrix_X_08 = table_08.X_scen_3;
        matrix_Y_08 = table_08.Y_scen_3;
        matrix_Z_08 = table_08.Z_scen_3_avg;
        
        matrix_X_16 = table_16.X_scen_3;
        matrix_Y_16 = table_16.Y_scen_3;
        matrix_Z_16 = table_16.Z_scen_3_avg;
        
        matrix_X_32 = table_32.X_scen_3;
        matrix_Y_32 = table_32.Y_scen_3;
        matrix_Z_32 = table_32.Z_scen_3_avg;
        
        matrix_X_64 = table_64.X_scen_3;
        matrix_Y_64 = table_64.Y_scen_3;
        matrix_Z_64 = table_64.Z_scen_3_avg;
        
    else
        
        disp('Error in defining the eruptive scenario')
        
    end

    % Heights from the nomograms
    Z_estimated_08 = interp2(matrix_X_08, matrix_Y_08, matrix_Z_08, x_down, y_down, 'spline');
    Z_estimated_16 = interp2(matrix_X_16, matrix_Y_16, matrix_Z_16, x_down, y_down, 'spline');
    Z_estimated_32 = interp2(matrix_X_32, matrix_Y_32, matrix_Z_32, x_down, y_down, 'spline');
    Z_estimated_64 = interp2(matrix_X_64, matrix_Y_64, matrix_Z_64, x_down, y_down, 'spline');
     
    vector_heights_08 = Z_estimated_08*ones(1,length(vector_density_08));
    vector_heights_16 = Z_estimated_16*ones(1,length(vector_density_16));
    vector_heights_32 = Z_estimated_32*ones(1,length(vector_density_32));
    vector_heights_64 = Z_estimated_64*ones(1,length(vector_density_64));
    
    all_X_gross = [vector_density_08 vector_density_16 vector_density_32 vector_density_64];
    all_Y_gross = [vector_size_08 vector_size_16 vector_size_32 vector_size_64];
    all_Z_gross = [vector_heights_08 vector_heights_16 vector_heights_32 vector_heights_64];
    
    % Find negative heights if present
    idxs_good = find(all_Z_gross>0);
    
    if length(idxs_good)~=length(all_X_gross)
        disp('Negative heights are present!')
        disp('(Please check if the combination of downwind and crosswind ranges for the interested size');
        disp('is outside the region of computed heights in the nomogram of the paper)')
    end
    
%     all_X = all_X_gross(idxs_good);
%     all_Y = all_Y_gross(idxs_good);
%     all_Z = all_Z_gross(idxs_good);
    
    all_X = [vector_density_08 vector_density_16 vector_density_32 vector_density_64];
    all_Y = [vector_size_08 vector_size_16 vector_size_32 vector_size_64];
    all_Z = [vector_heights_08 vector_heights_16 vector_heights_32 vector_heights_64];
    
    z_out_2 = griddata(all_X, all_Y, all_Z, density, pSizecm, 'linear');
    
    if isnan(z_out_2) == 1
        disp('Linear interpolation with griddata not possible, switch to v4 option')
        disp('Interpolation not accurate! Estimated height cannot be reliable')
        z_out_2 = griddata(all_X, all_Y, all_Z, density, pSizecm, 'v4');
    end
    
    figure(2)
    plot3(all_X, all_Y, all_Z, '.')
    hold on
    plot3(density, pSizecm, z_out_2, 'or')
    xlabel('Density (kg/m^3)')
    ylabel('Size (cm)')
    zlabel('Height (m)')

end