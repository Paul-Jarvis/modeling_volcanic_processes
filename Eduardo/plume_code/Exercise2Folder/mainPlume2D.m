function  plume_out = mainPlume2D(dem)

    close all;

    % Read the input file in the folder
    disp('Reading of the input file...');
    cd tools
    inp_var = readInputFile;
    cd ..
    
    inp_var.tempGrad_t = inp_var.tempGrad_t/1000;                          % Conversion of the temperature gradient: [K]->[K/m]
    inp_var.tempGrad_s = inp_var.tempGrad_s/1000;                          % Conversion of the temperature gradient: [K]->[K/m]

    cd tools
    wind_direction = changeWindAngle(inp_var.wind_azimuth);                % Wind azimuth
    inp_var.wind_direction = wind_direction;
    cd ..
                  
    % Evaluate the wind profile (from 0 to zmax)...
    tic
    disp('Preparing atmospheric data...');
    cd tools    
    [table_wind] = windManager(inp_var);                                   % Output: [m], [m/s]
    [table_temperature] = temperatureManager(inp_var);                     % Output: [m], [K]
    atmo_var.table_wind = table_wind;
    atmo_var.table_temperature = table_temperature;
    cd ..
    toc
    
    % Run plume model...
    tic
    cd plume
    disp('Running plume model: start')   
    [plume_out] = plumeModelMain(inp_var, atmo_var); 
    disp('Running plume model: end')
    cd ..
    toc
   
    % Plotting part
    cd tools
    plot3D2DDem(plume_out, dem)
    cd ..
    
    clc
    fprintf('Mass Eruption Rate @ the vent [kg/s]: %5.1e \n', plume_out.m(1));
    fprintf('\n');
    fprintf('Vent altitude a.s.l. [m]: %5.1f \n', plume_out.z(1));
    fprintf('\n');
    fprintf('Plume height a.s.l. [m]: %5.1f \n', plume_out.col_Height);
    fprintf('\n');
    fprintf('Height of the erupted column [m] (i.e. plume height-vent altitude): %5.1f \n', plume_out.col_Height-plume_out.z(1));
    fprintf('\n');
    fprintf('Height of the NBL a.s.l. [m]: %5.1f \n', plume_out.nbl_Height);
    fprintf('\n');
    fprintf('Radius of the plume at the vent [m]: %5.1f \n', plume_out.r(1));
    fprintf('\n');
    fprintf('Radius of the plume at the NBL [m]: %5.1f \n', plume_out.nbl_r);
    fprintf('\n');
    fprintf('Radius of the plume at the maximum height [m]: %8.1f \n', plume_out.r(end));
    fprintf('\n');
    fprintf('Exit velocity of the plume at the vent [m/s]: %5.1f \n', plume_out.vertical_velox(1));
    fprintf('\n');
    fprintf('Exit velocity of the plume at the NBL [m/s]: %5.1f \n', plume_out.nbl_u);
%     fprintf('\n');
%     fprintf('Parameter PI for weak plumes: %3.1f \n', 0.015*(plume_out.col_Height-plume_out.z(1))/(1.8*mean(plume_out.ambientwind)) * (0.1/0.5)^2 );
    
end

