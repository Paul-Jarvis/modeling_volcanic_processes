function outputTable = forwardMainGBF
    
    % Main function for the inversion of the balistic code.

    % Some notes on how it works the run:
    % - "gbf_FULL.jar" must always be present because this is the name of the
    % actual "jar" file we use
    % - "config/VulcanoDragDefault.conf" is the single configuration file
    % we run. 
    % 
    % INPUT: all the inputs are contained inside the configuration file
    % within the folder "config"
    
    % OUTPUT: the output file is inside the folder "results" (that YOU must
    % create) and the actual file is named "listImpacts.dat". Inside this
    % file each bomb has a row.


    close all;
    data = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT PARAMETERS FOR THE EXERCISE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    data.N_objects = 20000;                                                   % Number of objects released (#)
    data.density_mean = 1500;                                              % Mean value of clast density (kg/m^3)
    data.density_std = 1000;                                                  % Std value of clast density distribution (kg/m^3)
    data.phi_mean = -8;                                                    % Mean size of the clast (in phi scale)
    data.phi_std = 2;                                                      % Std size of the clast distribution (in phi scale)
    data.velox_mean = 250;                                                 % Mean velocity of the clast (m/s)
    data.velox_std = 100;                                                    % Std size of the clast velocity (m/s)
    data.angle_tilt = 0;                                                   % Tilt angle (degree)
    data.angle_spread = 90;                                                % Spread angle around the tilt
    data.angle_azimuth = 0;                                                % Azimuth tilt (degrees, 0 -> north; 90 -> east)
    data.wind_intensity = 0;                                               % Wind intensity (m/s)
    data.wind_direction = 0;                                               % Wind direction (degrees, 0 -> north; 90 -> east)
    data.X_vent = 496634;                                                  % Eastern coordinate of the vent
    data.Y_vent = 4250706;                                                 % Northern coordinate of the vent
    data.reduced_drag = 200;                                               % Reduced drag radius
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                               SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    data.filename_dem = 'dem/dem_10m.txt';                                 % Name of the DEM file
    data.filename_output_raw = 'results/listImpacts.dat';                  % Name of the output file
    data.filename_output_clean = 'results/listImpactsClean.dat';           % Name of the output file (cleaned)
    
    % Reading the DEM file
    [data.dem_X_final, data.dem_Y_final, data.dem_Z_final] = readDEM(data.filename_dem);
    data.Z_vent = interp2(data.dem_X_final,data.dem_Y_final,data.dem_Z_final,data.X_vent,data.Y_vent)
        
    tic
        
    fprintf('Particles released %7i\n', data.N_objects);
    fprintf('running...\n');
        
    % Create single input file for the i-esim iteration
    createInputFile(data);
        
 
    % Run the single simulation
    % N.B.: We run the java "jar" file using the Matlab function "system"
%     [status result] = system('java -jar gbf_FULL_old.jar config/input.conf 2')
    [status result] = system('java -jar gbf_FULL_lastversion.jar config/input.conf 2');
        
    fprintf('single simulation done!\n');
    toc
        
    % We clean the output file creating a new file without commas
    replaceCommas(data.filename_output_raw, data.filename_output_clean);

    % We read the output file from the clean output one
    data.output_all = readOutputFile(data);
        
    tabella_all = data.output_all;
    contatore_noia = 1;
    
    for i=1:data.N_objects
        
        x_clast = tabella_all(i,1);
        y_clast = tabella_all(i,2);
        
        distance = sqrt((data.X_vent-x_clast)^2+(data.Y_vent-y_clast)^2);

%         if  (distance<1.3*1000 && distance>0.68*1000 && x_clast<data.X_vent && y_clast>0.99*data.Y_vent)
%         if  (distance<0.5*1000 && distance>0.2*1000 && x_clast<data.X_vent && y_clast>0.99*data.Y_vent)
%             if  (distance<0.3*1000)
           if 2>1 
            tabella_small(contatore_noia,:) = tabella_all(i,:);
            tabella_small(contatore_noia,11) = sqrt( 2 * tabella_small(contatore_noia,6) / tabella_small(contatore_noia,4) );
            tabella_small(contatore_noia,7) = tabella_small(contatore_noia,4)/(pi/6*tabella_small(contatore_noia,5)^3);
            contatore_noia = contatore_noia+1;
        
        end
        
    end
    
    data.output = tabella_small;
    % Plotting the results of the single simulation
    plotResults(data);
    
    
    outputTable = data.output;

end

function createInputFile(data)

    fid = fopen('config/input.conf','wt');
    
    seed_min = 1;
    seed_max = 100000;
    seed = round((seed_max-seed_min).*rand(1,1) + seed_min);
    
    % Random seed (to not be modified)
    fprintf(fid, 'rng {\n');
    fprintf(fid, '\tseed = %i\n', seed);
    fprintf(fid, '}\n');
    
    % Terrain (to not be modified?)
    fprintf(fid, 'terrain {\n');
    fprintf(fid, '\tdemFile = "%s"\n', data.filename_dem);
    fprintf(fid, 'vent {\n');
    fprintf(fid, '\tE = %i\n', data.X_vent);
    fprintf(fid, '\tN = %i\n', data.Y_vent);
    fprintf(fid, '\taltitude = %f\n', data.Z_vent+5);
    fprintf(fid, '\t}\n');
    fprintf(fid, '}\n');
    
    % Source conditions
    fprintf(fid, 'source {\n');
    fprintf(fid, '\tdensAvg = %f\n', data.density_mean);
    fprintf(fid, '\tdensStd = %f\n', data.density_std);
    fprintf(fid, '\tphiAvg = %f\n', data.phi_mean);
    fprintf(fid, '\tphiStd = %f\n', data.phi_std);
    fprintf(fid, '\tvelocityAvg = %f\n', data.velox_mean);
    fprintf(fid, '\tvelocityStd = %f\n', data.velox_std);
    fprintf(fid, '\tspread = %f\n', data.angle_spread);
    fprintf(fid, '\ttilt = %f\n', data.angle_tilt);
    fprintf(fid, '\tazimuth = %i\n', data.angle_azimuth);
    fprintf(fid, '}\n');
    
    % Wind conditions
    fprintf(fid, 'wind {\n');
    fprintf(fid, '\tspeed = %f\n', data.wind_intensity);
    fprintf(fid, '\tdirection = %f\n', data.wind_direction);
    fprintf(fid, '}\n');
    
    % Drag conditions
    fprintf(fid, 'drag {\n');
    fprintf(fid, '\ttimeStep = 0.01\n');
    fprintf(fid, '\tpressure0 = 1.01325e5\n');
    fprintf(fid, '\ttemp0 = 298\n');
    fprintf(fid, '\tthermalLapse = -6.5e-3\n');
    fprintf(fid, '\treducedDragRadius = %f\n', data.reduced_drag);
    fprintf(fid, '}\n');
    
    % Number of particles and output file name
    fprintf(fid, 'experiment {\n');
    fprintf(fid, '\tsize = %i\n', data.N_objects);
    fprintf(fid, '\toutputFile = "%s"\n', data.filename_output_raw);
    fprintf(fid, '}\n');
    
%     fprintf(fid, '%s\n%s\n%s', string1, string2, string3);
    fclose(fid);

end

function tabella = readOutputFile(data)

    fid = fopen(data.filename_output_clean);
    out = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f', 'headerlines', 1);
%     out = textscan(fid,'%f %f %f %f %f %f %f %f %f', 'headerlines', 1);
    fclose(fid);
    
    tabella = cell2mat(out);
    
end

function [X_original, Y_original, Z_original, X_final, Y_final, Z_final] = readDEM(filename)
    
    fid_1 = fopen(filename);
    head_dem = textscan(fid_1, '%s', 12);
    fclose(fid_1);
    
    xllcorner   = round(str2double(head_dem{1}{6}));
    yllcorner   = round(str2double(head_dem{1}{8}));
    cellsize    = str2double(head_dem{1}{10});
    NODATA      = str2double(head_dem{1}{12});
    
    Z           = dlmread(filename,' ',6,0);
    Z           = Z(:,1:end-1);
    Z(Z==NODATA)= nan;
%     Z(Z==NODATA)= -9999;
    [X, Y]      = meshgrid(xllcorner:cellsize:xllcorner+(size(Z,2)-1)*cellsize,...
                           yllcorner:cellsize:yllcorner+(size(Z,1)-1)*cellsize);

%     X_final = (X-X(1))*10 + xllcorner;
%     Y_final = (Y-Y(1))*10 + yllcorner;
%     Z_final = flipud(Z);
    
    X_original = X;
    Y_original = Y;
    Z_original = flipud(Z);
      
end

function replaceCommas(FileNameOld, NewFileName)
    
    Data = fileread(FileNameOld);
    Data = strrep(Data, ',', '.');
    FID = fopen(NewFileName, 'w');
    fwrite(FID, Data, 'char');
    fclose(FID);

end

function plotResults(data)

    % Not magnified figure
    
    figure(1); 
    
    surf(data.dem_X_final,data.dem_Y_final,data.dem_Z_final); 
    shading interp;  
%     colormap(gray);
    
    hold on
    
    plot3(data.X_vent, data.Y_vent, data.Z_vent, 'pb');
    
    hold on
    
    plot3(data.output(:,1), data.output(:,2), data.output(:,3)+10, '.r');
    
    
    view([0 90])

%     % Magnified figure
%     
%     figure(2); 
%     
%     surf(data.dem_X_final,data.dem_Y_final,data.dem_Z_final); 
%     shading interp;  
%     colormap(gray);
%     
%     hold on
%     
%     plot3(data.X_vent, data.Y_vent, data.Z_vent, 'pb')
%     
%     hold on
%     plot3(data.output_iteration(1).table(:,1), data.output_iteration(1).table(:,2), data.output_iteration(1).table(:,3)+10, '.r')
%     
%     hold on
%     
%     plot3(data.table_measured(:,1), data.table_measured(:,2),data.table_measured(:,9)+10, '.c');
%     
%     xlim([4.963e5 4.97e5])
%     ylim([4.2502e6 4.2512e6])
%     view([0 90])

    
end

