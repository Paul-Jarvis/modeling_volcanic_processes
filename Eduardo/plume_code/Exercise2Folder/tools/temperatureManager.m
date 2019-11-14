function table_temperature = temperatureManager(inp_var)

    % This function provides the shape of the temperature profile
    % H1 = tropopause height
    % H2 = stratopause height
    % z = total height from 0 to 100 km (meters)
    % theta_a0 = temperature at the vent (K)
    % temp_Grad_Tropo = temperature gradient in the troposphere (K)
    % temp_Grad_Strato = temperature gradient in the troposphere (K)
    
    theta_a0 = inp_var.theta_a0;                                           % Local variables
    temp_Grad_Tropo = inp_var.tempGrad_t;
    temp_Grad_Strato = inp_var.tempGrad_s;
    
    file_temp_profile = inp_var.file_temperature_profile;
    file_run_folder = inp_var.file_run_folder;                             % Local variables
    file_run_extention = inp_var.file_run_extention;
    
    
    z_0 = inp_var.z_0;
    z_max = inp_var.z_max;
    N_points = inp_var.N_points;
    H1 = inp_var.H1;
    H2 = inp_var.H2;
    z_temp = linspace(z_0, z_max, N_points);
    
%     table_temperature = Woods1988(theta_a0, temp_Grad_Tropo, temp_Grad_Strato, z_temp, H1, H2);
    table_temperature = readTempFromFile(file_run_folder, file_run_extention, file_temp_profile, z_temp);
    
end

% Temperature profile defined by (Woods, 1988)
function table_temperature = Woods1988(theta_a0, temp_Grad_Tropo, temp_Grad_Strato, z, H1, H2)

    for i=1:length(z)
        table_temperature(i,1) = z(i);
        if (z(i) <= H1)
            table_temperature(i,2) = theta_a0 + temp_Grad_Tropo*z(i);
        elseif (z(i)>H1) && (z(i)<=H2)
            table_temperature(i,2) = theta_a0 + temp_Grad_Tropo*H1;
        else
            table_temperature(i,2) = theta_a0 + temp_Grad_Tropo*H1 + temp_Grad_Strato*(z(i)-H2);
        end
    end

end

function table_temperature = readTempFromFile(file_run_folder, file_run_extention, file_temp_profile, z_temp)
% Problem with the number of tabs in the input file: now it is set to 1 and
% 4. But it must be canghed in 1 and 2 (so just a single tab between data)
    tempprofilepath = strcat('../','runs','/',file_run_folder,'/',file_temp_profile,'.',file_run_extention)
    temp_data = importdata(tempprofilepath); % [km , Celsius]

    % Conversions
    temp_data(:,1) = temp_data(:,1) * 1000; % Conversion from km to m of the z-coordinate
    temp_data(:,4) = temp_data(:,4) + 273.15; % Conversion from km to m of the z-coordinate
    
    table_temperature = zeros(length(z_temp),2);
    
    for i=1:length(z_temp)
        table_temperature(i,1) = z_temp(i);
        table_temperature(i,2) = interp1(temp_data(:,1),temp_data(:,4),z_temp(i),'linear','extrap'); % Wind at the baseHeight
    end

end

