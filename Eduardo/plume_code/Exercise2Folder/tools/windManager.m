function table_wind = windManager(inp_var)

    % This function provides the shape of the wind
    % file_run_folder = file in the folder "runs" where are contained all
    % input files
    % H1 = tropopause height
    % H2 = stratopause height
    % z : in meters
    
    % NB: input file must be in [km , m/s]
    
    file_run_folder = inp_var.file_run_folder;                             % Local variables
    file_run_extention = inp_var.file_run_extention;
    file_wind_profile = inp_var.file_wind_profile;
    z_0 = inp_var.z_0;
    z_max = inp_var.z_max;
    N_points =inp_var.N_points;
    H1 = inp_var.H1;
    H2 = inp_var.H2;
    
    z_wind = linspace(z_0,z_max,N_points);
    % table_wind = jagged(z_wind, H1, H2); % Jagged profile
    table_wind = readWindFromFile(file_run_folder, file_run_extention, file_wind_profile, z_wind); 
    
    % Test check for NaN values
    [row, col] = find(isnan(table_wind));
    if (isempty(row)~=1) || (isempty(col)~=1)
        disp('NaN detected in wind input file');
    end
    
end

function table_wind = readWindFromFile(file_run_folder, file_run_extention, file_wind_profile, z_wind)
    windprofilepath = strcat('../','runs','/',file_run_folder,'/',file_wind_profile,'.',file_run_extention);
    wind_data = importdata(windprofilepath); % [km , m/s]
    wind_data(:,1) = wind_data(:,1) * 1000; % Conversion from km to m
    table_wind = zeros(length(z_wind),2);
    
    for i=1:length(z_wind)
        table_wind(i,1) = z_wind(i);
        table_wind(i,2) = interp1(wind_data(:,1),wind_data(:,2),z_wind(i),'linear','extrap'); % Wind at the baseHeight
    end

end

% Jagged profile, from Bonadonna and Phillips 2003
function table_wind = jagged(z_wind, H1, H2)
    Vmax = 20;  
    for i=1:length(z_wind)
        table_wind(i,1) = z_wind(i);
        
        if z_wind(i) <= H1
            table_wind(i,2) = Vmax*z_wind(i)/H1;
        elseif (z_wind(i)> H1) && (z_wind(i)<=H2)
            table_wind(i,2) = -0.9*Vmax*(z_wind(i)-H1)/(H2-H1) + Vmax;
        else
            table_wind(i,2) = 0.1*Vmax;
        end
                  
    end
    
end

