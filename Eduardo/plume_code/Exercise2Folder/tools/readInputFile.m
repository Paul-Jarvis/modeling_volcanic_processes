 function [inp_var] = readInputFile

    fid = fopen('../input.inp');
   
%   Atmospheric conditions
    H1 = textscan(fid,'%s\n','CommentStyle','%');                          % Tropopause height [m]
    H1 = str2double(H1{:});
    inp_var.H1 = H1;
    
    H2 = textscan(fid,'%s\n','CommentStyle','%');                          % Stratopause height [m]
    H2 = str2double(H2{:});
    inp_var.H2 = H2;
    
    theta_a0 = textscan(fid,'%s\n','CommentStyle','%');                    % Atmospheric temperature at the vent [K]
    theta_a0 = str2double(theta_a0{:});
    inp_var.theta_a0 = theta_a0;
    
    tempGrad_t = textscan(fid,'%s\n','CommentStyle','%');                  % Temperature gradient troposphere [K/1000 m]
    tempGrad_t = str2double(tempGrad_t{:});
    inp_var.tempGrad_t = tempGrad_t;
    
    tempGrad_s = textscan(fid,'%s\n','CommentStyle','%');                  % Temperature gradient stratosphere [K/1000 m]
    tempGrad_s = str2double(tempGrad_s{:});
    inp_var.tempGrad_s = tempGrad_s;
    
    alpha = textscan(fid,'%s\n','CommentStyle','%');                       % Radial entrainment
    alpha = str2double(alpha{:});
    inp_var.alpha = alpha;
    
    beta = textscan(fid,'%s\n','CommentStyle','%');                        % Wind entrainment 
    beta = str2double(beta{:});
    inp_var.beta = beta;
    
    r_h = textscan(fid,'%s\n','CommentStyle','%');                         % Relative humidity [%]
    r_h = str2double(r_h{:});
    inp_var.r_h = r_h;
    
    wind_azimuth = textscan(fid,'%s\n','CommentStyle','%');                % Azimuth wind angle
    wind_azimuth = str2double(wind_azimuth{:});
    inp_var.wind_azimuth = wind_azimuth;
    
    P_0 = textscan(fid,'%s\n','CommentStyle','%');                         % Atmospheric pressure at sea level [Pa]
    P_0 = str2double(P_0{:});
    inp_var.P_0 = P_0;
    
    z_0 = textscan(fid,'%s\n','CommentStyle','%');                         % Atmosphere division profile starting (=0): [m]
    z_0 = str2double(z_0{:});
    inp_var.z_0 = z_0;
    
    z_max = textscan(fid,'%s\n','CommentStyle','%');                       % Maximum cautelative height for atmosphere (= higher than plume): [m]
    z_max = str2double(z_max{:});
    inp_var.z_max = z_max;
    
    N_points = textscan(fid,'%s\n','CommentStyle','%');                    % Number of points to divide the entire vertical atmosphere: #
    N_points = str2double(N_points{:});
    inp_var.N_points = N_points;
    
    viscAir = textscan(fid,'%s\n','CommentStyle','%');                     % Air viscosity: [Pa s]
    viscAir = str2double(viscAir{:});
    inp_var.viscAir = viscAir;
    
    viscWater = textscan(fid,'%s\n','CommentStyle','%');                   % Water viscosity: [Pa s]
    viscWater = str2double(viscWater{:});
    inp_var.viscWater = viscWater;
    
    rho_air = textscan(fid,'%s\n','CommentStyle','%');                     % Air density at the sea level for aggregates porosity: [kg/m^3]
    rho_air = str2double(rho_air{:});
    inp_var.rho_air = rho_air;
    
    B_V_freq_tropo = textscan(fid,'%s\n','CommentStyle','%');              % Brunt-Vaisala frequency f the troposphere: [1/s]
    B_V_freq_tropo = str2double(B_V_freq_tropo{:});
    inp_var.B_V_freq_tropo = B_V_freq_tropo;
    
    B_V_freq_strato = textscan(fid,'%s\n','CommentStyle','%');             % Brunt-Vaisala frequency f the stratosphere: [1/s]
    B_V_freq_strato = str2double(B_V_freq_strato{:});
    inp_var.B_V_freq_strato = B_V_freq_strato;
    
%    Eruption initial conditions
    vent_Height = textscan(fid,'%s\n','CommentStyle','%');                 % Vent height a.s.l [m]
    vent_Height = str2double(vent_Height{:});
    inp_var.vent_Height = vent_Height;
    
    u_0 = textscan(fid,'%s\n','CommentStyle','%');                         % Plume velocity [m/s]
    u_0 = str2double(u_0{:});
    inp_var.u_0 = u_0;
    
    theta_0 = textscan(fid,'%s\n','CommentStyle','%');                     % Plume temperature [K]
    theta_0 = str2double(theta_0{:});
    inp_var.theta_0 = theta_0;
    
    n_0_s = textscan(fid,'%s\n','CommentStyle','%');                       % Mass fraction solid phase
    n_0_s = str2double(n_0_s{:});
    inp_var.n_0_s = n_0_s;
    
    n_0_l = textscan(fid,'%s\n','CommentStyle','%');                       % Mass fraction liquid phase
    n_0_l = str2double(n_0_l{:});
    inp_var.n_0_l = n_0_l;
    
    n_0_d = textscan(fid,'%s\n','CommentStyle','%');                       % Mass fraction dry phase
    n_0_d = str2double(n_0_d{:});
    inp_var.n_0_d = n_0_d;
    
    n_0_v = 1 - n_0_s - n_0_l - n_0_d;                                     % Mass fraction water vapor
    inp_var.n_0_v = n_0_v;
    
    lat_volcano = textscan(fid,'%s\n','CommentStyle','%');                 % Vent latitude
    lat_volcano = str2double(lat_volcano{:});
    inp_var.lat_volcano = lat_volcano;
    
    lon_volcano = textscan(fid,'%s\n','CommentStyle','%');                 % Vent longitude
    lon_volcano = str2double(lon_volcano{:});
    inp_var.lon_volcano = lon_volcano;
    
    diffusivity = textscan(fid,'%s\n','CommentStyle','%');                 % Diffusivity in the plume
    diffusivity = str2double(diffusivity{:});
    inp_var.diffusivity = diffusivity;
    
    rho_l = textscan(fid,'%s\n','CommentStyle','%');                       % Density of liquid water in the plume [kg/m^3]
    rho_l = str2double(rho_l{:});
    inp_var.rho_l = rho_l;
    
    rho_s = textscan(fid,'%s\n','CommentStyle','%');                       % Density of solid particles in the plume [kg/m^3]
    rho_s = str2double(rho_s{:});
    inp_var.rho_s = rho_s;
    
    mm_0 = textscan(fid,'%s\n','CommentStyle','%');                        % Mass flow rate at the vent [Kg/s]
    mm_0 = str2double(mm_0{:});
    inp_var.mm_0 = mm_0;
    
    radius_vent = textscan(fid,'%s\n','CommentStyle','%');                 % Plume radius: 0 = direct mass flow rate; otherwise indirect mass flow rate 
    radius_vent = str2double(radius_vent{:});
    inp_var.radius_vent = radius_vent;
    
    omega = textscan(fid,'%s\n','CommentStyle','%');                       % 0, no condensation, 1/(2.8*3600) moderate, 1/(1.7*60) rapid omega matters
    omega = str2double(omega{:});
    inp_var.omega = omega;
    
    tmass = textscan(fid,'%s\n','CommentStyle','%');                       % Total mass erupted [Kg]
    tmass = str2double(tmass{:});
    inp_var.tmass = tmass;
    
%   Sedimentation variables and single run folder
    file_run_folder = textscan(fid,'%s\n','CommentStyle','%');             % Name of the folder for the single run
    file_run_folder = num2str(cell2mat(file_run_folder{:}));
    inp_var.file_run_folder = file_run_folder;
    
    file_wind_profile = textscan(fid,'%s\n','CommentStyle','%');           % Name of the wind profile
    file_wind_profile = num2str(cell2mat(file_wind_profile{:}));
    inp_var.file_wind_profile = file_wind_profile;
    
    file_temperature_profile = textscan(fid,'%s\n','CommentStyle','%');    % Name of the temperature profile
    file_temperature_profile = num2str(cell2mat(file_temperature_profile{:}));
    inp_var.file_temperature_profile = file_temperature_profile;
    
    file_run_extention = textscan(fid,'%s\n','CommentStyle','%');          % Extention of files in the "runs" folder
    file_run_extention = num2str(cell2mat(file_run_extention{:}));
    inp_var.file_run_extention = file_run_extention;
    
    sedimentation_model = textscan(fid,'%s\n','CommentStyle','%');         % Model: strong [0]; weak [1]
    sedimentation_model = str2double(sedimentation_model{:});
    inp_var.sedimentation_model = sedimentation_model;
    
    sedimentation_evaluation = textscan(fid,'%s\n','CommentStyle','%');    % Evaluation: base [0]; ground [1]; both [2]
    sedimentation_evaluation = str2double(sedimentation_evaluation{:});
    inp_var.sedimentation_evaluation = sedimentation_evaluation;
    
    eccentricity = textscan(fid,'%s\n','CommentStyle','%');                % Eccentricity value for strong model
    eccentricity = str2double(eccentricity{:});
    inp_var.eccentricity = eccentricity;
    
    shapeF = textscan(fid,'%s\n','CommentStyle','%');                      % Shape factor for a gravitationally spreading plume (Lambda in the paper)
    shapeF = str2double(shapeF{:});
    inp_var.shapeF = shapeF;
    
    Ngrid = textscan(fid,'%s\n','CommentStyle','%');                       % N points to divide the atmosphere within the basecurrent: #
    Ngrid = str2double(Ngrid{:});
    inp_var.Ngrid = Ngrid;
    
    map_choice = textscan(fid,'%s\n','CommentStyle','%');                  % Map choice: no sedimentation map [0]; yes sedimentation map [1]
    map_choice = str2double(map_choice{:});
    inp_var.map_choice = map_choice;
    
%   Thermodynamic parameters and constants
    C_d = textscan(fid,'%s\n','CommentStyle','%');                         % Specific heat capacity: dry air
    C_d = str2double(C_d{:});
    inp_var.C_d = C_d;
    
    C_v = textscan(fid,'%s\n','CommentStyle','%');                         % Specific heat capacity: water vapour
    C_v = str2double(C_v{:});
    inp_var.C_v = C_v;
    
    C_l = textscan(fid,'%s\n','CommentStyle','%');                         % Specific heat capacity: liquid water
    C_l = str2double(C_l{:});
    inp_var.C_l = C_l;
    
    C_s = textscan(fid,'%s\n','CommentStyle','%');                         % Specific heat capacity: solid fraction
    C_s = str2double(C_s{:});
    inp_var.C_s = C_s;
    
    L = textscan(fid,'%s\n','CommentStyle','%');                           % Latent heat vaporization
    L = str2double(L{:});
    inp_var.L = L;
    
    R_d = textscan(fid,'%s\n','CommentStyle','%');                         % Specific gas constant: dry air
    R_d = str2double(R_d{:});
    inp_var.R_d = R_d;
    
    R_v = textscan(fid,'%s\n','CommentStyle','%');                         % Specific gas constant: water vapour
    R_v = str2double(R_v{:});
    inp_var.R_v = R_v;
    
    g = textscan(fid,'%s\n','CommentStyle','%');                           % Gravitational acceleration: (m/s^2)
    g = str2double(g{:});
    inp_var.g = g;
    
    k_B = textscan(fid,'%s\n','CommentStyle','%');                         % Boltzmann constant: (J/K)
    k_B = str2double(k_B{:});
    inp_var.k_B = k_B;
           
    % Aggregation variables
    file_GS = textscan(fid,'%s\n','CommentStyle','%');                     % Name of the wind profile
    file_GS = num2str(cell2mat(file_GS{:}));
    inp_var.file_GS = file_GS;
        
    N_bins = textscan(fid,'%s\n','CommentStyle','%');           		   % Number of bins
    N_bins = str2double(N_bins{:});
    inp_var.N_bins = N_bins;
    
    phi_min = textscan(fid,'%s\n','CommentStyle','%');                     % Minimum phi size
    phi_min = str2double(phi_min{:});
    inp_var.phi_min = phi_min;
    
    phi_max = textscan(fid,'%s\n','CommentStyle','%');                     % Maximum phi size
    phi_max = str2double(phi_max{:});
    inp_var.phi_max = phi_max;
    
    rho_interstitial = textscan(fid,'%s\n','CommentStyle','%');           		   % Diameter of the monomer: [micron]
    rho_interstitial = str2double(rho_interstitial{:});
    inp_var.rho_interstitial = rho_interstitial;
            
    KF = textscan(fid,'%s\n','CommentStyle','%');                          % Prefractal constant: [#]
    KF = str2double(KF{:});
    inp_var.KF = KF;
    
    mass_smolu_box = textscan(fid,'%s\n','CommentStyle','%');      		   % Mass contained in the box for Smoluchowski: [kg]
    mass_smolu_box = str2double(mass_smolu_box{:});
    inp_var.mass_smolu_box = mass_smolu_box;
    
    gamma_shear = textscan(fid,'%s\n','CommentStyle','%');                 % Fluid laminar shear: [#]
    gamma_shear = str2double(gamma_shear{:});
    inp_var.gamma_shear = gamma_shear;
    
    DF = textscan(fid,'%s\n','CommentStyle','%');                          % Fractal dimension: [#]
    DF = str2double(DF{:});
    inp_var.DF = DF;
        
    epsilon = textscan(fid,'%s\n','CommentStyle','%');                     % Dissipation kinetic energy: [m^2/s^3]
    epsilon = str2double(epsilon{:});
    inp_var.epsilon = epsilon;
    
    stick_eff_ice = textscan(fid,'%s\n','CommentStyle','%');               % Sticking efficiency for ice: [#]
    stick_eff_ice = str2double(stick_eff_ice{:});
    inp_var.stick_eff_ice = stick_eff_ice;
        
    N_times = textscan(fid,'%s\n','CommentStyle','%');                     % Number of times 
    N_times = str2double(N_times{:});
    inp_var.N_times = N_times;
    fclose(fid);
end