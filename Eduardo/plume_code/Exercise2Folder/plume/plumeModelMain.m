function [plume_out] = plumeModelMain(inp_var, atmo_var)
    
    g = inp_var.g;                  C_d = inp_var.C_d;                     % Local variables: from input file
    C_v = inp_var.C_v;              C_l = inp_var.C_l;
    C_s = inp_var.C_s;              R_d = inp_var.R_d;
    R_v = inp_var.R_v;              L = inp_var.L;
    rho_l = inp_var.rho_l;          rho_s = inp_var.rho_s;
    vent_Height = inp_var.vent_Height;    H1 = inp_var.H1;
    H2 = inp_var.H2;                theta_a0 = inp_var.theta_a0 ;
    P_0 = inp_var.P_0;              tempGrad_t = inp_var.tempGrad_t;
    tempGrad_s = inp_var.tempGrad_s;    alpha = inp_var.alpha;
    beta = inp_var.beta;            u_0 = inp_var.u_0;
    theta_0 = inp_var.theta_0;      n_0_s = inp_var.n_0_s;
    n_0_d = inp_var.n_0_d;          n_0_l = inp_var.n_0_l;
    n_0_v = inp_var.n_0_v;          omega = inp_var.omega;
    mm_0 = inp_var.mm_0;            r_h = inp_var.r_h;
    z_max = inp_var.z_max;          radius_vent = inp_var.radius_vent;
    viscAir = inp_var.viscAir;      mu_kyn_air = viscAir/inp_var.rho_air;   
    k_B = inp_var.k_B;   
        
    table_wind = atmo_var.table_wind;                                      % Local variables: form atmosphere
    table_temperature = atmo_var.table_temperature;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                %%%%% Initial conditions for plume %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Gas phase density: vapor+dry
    R_g0 = n_0_v*R_v/(n_0_v+n_0_d) + (n_0_d/(n_0_v+n_0_d))*R_d;            % Gas constant of the mixture of gases: vapor+dry air
    rho_g0 = P_0/R_g0/theta_0;                                             % Density of the gas phase inside the plume
    rho_B0 = ((n_0_v+n_0_d)/rho_g0 + n_0_l/rho_l + n_0_s/rho_s)^(-1)       % Bulk density
        
radius_vent
    
    % Evaluation of the mass flow rate: if radius_vent = 0, direct evaluation of the MFR
    if (radius_vent == 0)
        disp('Mass flux rate from user (SEE INPUT FILE)')
        m_0 = mm_0;
    else
        disp('Mass flux rate estimated (IGNORING INPUT FILE)')
%         radius_vent = (mm_0/(pi*rho_B0*u_0))^(1/2)
        mm_0 = pi*rho_B0*u_0*radius_vent^2;
        m_0 = mm_0;
    end    
    
    disp('Initial mass flux rate')
    m_0
    
    % Initial flux rates at the vent (kg/s)
    m_0 = m_0/pi; 
    m_v0 = n_0_v*m_0;                   
    m_l0 = n_0_l*m_0;
    m_s0 = n_0_s*m_0;                   
    m_d0 = n_0_d*m_0;     
    psi_0 = m_0*u_0;                    
    C_B0 = (m_d0*C_d + m_v0*C_v + m_l0*C_l + m_s0*C_s)/m_0;
    Q_0 = m_0*C_B0*theta_0;
    angle_0 = 0.5*pi;
    x_0 = 0;      
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% Initial conditions %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    IC = [x_0 vent_Height m_d0 m_v0 m_l0 m_s0 psi_0 angle_0 Q_0 P_0];      % initial conditions of the Odes
    Sspan = [vent_Height:50:z_max];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                %%%%% solve system of ODE's %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    % Here we set the error tolerances for the ODE and the stop point for
    % velocity = 0 m/s (see stopPlume.m file) 'NonNegative', 6:6+n_bins-1, 
    options = odeset('RelTol',1e-5,'AbsTol',1e-5,... 
                     'Events',@(z,y) stopHybrid(y, ...
                                                C_d, ...
                                                C_v, ...
                                                C_l, ...
                                                C_s, ...
                                                R_d, ...
                                                R_v, ...
                                                r_h, ...
                                                atmo_var, ...
                                                rho_s, ...
                                                rho_l));                                            
    disp('ODE time');
    tic
    
    % Here we solve the ODE
    [S,Y,TE,YE,IE] =  ode45(@(s,y) odeHybrid(y, ...
                                             g, ...
                                             C_d, ...
                                             C_v, ...
                                             C_l, ...
                                             C_s, ...
                                             R_d, ...
                                             R_v, ...
                                             L,   ...
                                             H1,    ...
                                             H2,    ...
                                             theta_a0, ...
                                             tempGrad_t, ...
                                             tempGrad_s, ...
                                             alpha,      ...
                                             beta,       ...
                                             omega,      ...
                                             r_h, ...
                                             atmo_var, ...
                                             viscAir, ...
                                             mu_kyn_air, ...
                                             k_B, ...
                                             rho_s, ...
                                             rho_l), Sspan, IC, options);
    toc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %%%%% Quantities evaluated for each time %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Order of the solutions
    % [x_0 vent_Height m_d0 m_v0 m_l0 m_s0_frac psi_0 angle_0 Q_0 P_0];
    
    x         = Y(:,1);                                                    % Y(1) --> x coordinate of the plume
    z         = Y(:,2);                                                    % Y(2) --> z coordinate of the plume
    m_d       = pi*Y(:,3);                                                 % Y(3) --> Variation of the mass flow rate: dry air
    m_v       = pi*Y(:,4);                                                 % Y(4) --> Variation of the mass flow rate: water vapor
    m_l       = pi*Y(:,5);                                                 % Y(5) --> Variation of the mass flow rate: liquid water
    m_s       = pi*Y(:,6);                                                 % Y(6) --> Variation of the momentum
    psi       = pi*Y(:,7);                                                 % Y(7) --> Variation of the momentum
    angle     = Y(:,8);                                                    % Y(8) --> Variation of the angle of the plume
    Q         = pi*Y(:,9);                                                 % Y(9) --> Variation of the heat flow rate
    P         = Y(:,10);                                                   % Y(10) --> Variation of the atmospheric pressure with altitude
    
    height_t  = z(end)
    eps = R_d/R_v;
    m = m_d + m_v + m_l + m_s;                                             % Mass flow rate (kg/s)
    n_d = m_d./m;
    n_v = m_v./m;
    n_l = m_l./m;
    n_s = m_s./m;
    u = psi./m;
    C_B = n_d.*C_d + n_v.*C_v + n_l.*C_l + n_s.*C_s;                       % Heat capacity of the plume
    theta = (1./C_B).*(Q./m);
    R_g = n_v.*R_v./(n_v+n_d) + (n_d./(n_v+n_d)).*R_d;                     % Gas constant of the mixture of gases: vapor+dry air
    rho_g = P./R_g./theta;                                                 % Density of the gas phase inside the plume
    rho_B = ((n_v+n_d)./rho_g + n_l./rho_l + n_s./rho_s).^(-1);            % Bulk density of the plume
    r = (m./(pi*rho_B.*u)).^(1/2);
    
    Q_V = pi*m./rho_B;                                                     % Volume flow rate (m^3/s)
    density_aggr = rho_B.*n_s;
    phi_s = m_s./(rho_s.*u.*r.^2);
     
    % Evaluation of atmospheric values: wind
    
    ambT = zeros(size(S));  
    Va = zeros(size(S));
    
    for i=1:length(ambT)
        ambT(i) = interp1(table_temperature(:,1),table_temperature(:,2),z(i)); 
        Va(i) = interp1(table_wind(:,1),table_wind(:,2),z(i));                 
    end
    
    % Evaluation of atmospheric values: humidity
    
    cd ../tools
        [w_s, w_a, eps] = humidityManager(R_d, R_v, ambT, P, r_h);
    cd ../plume
    rho_aB    = P./(R_v.*ambT).*(1+w_a)./(w_a+eps);                        % Environmental density
        


%     % Preparing the output

    plume_out.s = S;
    plume_out.col_Height = z(end);                                         % Maximum height of the plume (beyond NBL): [m]
    
    if length(TE)>1
        plume_out.nbl_Height = TE(2);                                          % Height of the NBL: [m]
    else
        plume_out.nbl_Height = TE(1);                                          % Height of the NBL: [m]
    end
    plume_out.nbl_r = interp1(z, r, plume_out.nbl_Height); 
%     plume_out.volFlowRate = Q_V_nbl;                                       % Volume flow rate at the NBL: [m^3/s]
    plume_out.u = u;                                                       % Velocity of the plume at different heights: [m/s]
    plume_out.m = m;
    plume_out.angle = angle;                                               % Angle of the plume(?)
    plume_out.z = z;                                                       % Total heights reached by the plume in the ODE: [m]
    plume_out.x = x;                                                       % Total heights reached by the plume in the ODE: [m]
%     plume_out.z_partial = z_verti_plot;                                    % Partial heights starting from idx_start: [m]
%     plume_out.rho_B = rho_B;                                               % Density of the plume: [kg/m^3]
%     plume_out.rho_aB = rho_aB;                                             % Density of the atmosphere: [kg/m^3]
%     plume_out.time_seconds = time_seconds;                                 % Different times at which the plume reaches different heights: [s]
%     plume_out.mass_dry = m_d;                                              % Mass of dry air contained in the plume: [kg]
%     plume_out.mass_vapor = m_v;                                            % Mass of the vapor contained in the plume: [kg]
%     plume_out.mass_liquid = m_l;                                           % Mass of the liquid fraction contained in the plume: [kg]
%     plume_out.mass_solid = m_s;                                            % Mass of the solid fraction contained in the plume: [kg]
    plume_out.vertical_velox = u.*sin(angle);                           % Vertical component of the plume velocity starting from z = idx_start: [m/s]
    plume_out.nbl_u = interp1(z, plume_out.vertical_velox, plume_out.nbl_Height); 
%     plume_out.temperature = theta;                                         % Temperature of the plume starting from z = idx_start: [K] 
    plume_out.r = r;                                                       % Radius of the plume starting from z = idx_start: [K]
    plume_out.density = rho_B;                                             % Density of the plume starting from z = idx_start: [kg/m^3]
    plume_out.ambientdensity = rho_aB;                                     % Density of the plume starting from z = idx_start: [kg/m^3]
    plume_out.ambientwind = Va;                                            % Wind velocity starting from z = idx_start: [m/s]
    plume_out.ambientT = ambT;                                             % Ambient temperature starting from z = idx_start: [K]
%     plume_out.m_d = m_d;                                                   % Mass of dry air inside the plume at different z-heights: [Kg/s]
%     plume_out.m_v = m_v;                                                   % Mass of water vapor inside the plume at ALL different z-heights: [Kg/s]
%     plume_out.m_l = m_l;                                                   % Mass of liquid water inside the plume at ALL different z-heights: [Kg/s]
%     plume_out.m_s = m_s;                                                   % Mass of soli fraction inside the plume at ALL different z-heights: [Kg/s]
    plume_out.theta = theta;                                               % Plume temperature at ALL different z-heights: [K]

    plume_out.n_gas = n_d+n_v;
    plume_out.n_s = n_s;
end