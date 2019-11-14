function dy_ds = odeHybrid(y, ...
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
                           rho_l)
                                   
    table_wind = atmo_var.table_wind;                                      % Local variables
    table_temperature = atmo_var.table_temperature;
    
    x         = y(1);
    z         = y(2);
    m_d       = y(3);                                                      % Mass flux rate dry air
    m_v       = y(4);                                                      % Mass flux rate water vapor
    m_l       = y(5);                                                      % Mass flux rate liquid water
    m_s       = y(6);                                                      % Mass flux rates solid 
    psi       = y(7);                                                      % Momentum equation
    angle     = y(8);                                                      % Angle equation
    Q         = y(9);                                                      % Heat equation
    P         = y(10);                                                     % Atmospheric pressure
    
    % _____________________________________________________________________    
    % Get the wind value at a given point z: in this step z is a scalar
    Va = interp1(table_wind(:,1),table_wind(:,2),z);
    
    % Get the atmosphere temperature at a given point z
    theta_a = interp1(table_temperature(:,1),table_temperature(:,2),z);
    
    % Get the humidity value at a given temperature value
    cd ../tools
        [w_s, w_a, eps] = humidityManager(R_d, R_v, theta_a, P, r_h);
    cd ../plume
    % _____________________________________________________________________
    
    % set up ODE's
    rho_aB    = P/(R_v*theta_a)*(1+w_a)/(w_a+eps);
    rhophi_av = P/(R_v*theta_a)*w_a/(w_a+eps);                             % Atmosphere: volume fraction of water vapor * rho_atmosphere
    rhophi_ad = P/(R_v*theta_a)*1/(w_a+eps);                               % Atmosphere: volume fraction of dry air * rho_atmosphere
    C_aB      = (C_d + w_a*C_v)/(1+w_a);                                   % Atmosphere: heat capacity
    m         = m_d + m_v + m_l + m_s;                                     % Total mass at each step
    n_s       = m_s/m;                                                     % Solid mass fraction
    n_d       = m_d/m;                                                     % Dry air mass fraction
    n_v       = m_v/m;                                                     % Water vapor mass fraction
    n_l       = m_l/m;                                                     % Liquid water mass fraction
    u         = psi/m;                                                     % Velocity of the plume
    C_B       = (m_d*C_d + m_v*C_v + m_l*C_l + m_s*C_s)/m;                 % Heat capacity of the plume
    theta     = (1/C_B)*(Q/m);                                             % Temperature
    R_g       = n_v*R_v/(n_v+n_d) + (n_d/(n_v+n_d))*R_d;                   % Gas constant of the mixture of gases: vapor+dry air
    rho_g     = P/R_g/theta;                                               % Density of the gas phase inside the plume
    rho_B = ((n_v+n_d)/rho_g + n_l/rho_l + n_s/rho_s)^(-1);                % Bulk density of the plume (Eduardo)
    r         = (m/(rho_B*u))^(1/2);                                       % Radius of the plume
    rhophi_v  = (P/R_v/theta)*m_v/(m_v+eps*m_d);
   
    % ENTRAINMENT ASSUMPTION
    u_eps = alpha*abs(u-Va*cos(angle)) + beta*abs(Va*sin(angle));          % Entrainment velocity
   
         
    % ODE EQUATIONS
    dx_ds     = cos(angle);                                                % equation S26
    dz_ds     = sin(angle);                                                % equation S27
    dm_d_ds   = 2*u_eps*r*rhophi_ad;                                       % equation S1
    dm_v_ds   = 2*u_eps*r*rhophi_av - omega*rhophi_v*r^2;                  % equation S2
    dm_l_ds   = omega*rhophi_v*r^2;                                        % equation S3
    dm_s_ds   = 0;                                                  
    dm_ds     = dm_d_ds + dm_v_ds + dm_l_ds + dm_s_ds;
    dpsi_ds   = g*(rho_aB - rho_B)*r^2*sin(angle) + Va*cos(angle)*dm_ds;   % equation S5
    dangle_ds = (1/psi)*(g*(rho_aB - rho_B)*r^2*cos(angle) - Va*sin(angle)*dm_ds); % equation S6
    dQ_ds     = C_aB*theta_a*dm_ds -rho_B*u*r^2*g*sin(angle) + L*dm_l_ds;  % equation S7
    dP_ds     = -rho_aB*g*sin(angle);                                      % equation S15
    dy_ds    = zeros(length(y),1);
    dy_ds(1) = dx_ds;                                                      % Y(1) --> x coordinate of the plume
    dy_ds(2) = dz_ds;                                                      % Y(2) --> z coordinate of the plume
    dy_ds(3) = dm_d_ds;                                                    % Y(3) --> Variation of the mass flow rate: dry air
    dy_ds(4) = dm_v_ds;                                                    % Y(4) --> Variation of the mass flow rate: water vapor
    dy_ds(5) = dm_l_ds;                                                    % Y(5) --> Variation of the mass flow rate: liquid water
    dy_ds(6) = dm_s_ds;                                                    % Y(6) --> Variation of solid mass 
    dy_ds(7) = dpsi_ds;                                                    % Y(...) --> Variation of the momentum
    dy_ds(8) = dangle_ds;                                                  % Y(...) --> Variation of the angle of the plume
    dy_ds(9) = dQ_ds;                                                      % Y(...) --> Variation of the heat flow rate
    dy_ds(10) = dP_ds;                                                     % Y(...) --> Variation of the atmospheric pressure with altitude

end
