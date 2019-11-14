function [value, isterminal, direction] = stopHybrid(y, ...
                                                     C_d, ...
                                                     C_v, ...
                                                     C_l, ...
                                                     C_s, ...
                                                     R_d, ...
                                                     R_v, ...
                                                     r_h, ...
                                                     atmo_var, ...
                                                     rho_s, ...
                                                     rho_l)

                     
    z         = y(2);
    m_d       = y(3);
    m_v       = y(4);
    m_l       = y(5);
    m_s       = y(6);                                 
    psi       = y(7);
    angle     = y(8);
    Q         = y(9);
    P         = y(10);

    table_temperature = atmo_var.table_temperature;                        
    
    % get the temperature at the height z [K]
    theta_a = interp1(table_temperature(:,1),table_temperature(:,2),z); 

    % humidity
    cd ../tools
        [w_s, w_a, eps] = humidityManager(R_d, R_v, theta_a, P, r_h);
    cd ../plume

    % set up ODE's
    rho_aB = P/(R_v*theta_a)*(1+w_a)/(w_a+eps);
    m = m_d + m_v + m_l + m_s;
    n_s = m_s/m;                                                           % Solid mass fraction
    n_d = m_d/m;                                                           % Dry air mass fraction
    n_v = m_v/m;                                                           % Water vapor mass fraction
    n_l = m_l/m;                                                           % Liquid water mass fraction
    u = psi/m;
    C_B = (m_d*C_d + m_v*C_v + m_l*C_l + m_s*C_s)/m;
    theta = (1/C_B)*(Q/m);
    R_g = n_v*R_v/(n_v+n_d) + (n_d/(n_v+n_d))*R_d;                         % Gas constant of the mixture of gases: vapor+dry air
    rho_g = P/R_g/theta;                                                   % Density of the gas phase inside the plume
    rho_B = ((n_v+n_d)/rho_g + n_l/rho_l + n_s/rho_s)^(-1);                % Bulk density of the plume (Eduardo)

    % Locate the height when velocity passes through zero
    value = [u*sin(angle)-0.01; rho_B-rho_aB];                             % Detect velocity = 0 , the 0.001 helps for stability; detecting neutral buoyancy level
    isterminal = [1; 0];                                                   % Stop the integration
    direction = [0; 0];                                                    % detect all zeros (default)