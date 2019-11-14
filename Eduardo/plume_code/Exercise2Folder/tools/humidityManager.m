function [w_s, w_a, eps] = humidityManager(R_d, R_v, theta_a, P, r_h)
    % This function provides the shape of the humidiy values: it works with
    % scalar or vectors

    [w_s, w_a, eps] = Degrutyer2013(R_d, R_v, theta_a, P, r_h);
end

% Humidity value (Degrutyer, 2013)
function [w_s, w_a, eps] = Degrutyer2013(R_d, R_v, theta_a, P, r_h)

    eps = R_d/R_v;
    
    % Tetens' formula for estimation of the saturation vapor pressure
    es = 100.*6.112*exp(17.67*(theta_a-273.15)./(theta_a +243.5-273.15));  % Saturation vapor pressure
    
    % Mass ratio of water vapor to air under saturated conditions
    w_s = 1./eps .* es./(P-es);                                            % Mixing ratio (saturated)
    
    % Mass ratio of water vapor to dry air in the atmosphere
    w_a = r_h.*w_s;                                                        % Mixing ratio (not saturated)
    
end
