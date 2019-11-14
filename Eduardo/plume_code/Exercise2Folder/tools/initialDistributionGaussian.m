function [phi, dWdPhi, D, D_centres, delta_D, dWdPhi_centres, phi_centres, ...
            volume_particle, normalisation_bin, delta_Phi, m_centres, m_extremes, vol_centres, vol_extremes] = ... 
                                                                           initialDistributionGaussian(n_bins, phi_min, phi_max, a1, b1, c1, rho)
        
    % Create phi vector
    phi = linspace(phi_min, phi_max, n_bins+1);
    dWdPhi = a1.*exp(-((phi-b1)./c1).^2);
    D = 1e-3.*2.^(-phi);                                                   % [m]
	   
    for i = 1:n_bins
        dWdPhi_average(i) = (dWdPhi(i+1)+dWdPhi(i))/2;
        phi_centres(i) = (phi(i+1)+phi(i))/2;
        delta_Phi(i) = abs(phi(i+1)-phi(i));
        normalisation_bin(i) = delta_Phi(i);
        iesim_term(i) = dWdPhi_average(i).*delta_Phi(i);
        D_centres(i) = 1e-3*2^(-phi_centres(i));                           % [m]
        delta_D(i) = abs(D(i+1)-D(i));
        volume_particle(i) = pi/6.0*D_centres(i)^3;                        % Volume associated with particle
        
    end
    
    % Mass associated with a single particle
%     m_centres = rho.*pi./6.*10.*2.^(-3.*phi_centres);
%     m_extremes = rho.*pi./6.*10.*2.^(-3.*phi);
    m_centres = rho.*pi./6.*(D_centres.^3);
    m_extremes = rho.*pi./6.*(D.^3);
    vol_centres = pi./6.*(D_centres.^3);
    vol_extremes = pi./6.*(D.^3);
    
    normalisation_value = 1/sum(iesim_term);
	dWdPhi = normalisation_value.*(a1.*exp(-((phi(i)))-b1)./c1).^2;
    dWdPhi_centres = dWdPhi_average.*normalisation_value;
    	
end