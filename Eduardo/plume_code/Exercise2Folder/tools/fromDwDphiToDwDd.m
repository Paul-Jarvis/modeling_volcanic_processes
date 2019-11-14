function [mass_fraction] = fromDwDphiToDwDd(n_bins, M_tot, D_centres, delta_D, dWdPhi_centres, volume_particle, rho)

    % Function to transform from dw/dphi to dw/dD

    % dWdD distribution
	for i=1:n_bins
        jacobian1(i) = 1/D_centres(i)/log(2);
		dWdD_centres1(i) = dWdPhi_centres(i)*jacobian1(i);
		mass_fraction(i) = dWdD_centres1(i)*delta_D(i);
    end
       
end