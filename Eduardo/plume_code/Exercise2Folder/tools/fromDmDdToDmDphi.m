function [matrix_dwdphi] = fromDmDdToDmDphi(matrix_in, phi_centres, delta_D)

% This function first evaluates dw/dD, then it computes dw/dphi
% Matrix_in represents the mass fractions for different times and different sizes
% (row = times; columns = sizes)

    % dD/dphi
    jacobian2 = ((1e-3).*2.^(-phi_centres).*log(2));
    [times_aggr, bins_aggr] = size(matrix_in);
                                                        
    % Transformation from wt_dD a wt_dphi 
	for i=1:times_aggr
%         fraction = matrix_in(i,:)./sum(matrix_in(i,:));
        fraction = matrix_in(i,:);
        matrix_dwdphi(i,:) = fraction./delta_D.*jacobian2;
    end    
       		
end