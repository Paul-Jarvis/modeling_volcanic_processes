function [matrix, memoria] = preliminaryMatrix(n_bin, vol, vol_extremes)

    matrix = zeros(n_bin, n_bin, n_bin);
    idx = 1;
    
    for i=1:n_bin
        for j=1:n_bin
            for k=1:n_bin
                somma = vol(j)+vol(k);
                if ((somma>=vol_extremes(i)) && (somma<vol_extremes(i+1)))
                    matrix(j,k,i) = 1;
                    memoria(idx) = i;
                    idx = idx+1;
                else
                    matrix(j,k,i) = 0;
                end
            end
        end
    end

end