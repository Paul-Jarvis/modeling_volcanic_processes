function wf_plot = getWf2D(m_s_bins, m, n_s, N_m, N_V, X, Y, S)

    dim = size(m_s_bins);
    N_times = dim(1);    
    wf_plot = zeros(N_times,N_m,N_V);
    
    for k=1:N_times
%         mat = fromRowToMatrix(N_p_time(k,1:end), N_m, N_V);
        m_s_bins_mat = vec2mat(m_s_bins(k,:),N_V);
        
        for i = 1:N_m
            for j=1:N_V
%                 data_plot(k,i,j) = mat(i,j)./(dx_discr(i).*dy_discr(j));
                wf_plot(k,i,j) = m_s_bins_mat(i,j)/(m(k)*n_s(k));
            end
        end
        
%         sum(sum(wf_plot(k,:,:)))
%         input('dddd')
    end
      
    
end

