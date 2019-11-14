function testDraw
    close
    data_plot_aggr = getDataAggregation;
    
    % t=0
    times = data_plot_aggr(40).times;
    sizes = data_plot_aggr(40).sizes;
    dwdphi = data_plot_aggr(40).dwdphi;
    norm = data_plot_aggr(40).normal
    %stocazzo = sum(dwdphi./norm)
    
    for i=1:length(sizes)
        pezzo(i) = dwdphi(i)*norm(i);
    end
    
%     pezzo
     sticazzi = sum(pezzo)
    
    %semilogx(sizes, dwdphi,'-')
    plot(sizes, dwdphi,'-')
    
end