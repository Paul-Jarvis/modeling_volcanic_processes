function diff = aggregationQuantify(plume_out)
    % xi are the mass fractions as a matrix: rows = different heights;
    % columns = different sizes
    sizes = plume_out.aggregation.sizes;
    dwdphi = plume_out.aggregation.dwdphi(:,:);
    z = plume_out.z;
    
    % In this way you plot directly the mass fraction
    for i=1:length(z)
        xi(i,:) = dwdphi(i,:).*plume_out.delta_Phi;
    end
        
    
    diff = (xi(end,:)-xi(1,:))./xi(1,:)*100;
    
%     figure(98)
%     plot(sizes, diff, 'o-')
%     xlabel('phi');
%     ylabel('Mass fraction difference (%)');

end