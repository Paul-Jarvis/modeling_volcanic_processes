function data_plot_aggr = getDataAggregation
    % Output data are inside "aggregation.out" by default. This file is
    % inside the folder "aggregation"
    
    filename = '../aggregation/aggregation.out';
    
    fid = fopen(filename);                                                % Open file
    idx = 1;
    
    while (~feof(fid)) 
        riga = fscanf(fid,'%f %i',2);                                      % Number of size classes: [#]
        
        if isempty(riga) == 1
            break;
        end
        
        time = riga(1);
        bins = riga(2);
        
        for i = 1:bins
            riga = fscanf(fid,'%f %f %f',3);                                  % Particle sizes and wt% (vector Nclass x 2): [mm], [%]
            data(i,1) = riga(1);
            data(i,2) = riga(2);
            data(i,3) = riga(3);
        end
        
        data_plot_aggr(idx).times = time;
        data_plot_aggr(idx).sizes = data(:,1);
        data_plot_aggr(idx).dwdphi = data(:,2);
        data_plot_aggr(idx).normal = data(:,3);
        
        idx = idx+1;
    end
    
end