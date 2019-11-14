function [phi, phi_centres, D, D_centres, m_centres, m_extremes, vol_centres, vol_extremes, dWdPhi_centres, delta_Phi, delta_D] = createGrainSize(inp_var)

% NB. Mass fraction must be from the largest size to the smallest, i.e. from
% negative to positive phi

    file_grainsize = inp_var.file_GS
    file_run_extention = inp_var.file_run_extention;
    file_run_folder = inp_var.file_run_folder;
    file_grainsize = strcat('../','runs','/',file_run_folder,'/',file_grainsize,'.',file_run_extention)
    
    fileID = fopen(file_grainsize,'r');
    GSD = fscanf(fileID,'%f %f', [2 inf]);
    fclose(fileID);
    GSD = GSD';
   
    delta_phi = abs(GSD(1,1) - GSD(2,1));
    vec = GSD(:,1);
    phi = zeros(1,length(vec)+1);
    phi(1) = min(GSD(1,:))-delta_phi;
    phi(2:end) = vec;
    D = 1e-3.*2.^(-phi); 
    
    for i=1:length(GSD(:,2))        
        phi_centres(i) = (phi(i+1)+phi(i))/2;
        wf(i) = GSD(i,2)/100;
        delta_Phi(i) = abs(phi(i+1)-phi(i));
        delta_D(i) = abs(D(i+1)-D(i));
        dWdPhi_centres(i) = wf(i)/delta_Phi(i);
    end
   
    D_centres = 1e-3.*2.^(-phi_centres);
    m_centres = inp_var.rho_s.*pi./6.*(D_centres.^3);
    m_extremes = inp_var.rho_s.*pi./6.*(D.^3);
    vol_centres = pi./6.*(D_centres.^3);
    vol_extremes = pi./6.*(D.^3);

    % Invert all elements from phi negative-->positive to
    % positive-->negative
    
    phi = phi(end:-1:1);
    phi_centres =  phi_centres(end:-1:1)
    D =  D(end:-1:1); 
    D_centres =  D_centres(end:-1:1); 
    m_centres =  m_centres(end:-1:1); 
    m_extremes = m_extremes(end:-1:1); 
    vol_centres = vol_centres(end:-1:1);
    vol_extremes = vol_extremes(end:-1:1);
    dWdPhi_centres = dWdPhi_centres(end:-1:1); 
    delta_Phi = delta_Phi(end:-1:1);
    delta_D = delta_D(end:-1:1);
    
    
end