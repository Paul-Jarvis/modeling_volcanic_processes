function [X, Y, dx_discr, dy_discr, x_axis_discr, y_axis_discr, x_discr, ...
          y_discr, phi_extremes, phi_c, D_extremes, D_c, wf, P, IC_m_sij, volumes_c_phi] = createBidimensionalGrid(inp_var, rho_B0, n_0_s, m_s0, Chi)

    phi_infinite = -10;                                                    % Largest phi admitted for the volume
    
    % ------ We read the GS file: phi extremes, wf !!! ----- %
    % Remember: the grain-size used as input is formatted in the following
    % way:
    % - phi_extreme 1 wf1
    % - phi_extreme 2 wf2
    % ...
    % - phi_extreme N -999
    % This means that "phis" are the extreme of the bins!!!
    
    file_grainsize = inp_var.file_GS
    file_run_extention = inp_var.file_run_extention;
    file_run_folder = inp_var.file_run_folder;
    file_grainsize = strcat('../','runs','/',file_run_folder,'/',file_grainsize,'.',file_run_extention)
    fileID = fopen(file_grainsize,'r');
    GSD = fscanf(fileID,'%f %f', [2 inf]);
    fclose(fileID);
    GSD = GSD';
    phi_extremes = GSD(:,1);
    D_extremes = 1e-3.*2.^(-phi_extremes);
    volumes_extremes_ini = pi./6.*D_extremes.^3;
    mass_extremes = inp_var.rho_s*volumes_extremes_ini;
    wf = GSD(:,2)/100;                                                     % WF normalized between 0 and 1
    wf = wf(1:end-1);   
    
    % ----- Bidimensional scheme size: N_m * N_V ----- %
    
    N_discr_mass = length(wf);                                             % Phis are 1 more than wfs in the input file
    N_discr_volumes = N_discr_mass+0;
    
    % ----- Vector of masses definition ----- %
        
    for i = 1:N_discr_mass
        delta_phi(i) = abs(phi_extremes(i+1)-phi_extremes(i));
        phi_c(i) = (phi_extremes(i+1)+phi_extremes(i))/2;
        D_c(i) = 1e-3.*2.^(-phi_c(i));                                        % [m]
        delta_D(i) = abs(D_extremes(i+1)-D_extremes(i));
        volumes_ini_c(i) = pi./6*D_c(i)^3;
        mass_c(i) = inp_var.rho_s*volumes_ini_c(i);
        delta_mass(i) = abs(mass_extremes(i+1)-mass_extremes(i));
    end
    
    % ----- Vector of volumes definition ----- %
    
    D_infinite = 1e-3*2^(-phi_infinite);                                   % [m]
    volume_infinite = pi/6*D_infinite^3;
    volumes_beyond_extremes = logspace(log10(volumes_extremes_ini(end)), log10(volume_infinite), N_discr_volumes-N_discr_mass+1);
    volumes_extremes_ini = volumes_extremes_ini';
    volumes_beyond_extremes = volumes_beyond_extremes;
    
    volumes_extremes = [volumes_extremes_ini volumes_beyond_extremes(2:end)];
    
    
    for i=1:N_discr_volumes
        
        if (i<=N_discr_mass)
            volumes_c(i) = volumes_ini_c(i);
            delta_volumes(i) = abs(volumes_extremes(i+1)-volumes_extremes(i));
        else
            volumes_c(i) = (volumes_extremes(i+1)+volumes_extremes(i))/2;
            delta_volumes(i) = abs(volumes_extremes(i+1)-volumes_extremes(i)); 
        end
        
%         if (i<=N_discr_mass)
%             volumes_c(i) = volumes_ini_c(i);
%             delta_volumes(i) = abs(volumes_extremes_ini(i+1)-volumes_extremes_ini(i));
%         else
%             volumes_c(i) = (volumes_beyond_extremes(i+1)+volumes_beyond_extremes(i))/2;
%             delta_volumes(i) = abs(volumes_beyond_extremes(i+1)-volumes_beyond_extremes(i)); 
%         end
        
    end
    d_vol_c = (6.*volumes_c/pi).^(1./3);
    volumes_c_phi = -log2(d_vol_c./1e-3);
       
    % ----- Bidimensional grid creation ----- %
    
    [X, Y, dx_discr, dy_discr, x_axis_discr, y_axis_discr, x_discr, y_discr] = grid(mass_c, volumes_c, delta_mass, delta_volumes, mass_extremes, volumes_extremes);
    
    % ----- Points generation ----- %
    
    [P, IC_m_sij] = createPoint(X, Y, N_discr_mass, N_discr_volumes, wf, rho_B0, n_0_s, m_s0);
    IC_m_sij = IC_m_sij';
end

function [X, Y, dx_discr, dy_discr, x_axis_discr, y_axis_discr, x_discr, y_discr] = grid(mass_c, volumes_c, delta_mass, delta_volumes, mass_extremes, volumes_extremes)
  
    
    % Creation of the bidimensional grid
    x_discr = mass_c;
    y_discr = volumes_c;
    N_discr_x = length(x_discr);
    N_discr_y = length(y_discr);
    x_axis_discr = mass_extremes;
    y_axis_discr = volumes_extremes;
           
    dx_discr = delta_mass;
    dy_discr = delta_volumes;
    
    [X, Y] = meshgrid(x_discr,y_discr);
    
    X = X';
    Y = Y';    

end


function [P, IC_m_sij] = createPoint(m_i, V_j, N_i, N_j, wf, rho_B0, n_0_s, m_s0)

% Here we create the point structure. 

% Each point has the following structure:
% - i = mass coordinate (from 1 to N_i)
% - j = volume coordinate (from 1 to N_j)
% - P(i,j).m = mass related to the cell
% - P(i,j).V = volume related to the cell
% - P(i,j).rho = density related to the cell
% - P(i,j).N_p = number of particles per unit volume in each cell

    % Initialization of the structure
    
    P(N_i,N_j).m = 0;
    P(N_i,N_j).rho = 0;
    P(N_i,N_j).wf = 0;
    P(N_i,N_j).V = 0;
    P(N_i,N_j).N_p = 0;
    P(N_i,N_j).true = 0;
    P(N_i,N_j).m_s_ij = 0;
    IC_m_sij = zeros(N_i*N_j,1);
    
    contatore = 1;
    
    % We initialize weight fractions on the diagonal
    
    for i = 1:N_i
        for j = 1:N_j
            
            if (j==i)
                P(i,j).m = m_i(i,j);
                P(i,j).V = V_j(i,j);
                P(i,j).rho = P(i,j).m/P(i,j).V;
                P(i,j).N_p = rho_B0*n_0_s*wf(i)/m_i(i);
                P(i,j).wf = wf(i);
                P(i,j).m_s_ij = m_s0*wf(i);
                
            else
                P(i,j).m = m_i(i,j);
                P(i,j).V = V_j(i,j);
                P(i,j).rho = P(i,j).m/P(i,j).V;
                P(i,j).N_p = 0;
                P(i,j).wf = 0;
                P(i,j).m_s_ij = 0;
                
            end
            IC_m_sij(contatore) = P(i,j).m_s_ij;
%             IC_m_sij(contatore) = P(i,j).wf;
            contatore = contatore+1;
        end
    end
    
end