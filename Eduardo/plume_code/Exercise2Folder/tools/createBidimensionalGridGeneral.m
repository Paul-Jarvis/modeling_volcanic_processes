function [Points] = createBidimensionalGridGeneral(inp_var, rho_B0, n_0_s, m_s0)

    
    % ------ We read the GS file: phi extremes, wf !!! ----- %
    % Remember: the grain-size used as input is formatted in the following
    % way:
    % - phi_extreme 1 wf1
    % - phi_extreme 2 wf2
    % ...
    % - phi_extreme N -999
    % This means that "phis" are the extreme of the bins!!!
    
%     file_grainsize = inp_var.file_GS
%     file_run_extention = inp_var.file_run_extention;
%     file_run_folder = inp_var.file_run_folder;
%     file_grainsize = strcat('../','runs','/',file_run_folder,'/',file_grainsize,'.',file_run_extention)
%     fileID = fopen(file_grainsize,'r');
%     GSD = fscanf(fileID,'%f %f', [2 inf]);
%     fclose(fileID);
%     GSD = GSD';
%     phi_extremes = GSD(:,1);
%     D_extremes = 1e-3.*2.^(-phi_extremes);
%     volumes_extremes_ini = pi./6.*D_extremes.^3;
%     mass_extremes = inp_var.rho_s*volumes_extremes_ini;
%     wf = GSD(:,2)/100;                                                     % WF normalized between 0 and 1
%     wf = wf(1:end-1);   
    
    % ----- Bidimensional scheme size: n_diagonal + n_not_diagonal ----- %
    phi(1) = 5;
    phi(2) = -5;
    n_diagonal = 100;                                                       % Phis are 1 more than wfs in the input file
    n_not_diagonal = 15; 
    n_dense = round(8/9*n_not_diagonal);                                   % Dense part of the rectangular grid
    dense_factor = 1.5;                                                    % The higher the less dense is the grid around the diagonal
    
    d_min = 1e-3*2^-(phi(1));                                              % Minimum diameter (m)
    d_max = 1e-3*2^-(phi(2));                                              % Maximum diameter (m)
        
    V_min = pi/6*d_min^3;
    V_max = pi/6*d_max^3;
 
    m_min = inp_var.rho_s*V_min;
    m_max = inp_var.rho_s*V_max;
 
    m_min_log = log10(m_min);
    m_max_log = log10(m_max);
    V_min_log = log10(V_min);
    V_max_log = log10(V_max);
    
    % Points on the diagonal
    m_vec_axis_log = linspace(m_min_log, m_max_log, n_diagonal);               % "x" coordinates of the diagonal points
    V_vec_axis_log = linspace(V_min_log, V_max_log, n_diagonal);               % "y" coordinates of the diagonal points
    
    % Dense points of the triangular grid: vec_good refers to points on the
    % triangular grid
    sufficient = 0;
    counter = 1;
    m_vec_good_log = [];                                                       % "x" coordinates of the triangular points                                    
    V_vec_good_log = [];                                                       % "y" coordinates of the triangular points
         
    while (sufficient<1)
        
        m_vec_log = m_min_log + rand(1,1)*(m_max_log-m_min_log);
        V_vec_log = V_min_log + rand(1,1)*(V_max_log-V_min_log);
        
        if ( (V_vec_log >= (m_vec_log-log10(inp_var.rho_s))) && (V_vec_log <=  (m_vec_log-log10(inp_var.rho_s)+dense_factor)) )
            m_vec_good_log(counter) = m_vec_log;
            V_vec_good_log(counter) = V_vec_log;
            counter = counter + 1;
        end
        
        if (length(m_vec_good_log) >= n_dense)
            sufficient = 1;
        end
    end
    
    % Sparse points of the triangular grid
    sufficient = 0;
        
    while (sufficient<1)
        
        m_vec_log = m_min_log + rand(1,1)*(m_max_log-m_min_log);
        V_vec_log = V_min_log + rand(1,1)*(V_max_log-V_min_log);
        
        if (V_vec_log >= (m_vec_log-log10(inp_var.rho_s)+dense_factor))
            m_vec_good_log(counter) = m_vec_log;
            V_vec_good_log(counter) = V_vec_log;
            counter = counter + 1;
        end
        
        if (length(m_vec_good_log)>=n_not_diagonal)
            sufficient = 1;
        end
    end
    
    % I recall points on the diagonal and in the grid 
    points_diagonal_log(:,1) = m_vec_axis_log;
    points_diagonal_log(:,2) = V_vec_axis_log;
    points_triangle_log(:,1) = m_vec_good_log;
    points_triangle_log(:,2) = V_vec_good_log;
    
    Points = createGrid(rho_B0, n_0_s, m_s0, points_diagonal_log, points_triangle_log);
       
%     [P, IC_m_sij] = createPoint(X, Y, N_discr_mass, N_discr_volumes, wf, rho_B0, n_0_s, m_s0);
%     IC_m_sij = IC_m_sij';
end

function data = createGrid(rho_B0, n_0_s, m_s0, points_diagonal, points_triangle)

    % points_diagonal and points_triangle are in logarithmic scale

    dim = size(points_diagonal);
    rows_diagonal = dim(1);
    columns_diagonal = dim(2);
    
    dim = size(points_triangle);
    rows_triangle = dim(1);
    columns_triangle = dim(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% Geometric creation of cells %%%%%%%%%%%%%%%%%%%%%%
    
    % Squares on the diagonal
    for i=1:length(points_diagonal(:,1))
        
        % We search for the latus of the box (useless)
        point_net_col_1 = points_diagonal(:,1);
        point_net_col_2 = points_diagonal(:,2);
        point_net_col_1([i]) = [];
        point_net_col_2([i]) = [];
        points_net(:,1) = point_net_col_1;
        points_net(:,2) = point_net_col_2;
        [k, d] = dsearchn(points_net, points_diagonal(i,:));
        latus = d/sqrt(2);
        
        % We define x and y coordinates of the vertixes of each square 
        v1(i,1) = points_diagonal(i,1)+latus/2;
        v1(i,2) = points_diagonal(i,2)+latus/2;
        v2(i,1) = points_diagonal(i,1)+latus/2;
        v2(i,2) = points_diagonal(i,2)-latus/2;
        v3(i,1) = points_diagonal(i,1)-latus/2;
        v3(i,2) = points_diagonal(i,2)-latus/2;
        v4(i,1) = points_diagonal(i,1)-latus/2;
        v4(i,2) = points_diagonal(i,2)+latus/2;
        vert_cell = [v1(i,1) v1(i,2); ...
                     v2(i,1) v2(i,2); ...
                     v3(i,1) v3(i,2); ...
                     v4(i,1) v4(i,2)];                                     % Single box coordinates
        
        % We define the polygon once the vertixes are known         
        KVert = convhulln(vert_cell); 
        cell_diag(i).table = vert_cell;                                    % Each square property
        cell_diag(i).KVert = KVert;
        
        % We take two points for each square in order to define the
        % boundaries for the internal triangulation (We need it to define
        % the constrain of the triangular sub-division of the space)
        
        external_poly((i*2)-1,1) = v4(i,1);
        external_poly((i*2)-1,2) = v4(i,2);
        external_poly(i*2,1) = v1(i,1);
        external_poly(i*2,2) = v1(i,2);                                    % Points to use later to constrain triangulation
    end
    
    % We add two points to close the polygon for the constrain
    external_poly(end+1,1) = external_poly(end,1);
    external_poly(end,2) = external_poly(end-1,2)+latus;                   % + latus
    external_poly(end+1,1) = external_poly(1,1);
    external_poly(end,2) = external_poly(end-1,2);                         % last point
    
    
    %%%%%%%%%%%%%%%% Net good points for triangulation %%%%%%%%%%%%%%%%%%%%
    
    % We define the points for the triangulation: boundary points+random
    % points. But we have the problem that some random points can fall
    % inside the boxes on the diagonal. These points must be avoided.
    
    % First points are the boundaries: from 1 to fin
    points_triangle_more(:,1) = external_poly(:,1);
    points_triangle_more(:,2) = external_poly(:,2);
    fin = length(points_triangle_more(:,1));
    points_triangle_more(fin+1:fin+rows_triangle,1) = points_triangle(:,1);
    points_triangle_more(fin+1:fin+rows_triangle,2) = points_triangle(:,2);
    
    % Points_triangle_more contains all the points on the triangular grid: boundary + random
    for i=1:fin
        constrains(i,1)  = i;
        constrains(i,2)  = i+1;
    end
    
    constrains(end,1) = fin;
    constrains(end,2) = 1;
    
    % Now we create our bidimensional grid: we will use this grid to define
    % the fixed pivot. We have still the problem that some points are
    % inside the boxes on the diagonal
    DT = delaunayTriangulation(points_triangle_more, constrains);                      % Delaunay triangulation   
    matrix_vertixes = DT.ConnectivityList;                                 % Matrix: n_triangles * 3
    dim = size(matrix_vertixes);
    n_of_triangles = dim(1);                                               % Number of triangles
    IO = isInterior(DT);
    data.DT = DT;
    data.DT_IO = IO;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%    Actual grid used in the code    %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Final data: values along the diagonal
    
    wf = rand(1,rows_diagonal); 
    S = sum(wf); 
    wf = wf/S; % renormalisation
    
    for i=1:rows_diagonal
        data.log_m(i) = points_diagonal(i,1);
        data.log_V(i) = points_diagonal(i,2);
        data.m(i) = 10^(data.log_m(i));
        data.V(i) = 10^(data.log_V(i));
        data.area(i) = polyarea(cell_diag(i).table(:,1),cell_diag(i).table(:,2));   
        data.birth(i).id1 = [];
        data.birth(i).id2 = [];
        data.birth(i).value = [];
%         data.wf(i) = 0;
        data.wf(i) = wf(i);
        data.N0(i) = rho_B0*n_0_s*data.wf(i)/data.m(i);
        data.N(i) = 0;
        data.m_s_ij(i) = m_s0*data.wf(i);
        data.vertixes(i).table = cell_diag(i).table;
        data.diagonal = i;
    end
    
    data.m_min = data.m(1);
    data.m_max = data.m(end);
    data.V_min = data.V(1);
    data.V_max = data.V(end);
    
    % Centroid of each triangle
    counter = rows_diagonal+1;
    internal_counter = 1;
    
    for i=1:length(IO)
        if IO(i)>0
            id_vertixes = DT.ConnectivityList(i,:);
            v1 = points_triangle_more(id_vertixes(1),:);
            v2 = points_triangle_more(id_vertixes(2),:);
            v3 = points_triangle_more(id_vertixes(3),:);
            x = [v1(1,1) v2(1,1) v3(1,1)];
            y = [v1(1,2) v2(1,2) v3(1,2)];
            table(:,1) = x;
            table(:,2) = y;
                        
            data.log_m(counter) = mean(x);
            data.log_V(counter) = mean(y);
            data.m(counter) = 10^(data.log_m(counter));
            data.V(counter) = 10^(data.log_V(counter));
            data.area(counter) = polyarea(x,y);   
            data.birth(counter).id1 = [];
            data.birth(counter).id2 = [];
            data.birth(counter).value = [];   
            data.wf(counter) = 0;
            data.m_s_ij(counter) = m_s0*data.wf(counter);
            data.N0(counter) = rho_B0*n_0_s*data.wf(counter)/data.m(counter);
            data.N(counter) = 0;
            data.vertixes(counter).table = table;
            data.triangles = internal_counter;
            counter = counter+1;
            internal_counter = internal_counter+1;
        else
            id_vertixes = DT.ConnectivityList(i,:);
            v1 = points_triangle_more(id_vertixes(1),:);
            v2 = points_triangle_more(id_vertixes(2),:);
            v3 = points_triangle_more(id_vertixes(3),:);
            x = [v1(1,1) v2(1,1) v3(1,1)];
            y = [v1(1,2) v2(1,2) v3(1,2)];
        end
    end
    
    data.n_bins = data.diagonal+data.triangles;
    data.n_bins_diagonal = data.diagonal;
    data.n_bins_triangles = data.triangles;
        
%     data.wf(3) = 1;
%     data.m_s_ij(3) = m_s0*data.wf(3);
  
    tab(:,1) = data.log_m;
    tab(:,2) = data.log_V;
    
    DT2 = delaunayTriangulation(tab);
    data.DT2 = DT2;
    
end