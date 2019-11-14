function data = splittingGeneralMatrix(data)
    
    DT2 = data.DT2;
    matrix_triangle_vertixes = DT2.Points;
    matrix_id_triangles = DT2.ConnectivityList;
    dim = size(matrix_id_triangles);
    n_of_triangles = dim(1);
    
    point_sum = zeros(1,2);

    % i-esim cell
    for i=1:data.n_bins
        
        % j-esim cell
        for j=1:data.n_bins
            
            % Relation on masses
            m_sum_log = log10( data.m(i) + data.m(j) );
            
            % Relation on volumes
            Chi = 1;
%             V_sum_log = log10( data.V(i) + data.V(j) + Chi*(data.V(i)+data.V(j)));
            V_sum_log = log10( data.V(i) + data.V(j) +0.00000000001*(data.V(i) + data.V(j)));
                              
            point_sum(1,1) = m_sum_log;
            point_sum(1,2) = V_sum_log;
            
            ti = pointLocation(DT2, point_sum);
            
            if isnan(ti)<1
                
                tria_vertixes = matrix_id_triangles(ti,:);                     % Here we find the three vertices of the triangle containing the point_sum
                
                vertex_1 = tria_vertixes(1);                               % First vertex
                vertex_2 = tria_vertixes(2);                               % Second vertex
                vertex_3 = tria_vertixes(3);                               % Third vertex
                                
                m1 = 10^(matrix_triangle_vertixes(vertex_1,1));
                V1 = 10^(matrix_triangle_vertixes(vertex_1,2));
                
                m2 = 10^(matrix_triangle_vertixes(vertex_2,1));
                V2 = 10^(matrix_triangle_vertixes(vertex_2,2));
                
                m3 = 10^(matrix_triangle_vertixes(vertex_3,1));
                V3 = 10^(matrix_triangle_vertixes(vertex_3,2));
                
                m_sum = 10^(m_sum_log);
                V_sum = 10^(V_sum_log);
                
                den = (V1*m2 - V2*m1 - V1*m3 + V3*m1 + V2*m3 - V3*m2);
                
                value1 = (V2*m3 - V3*m2 - V2*m_sum + V_sum*m2 + V3*m_sum - V_sum*m3) / den;
                value2 =  -(V1*m3 - V3*m1 - V1*m_sum + V_sum*m1 + V3*m_sum - V_sum*m3) / den;
                value3 = (V1*m2 - V2*m1 - V1*m_sum + V_sum*m1 + V2*m_sum - V_sum*m2) / den;
                
                if (abs(value1)<1e-3)
                    value1 = 0;
                end
                   
                if (abs(value2)<1e-3)
                    value2 = 0;
                end
                
                if (abs(value3)<1e-3)
                    value3 = 0;
                end
                data.birth(vertex_1).id1 = [data.birth(vertex_1).id1; i];
                data.birth(vertex_1).id2 = [data.birth(vertex_1).id2; j];
                data.birth(vertex_1).value = [data.birth(vertex_1).value; value1]; 
                
                data.birth(vertex_2).id1 = [data.birth(vertex_2).id1; i];
                data.birth(vertex_2).id2 = [data.birth(vertex_2).id2; j];
                data.birth(vertex_2).value = [data.birth(vertex_2).value; value2]; 
                
                data.birth(vertex_3).id1 = [data.birth(vertex_3).id1; i];
                data.birth(vertex_3).id2 = [data.birth(vertex_3).id2; j];
                data.birth(vertex_3).value = [data.birth(vertex_3).value; value3]; 
                 
            end
            
        end
        
    end
    
end
