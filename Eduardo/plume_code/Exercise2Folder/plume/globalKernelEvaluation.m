function [kern_out, v_ij] = globalKernelEvaluation(dlist, V_settl, epsilon, mu_kyn_air, k_B, T_0, g, boolean_vector, rho_p, rho_f, V_plume)
    
    idx = length(dlist);
    gamma_shear_turb = sqrt(1.3*epsilon/mu_kyn_air);
    kern_out = zeros(idx,idx);   
    cost = 2*k_B*T_0/(3*mu_kyn_air);     
    gamma_shear_turb = sqrt(1.3*epsilon);
    V_plume = 0.001*V_plume;
    
    for i=1:idx
        rel_time_i = rho_p(i)*dlist(i)^2/(18*mu_kyn_air*rho_f);
                
        for j=i:idx
            rel_time_j = rho_p(j)*dlist(j)^2/(18*mu_kyn_air*rho_f);
            
            % ________________ FALL 3D ____________________________________
            brownian_term = cost * (dlist(i) + dlist(j))^2 / (dlist(i)*dlist(j));                  
            sedimentation_term = pi/4*(dlist(i)+dlist(j))^2*abs(V_settl(i)-V_settl(j));   
            shear_turbulent = gamma_shear_turb/6 * (dlist(i)+dlist(j))^3;
           
            v_ij(i,j) = (brownian_term+sedimentation_term+shear_turbulent)*4./(pi*(dlist(i)+dlist(j))^2);
            turbulent_inertial = pi*epsilon^(3/4) / (g*mu_kyn_air^(1/4)) * (dlist(i)/2+dlist(j)/2)^2 * abs(v_ij(i,j));         
%             turbulent_inertial = pi*epsilon^(3/4) / (g*mu_kyn_air^(1/4)) * (dlist(i)/2+dlist(j)/2)^2 * abs(V_settl(i)-V_settl(j));         
%             turbulent_inertial = epsilon^(3/4)/(mu_kyn_air^(1/4)) * (dlist(i)/2+dlist(j)/2)^2 * abs(rel_time_i-rel_time_j);
            


%             % _____________ SEDIMENTATION AS A PLUME FUNCTION ___________
%             % Pay attention to the right "boolean vector" (must be defined
%             % in odeHybrid
%             
%             if (boolean_vector(i) < 1) && (boolean_vector(j) < 1)
%                 % No gravitational collection
%                 brownian_term = cost * (dlist(i) + dlist(j))^2 / (dlist(i)*dlist(j));                  
%                 sedimentation_term = pi/4*(dlist(i)+dlist(j))^2*abs(V_settl(i)-V_settl(j));   
%                 shear_turbulent = gamma_shear_turb/6 * (dlist(i)+dlist(j))^3;
%            
%                 v_ij(i,j) = (brownian_term+sedimentation_term+shear_turbulent)*4./(pi*(dlist(i)+dlist(j))^2);
%                 turbulent_inertial = pi*epsilon^(3/4) / (g*mu_kyn_air^(1/4)) * (dlist(i)/2+dlist(j)/2)^2 * abs(v_ij(i,j));         
%                 sedimentation_term = 0;
%             else
%                 % Yes gravitational collection
%                 brownian_term = cost * (dlist(i) + dlist(j))^2 / (dlist(i)*dlist(j));                  
%                 sedimentation_term = pi/4*(dlist(i)+dlist(j))^2*abs(V_settl(i)-V_settl(j));   
%                 shear_turbulent = gamma_shear_turb/6 * (dlist(i)+dlist(j))^3;
%            
%                 v_ij(i,j) = (brownian_term+sedimentation_term+shear_turbulent)*4./(pi*(dlist(i)+dlist(j))^2);
%                 turbulent_inertial = pi*epsilon^(3/4) / (g*mu_kyn_air^(1/4)) * (dlist(i)/2+dlist(j)/2)^2 * abs(v_ij(i,j));         
%             end
            
            
            
            % ________________ KINETIC KERNEL _____________________________
            % Pay attention to the right "boolean vector" (must be defined
            % in odeHybrid
            
%             if (boolean_vector(i) < 1) && ((boolean_vector(j) < 1))
%                 % Saffman-Turner theory
%                 
%                 brownian_term = cost * (dlist(i) + dlist(j))^2 / (dlist(i)*dlist(j));                  
%                 sedimentation_term = pi/4*(dlist(i)+dlist(j))^2*abs(V_settl(i)-V_settl(j));   
%                 shear_turbulent = gamma_shear_turb/6 * (dlist(i)+dlist(j))^3;
%                 
%                 % Friedlander
% %                 turbulent_inertial = epsilon^(3/4)/(mu_kyn_air^(1/4)) * (dlist(i)/2+dlist(j)/2)^2 * abs(rel_time_i-rel_time_j);
% 
%                 % Jacobson
%                 turbulent_inertial = pi*epsilon^(3/4) / (g*mu_kyn_air^(1/4)) * (dlist(i)/2+dlist(j)/2)^2 * abs(V_settl(i)-V_settl(j));         
%                 kinetic_term = 0;
% 
%             elseif (boolean_vector(i) > 1) && ((boolean_vector(j) < 1))
%                 % Abrahamson theory
%                 brownian_term = cost * (dlist(i) + dlist(j))^2 / (dlist(i)*dlist(j));         
%                 sedimentation_term = pi/4*(dlist(i)+dlist(j))^2*abs(V_settl(i)-V_settl(j));   
%                 shear_turbulent = 0;
%                 turbulent_inertial = 0;
%                 kinetic_term = 2^(3/2)/4 * sqrt(pi) * V_plume * (dlist(i)+dlist(j))^2;
%             
%             elseif (boolean_vector(i) < 1) && ((boolean_vector(j) > 1))
%                 % Abrahamson theory
%                 brownian_term = cost * (dlist(i) + dlist(j))^2 / (dlist(i)*dlist(j));         
%                 sedimentation_term = pi/4*(dlist(i)+dlist(j))^2*abs(V_settl(i)-V_settl(j));   
%                 shear_turbulent = 0;
%                 turbulent_inertial = 0;
%                 kinetic_term = 2^(3/2)/4 * sqrt(pi) * V_plume * (dlist(i)+dlist(j))^2;    
%                 
%             else
%                 % Abrahamson theory
%                 brownian_term = cost * (dlist(i) + dlist(j))^2 / (dlist(i)*dlist(j));         
%                 sedimentation_term = pi/4*(dlist(i)+dlist(j))^2*abs(V_settl(i)-V_settl(j));   
%                 shear_turbulent = 0;
%                 turbulent_inertial = 0;
%                 kinetic_term = 2^(3/2)/4 * sqrt(pi) * 1.44 * V_plume * (dlist(i)+dlist(j))^2;   % 1.44 is sqrt(2)
%                 
%             end
            
            kern_out(i,j) =  (brownian_term +  turbulent_inertial + shear_turbulent + sedimentation_term);
           
        end
    end
    kern_out = kern_out+triu(kern_out,1)';
    
end


%             brownian_term = cost * (dlist(i) + dlist(j))^2 / (dlist(i)*dlist(j));                  
%             sedimentation_term = pi/4*(dlist(i)+dlist(j))^2*abs(V_settl(i)-V_settl(j));   
%             shear_laminar = gamma_shear/6 * cost_shape^3 * (dlist(i)+dlist(j))^3;       
%             shear_turbulent = gamma_shear_turb/6 * (dlist(i)+dlist(j))^3;
%             
%             v_ij(i,j) = (brownian_term+sedimentation_term+shear_turbulent)*4./(pi*(dlist(i)+dlist(j))^2);
%             
%             turbulent_inertial = pi*epsilon^(3/4) / (g*mu_kyn_air^(1/4)) * (dlist(i)/2+dlist(j)/2)^2 * abs(v_ij(i,j));         
%             
%             %shear_total = max(shear_laminar,shear_turbulent);               
%             shear_total = shear_turbulent;               
%             
%             kern_out(i,j) =  (brownian_term + sedimentation_term + turbulent_inertial + shear_total);
