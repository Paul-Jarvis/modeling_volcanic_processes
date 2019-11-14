function stick_eff = sticking_eff(dlist, v_particles, mu_l, rho_particle, m_k, choice, flag, v_ij)

    if (choice == 0)
        % Ice sticking
        stick_eff = stickingIceWet(dlist, v_particles, mu_l, rho_particle, m_k, flag);
%         disp('Ice')
    else
        % Wet sticking
        stick_eff = stickingLiquidWet(dlist, v_particles, mu_l, rho_particle, v_ij);
%         disp('Humid')
    end
     
end

% _________________________________________________________________________

function stick_eff = stickingLiquidWet(dlist, v_particles, mu_l, rho_particle, v_ij)

    St_cr = 1.3;
    q = 0.8;
    viscl = 5.4e-4;
    stick_eff = zeros(length(dlist),length(dlist));
    
    for i=1:length(dlist)
        v_i(i) = v_particles(i);
        for j=i:length(dlist)
            v_j = v_particles(j);
%             St =  8.0*rho_particle*abs(v_i(i)-v_j)*dlist(i)*dlist(j)/(9.0*viscl*(dlist(i)+dlist(j)));
            St =  8.0*rho_particle*abs(v_ij(i,j))*dlist(i)*dlist(j)/(9.0*viscl*(dlist(i)+dlist(j)));
            stick_eff(i,j) = 1/(1+(St/St_cr)^q);
        end
    end
    
    stick_eff = stick_eff+triu(stick_eff,1)';
    stick_eff = avg_values_diagonal(stick_eff);

end


function stick_eff = stickingIceWet(dlist, v_particles, mu_l, rho_particle, m_k, flag)

% "Flag" is related to the position beyond which the sticking must decrease
     
    a=1;
    b=flag;
    c=2;                                                                   % The higher the value, the smoother the decrease
       
    for i=1:length(dlist)
        for j=i:length(dlist)
            stick_eff(i,j) = 0.09*(a-1./(1+exp(-(j-b)./c)));
        end
    end
%     b = 0.47;
%     B = 6.1e-4;    
%     for i=1:length(dlist)
%         v_i(i) = v_particles(i);
%         for j=i:length(dlist)
%             v_j = v_particles(j);
% %             stick_eff(i,j) = B * (m_k(i)^b + m_k(j)^b) / ( pi/4 * (dlist(i)+dlist(j))^2 * abs(v_i(i)-v_j));
%             stick_eff(i,j) = ;
%         end
%     end
    
    stick_eff = stick_eff+triu(stick_eff,1)';
    stick_eff = avg_values_diagonal(stick_eff);
    
end


function matrix = avg_values_diagonal(matrix)
    
    [rows, columns] = size(matrix);
    
    for i=2:rows-1
        matrix(i,i) = mean([matrix(i,i-1),matrix(i,i+1)]); 
    end
    matrix(1,1) = matrix(1,2);
    matrix(end,end) = matrix(end,end-1);
end