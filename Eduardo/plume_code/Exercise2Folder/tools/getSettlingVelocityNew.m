function [v] = getSettlingVelocityNew(rho_p, rho_a, D, g, viscAir)
% Remember that viscAir = inp_var.viscAir
         
    ReCrit1 = (2.4)^2;
    ReCrit2 = (10/0.43)^2;
    dCrit1 = (ReCrit1.*18.*viscAir.^2./(g.*rho_a.*(rho_p-rho_a))).^(1/3);
    dCrit2 = (ReCrit2.*viscAir./(rho_a.*((4*g.^2.*(rho_p-rho_a).^2)./(225.*rho_a.*viscAir)).^(1/3))).^(1/2);
    
    if (D < dCrit1)
        v = ((rho_p-rho_a)*(g*D^2))/(18*viscAir);
    elseif (D > dCrit1) && (D < dCrit2)
        v = D*((4*g^2*(rho_p-rho_a)^2)/(225*rho_a*viscAir))^(1/3);
    else
        v = ((3.1*g*D)*(rho_p-rho_a)/rho_a)^(1/2);
    end
%     
%     
%     % Bagheri's version for drag coefficient
%         
%     for i = 1:length(rho)
%         rho_particle = pDense;
%         rho_air = rho(i);
%         d = pSize;
%         
%         % Formula di Arastoopour
%         Cd = 1; atol = 1e-20; rtol = 1e-7;
%         V_settl = sqrt(4*g*(rho_particle-rho_air)*d/(3*Cd*rho_air));
%         V_old = V_settl;
%         semaforo = 0;    iterazione = 0;    iterazione_MAX = 10000;
%         
%         while (semaforo < 1)
%             iterazione = iterazione+1;
%             Reynolds = (V_settl*d/mu_kyn_air);
%             if (Reynolds <= 1000)
%                 Cd = 24/Reynolds*(1+0.15*Reynolds^(0.687));
%             else
%                 Cd = 0.44;
%             end
%             V_settl = sqrt(4*g*(rho_particle-rho_air)*d/(3*Cd*rho_air));
%     
%             if(abs(V_settl-V_old)<=1e-6)
%                 semaforo = 1;
%             elseif(iterazione > iterazione_MAX)
%                 semaforo = 1;
%             else
%                 V_old = V_settl;
%             end    
%         end
%     end % Fine for


end