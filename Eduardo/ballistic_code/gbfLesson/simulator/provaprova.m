function provaprova
    
close all;

    t = linspace(0,50,40000);

    f1 = sin( (pi/12) .* t);
    f2 = sin( (pi/14) .*t + 3);
    f3 = sin( (pi/16) .*t);

    
    figure(1)
    plot(t,f1,'-');  
    hold on; 
    plot(t,f3,'-g');
    
    figure(2)
    plot(t,f1+f2+f3,'-k');
   
    
end