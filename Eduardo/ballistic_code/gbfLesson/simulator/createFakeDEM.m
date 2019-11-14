function createFakeDEM

    filename_dem = 'dem/dem_10m.txt';
    [X_original, Y_original, Z_original, X_final, Y_final, Z_final] = readDEM(filename_dem);


end

function [X_original, Y_original, Z_original, X_final, Y_final, Z_final] = readDEM(filename)
    
    x_low = 496474;  
    x_up = 496854;
    y_low = 4250466;
    y_up =  4250846;
      
    x_vent = [496634];
    y_vent = [4250706];
    
    z_plateau = 218;
    z_pleteau_new = 140;
    z_conduit = 130;
    radius = 30;

    fid_1 = fopen(filename);
    head_dem = textscan(fid_1, '%s', 12);
    fclose(fid_1);
    head_dem{1}
   
    xllcorner   = round(str2double(head_dem{1}{6}));
    yllcorner   = round(str2double(head_dem{1}{8}));
    cellsize    = str2double(head_dem{1}{10});
    NODATA      = str2double(head_dem{1}{12});
    
    Z           = dlmread(filename,' ',6,0);
    Z           = Z(:,1:end-1);
    Z(Z==NODATA)= nan;
    [X, Y]      = meshgrid(xllcorner:cellsize:xllcorner+(size(Z,2)-1)*cellsize,...
                           yllcorner:cellsize:yllcorner+(size(Z,1)-1)*cellsize);


    X_original = X;
    Y_original = Y;
    Z_original = flipud(Z);

       
    % Here we modify the topography in order to have the original crater
    dim = size(Z_original);
    n_rows = dim(1)
    n_columns = dim(2)

    X_final = X_original;
    Y_final = Y_original;
    Z_final = Z_original;
    
            
    for i_row = 1:n_rows    
        i_row
        for i_col = 1:n_columns     
                
            % If over the square
            
            if (X_original(i_row,i_col)>x_low && X_original(i_row,i_col)<x_up)
                
                if (Y_original(i_row,i_col)>y_low && Y_original(i_row,i_col)<y_up )
                    
                    if Z_original(i_row, i_col) < z_plateau
                        
                        Z_original(i_row, i_col) = z_pleteau_new;
                        
                        d = sqrt( (X_original(i_row,i_col)-x_vent)^2 + (Y_original(i_row,i_col)-y_vent)^2);
                        
%                         if d <=radius 
%                             Z_original(i_row, i_col) = z_conduit;
%                         end
                    
                    end
                    
                    
                end
                
            end
            
            
        end
    end
    
    figure(12121212)
    surf(X_original,Y_original,Z_original)
    
    Z_original = flipud(Z_original);
    
    fid = fopen('dem/VulcanoFake2mFlat.txt','wt');
    
    fprintf(fid, '%s  %s\n', head_dem{1}{1}, head_dem{1}{2});
    fprintf(fid, '%s  %s\n', head_dem{1}{3}, head_dem{1}{4});
    fprintf(fid, '%s  %s\n', head_dem{1}{5}, head_dem{1}{6});
    fprintf(fid, '%s  %s\n', head_dem{1}{7}, head_dem{1}{8});
    fprintf(fid, '%s  %s\n', head_dem{1}{9}, head_dem{1}{10});
    fprintf(fid, '%s  %s\n', head_dem{1}{11}, head_dem{1}{12});
    
    for i_row = 1:n_rows    
        i_row
        for i_col = 1:n_columns  
    
            % Write txt file
                       
            if isnan(Z_original(i_row, i_col))>0
%                 Z_original(i_row, i_col) = 0;
                Z_original(i_row, i_col) = -9999;
            end
            
            if (i_col<n_columns)
                fprintf(fid, '%i ', Z_original(i_row, i_col));
            else
                fprintf(fid, '%i\n', Z_original(i_row, i_col));
            end
            
        end
    end
      
     fclose(fid);
    
     
     
     
     shading interp;  
end