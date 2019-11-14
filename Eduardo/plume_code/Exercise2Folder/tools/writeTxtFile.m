function writeTxtFile(inp_var, plume_out)

    sizes = plume_out.aggregation.sizes; 
    mf = plume_out.mass_fractions;
    [rows, columns] = size(mf);
        
    
    file_out = 'out';
    file_run_folder = inp_var.file_run_folder;                             % Local variables
    file_run_extention = inp_var.file_run_extention;
    
    output_path = strcat('../','runs','/',file_run_folder,'/',file_out,'.',file_run_extention);
       
    fileID = fopen(output_path,'w');
    for i=1:length(sizes)-1
        fprintf(fileID, '%3.3f ', sizes(i));
    end
    
    fprintf(fileID,'%3.3f\r\n', sizes(end));
    
    for i=1:length(plume_out.z)
        
        for j=1:columns-1
            fprintf(fileID, '%3.3f ', mf(i,j));     
        end
        fprintf(fileID,'%3.3f\r\n', mf(i,end));
        
    end
%     
%     fprintf(fileID,'%i\r\n', count); 
% 
    fclose(fileID);
end