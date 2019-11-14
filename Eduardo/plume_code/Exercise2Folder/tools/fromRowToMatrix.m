function mat = fromRowToMatrix(vec, N_m, N_V)
    % Vec is a column vector (N_m*N_V,1)
    % First N_m values are (1:x,1)
    % Second N_m values are (1:x,2)

   mat = zeros(N_m,N_V);
   
   for i=1:N_m
%        mat(1:N_m,i) = vec( (i-1)*N_m+1 : (i-1)*N_m+N_m );
       mat(i,1:N_V) = vec( (i-1)*N_V+1 : (i-1)*N_V+N_V );
   end

end