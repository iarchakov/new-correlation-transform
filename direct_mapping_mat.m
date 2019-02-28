% Programming implementation of the new method of unconstrained 
% transformation for correlation matrices suggested 
% in Archakov and Hansen (2018) 
%
% Direct mapping from a non-singular correlation matrix "C"
% to a unique corresponding real vector "gamma"
% ------------------------------------------------------------------------



function gamma = direct_mapping_mat(C)
    
    gamma = [];
    
    % Check if input vector is of suitable dimensions and 
    % tolerance value belongs to a proper interval
    if (size(C,1) == size(C,2)) && isequal(diag(C),ones(size(C,1),1))...
            && all(eig(C) > 0) && isequal(C,C')
   
       
        % Apply matrix log-transformation to C and get 
        % off-diagonal elements
        A = logm(C);
        gamma = A(logical(tril(ones(size(C,1)),-1)));
        
    else
        fprintf('Error : input is of wrong format');
    end