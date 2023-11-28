function S = make_smooth_mat(L)
% MAKE_SMOOTH_MAT creates the sparse matrix S used to smooth with 
% a length 5 moving average defined as (1/13)*[1 3 5 3 1].  
% Smooths at ends following Ayesha's implementation in smoothst.m file.  
%
% S = SMOOTH_MAT_QUARTER(L) returns the sparse matrix S used 
%   to smooth. Size of S is L x L. Smooths at ends. 
%
% L: Length of column vector to smooth 

% Tiffany Anderson
% Created: 3/27/13



if L==1 % One paremeter array (scalar)
    
    S = sparse(1);

elseif L==2 % Two paremeter array

    S = sparse([5/8 3/8;...
                3/8 5/8]);
    
elseif L==3 % Three parameter array
    
    S = sparse([5/9 3/9 1/9;...
                3/11 5/11 3/11;...
                1/9 3/9 5/9]);
            
elseif L==4 % Four parameter array
    
    S = sparse([5/9 3/9 1/9 0;...
                3/12 5/12 3/12 1/12;...
                1/12 3/12 5/12 3/12;...
                0 1/9 3/9 5/9]);

else % Five or more parameter array
    
    row1 = [5/9 3/9 1/9];
    row2 = [3/12 5/12 3/12 1/12];
    rowL2 = fliplr(row2);
    rowL = fliplr(row1);
    
    smooth_vec = [1/13 3/13 5/13 3/13 1/13]; 

    % define arrays to hold sparse matrix indices and value
    i = [ones(3,1);...                            % first row
         2*ones(4,1);...                          % second row
         reshape(repmat(3:L-2,5,1),(L-4)*5,1);... % third row to row L-2 
         (L-1)*ones(4,1);...                      % row L-1 (second to last row)
         L*ones(3,1)];                            % row L (last row)
    
    j = [(1:3)';...                            % first row
         (1:4)';...                            % second row
         reshape([(1:L-4);(1:L-4)+1;(1:L-4)+2;(1:L-4)+3;(1:L-4)+4],...
                 (L-4)*5,1);...             % third row to row L-2 
         (L-3:L)';...                          % row L-1 (second to last row)
         (L-2:L)'];                            % row L (last row)
     
    s = [row1';row2';repmat(smooth_vec',L-4,1);rowL2';rowL'];
    
    S = sparse(i,j,s);
    
end















