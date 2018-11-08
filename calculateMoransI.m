function morans_I = calculateMoransI(pattern)
% This function calculates the basic version of Moran's I for a supplied
% pattern. Moran's I is a measure of spatial autocorrelation, here
% calculated using only nearest neighbours (von Neumann neighbourhood)

% First, read out the dimensions of the pattern
[Ny, Nx] = size(pattern);

datalist = pattern(:);

% Calculate pattern density (this is the 'mean' of a 0/1 pattern)
density = sum(pattern(:)) / numel(pattern);

% Store the transpose of the pattern just because MATLAB doesn't otherwise
% allow the notation we want
patternT = pattern';

% Now, set up the mesh structure to allow for one big vectorised
% calculation

% Create a big matrix specifying which sites are real connections
weight_matrix = [  [ ones(2,4); zeros(2,4) ], [ ones(3,2*(Nx-2)+2*(Ny-2)); zeros(1,2*(Nx-2)+2*(Ny-2)) ], ones(4,(Nx-2)*(Ny-2)) ];
                          % corners                                  % sides                                   % interior            

% Create a matrix of all the interior points. This just makes the following
% step easier
interior_matrix = (2:Nx-1)' + (Nx:Nx:(Ny-2)*Nx);
                       
                       
% Create a matrix that lists all neighbours of each point
neigh_matrix = [  [ 2; Nx+1; 1; 1 ], [ Ny*Nx-1; (Ny-1)*Nx; 1; 1 ], [ Nx-1; 2*Nx; 1; 1 ], [ (Ny-1)*Nx+2; (Ny-2)*Nx+1; 1; 1 ], ...         % Corners
                  [ 1:Nx-2; 3:Nx; Nx+2:2*Nx-1; ones(1,Nx-2) ],  ...                                                                      % Bottom
                  [ (Ny-1)*Nx + (1:Nx-2); (Ny-1)*Nx + (3:Nx); (Ny-2)*Nx + (2:Nx-1); ones(1,Nx-2) ], ...                                  % Top
                  [ 1:Nx:(Ny-3)*Nx+1; 2*Nx+1:Nx:(Ny-1)*Nx+1; Nx+2:Nx:(Ny-2)*Nx+2; ones(1,Ny-2) ], ...                                    % Left
                  [ Nx:Nx:(Ny-2)*Nx; 3*Nx:Nx:Ny*Nx; 2*Nx-1:Nx:(Ny-1)*Nx-1; ones(1,Ny-2) ], ...                                           % Right
                  [ interior_matrix(:)' + 1; interior_matrix(:)' - 1; interior_matrix(:)' + Nx; interior_matrix(:)' - Nx] ];             % Interior
              
% Rearrange the matrices so that the columns match to the sites in order
order = [1, Ny*Nx, Nx, (Ny-1)*Nx+1, 2:Nx-1, (Ny-1)*Nx+(2:Nx-1), Nx+1:Nx:(Ny-2)*Nx+1, 2*Nx:Nx:(Ny-1)*Nx, interior_matrix(:)'];
weight_matrix(:,order) = weight_matrix;
neigh_matrix(:,order) = neigh_matrix;                
         
% Calculate Moran's I in a vectorised fashion using the created matrices
morans_I = (Ny*Nx) / sum(weight_matrix(:)) * sum( (patternT(:) - density)' .* sum(weight_matrix .* (patternT(neigh_matrix) - density), 1) ) / sum( (patternT(:) - density).^2 );

end

