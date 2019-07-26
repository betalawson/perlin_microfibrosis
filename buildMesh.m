function mesh = buildMesh(Nx, Ny, pixel_width)
% This function simply builds a basic struct that contains the mesh
% information, to be used by the other functions

% Create a set of points that fall in centres of pixels
xv = linspace(pixel_width/2, pixel_width*(Nx-1/2), Nx);
yv = linspace(pixel_width/2, pixel_width*(Ny-1/2), Ny);
[Y,X] = meshgrid(yv,xv);
points = [X(:) Y(:)];

% Store in mesh struct
mesh.points = points;
mesh.Nx = Nx;
mesh.Ny = Ny;

end

