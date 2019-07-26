function createChasteMeshfiles(filename, theta, fibmap)
% This function creates a set of meshfiles for use in Chaste, from a
% provided 2D fibrosis map. Chaste uses triangular elements, so each
% non-obstructed site is divided into two triangular elements. Obstructed
% sites are non-conducting and have no associated elements - boundary
% conditions are constructed to reflect this. Nodepoints that are
% completely blocked in by fibrosis are also removed. The stimulus site is
% positioned in a way that corresponds with the setup in embedfibmap.m
%
% THIS FILE ASSUMES THAT THERE IS NO FIBROSIS ON DOMAIN EDGES
%
%   INPUTS
% -----------
%
% filename:        the filename for the meshfiles created
% theta:           the orientation of fibres (used in anisotropic conduction - .ortho file)
% fibmap:          a binary matrix with 1's specifying locations occupied by fibrosis


% Define the physical length of each pixel
pixel_width = 1/136;        % (mm)

% Define the radius of the stimulus site
stim_rad = 0.2;             % (mm)

% Define the position of the stimulus site (as proportions of total domain)
stim_X = 0.1;
stim_Y = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Read out the size of the fibrosis map
[Ny, Nx] = size(fibmap);

% Calculate the associated lengths in physical distances - they are also
% converted to centimetres, as used by Chaste
Lx = Nx * pixel_width * 0.1;
Ly = Ny * pixel_width * 0.1;

% Convert the stimulus position into real location
stim_X = Lx * stim_X;
stim_Y = Ly * stim_Y;

% Also convert the stimulus radius to cm
stim_rad = stim_rad * 0.1;

% The fibrosis map is flipped, because fibrosis maps are images (with (0,0)
% in top left), which then get put into  traditional (x,y) space for space 
% (with (0,0) in bottom left)
fibmap = flipud(fibmap);

% Create a list of blocked sites. Because the (:) operator works
% columnwise, a transpose is used. 
blocklist = fibmap'; blocklist = logical(blocklist(:));

% Now populate the problem domain with points
x = linspace(0, Lx, Nx+1);
y = linspace(0, Ly, Ny+1);
[X,Y] = meshgrid(x,y);
XT = X'; YT = Y';   % Transpose so linear indexing will work
points = [XT(:), YT(:)];

% The total number of points and elements is stored for convenience
Np = (Nx+1)*(Ny+1);

%%% NODES

% Create four vectors that together allow us to search for locations
% which are composed of 2x2 blocks of fibrosis (the centre node in
% these needs to be removed as it is inaccessible)
b1 = fibmap(1:Ny-1, 1:Nx-1); b1 = b1(:);
b2 = fibmap(1:Ny-1, 2:Nx); b2 = b2(:);
b3 = fibmap(2:Ny, 1:Nx-1); b3 = b3(:);
b4 = fibmap(2:Ny, 2:Nx); b4 = b4(:);

% When all four blocks around a node are fibrotic, that node is marked
% for removal
remove_squares = find(all([b1, b2, b3, b4], 2));

% The location is specified as the lower left square, and this must be
% converted then to the node number
[remove_js, remove_is] = ind2sub([Ny-1 Nx-1], remove_squares);
remove_nodes = (remove_js - 1) * (Nx+1) + remove_is + Nx+1 + 1;

% Removed nodes are sorted so the reordering stuff below works
remove_nodes = sort(remove_nodes, 'ascend');

% Now remove these nodes
points(remove_nodes,:) = [];

% Also insert these dead nodes as NaNs into a node reference list
% (pushing other elements to the right and thus creating a new node
% numbering
node_IDs = (1:Np)';

% Loop over NaNs to be inserted
for k = 1:length(remove_nodes)
    remove_node = remove_nodes(k);
    % Insert a NaN into the vector
    node_IDs = [node_IDs(1:remove_node - 1); NaN; node_IDs(remove_node:end)];
end

% Trim the vector, cutting off unneeded numbers
node_IDs = node_IDs(1:Np);

% Also, update number of nodes
Np = Np - length(remove_nodes);

% Finally, mark which nodes are stimulus locations or not, using the
% defined location of stimulus site and its radius
snode = sqrt((points(:,1) - stim_X).^2 + (points(:,2) - stim_Y).^2) < stim_rad;

%%% ELEMENTS

% First, the edge structure is set up. This is done first because after
% nodes are removed, the simple mappings to "node to east", "node to north"
% and "node to northeast" are broken.

% The bottom left corners for each element are used as a starting point.
% First, we grab out all bottom left corners, which are all sites except
% the top row and rightmost column of nodes

base_locs = (1:Nx)' + (0:Ny-1) * (Nx+1);
base_locs = base_locs(:);
base_locs = base_locs(~blocklist);

% Now, the node numbering information in node_IDs is used to give the
% correct numbers to these
c00 = node_IDs(base_locs);
c10 = node_IDs(base_locs + 1);    % Add 1 to move right
c01 = node_IDs(base_locs + Nx+1); % Add Nx+1 to move up
c11 = node_IDs(base_locs + Nx+2); % Add Nx+1 and 1 to move up and right

% Also, the number of elements, Ne is stored here for convenience. There
% are two triangles for each square mesh location
Ne = 2 * length(c00);

% Create a matrix of the triangular elements. For now, place all upper
% triangles together, then all lower triangles together
elems = [ [(1:2:Ne)', c00, c01, c11],  [(2:2:Ne)', c00, c10, c11] ];
% Now, we want to read across rows, but linear indexing works with columns,
% so transpose
elems = elems';
% Read down the columns, which will now alternate between the first and
% second sets of elements
elems = elems(:);
% Now convert this back into the appropriate shape
elems = reshape(elems, Ne, 4)';


%%% EDGES

% Only boundary edges are specified. First, initialise the edge nodes
% matrix, which will store beginning and end nodes for all boundary edges
edge_nodes = [];

% We also must designate which nodes are boundary nodes. First we
% initialise a vector of all 0's for each node, that will then be updated
% as we determine boundary edges
bnode = zeros(Np,1);

% Now work through the four edges of the domain. For now this assumes all
% boundary sites are non-fibrotic

% Bottom edge
start_nodes = node_IDs(1:Nx);
end_nodes = node_IDs(2:Nx+1);
edge_nodes = [edge_nodes; [start_nodes, end_nodes]];
bnode(start_nodes) = 1;
bnode(end_nodes) = 1;

% Top edge
start_nodes = node_IDs(Ny*(Nx+1)+(1:Nx));
end_nodes = node_IDs(Ny*(Nx+1)+(2:Nx+1));
edge_nodes = [edge_nodes; [start_nodes, end_nodes]];
bnode(start_nodes) = 1;
bnode(end_nodes) = 1;

% Left edge
start_nodes = node_IDs((0:Ny-1)*(Nx+1)+1);
end_nodes = node_IDs((1:Ny)*(Nx+1)+1);
edge_nodes = [edge_nodes; [start_nodes, end_nodes]];
bnode(start_nodes) = 1;
bnode(end_nodes) = 1;

% Right edge
start_nodes = node_IDs((0:Ny-1)*(Nx+1)+Nx+1);
end_nodes = node_IDs((1:Ny)*(Nx+1)+Nx+1);
edge_nodes = [edge_nodes; [start_nodes, end_nodes]];
bnode(start_nodes) = 1;
bnode(end_nodes) = 1;

% Now, all edges of a non-fibrotic region that border with a fibrotic
% region create additional boundary edges. The approach used is similar to
% the check for 2x2 fibrotic squares above.

% Will check to right, so grab out fibrotic status of all sites that have a
% site to the right (no rightmost edge), and then those fibrotic sites to
% the right
b1 = fibmap(1:Ny, 1:Nx-1); b1 = b1(:);
b2 = fibmap(1:Ny, 2:Nx); b2 = b2(:);

% Now any boundary between these two sites where one is fibrotic and one is
% non-fibrotic will create a boundary edge
lr_bound = find(all([~b1, b2], 2) | all([b1, ~b2], 2));

% Associated edge nodes with right boundary are the right and right-up node
% First, find (i,j) locations of these sites
[js, is] = ind2sub([Ny Nx-1], lr_bound);

% Convert the is and js into node locations - right and right-up. node_IDs
% vector is used to apply the renumbering
start_nodes = node_IDs((js-1) * (Nx+1) + is + 1);
end_nodes = node_IDs((js-1) * (Nx+1) + is + Nx+1 + 1);

% Add these to the list of edge nodes, and make sure they are marked as
% boundary nodes
edge_nodes = [edge_nodes; [start_nodes, end_nodes]];
bnode(start_nodes) = 1;
bnode(end_nodes) = 1;


% Now repeat this process for checking up

% Grab out all valid sites
b1 = fibmap(1:Ny-1, 1:Nx); b1 = b1(:);
b2 = fibmap(2:Ny, 1:Nx); b2 = b2(:);
% Find those with a mismatch (creating a boundary)
ud_bound = find(all([~b1, b2], 2) | all([b1, ~b2], 2));
% Read out (i,j) locations
[js, is] = ind2sub([Ny-1 Nx], ud_bound);
% Convert these into node locations (up, up-right)
start_nodes = node_IDs((js-1) * (Nx+1) + is + Nx+1);
end_nodes = node_IDs((js-1) * (Nx+1) + is + Nx+1 + 1);
% Add these to the list of edge nodes, and make sure they are marked as
% boundary nodes
edge_nodes = [edge_nodes; [start_nodes, end_nodes]];
bnode(start_nodes) = 1;
bnode(end_nodes) = 1;

% Read out number of boundary edges
Nb = size(edge_nodes,1);

%%% ORTHO

% Determine the vector components for fibre and transverse directions
f_x = cos(theta);
f_y = sin(theta);
t_x = -sin(theta);
t_y = cos(theta);

% Make these into a vector
ortho_vec = [f_x; f_y; t_x; t_y];


% With everything now dealt with, perform the file writing

% Initialise files to write in
nodefile = fopen([filename, '.node'],'w');
elefile = fopen([filename, '.ele'],'w');
edgefile = fopen([filename, '.edge'],'w');
orthofile = fopen([filename, '.ortho'],'w');

% Node file begins with:
%    #nodes    #dims    #extra_vars    #boundary_markers
fprintf(nodefile, '%d\t%d\t%d\t%d\n', Np, 2, 2, 1);
% Then, node data is:
%    ID     x-coord     y-coord    var1(stimulus?)    var2(boundary?)
fprintf(nodefile, '%d\t%g\t%g\t1\t%d\t%d\n', [(1:Np)', points(:,1), points(:,2), snode, bnode]');

% Elements file begins with:
%    #elements   #nodes-per-element   #extra_vars (here zero)
fprintf(elefile, '%d\t3\t0\n', Ne);
% Element data is:
%    ID    node1     node2    node3
fprintf(elefile, '%d\t%d\t%d\t%d\n', elems');

% Edges file begins with:
%    #edges     #boundary_specifiers
fprintf(edgefile, '%d\t1\n', Nb);
% Edge data is:
%    ID    node1      node2    boundary?
fprintf(edgefile, '%d\t%d\t%d\t1\n', [(1:Nb)', edge_nodes]' );

% Ortho file begins with:
%     #elements
fprintf(orthofile, '%d\n', Ne);
% Ortho data is:
%     vec1x    vec1y    vec2x    vec2y
fprintf(orthofile, '%f\t%f\t%f\t%f\n', repmat(ortho_vec, 1, Ne) );


% Before finishing, write to each file a comment, just to prevent any empty
% endlines down the bottom from possibly screwing things up
fprintf(orthofile, '# Generated by createChasteMeshfiles.m');
fprintf(nodefile, '# Generated by createChasteMeshfiles.m');
fprintf(elefile, '# Generated by createChasteMeshfiles.m');
fprintf(edgefile, '# Generated by createChasteMeshfiles.m');

% Finally, close all files
fclose all;