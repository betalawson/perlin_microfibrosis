function fibmesh = embedPattern(pattern)
% This function takes a fibrotic pattern and embeds it into a larger region
% of non-fibrotic tissue, with a region-growing algorithm used to create a
% less rigidly rectangular shape in which the pattern is placed. The
% automatic sizing for the full mesh uses a setup where the pattern is
% stimulated from the side, with wider margins on left/right than
% above/below. This design was created assuming a rectangular pattern that
% is taller than it is wide

% Define the seed spacing (smaller number places more seeds, leads to more
% fragmented patterns)
base_seed_spacing = 25;

% Read out the dimensions of the pattern
[Py,Px] = size(pattern);

% Use these to automatically calculate the size of the full mesh
Nx = round( Px * (3/2 + 1 + 1) );      % 3/2 the vertical size on the stimulus side, same width as pattern on the far side
Ny = round( Py * (5/16 + 1 + 5/16) );  % 5/16 the vertical size of the pattern for the top and bottom margins

% Initialise the pattern
fibmesh = zeros(Ny, Nx);

% Define where the pattern placement 'box' starts and ends
box_start_x = round(3/2*Px);
box_start_y = round(5/16*Py);

% Use the region-growing algorithm to create a rougher shape than a pure
% rectangle in which to place the pattern

% Calculate appropriate seed spacing distances in the x and y directions
seed_spacing_x = base_seed_spacing * sqrt(Px/Py);
seed_spacing_y = base_seed_spacing * sqrt(Py/Px);

% Define the margins in which no seeds are placed
marg_x = 0.2 * Px * sqrt(Py/Px);
marg_y = 0.2 * Py * sqrt(Px/Py);

% From these, calculate the "true placement lengths"
place_x = Px - 2 * marg_x;
place_y = Py - 2 * marg_y;

% Calculate the number of seeds to place in the x and y directions
n_seeds_x = round( place_x / seed_spacing_x);
n_seeds_y = round( place_y / seed_spacing_y);

% Define the locations of the seeds
x_seeds = box_start_x + ( ( marg_x + (place_x/(2*n_seeds_x)) ):(place_x/n_seeds_x):(place_x + marg_x) );
y_seeds = box_start_y + ( ( marg_y + (place_y/(2*n_seeds_y)) ):(place_y/n_seeds_y):(place_y + marg_y) );

% Round these locations to whole numbers
x_seeds = round(x_seeds);
y_seeds = round(y_seeds);

% Use meshgrid to create a full list of seeding locations
[X, Y] = meshgrid(x_seeds, y_seeds);
seed_locs = [X(:) Y(:)];

% Now run the region-growing algorithm to create a mask in which the
% pattern is placed
mask = regionGrower(Nx, Ny, seed_locs, 0.85 * Px * Py, 0);

% Paste the pattern into the placement box
fibmesh(box_start_y:box_start_y+Py-1, box_start_x:box_start_x+Px-1) = pattern;

% Apply the mask to destroy the perfect rectangular shape
fibmesh = fibmesh .* mask;

end

function occ = regionGrower(Nx, Ny, seed_locs, n_sites, diagonals)
% This function creates regions by growing out randomly from each of the
% seed locations specified as an input, onto a mesh of total size Nx x Ny.
% No code is included for boundary handling, so n_sites should be specified
% such that a boundary collision does not occur.

% Initialise occupancy matrix
occ = zeros(Ny*Nx,1);

% Initialise routine by placing seed objects
placed = 0;
choosesites = [];
for k = 1:size(seed_locs,1)
   
    ce = ( seed_locs(k,2) - 1) * Nx + seed_locs(k,1);
    occ(ce) = 1;
    choosesites = [choosesites, ce];
    placed = placed + 1;
    
end

% Define the directions that can be selected at each point. Default to
% straight directions only, but allow diagonals if specified
dirs = [+1, -1, +Nx, -Nx];                              
if nargin >= 5
    if diagonals
        dirs = [+1, -1, +Nx, -Nx, +Nx+1, +Nx-1, -Nx+1, -Nx-1];
    end
end

% Repeatedly place items until requested number placed
while placed < n_sites
   
    % Pick a random occupied site
    r = ceil(rand * length(choosesites));
    
    % Check if it's completely surrounded
    if all(occ(choosesites(r) + dirs))
        
        % If surrounded, make sure this site isn't chosen again (massively
        % improves efficiency)
        choosesites(choosesites == choosesites(r)) = [];
        
    else
    
        % Otherwise, pick a random direction and place in that direction if
        % it's not occupied
        r2 = ceil(rand * length(dirs));
    
        % Place an object here if this site is not occupied
        trial_site = choosesites(r) + dirs(r2);
        if occ(trial_site) == 0
            choosesites = [choosesites, trial_site];
            occ(trial_site) = 1;
            placed = placed + 1;            
        end    
    end
end

% Reshape back into correct form
occ = reshape(occ', Nx, Ny)';

end