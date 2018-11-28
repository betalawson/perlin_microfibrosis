function fibro_maps = FibrosisImageLoader(reduce, cut_compact)
% This function loads the source image of microfibrosis histology and
% converts the four different patterns into presence/absence maps via
% thresholding. The user specifies a single input, a boolean specifying
% whether or not to use 'reduced' patterns (representative sections of the
% full patterns). Maps are stored in a cell vector, with ordering
%        {interstitial, compact, diffuse, patchy}
%
% If a second argument is supplied, and set to true (1), then additional
% material will be cut from the compact pattern. Only used if reduce is
% also specified.

% First, load in the raw image
im = imread('micro_fibrosis.jpg');

% Define the threshold value in the green channel that determines fibrosis
% presence (decrease vaulue to increase fibrotic density - current value
% selected by experimentation)
green_threshold_val = 180;

% Define the size of the representative subsections (used if user specifies
% 'reduce' flag
subregion_width = 250;
subregion_height = 400;

% Define starting locations of the representative subsection for each image
% Each location is specified from the top-left corner
subregion_starts = [
    [100,120];
    [100,175];
    [50,125];
    [85,110]
    ];

% Create vectors of the amount of blue averaged along rows and columns
% (blue is best colour for picking out the starts and ends of each image
% in this case)
rows_B = mean(im(:,:,3),2);
cols_B = mean(im(:,:,3),1);

% Define cutoff values for where images start and end in terms of these
% average values across columns/rows (values selected by plotting the above
% vectors and inspecting the output)
blue_cutoff_rows = 140;
blue_cutoff_cols = 180;

% All images start and end at the same vertical location, so top/bottom
% image boundaries may be found outside of the below loop
vert_start = find( rows_B < blue_cutoff_rows, 1);
vert_finish = find( rows_B < blue_cutoff_rows, 1, 'last');

% Assume if not specified a flag for 'reduce', that no reduction is desired
if nargin < 1
    reduce = 0;
end

% Initialise the cell variable used to store the final patterns
fibro_maps = cell(1,4);

% Loop over the four images, thresholding each to get a presence/absence
% map for fibrosis
search_start = 1;
for k = 1:4
    
    % Find first point to the right of current location that goes below
    % the cutoff blue amount. Cutoff selected again by
    horz_start = search_start - 1 + find(cols_B(search_start:end) < blue_cutoff_cols, 1);
    % Find next point that goes above cutoff blue amount (subtract one so
    % this point is not included)
    horz_end = horz_start - 1 + find(cols_B(horz_start:end) >= blue_cutoff_cols, 1) - 1;
    
    % Now that image bounds are known
    image_G = squeeze(im(vert_start:vert_finish, horz_start:horz_end,2));
    
    % If using reduced images, cut out the representative sections (as
    % specified above) from the full images
    if reduce
        image_G = image_G( subregion_starts(k,2):subregion_starts(k,2)+subregion_height-1 , subregion_starts(k,1):subregion_starts(k,1)+subregion_width-1 );
    end
    
    % Convert this image into a presence/absence map for collagen by
    % thresholding the green channel as specified
    fibro_maps{k} = image_G < green_threshold_val;
    
    % Update beginning of search for next image
    search_start = horz_end+1;
    
end

% Remove the questionable parts of the compact pattern using a 'fill
% remove' command on each contiguous piece of fibrosis at the boundary
% section.
if reduce
    fibro_maps{2} = fill_remove(fibro_maps{2}, [124, 399]);
    fibro_maps{2} = fill_remove(fibro_maps{2}, [100, 390]);
    fibro_maps{2} = fill_remove(fibro_maps{2}, [100, 300]);
else
    fibro_maps{2} = fill_remove(fibro_maps{2}, [17, 538]);
    fibro_maps{2} = fill_remove(fibro_maps{2}, [55, 504]);
    fibro_maps{2} = fill_remove(fibro_maps{2}, [62, 496]);
    fibro_maps{2} = fill_remove(fibro_maps{2}, [200, 480]);
    fibro_maps{2} = fill_remove(fibro_maps{2}, [190, 575]);
    fibro_maps{2} = fill_remove(fibro_maps{2}, [222, 573]);
    fibro_maps{2} = fill_remove(fibro_maps{2}, [210, 590]);
    fibro_maps{2} = fill_remove(fibro_maps{2}, [217, 584]);
end

% Remove the top part of the image for compact, if using reduced patterns
% and specifically requested
if nargin > 1
    if reduce && cut_compact
        fibro_maps{2} = fill_remove(fibro_maps{2}, [152, 8]);
        fibro_maps{2} = fill_remove(fibro_maps{2}, [193, 18]);
        fibro_maps{2} = fill_remove(fibro_maps{2}, [110, 30]);
        fibro_maps{2} = fill_remove(fibro_maps{2}, [96, 12]);
        fibro_maps{2} = fill_remove(fibro_maps{2}, [142, 30]);
        fibro_maps{2} = fill_remove(fibro_maps{2}, [226, 30]);
        fibro_maps{2} = fill_remove(fibro_maps{2}, [248, 34]);
        fibro_maps{2} = fill_remove(fibro_maps{2}, [231, 18]);
    end
end


% Visualise the four patterns to confirm everything's working okay
% subplot(1,4,1);
% imagesc(fibro_maps{1});
% axis equal;
% subplot(1,4,2);
% imagesc(fibro_maps{2});
% axis equal;
% subplot(1,4,3);
% imagesc(fibro_maps{3});
% axis equal;
% subplot(1,4,4);
% imagesc(fibro_maps{4});
% axis equal;

end


function pattern = fill_remove(pattern, seed_loc)
% This function performs a fill remove operation that zeros all sites
% connected to the seed location in the given pattern. A Moore
% neighbourhood is used (diagonals count as connected)

pattern = remove_and_propagate(pattern, seed_loc(2), seed_loc(1) );

% Subfunction that removes the current site and spills into all filled
% neighbours. Currently creates 8 recursive calls and THEN checks how
% valid they are. Probably more efficient (but much uglier) to check
% before doing a recursive call
    function A = remove_and_propagate(A,i,j)
        
        % Ensure that processing only occurs within image boundaries, and
        % if current site is filled (if test fails before trying to check
        % A(i,j) when outside bounds, so this is safe
        if i >= 1 && j >= 1 && i <= size(A,1) && j <= size(A,2) && A(i,j) == 1
            
            % Remove this site
            A(i,j) = 0;
            
            % Propagate to all neighbours (Moore neighbourhood)
            A = remove_and_propagate(A,i-1,j);
            A = remove_and_propagate(A,i-1,j-1);
            A = remove_and_propagate(A,i,j-1);
            A = remove_and_propagate(A,i+1,j-1);
            A = remove_and_propagate(A,i+1,j);
            A = remove_and_propagate(A,i+1,j+1);
            A = remove_and_propagate(A,i,j+1);
            A = remove_and_propagate(A,i-1,j+1);
            
        end
        
    end


end