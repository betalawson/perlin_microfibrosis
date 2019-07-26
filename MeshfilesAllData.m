function MeshfilesAllData

% This function creates meshfiles for Chaste from all of the data
% visualised in the paper, including the with and without fibre cases for
% compact and diffuse fibrosis

% First, list the pattern filenames
pattern_names = {'interstitial', 'compact', 'diffuse', 'patchy'};
pattern_names_reduced = {'int', 'cmp', 'dif', 'pat'};
pattern_names_nofibres = {'compact', 'diffuse'};
pattern_names_nofibres_reduced = {'cmp', 'dif'};

% Specify how many patterns to use (set to same as figures)
n_select = 6;

%%% PATTERN MATCHING - FIBRES TURNED OFF WHERE NEEDED
for j = 1:length(pattern_names)
    
    % Load particle data
    load(['Results/',pattern_names{j},'2000_full'], 'particles');
    
    % Reduce to the set of only unique particles by using unique command
    [~, I] = unique(particles.thetas, 'rows');
    particles.vals = particles.vals(I);
    particles.Ds = particles.Ds(I);
    
    % Sort the unique particles by their discrepancies
    [particles.Ds, I] = sort(particles.Ds);
    particles.vals = particles.vals(I);
    particles.thetas = particles.thetas(I,:);
    
    % Select particles corresponding to quantiles (to give even spread of
    % discrepancy values)
    quantiles = round(linspace(1,length(particles.Ds),n_show));
   
    % Loop over patterns to make meshfiles for
    for k = 1:n_select
        
        % Embed pattern in a larger mesh
        fibmesh = embedPattern(particles.vals{quantiles(k)});
        % Read out fibre orientation
        theta = particles.thetas(quantiles(k),8);
        % Creat the Chaste meshfiles
        createChasteMeshfiles([pattern_names_reduced{k},'Gen',num2str(k)], theta, fibmesh);
        
    end
    
end


%%% PATTERN MATCHING - WITH FIBRES
for j = 1:length(pattern_names_nofibres)
    
    % Load particle data
    load(['Results/',pattern_names_nofibres{j},'2000_full_nofibres'], 'particles');
    
    % Reduce to the set of only unique particles by using unique command
    [~, I] = unique(particles.thetas, 'rows');
    particles.vals = particles.vals(I);
    particles.Ds = particles.Ds(I);
    
    % Sort the unique particles by their discrepancies
    [particles.Ds, I] = sort(particles.Ds);
    particles.vals = particles.vals(I);
    particles.thetas = particles.thetas(I,:);
    
    % Select particles corresponding to quantiles (to give even spread of
    % discrepancy values)
    quantiles = round(linspace(1,length(particles.Ds),n_show));
   
    % Loop over patterns to make meshfiles for
    for k = 1:n_select
        
        % Embed pattern in a larger mesh
        fibmesh = embedPattern(particles.vals{quantiles(k)});
        % Read out fibre orientation - in variable location 6 now because fibre variables are gone
        theta = particles.thetas(quantiles(k),6);
        % Creat the Chaste meshfiles
        createChasteMeshfiles([pattern_names_nofibres_reduced{k},'NFGen',num2str(k)], theta, fibmesh);
        
    end
    
end