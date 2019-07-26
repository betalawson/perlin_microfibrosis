function generateSeedData(N_seeds, N_freqs, seed)
% This function creates seed data and saves it as fibro_seedinfo.mat
% This function only needs to be called when seed data needs to be changed
% or (re)generated. Usage:
%
% generateSeedData( N_seeds, (N_offsets), (seed) );
%
% N_seeds:   number of unique seeds to generate
% N_freqs:   optional argument specifying the maximum number of octave
%            layers that will ever be required (default: 8)
% seed:      an integer that, if provided, will be used for seeding
%            MATLAB's random number generator


% If the user provided a seed, use that seed to generate the 
if nargin > 2
    rng(seed);
else
    rng('shuffle');
end

% If the user didn't specify the number of offsets to use, assume a decent
% safe number like eight
if nargin == 1
    N_freqs = 8;
end

% Create the requested number of permutation and offset tables
for k = 1:N_seeds
    
    % Permutation tables for this seed
    for j = 1:N_freqs
        permute_tables{k}(j,:) = int32(randperm(256) - 1);
    end
    % Offset table for this seed
    offset_tables{k} = rand(N_freqs, 2) - 0.5;
    
end

% Save results
save('fibro_seedinfo.mat', 'permute_tables', 'offset_tables');
    