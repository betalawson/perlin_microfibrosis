function [t, V, APD, AT, nodemap] = load_results(path, filename)
% This funtion handles the actual reading of the h5 files and the
% re-ordering of nodes. Because the re-ordering depends on information that
% is contained in the base results file, it is also used here for
% interpretation of the data from the snapshots file.
%
% If a filename is provided, this function assumes a standard structure in
% which the 

% Filenames assumed if not provided, otherwise initialised if given
if nargin == 1
    main_path = fullfile(path, 'results.h5');
    snapshots_path = fullfile(path, 'snapshots.h5');
    perm_path = fullfile(path, 'permutation.txt');
else
    main_path = fullfile(path, [filename,'.h5']);
    snapshots_path = fullfile(path, [filename,'_snapshots.h5']);
    perm_path = fullfile(path, [filename,'.txt']);
end

% Read out the time and voltage data
t = h5read(main_path, '/Data_Unlimited');
V = squeeze(h5read(main_path, ['/Data']));

% Read out the APD and activation time data
APD = squeeze(h5read(snapshots_path, '/APD'));
AT = squeeze(h5read(snapshots_path, '/Activation'));

% Initialise default nodemap
nodemap = (1:size(V, 1))-1;

% Update nodemap if provided in h5 file
if ~h5readatt(main_path, '/Data', 'IsDataComplete')
    nodemap = h5readatt(main_path, '/Data', 'NodeMap');
end

% If a permutation .txt file has been provided, use it to re-order the data
if exist(perm_path, 'file')
    
    % Open file and read out the permutation data
    file_obj=fopen(perm_path);
    perm=cell2mat(textscan(file_obj,'','headerlines',1,'delimiter',' ','collectoutput',1));
    fclose(file_obj);
    
    % Apply re-ordering, according to whether nodemap was found or not
    if ~isempty(nodemap)
        iperm = zeros(size(perm, 1), 1);
        iperm(perm(:, 2)+1) = 1:length(iperm);
        nodemap = iperm(nodemap+1)-1;
        [nodemap, k] = sort(nodemap);
        V = V(k, :);
        APD = APD(k,:);
        AT = AT(k,:);
    else
        V = V(perm(:, 2)+1, :);
        APD = APD(perm(:, 2)+1, :);
        AT = AT(perm(:, 2)+1, :);
    end
    
end