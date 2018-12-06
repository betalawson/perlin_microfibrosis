function [part_thetas, part_outputs, part_summaries, part_Ds] = performABCSMC_Fibro(N_parts, f_simulate, f_summaries, f_discrepancy, target_data, params_mins, params_maxs, scale_param, desired_D, visualise, varargin)

% This function performs SMC-ABC as laid out by Drovandi and Pettitt (2011)
% Sets of particles satisfying a set of intermediary "discrepancy cutoffs"
% are generated by successive resampling and mutation steps in the vein of
% traditional SMC. The function is designed to be modular, so accepts as
% input the functions that will be used to simulate the model, and to
% calculate the discrepancy with the target data. Further options may be 
% specified by the user within the body of the program.
%
% The 'transform' jumping distribution method was provided by Christopher
% Drovandi, as laid out in South et al. (2018)
%
%
% Usage:   [part_thetas, part_vals] = SMC_ABC(N_parts, f_simulate, f_discrep, target_data, theta_mins, theta_maxs, desired_D, visualise, f_visualise)
%
% INPUTS
% -------
%
% N_parts - The number of particles to generate. Corresponds to the
%           expected number of outputs satisfying the discrepancy criteria
%           at the end of the method
%
% f_simulate - The function used to simulate the model. This should
%              generate only one output (but a struct could be used to 
%              store varied information/multiple outputs). The input is a
%              single vector of parameters
%
% f_summaries - The function used to calculate the summary statistics for
%               given data. Takes as input a single object, which is the 
%               output of f_simulate. The target data should also be
%               specified in this format
%
% f_discrepancy - The function used to calculate the discrepancy between
%                 two sets of metrics. As a basis this would be the
%                 Euclidean or Mahalanobis distance, but is supplied
%                 separately so that extra problem-specific modifications
%                 can be made (e.g. taking into account periodicity of
%                 angular measures)
%
% target_data - The data, of equivalent type to that generated by
%               f_simulate, that the SMC_ABC routine seeks to generate
%               samples sufficiently close to
%               
% params_mins - The minimum values allowed for the parameters of the model
%
% params_maxs - The maximum values allowed for the parameters of the model
%
% scale_param - Vector of logicals of length equal to number of parameters.
%               Specifies which parameters are 'scale' parameters that are
%               to be transformed logarithmically so as to have a uniform
%               prior.
%
% desired_D - The discrepancy that the algorithm will attempt to achieve
%             (typical D values depend on application, so this is
%             externally specified by the user)
%
% visualise - A boolean flag specifying if the user wants to visualise the
%             output during the course of the algorithm
%
% (f_visualise) - A function (supplied by user because likely
%                 application-specific) that visualises the results.
%                 Optional argument, only used if visualise flag is set.
%                 Syntax must be
%                    visFunc( params, outputs, summaries, discrepancies )
%                 
% (options) - A struct containing modifications to the SMC-ABC algorithm's
%             default options/parameters. Information regarding these is
%             contained within the file. User may set:
%
%             options.jumping_type ('normal', 'transform')
%             options.resample_weighting ('equal', 'weighted')
%             options.metric_weighting ('mahalanobis', 'variance', 'none')
%             options.keep_fraction ([0,1])
%             options.max_MCMC_steps (positive integer)
%             options.verbose (0 - false, 1 - true)
%
% For optional outputs, (options) must come last (but can be supplied after
% 'visualise' flag if it is set to zero)
%
% OUTPUTS
% --------
%
% part_thetas - The locations in parameter space of the final particles
%
% part_vals - The model outputs corresponding to each particle


% Default options (can be changed here or adjusted on the fly by supplying options argument)
jumping_type = 'normal';  % 'normal' - use multivariate normal jumping distribution (with 'optimal' scaling of covariance matrix, 
                             %            which is defined by currrent particle locations)
                             % 'transform' - use a series of probability transforms to attain normal marginals, then fit a three 
                             %               component Gaussian mixture to the transformed data, and generate proposals from this 
                             %               Gaussian mixture

resample_weighting = 'equal';  % 'equal' - All kept particles are weighted equally for resampling steps
                               % 'weighted' - Particles are weighted according to their discrepancy values when resampling (taking 
                               %              worst kept particle as three sigma in a normal distribution, and then normalising
                               %              weights so they sum to one)
                               
metric_weighting = 'variance';  % 'mahalanobis' - Mahalanobis distance is used for calculations of discrepancy
                                   % 'variance' - Metrics are scaled according to their variance, but covariances are ignored
                                   % 'none' - Euclidean distance is used 
                                   % (All options still include rescaling of angle metrics with respect to eccentricity)
                               
keep_fraction = 0.5;         % The proportion of particles to keep during each resample step (recommended 0.5 as a starting value)

max_MCMC_steps = 200;         % The maximum number of 'jiggle' steps applied to all particles in the attempt to find unique locations
                            
verbose = 1;                 % Output information about particle uniqueness and discrepancy targets
                             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT HANDLING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over optional input arguments, and set them accordingly
f_visualise_provided = false;
options_provided = false;
for k = 1: nargin - 10    % 10 set inputs, only loop over extras
    
    % Check if this is a function-type input (in which case it will be f_visualise)
    if isa(varargin{k},'function_handle')
        f_visualise = varargin{k};
        f_visualise_provided = true;
    % Check if this is a struct-type input (in which case it will be options)
    elseif isa(varargin{k}, 'struct')
        options = varargin{k};
        options_provided = true;
    % Otherwise, give warning, but continue
    else
        fprintf('Warning: Optional input provided to performABCSMC_Fibro but not in an appropriate format \n'); 
    end
    
end

% Cannot visualise if no visualisation function provided
if visualise && ~f_visualise_provided
    fprintf('Warning: No visualisation function provided. Continuing without visualisation... \n');
    visualise = 0;
end

% If options were provided, overwrite the default options with the custom
% options that were provided
if isfield(options,'jumping_type')
    jumping_type = options.jumping_type;
end
if isfield(options,'resample_weighting')
    resample_weighting = options.resample_weighting;
end
if isfield(options,'metric_weighting')
    metric_weighting = options.metric_weighting;
end
if isfield(options,'keep_fraction')
    keep_fraction = options.keep_fraction;
end
if isfield(options,'max_MCMC_steps')
    max_MCMC_steps = options.max_MCMC_steps;
end
if isfield(options,'verbose')
    verbose = options.verbose;
end
							 
							 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE STARTS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the summaries for the target data
target_summaries = f_summaries(target_data);

% Calculate the ordinal location of the last particle kept for convenience
worst_keep = round(N_parts * keep_fraction);

% Read out the number of variables
N_theta = length(params_mins);

% Convert scaling parameters to log equivalents. The set of transformed
% parameters is termed 'theta' in this code
theta_mins = params_mins;
theta_mins(scale_param) = log(theta_mins(scale_param));
theta_maxs = params_maxs;
theta_maxs(scale_param) = log(theta_maxs(scale_param));

% Initialise particles using Latin Hypercube Sampling
part_lhs = lhsdesign(N_parts, N_theta);
part_thetas = theta_mins + part_lhs .* (theta_maxs - theta_mins);

% For simulation, convert scaling parameters back to true values
part_params = part_thetas;
part_params(scale_param) = exp(part_params(scale_param));

% Load seed data (generated by a separate file)
load('fibro_seedinfo.mat','permute_tables', 'offset_tables');
% Count number of seeds created
N_seeds = length(permute_tables);

% Generate the seed information for each particle, then simulate the model
% using this seed information, storing both its output and the associated
% summary statistics
parfor k = 1:N_parts
    
    % Create an individual simulator object for this praticle, using the
    % seed information
    seed_num = mod(k, N_seeds) + 1;     % Can only use as many seeds as were provided, so start looping if more particles are requested
    part_simulators{k} = @(params) f_simulate(params, permute_tables{seed_num}, offset_tables{seed_num});
    % Run the simulator and store summaries
    part_outputs{k} = part_simulators{k}(part_params(k,:));
    part_summaries(k,:) = f_summaries(part_outputs{k});
    
end


%%% Set up the weighting of metrics by
switch metric_weighting
    
    case {'None', 'none'}
        
        % Just use discrepancy function with no additional argument provided
        f_discrep = @(summaries) f_discrepancy(summaries, target_summaries);
        
    case {'Variance', 'variance'}
        
        % Calculate sample covariance matrix from the initial particles
        C_S = cov(part_summaries);
        % Use diagonals of this to extract only variances (no covariances)
        C_S = diag(diag(C_S));
        % Store the inverse, so that it only need be calculated once
        invC_S = inv(C_S);
        % Use discrepancy function with this inverse matrix included
        f_discrep = @(summaries) f_discrepancy(summaries, target_summaries, invC_S);
        
    case {'Mahalanobis', 'mahalanobis'}
    
        % Calculate sample covariance matrix from the initial particles
        C_S = cov(part_summaries);
        % Store the inverse, so that it only need be calculated once
        invC_S = inv(C_S);
        % Use discrepancy function with this inverse matrix included
        f_discrep = @(summaries) f_discrepancy(summaries, target_summaries, invC_S);
    
end
    

% Calculate discrepancies for each particle based on this distance measure
part_Ds = f_discrepancy(part_summaries, target_summaries, invC_S);

% Loop until stopping criteria is hit
looping = 1;
while looping
    
    % First, sort all the particles according to their discrepancy
    [part_Ds, ordering] = sort(part_Ds);
    part_outputs = part_outputs(ordering);
    part_thetas = part_thetas(ordering,:);
    part_summaries = part_summaries(ordering,:);
    
    % Now select the target discrepancy as the worst particle of the kept
    % fraction
    target_D = part_Ds(worst_keep);
    
    % Select which particles, from those that are kept, to resample the
    % discarded particles onto
    switch resample_weighting
        case 'equal'
            
            % Select particles at random from the kept particles for each
            % rejected particle
            selection = randi(worst_keep, N_parts-worst_keep, 1);
            
        case 'weighted'
            
            % Calculate weights for each particle. These are found by
            % assuming that the worst kept particle is three standard
            % deviations away from D = 0. Due to normalisation of weights,
            % which here takes place inside randsample, even if all
            % particles are far from D = 0 it will not be an issue.
            part_logweights = - 9 * part_Ds.^2 / part_Ds(worst_keep)^2;
            part_logweights = part_logweights + max(part_logweights);   % Normalise so largest weight is 1 (0 on log scale)
            part_weights = exp(part_logweights);
            
            % Now select particles from the kept particles, according to
            % their weights, to decide which particles to copy onto
            selection = randsample(worst_keep, N_parts-worst_keep, true, part_weights);
            
    end
    
    % Now perform the particle copy
    part_thetas(worst_keep+1:end,:) = part_thetas(selection,:);
    part_outputs(worst_keep+1:end) = part_outputs(selection);
    part_summaries(worst_keep+1:end,:) = part_summaries(selection,:);
    part_Ds(worst_keep+1:end) = part_Ds(selection);
    
    % Now the particles are 'jiggled' so that we return to having a set of
    % unique particles. This step can also be thought of as performing
    % local exploration, and is sometimes called 'mutation'. MCMC steps are
    % used to perform mutation, and an advantage of SMC approaches is that
    % the positions of all particles can be used to decide on a sensible
    % jumping distribution.
    
    % According to the jumping distribution, movement steps are processed
    % in different ways, so split up here
    switch jumping_type
        case 'normal'
            
            % Create a weighted covariance matrix to use in a multivariate
            % normal jumping distribution. 2.38^2/N is the 'optimal' factor
            % under some assumptions
            Ctheta = 2.38^2 / N_theta * cov(part_thetas);
            
            % Perform one MCMC step to evaluate how many are expected to be
            % required
            parfor k = 1:N_parts
                [part_thetas(k,:), part_outputs{k}, part_summaries(k,:), part_Ds(k), accepted(k)] = MVN_move(1, part_thetas(k,:), part_outputs{k}, part_summaries(k,:), part_Ds(k), Ctheta, theta_mins, theta_maxs, scale_param, target_D, part_simulators{k}, f_summaries, f_discrep); 
            end
            
            % Calculate the acceptance rate, and ensure that the case where
            % all or no particles are accepted does not result in a number
            % of MCMC steps that cannot be calculated
            est_accept_rate = mean(accepted);
            est_accept_rate(est_accept_rate == 0) = 1e-6;
            est_accept_rate(est_accept_rate == 1) = 1 - 1e-6;
            
            % Calculate the expected number of MCMC steps required from the
            % estimated acceptance rate. The number of steps cannot go
            % above a user-specified value
            R = ceil( log(0.05) / log(1 - est_accept_rate) );
            R = min([R,max_MCMC_steps]);
            
            % Perform the remaining MCMC steps
            parfor k = 1:N_parts
                [part_thetas(k,:), part_outputs{k}, part_summaries(k,:), part_Ds(k), accepted(k)] = MVN_move(R-1, part_thetas(k,:), part_outputs{k}, part_summaries(k,:), part_Ds(k), Ctheta, theta_mins, theta_maxs, scale_param, target_D, part_simulators{k}, f_summaries, f_discrep);
            end
            
        case 'transform'
            
            % Scale all particle thetas into [0, 1]
            part_Xs = (part_thetas - theta_mins) ./ (theta_maxs - theta_mins);
            
            % Now use a separate function to create a Gaussian mixture
            % model that becomes the jumping distribution, J (separate
            % because it gets messy)
            [J, part_logQs, Bmix_fit] = fitGMM(part_Xs);
            
            % Perform one MCMC step to evaluate how many are expected to be
            % required
            parfor k = 1:N_parts
                [part_thetas(k,:), part_outputs{k}, part_summaries(k,:), part_Ds(k), part_logQs(k), accepted(k)] = GMM_move(1, J, part_thetas(k,:), part_outputs{k}, part_summaries(k,:), part_Ds(k), part_logQs(k), N_theta, theta_mins, theta_maxs, scale_param, target_D, Bmix_fit, part_simulators{k}, f_summaries, f_discrep);
            end
            
            % Calculate the acceptance rate, and ensure that the case where
            % all or no particles are accepted does not result in a number
            % of MCMC steps that cannot be calculated
            est_accept_rate = mean(accepted);
            est_accept_rate(est_accept_rate == 0) = 1e-6;
            est_accept_rate(est_accept_rate == 1) = 1 - 1e-6;
            
            % Calculate the expected number of MCMC steps required from the
            % estimated acceptance rate. The number of steps cannot go
            % above a user-specified value
            R = ceil( log(0.05) / log(1 - est_accept_rate) );
            R = min([R,max_MCMC_steps]);
            
            % Perform the remaining MCMC steps
            parfor k = 1:N_parts
                [part_thetas(k,:), part_outputs{k}, part_summaries(k,:), part_Ds(k), part_logQs(k), accepted(k)] = GMM_move(R-1, J, part_thetas(k,:), part_outputs{k}, part_summaries(k,:), part_Ds(k), part_logQs(k), N_theta, theta_mins, theta_maxs, scale_param, target_D, Bmix_fit, part_simulators{k}, f_summaries, f_discrep);
            end
            
    end
    
    % Count the number of unique particles
    [unique_thetas, ~] = unique(part_thetas, 'rows');
    N_unique = size(unique_thetas,1);
    
    % Output information on discrepancy and uniqueness if flag is set
    if verbose
        fprintf('Current discrepancy target is %g, number of unique particles is %g \n', target_D, N_unique);
    end
    
    % Check to see if the desired discrepancy target has been reached, and
    % terminate the loop if so. The structure of the loop means that a
    % mutation step will happen after the final resample, so we expect a
    % full sample of N_parts particles all satisfying the desired
    % discrepancy constraint
    if target_D < desired_D
        looping = 0;
    end
    
    % Also check to see if degeneracy is occurring - this happens when the
    % MCMC steps fail to find new locations and thus remain copies of
    % previous particles. This also terminates the loop when it becomes too
    % severe, but with additional output of a warning
    if N_unique <= N_parts / 2
        fprintf('WARNING: SMC loop terminated due to particle degeneracy. Target not reached! \n');
        looping = 0;
    end
    
    % Visualise Particles if flag is set
    if visualise
        
        % Call provided plotting function
        f_visualise(part_thetas, part_outputs, part_summaries, part_Ds);
        drawnow;
        
    end
    
end

end


function [GMModel, part_logQs, Bmix_fits] = fitGMM(part_Xs)
% This function takes the set of input particle locations scaled to [0, 1]
% and fits a Gaussian mixture model to them, so that samples can be
% generated that match the distribution of the particles

% Each parameter is handled separately in finding the appropriate
% transforms to get the data how we want it
log_Bprob = 0;
U = zeros(size(part_Xs));
Z = zeros(size(part_Xs));
log_Zprob = zeros(size(part_Xs,1),1);
Bmix_fits = zeros(size(part_Xs,2),5);
for k = 1:size(part_Xs,2)
       
    % Find a best-fit Beta distribution to particle locations
    % in this dimension of parameter space
    B_fit = betafit(part_Xs(:,k));
    
    % Use this as the starting point for finding the best fit
    % of a mixture of two beta distributions to the data (found
    % by maximising the log likelihood)
    [Bmix_fit,~] = fminsearch( @(x) -betamix_loglike(x, part_Xs(:,k)), [0, log(B_fit), log(B_fit)],optimset('Display','off'));
    
    % Convert the parameters of the Beta distributions as used
    % above into the standard Beta distribution parameters
    Bmix_fits(k,1) = 1/(1+exp(-Bmix_fit(1))); Bmix_fits(k,2:5) = exp(Bmix_fit(2:5));
    
    % Now use the CDF of this mixture of beta distributions to
    % obtain uniformly distributed data
    U(:,k) = Bmix_fits(k,1) * betainc( part_Xs(:,k), Bmix_fits(k,2), Bmix_fits(k,3) ) + (1 - Bmix_fits(k,1)) * betainc( part_Xs(:,k), Bmix_fits(k,4), Bmix_fits(k,5) );
    
    % Take this uniformly distributed data and make it normally
    % distributed
    Z(:,k) = norminv(U(:,k));
    
    % Calculate the log probability of each particle in terms
    % of the beta mixtures in each dimension. These combine
    % multiplicatively as is typical for likelihood, which is
    % additively in terms of logs
    log_Bprob = log_Bprob + log( Bmix_fits(k,1) * betapdf( part_Xs(:,k), Bmix_fits(k,2), Bmix_fits(k,3)) + (1-Bmix_fits(k,1)) * betapdf(part_Xs(:,k), Bmix_fits(k,4), Bmix_fits(k,5)));
    
    % Calculate the log probability of each particle in terms
    % of the Gaussian transform. These combine additively for
    % the same reason as above
    log_Zprob = log_Zprob + log( normpdf(Z(:,k), 0, 1) );
    
end

% Now fit a Gaussian mixture model to the normal-transformed particle
% locations
try 
    GMModel = fitgmdist(Z,3,'RegularizationValue',0.01,'Replicates',5,'Options',optimset('MaxIter',1000));
catch
    fprintf('ERROR: Gaussian mixture model fitting has failed! Entering debug mode. \n');
    keyboard
end

% Calculate the log probabilities of each particle according to this
% transformed distribution (from pdf of the GMM itself modified by the
% appropriate factors for the Beta and normal transforms
part_logQs = log(pdf(GMModel,Z)) + log_Bprob - log_Zprob;

end




function f = betamix_loglike(theta,x)
% This (internal) function is used to calculate log likelihood for a
% mixture of two beta distributions, with parameters supplied by theta, for
% the particles specified by x

% Transform to the real parameters of the Beta distribution
% (this setup is used to contrain the optimiser search
theta(1) = 1/(1+exp(-theta(1)));
theta(2:5) = exp(theta(2:5));

% Calculate the log likelihood from the PDF
f = sum( log( theta(1) * betapdf( x, theta(2), theta(3)) + (1-theta(1)) * betapdf( x, theta(4), theta(5) ) ) );

end


function [theta, output, summaries, D, moved] = MVN_move(N_moves, theta, output, summaries, D, Ctheta, theta_mins, theta_maxs, scale_param, target_D, f_simulate, f_summaries, f_discrep)
% This function performs MCMC updates ('move steps' or 'mutations') for a
% single particle, using a multivariate normal jumping distribution

% Initialise variable 'moved' as false
moved = 0;

% Perform the requested number of move steps
for k = 1:N_moves
    
    % Pick a new location from the multivariate normal centred around the
    % requested point
    prop_theta = mvnrnd(theta, Ctheta);
    
    % Instant rejection if particle moves outside boundary, so only process
    % if new proposed location is inside prior space
    if all(prop_theta <= theta_maxs) && all(prop_theta >= theta_mins)
        
        % Convert scaling parameters back into true values for simulation
        prop_params = prop_theta;
        prop_params(scale_param) = exp(prop_params(scale_param));
    
        % Calculate the associated model output, and discrepancy
        prop_output = f_simulate(prop_params);
        prop_summaries = f_summaries(prop_output);
        prop_D = f_discrep(prop_summaries);
    
        % Now accept this move if it still satisfies the current discrepancy
        % constraint (ratio of jumping distribution probabilities is one)
        if prop_D <= target_D
        
            % Update all values for this particle
            theta = prop_theta;
            output = prop_output;
            D = prop_D;
            summaries = prop_summaries;
        
            % Switch 'moved' flag to true
            moved = 1;
            
        end
    end
    
end

end




function [theta, output, summaries, D, logQ, moved] = GMM_move(N_moves, J, theta, output, summaries, D, logQ, N_theta, theta_mins, theta_maxs, scale_param, target_D, Bmix_fit, f_simulate, f_summaries, f_discrep)
% This function performs MCMC updates ('move steps' or 'mutations') for a
% single particle, using a Gaussian mixture model jumping distribution
% (previously found using function fitGMM and supplied to this function as
% 'J')

% Initialise variable 'moved' as false
moved = 0;

% Perform the requested number of move steps
for k = 1:N_moves
    
    % Pick a new location from the multivariate normal centred around the
    % requested point
    prop_Z = random(J);
    
    % Scale this back to uniform
    prop_U = normcdf(prop_Z,0,1);
    
    % Calculate the effect of this transform (multiplication of each
    % marginal is summing in terms of logs)
    log_Zprob = sum(log(normpdf(prop_Z,0,1)));
    
    % Loop over each parameter for the more complicated stuff
    log_Bprob = 0;
    prop_B = zeros(1,N_theta);
    for j = 1:N_theta
        
        % Undo the beta mixture transform (will still give a value on [0,1]
        % This is achieved using the quantile function of the beta mixture
        prop_B(j) = quantile_fun_betamix(prop_U(j), Bmix_fit(j,:) );
        
        % Process this parameter's multiplicative contribution to the
        % probability transform factor (for logQ calculation)
        log_Bprob = log_Bprob + log( Bmix_fit(j,1) * betapdf( prop_B(j), Bmix_fit(j,2), Bmix_fit(j,3) ) + (1 - Bmix_fit(j,1)) * betapdf( prop_B(j), Bmix_fit(j,4), Bmix_fit(j,5) ) );
        
    end
    
    % Transform from [0,1] back to the correct theta value
    prop_theta = theta_mins + prop_B .* (theta_maxs - theta_mins);
       
    % Calculate the probability of selecting this point using jumping
    % distribution J.
    % Previous particle location is not used, so we have
    %              q( theta_new | theta_old ) = q ( theta_new )
    prop_logQ = log_Bprob + log(pdf(J,prop_Z)) - log_Zprob;
    
    % EARLY REJECTION STEP - Before simulating, rejection might already
    % occur due to Metropolis-Hastings accept/reject algorithm, so do this
    % first
    r = rand;
    if r < exp(logQ - prop_logQ)   % M-H step accepted, so continue
        
        % Convert scaling parameters back into true values for simulation
        prop_params = prop_theta;
        prop_params(scale_param) = exp(prop_params(scale_param));
                
        % Calculate the associated model output, and discrepancy
        prop_output = f_simulate(prop_params);
        prop_summaries = f_summaries(prop_output);
        prop_D = f_discrep(prop_summaries);
        
        % Now accept this move if it still satisfies the current discrepancy
        % constraint (ratio of jumping distribution probabilities is one)
        if prop_D <= target_D
            
            % Update all values for this particle
            theta = prop_theta;
            output = prop_output;
            D = prop_D;
            summaries = prop_summaries;
            
            % Switch 'moved' flag to true
            moved = 1;
            
        end
        
    end
    
end

end





function x = quantile_fun_betamix(Q,Bmix_fit)
% This function takes an input quantile, Q, and the parameters for a 
% mixture of two beta distributions, and outputs the location where this
% quantile occurs. This is found using the bisection method

% Specify tolerance and maximum iterations for finding the quantile
tol = 1e-8;
max_iters = 200;

% Initialise bisection window as the whole domain of beta mixture, [0,1]
x_low = 0;
x_high = 1;

% Loop until quantile location is found to sufficient accuracy
err = 1e10;

iters = 0;
while err > tol && iters < max_iters
    
    % Increment iteration counter
    iters = iters + 1;
    
    % Bisect the window and check its value
    x_check = (x_low + x_high) / 2;
    % Find the value of the beta mixture CDF at this point
    Q_check = Bmix_fit(1) * betainc( x_check, Bmix_fit(2), Bmix_fit(3) ) + (1 - Bmix_fit(1)) * betainc( x_check, Bmix_fit(4), Bmix_fit(5) );
        
    % Check error at current point
    err = abs(Q_check - Q);
        
    % Select new window according to Q_check and the desired Q to be found
    if Q > Q_check
        x_low = x_check;
    else
        x_high = x_check;
    end

end

if iters == max_iters
   fprintf('Warning, finding quantile of beta mixture did not converge, error was %g \n', err); 
end
% Return found value
x = x_check;

end