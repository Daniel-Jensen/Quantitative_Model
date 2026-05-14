%% GMAR_Fit_Main.m
% Goal:
% Fit a GMAR(1) process to match specified distributional targets (evaluated
% on the Pareto-tail-adjusted income grid), then discretize.
%
% OUTPUT FILES:
%   Px_GMAR.txt            - Transition matrix of fitted GMAR
%   x_vec.txt              - Income grid with Pareto tail
%   x_quantile.mat         - Average income by decile (2 variants: Pareto, No Pareto)
%   fitted_gmar_params.mat - Struct of fitted parameters and moment comparison
clear; clc; close all;
addpath(genpath(pwd));


%% ========================================================================
%  SECTION 1: USER-DEFINED DISTRIBUTIONAL TARGETS
%  ========================================================================
%  The default momentFun computes six moments:
%    mean(y), std(y), skew(y), P90-P10 diff, top-10% share, excess kurtosis
%
%  NOTE: All moments are evaluated on the PARETO-ADJUSTED income grid.
%        The optimizer fits GMAR parameters so that the post-Pareto
%        distribution matches the targets below.
%  ========================================================================


momentFun = @(X, dist, P) momentFun_default(X, dist, P);

%  Target values (column vector, same order as momentFun output)
target_values = [ ...
     1.0;    % 1. Mean income (normalized)
     1.39;   % 2. Std dev of income
     7;      % 3. Skewness of income
     4.2;    % 4. P90 / P10 ratio of income
     0.41;   % 5. Top-10% income share
     12.48;  % 6. Excess kurtosis of income
     0.39];  % 7. Gini coefficient

%  Labels for display (must match order above)
target_labels = { ...
    'Mean      '; ...
    'Std       '; ...
    'Skew      '; ...
    'P90/P10   '; ...
    'Top10% share '; ...
    'ExKurtosis'; ...
    'Gini      '};

%  Optimization weights
target_weights = [0; 1; 1; 1; 1; 0; 1];   % mean and kurtosis dropped, 


%% ========================================================================
%  SECTION 2: GMAR PARAMETER BOUNDS AND INITIAL GUESS
%  ========================================================================
%  The GMAR(1) process is:   y_t = mu*(1-rho) + rho*y_{t-1} + e_t
%  where e_t is drawn from a Gaussian mixture with K=2 components.
%  ========================================================================

mu = 0;
%               rho     p1      mu1     sig1   sig2
theta0   =  [ 0.935,  0.85,  0.017,  0.165,  0.535 ];  % initial guess
theta_lb =  [ 0.001,  0.51, -0.50,   0.000,  0.000 ];  % lower bounds
theta_ub =  [ 0.999,  0.99,  0.50,   5.000,  5.000 ];  % upper bounds


%% ========================================================================
%  SECTION 3: DISCRETIZATION SETTINGS
%  ========================================================================

n      = 19;      % number of discrete grid points
nmom   = 2;       % conditional moments matched per node (1, 2, 3, or 4)
method = 'even';  % grid method
space  = 3;       % half-width of grid in unconditional std dev units

% Pareto tail parameters  (same as Discretization_Main.m)
pareto_top_mass = 0.15;   % Pareto tail starts at top-X% of distribution
kappa           = 1.6;    % Pareto tail coefficient

% Quantile settings
numel_q = 10;             % number of quantiles for x_quantile output

% AR(1) robustness check (higher persistence variant)
rho_robust = 0.95;


%% ========================================================================
%  SECTION 4: OPTIMIZE GMAR PARAMETERS TO FIT TARGETS
%  ========================================================================

fprintf('=================================================================\n');
fprintf('  GMAR_Fit_Main: Fitting GMAR(1) process to distributional targets\n');
fprintf('  NOTE: Targets are evaluated on the Pareto-adjusted income grid.\n');
fprintf('=================================================================\n\n');

%  Verify that the initial guess evaluates without error
fprintf('Checking initial guess...\n');
try
    val0 = gmarObjective(theta0, mu, n, nmom, method, space, ...
                         pareto_top_mass, kappa, ...
                         momentFun, target_values, target_weights);
    fprintf('  Initial objective value: %.6f\n\n', val0);
catch ME
    error(['Initial guess failed with error:\n  %s\n' ...
           'Check your momentFun definition and GMAR parameters.'], ME.message);
end

%  Run bounded optimizer
fprintf('Running optimizer (fmincon, interior-point)...\n');
opts_fit = optimoptions('fmincon', ...
    'Display',                'iter-detailed', ...
    'TolFun',                 1e-10, ...
    'TolX',                   1e-10, ...
    'MaxFunctionEvaluations', 10000, ...
    'MaxIterations',          2000, ...
    'Algorithm',              'interior-point');

[theta_fit, fval, exitflag] = fmincon( ...
    @(theta) gmarObjective(theta, mu, n, nmom, method, space, ...
                           pareto_top_mass, kappa, ...
                           momentFun, target_values, target_weights), ...
    theta0, [], [], [], [], theta_lb, theta_ub, [], opts_fit);

if exitflag <= 0
    warning('Optimizer did not converge cleanly (exitflag = %d). Results may be suboptimal.', exitflag);
end

%  Unpack fitted parameters
rho_fit  = theta_fit(1);
p1_fit   = theta_fit(2);
mu1_fit  = theta_fit(3);
sig1_fit = theta_fit(4);
sig2_fit = theta_fit(5);

p2_fit   = 1 - p1_fit;
mu2_fit  = -(p1_fit / p2_fit) * mu1_fit;

pvec_fit   = [p1_fit,  p2_fit  ];
muvec_fit  = [mu1_fit, mu2_fit ];
sigvec_fit = [sig1_fit, sig2_fit];

%  Display fitted parameters
fprintf('\n=== Fitted GMAR Parameters ===\n');
fprintf('  rho  = %10.6f\n', rho_fit);
fprintf('  p1   = %10.6f     p2   = %10.6f\n', p1_fit,  p2_fit);
fprintf('  mu1  = %10.6f     mu2  = %10.6f\n', mu1_fit, mu2_fit);
fprintf('  sig1 = %10.6f     sig2 = %10.6f\n', sig1_fit, sig2_fit);
fprintf('  Optimizer exit flag: %d\n\n', exitflag);


%% ========================================================================
%  SECTION 5: DISCRETIZE THE FITTED GMAR PROCESS
%  ========================================================================

fprintf('Discretizing fitted GMAR process...\n');

[pi_mix, level_mix] = discreteGMAR(mu, rho_fit, pvec_fit, muvec_fit, ...
                                    sigvec_fit, n, nmom, method, space);

%  Stationary distribution (power iteration)
dist_mix     = stationaryDist(pi_mix, n);      % n×1 column vector
cum_dist_mix = cumsum(dist_mix);               % n×1 cumulative distribution

%  Pre-Pareto moments (diagnostic only — optimizer matched post-Pareto)
achieved_pre = momentFun(exp(level_mix(:)), dist_mix, pi_mix);

fprintf('\n=== Moment Fit: Pre-Pareto (diagnostic) ===\n');
fprintf('%-22s  %10s  %10s  %+10s\n', 'Moment', 'Target', 'Achieved', 'Error');
fprintf('%s\n', repmat('-', 1, 58));
for im = 1:length(target_labels)
    err = achieved_pre(im) - target_values(im);
    fprintf('%-22s  %10.5f  %10.5f  %+10.5f\n', ...
            target_labels{im}, target_values(im), achieved_pre(im), err);
end
fprintf('\n');

level_mix_log = level_mix(:);           % n×1 income grid (used for Pareto step)

%  Store: transition matrix (GMAR)
%  Write row-by-row so Px_GMAR.txt is exactly pi_mix (row i = P(from i -> all j))
fid = fopen('Outputs/Px_GMAR.txt', 'w');
for i_row = 1:n
    fprintf(fid, '%12.16f  ', pi_mix(i_row, :));
    fprintf(fid, '\n');
end
fclose(fid);

%  Store: income grid without Pareto
fid = fopen('Outputs/x_vec_nopareto.txt', 'w');
fprintf(fid, '%12.16f\n', exp(level_mix_log));
fclose(fid);


%% ========================================================================
%  SECTION 6: AVERAGE PRODUCTIVITY BY QUANTILE  (No Pareto)
%  ========================================================================

x_quantile = ones(numel_q, 2);
x_quantile(:, 2) = computeQuantileAverages(exp(level_mix_log), dist_mix, cum_dist_mix, numel_q);


%% ========================================================================
%  SECTION 7: PARETO TAIL ADJUSTMENT
%  ========================================================================
%  Replaces the top part of the grid with a Pareto distribution,
%  following Hubmer, Krusell, and Smith.
%  (Mirrors the Pareto section of Discretization_Main.m exactly.)

xgrid_pareto = applyParetoTail(level_mix_log, dist_mix, pareto_top_mass, kappa);
xgrid_pareto = xgrid_pareto / sum(dist_mix .* xgrid_pareto);   % normalize mean to 1

%  Store: Pareto-adjusted grid
fid = fopen('Outputs/x_vec.txt', 'w');
fprintf(fid, '%12.16f\n', xgrid_pareto);
fclose(fid);

%  Quantile averages with Pareto
x_quantile(:, 1) = computeQuantileAverages(xgrid_pareto, dist_mix, cum_dist_mix, numel_q);

%  Recover true average productivities (multiply out the numel_q normalization)
x_quantile = x_quantile * numel_q;

%  Report post-Pareto moments (these are what the optimizer targeted)
achieved_par = momentFun(xgrid_pareto, dist_mix, pi_mix);

fprintf('\n=== Moment Fit: Post-Pareto (optimization target) ===\n');
fprintf('%-22s  %10s  %10s  %+10s\n', 'Moment', 'Target', 'Achieved', 'Error');
fprintf('%s\n', repmat('-', 1, 58));
for im = 1:length(target_labels)
    err_par = achieved_par(im) - target_values(im);
    fprintf('%-22s  %10.5f  %10.5f  %+10.5f\n', ...
            target_labels{im}, target_values(im), achieved_par(im), err_par);
end
fprintf('\n');



%% ========================================================================
%  SECTION 9: SAVE RESULTS
%  ========================================================================

%  Quantile averages (2 columns: col 1 = Pareto, col 2 = No Pareto)
save('Outputs/x_quantile.mat', 'x_quantile');

%  Fitted parameters and moment comparison (post-Pareto)
fitted_params = struct( ...
    'mu',    mu,       'rho',  rho_fit,  ...
    'p1',    p1_fit,   'p2',   p2_fit,   ...
    'mu1',   mu1_fit,  'mu2',  mu2_fit,  ...
    'sig1',  sig1_fit, 'sig2', sig2_fit, ...
    'n',     n,        'nmom', nmom,     ...
    'method', method,  'space', space);

save('Outputs/fitted_gmar_params.mat', 'fitted_params', 'theta_fit', ...
     'target_values', 'target_labels', 'target_weights', 'achieved_par');

fprintf('=== Results saved to Outputs/ ===\n');
fprintf('  Transition matrix:         Outputs/Px_GMAR.txt\n');
fprintf('  Income grid (no Pareto):   Outputs/x_vec_nopareto.txt\n');
fprintf('  Income grid (Pareto):      Outputs/x_vec.txt\n');
fprintf('  Decile averages:           Outputs/x_quantile.mat\n');
fprintf('  Fitted parameters:         Outputs/fitted_gmar_params.mat\n');
fprintf('Done.\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

% -------------------------------------------------------------------------
% momentFun_default
%   Default moment function called via the momentFun handle in Section 1.
%   Returns 7 moments:
%     [mean(y); std(y); skew(y); P90/P10 ratio; top-10% share; excess kurtosis; Gini]
%
%   Inputs:
%     X    - (n×1) income grid (levels, sorted ascending)
%     dist - (n×1) stationary probability weights (sum to 1)
%     P    - (n×n) transition matrix, P(i,j) = Prob(state j | state i)
%   Output:
%     moments - (7×1) column vector of computed moment values
% -------------------------------------------------------------------------
function moments = momentFun_default(X, dist, ~)

    X    = X(:);      % ensure column vectors
    dist = dist(:);

    % --- Income moments ---
    mu_x   = sum(dist .* X);
    dev    = X - mu_x;
    var_x  = sum(dist .* dev.^2);
    sig_x  = sqrt(var_x);
    skew_x = sum(dist .* (dev ./ sig_x).^3);
    kurt_x = sum(dist .* (dev ./ sig_x).^4) - 3;   % excess kurtosis (0 for Gaussian)

    % --- P90 / P10 income ratio (linearly interpolated) ---
    %   Grid is sorted ascending; interpolation avoids step-function error
    %   which is especially large in the sparse Pareto tail region.
    cumD = cumsum(dist);

    p10 = interp1([0; cumD], [X(1); X], 0.10, 'linear', 'extrap');
    p90 = interp1([0; cumD], [X(1); X], 0.90, 'linear', 'extrap');
    p90p10 = p90 / p10;

    % --- Top-10% income share of total income ---
    %   share = E[y * 1{top 10%}] / E[y]
    %   cut90 used as grid index to sum all mass above the 90th percentile
    cut90 = find(cumD >= 0.90, 1, 'first');
    if isempty(cut90), cut90 = length(X); end
    share_top10 = sum(dist(cut90:end) .* X(cut90:end)) / mu_x;

    % --- Gini coefficient (Brown formula on sorted grid) ---
    L      = cumsum(dist .* X) / mu_x;   % Lorenz curve ordinates
    L_prev = [0; L(1:end-1)];
    gini   = 1 - sum(dist .* (L + L_prev));

    moments = [mu_x; sig_x; skew_x; p90p10; share_top10; kurt_x; gini];

end


% -------------------------------------------------------------------------
% gmarObjective
%   Objective function for the GMAR parameter optimizer.
%   Returns the weighted sum of squared moment errors, evaluated on the
%   Pareto-adjusted income grid.
%
%   Inputs:
%     theta              - [rho, p1, mu1, sig1, sig2]
%     mu                 - fixed unconditional mean
%     n, nmom, method, space - discretization settings
%     pareto_top_mass    - Pareto tail mass threshold
%     kappa              - Pareto tail exponent
%     momentFun          - function handle @(X,dist,P) → [L×1]
%     target_values      - [L×1] target moments
%     target_weights     - [L×1] non-negative optimization weights
%   Output:
%     obj - scalar objective value (lower = better fit)
% -------------------------------------------------------------------------
function obj = gmarObjective(theta, mu, n, nmom, method, space, ...
                              pareto_top_mass, kappa, ...
                              momentFun, target_values, target_weights)

    % Unpack parameters
    rho  = theta(1);
    p1   = theta(2);
    mu1  = theta(3);
    sig1 = theta(4);
    sig2 = theta(5);

    p2   = 1 - p1;
    mu2  = -(p1 / p2) * mu1;   % zero-mean mixture shock constraint

    pvec   = [p1,  p2  ];
    muvec  = [mu1, mu2 ];
    sigvec = [sig1, sig2];

    try
        % Suppress discreteGMAR's internal warnings during optimization
        ws = warning('off', 'all');

        [P, X] = discreteGMAR(mu, rho, pvec, muvec, sigvec, n, nmom, method, space);

        warning(ws);    % restore warning state

        % Stationary distribution (power iteration, fast convergence)
        x_stat = ones(1, n) / n;
        for iter = 1:5000
            xnew = x_stat * P;
            if max(abs(xnew - x_stat)) < 1e-12
                break;
            end
            x_stat = xnew;
        end
        dist = (x_stat / sum(x_stat))';   % n×1 column vector

        % Apply Pareto tail, normalize to mean = 1, evaluate moments
        xgrid_par = applyParetoTail(X(:), dist, pareto_top_mass, kappa);
        xgrid_par = xgrid_par / sum(dist .* xgrid_par);

        % Compute moments and weighted squared error
        achieved = momentFun(xgrid_par, dist, P);
        errors   = (achieved - target_values) .* target_weights;
        obj      = sum(errors.^2);

    catch
        % Return a large penalty if GMAR discretization fails at this theta
        warning(ws);
        obj = 1e10;
    end

end


% -------------------------------------------------------------------------
% applyParetoTail
%   Replaces the top portion of the income grid with Pareto quantile
%   function values, following Hubmer, Krusell, and Smith.
%
%   Inputs:
%     X_log           - (n×1) log-income grid
%     dist            - (n×1) stationary distribution
%     pareto_top_mass - fraction of mass in the Pareto tail
%     kappa           - Pareto tail exponent
%   Output:
%     xgrid_pareto - (n×1) income grid (levels) with Pareto tail applied
% -------------------------------------------------------------------------
function xgrid_pareto = applyParetoTail(X_log, dist, pareto_top_mass, kappa)

    n     = length(X_log);
    X_log = X_log(:);
    dist  = dist(:);

    cum_dist_rev = cumsum(dist, 'reverse');
    cutoff = find(cum_dist_rev < pareto_top_mass, 1, 'first');
    if isempty(cutoff)
        xgrid_pareto = exp(X_log);
        return;
    end
    x_m = 1 - cum_dist_rev(cutoff);

    fx      = ones(n, 1);
    arg_par = ones(n, 1);
    for ii = cutoff:n
        fx(ii)      = 1 - cum_dist_rev(ii);
        arg_par(ii) = (fx(ii) - x_m) / (1 - x_m);
    end

    x_par = ones(n, 1);
    for ii = cutoff:n
        x_par(ii) = exp(X_log(cutoff)) / ((1 - arg_par(ii))^(1/kappa));
    end

    level_exp    = exp(X_log);
    xgrid_pareto = [level_exp(1:cutoff); x_par(cutoff+1:n)];

end


% -------------------------------------------------------------------------
% stationaryDist
%   Computes the stationary distribution of a Markov chain by power
%   iteration.
%
%   Output:
%     dist - (n×1) stationary probability vector
% -------------------------------------------------------------------------
function dist = stationaryDist(P, n)

    x   = ones(1, n) / n;   % uniform initialization
    err = 1;
    while err > 1e-12
        xnew = x * P;
        err  = max(abs(xnew - x));
        x    = xnew;
    end
    dist = (x / sum(x))';   % normalize and return as column vector

end


% -------------------------------------------------------------------------
% computeQuantileAverages
%   Computes average income within each of nq quantiles, following the
%   exact same logic as Discretization_Main.m.  Returns values normalized
%   by nq; multiply by nq to obtain true conditional expectations.
%
%   Output:
%     xq - (nq×1) probability-weighted average income per quantile,
%          normalized so that xq * nq = true conditional expectation
% -------------------------------------------------------------------------
function xq = computeQuantileAverages(grid, dist, cum_dist, nq)

    grid     = grid(:);
    dist     = dist(:);
    cum_dist = cum_dist(:);
    n        = length(grid);

    xq = ones(nq, 1);

    % --- First quantile ---
    iq = 1;
    a  = 1;
    b  = find(cum_dist > iq/nq, 1, 'first');
    if isempty(b), b = n; end

    if b > 1
        cum_prev = cum_dist(b-1);
    else
        cum_prev = 0;
    end
    xq(1) = sum(dist(a:b-1) .* grid(a:b-1)) + grid(b) * (iq/nq - cum_prev);
    a = b;

    % --- Middle quantiles ---
    for iq = 2:nq-1
        b = find(cum_dist > iq/nq, 1, 'first');
        if isempty(b), b = n; end

        if b == a
            xq(iq) = grid(b) / nq;
        else
            xq(iq) = sum(dist(a+1:b-1) .* grid(a+1:b-1)) + ...
                     (cum_dist(a) - (iq-1)/nq) * grid(a)  + ...
                     grid(b) * (iq/nq - cum_dist(b-1));
        end
        a = b;
    end

    % --- Last quantile ---
    xq(end) = sum(dist(a+1:end) .* grid(a+1:end)) + ...
              (cum_dist(a) - (nq-1)/nq) * grid(a);

end
