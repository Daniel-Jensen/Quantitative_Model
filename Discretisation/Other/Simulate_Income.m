%% Simulate_Income.m
% Simulates ONE individual's income path over N periods using the
% discretized Markov chain (Px_GMAR, x_vec) and analyses persistence.
%
% No Statistics / Econometrics Toolbox required.

clear; clc;
addpath(genpath(pwd));


%% ========================================================================
%  SETTINGS
%  ========================================================================

N   = 100000;   % simulation length
rng(42);        % reproducibility


%% ========================================================================
%  LOAD INPUTS
%  ========================================================================

load('Outputs/Px_GMAR.mat');           % pi_mix  (n×n transition matrix)
x_vec = load('Outputs/x_vec.txt');     % n×1 income grid (levels, mean = 1)
n     = length(x_vec);


%% ========================================================================
%  STATIONARY DISTRIBUTION  (power iteration)
%  ========================================================================

dist = ones(1, n) / n;
for iter = 1:20000
    dist_new = dist * pi_mix;
    if max(abs(dist_new - dist)) < 1e-14, break; end
    dist = dist_new;
end
dist = dist(:);                        % n×1


%% ========================================================================
%  SIMULATE INCOME PATH
%  ========================================================================

% Pre-compute cumulative row CDFs once
cumP = cumsum(pi_mix, 2);              % n×n

state = zeros(N, 1);

% Initial state drawn from stationary distribution
cumD0       = cumsum(dist);
state(1)    = find(cumD0 >= rand(), 1, 'first');

% Draw path
for t = 2:N
    s       = state(t-1);
    state(t) = find(cumP(s, :) >= rand(), 1, 'first');
end

income  = x_vec(state);               % N×1 income levels
log_inc = log(income);                % N×1 log income


%% ========================================================================
%  AUTOCORRELATIONS  (manual — no toolbox)
%  ========================================================================

max_lag = 10;

x_dm = log_inc - mean(log_inc);
var0 = mean(x_dm.^2);

ac = zeros(max_lag, 1);
for k = 1:max_lag
    ac(k) = mean(x_dm(1:end-k) .* x_dm(k+1:end)) / var0;
end


%% ========================================================================
%  OLS AR(1) ON LOG INCOME
%  ========================================================================

Y    = log_inc(2:end);
X    = [ones(N-1, 1), log_inc(1:end-1)];
beta = (X' * X) \ (X' * Y);

resid     = Y - X * beta;
R2        = 1 - mean(resid.^2) / mean((Y - mean(Y)).^2);
rho_ols   = beta(2);
half_life = log(0.5) / log(abs(rho_ols));


%% ========================================================================
%  DECILE TRANSITION MATRIX  (no toolbox percentiles)
%  ========================================================================

n_dec   = 10;
inc_sort = sort(income);
edges   = zeros(n_dec + 1, 1);
edges(1)       = 0;
edges(n_dec+1) = Inf;
for d = 1:n_dec-1
    idx       = round(d/n_dec * N);
    edges(d+1) = inc_sort(max(idx,1));
end

decile = zeros(N, 1);
for d = 1:n_dec
    mask          = income >= edges(d) & income < edges(d+1);
    decile(mask)  = d;
end

trans_count = zeros(n_dec, n_dec);
for t = 2:N
    trans_count(decile(t-1), decile(t)) = ...
        trans_count(decile(t-1), decile(t)) + 1;
end
row_sums  = sum(trans_count, 2);
trans_prob = trans_count ./ row_sums;   % row-stochastic


%% ========================================================================
%  PRINT RESULTS
%  ========================================================================

fprintf('=================================================================\n');
fprintf('  INCOME PERSISTENCE  —  simulated individual (N = %d)\n', N);
fprintf('=================================================================\n\n');

fprintf('--- Autocorrelations of log income ---\n');
fprintf('  Lag   AC\n');
for k = 1:max_lag
    fprintf('   %2d   %7.4f\n', k, ac(k));
end

fprintf('\n--- OLS AR(1) on log income ---\n');
fprintf('  rho_OLS    = %7.4f\n', rho_ols);
fprintf('  R-squared  = %7.4f\n', R2);
fprintf('  Half-life  = %7.2f  periods\n', half_life);

fprintf('\n--- Decile transition matrix (row = from, col = to, entries %%) ---\n');
fprintf('     ');
for d = 1:n_dec, fprintf(' D%d ', d); end
fprintf('\n');
for d = 1:n_dec
    fprintf('D%d   ', d);
    for j = 1:n_dec
        fprintf('%4.1f ', trans_prob(d,j)*100);
    end
    fprintf('\n');
end

fprintf('\n--- Diagonal: probability of staying in same decile ---\n');
for d = 1:n_dec
    fprintf('  Decile %2d :  %5.1f%%\n', d, trans_prob(d,d)*100);
end
fprintf('\n  Average stay-probability: %.1f%%\n\n', mean(diag(trans_prob))*100);

fprintf('Done.\n');
