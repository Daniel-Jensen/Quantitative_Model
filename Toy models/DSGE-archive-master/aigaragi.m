%% =========================================================
%  AIYAGARI (1994) — Minimal Implementation
%  VFI with continuous savings choice (interpolation)
%  to avoid staircase artifacts in policy functions.
%
%  Variables you can tune at the top:
%    n_a   : number of asset grid points
%    n_z   : number of income states
% =========================================================

clear; clc; close all;

%% ---- 1. PARAMETERS -------------------------------------
beta   = 0.96;    % discount factor
sigma  = 2.0;     % CRRA curvature
alpha  = 0.36;    % capital share
delta  = 0.08;    % depreciation
rho    = 0.90;    % persistence of log-income shock
sig_e  = 0.20;    % std of innovation

% Grid sizes (change these freely)
n_a    = 100;     % asset grid points
n_z    = 50;       % income states

%% ---- 2. INCOME PROCESS (Tauchen 1986) ------------------
[z_grid, Pi] = tauchen(n_z, 0, rho, sig_e, 3);
z_grid = exp(z_grid);
pi_stat = (Pi^1000) * (ones(n_z,1)/n_z);   % stationary distribution
z_grid  = z_grid / (z_grid' * pi_stat);     % normalise mean to 1
%% ---- 3. ASSET GRID -------------------------------------
a_min  = -0.5;    % borrowing limit (can be tightened)
a_max  = 30;
% Denser near borrowing constraint
a_grid = a_min + (a_max - a_min) * ...
         ((0:n_a-1)' / (n_a-1)).^1.4;

%% ---- 4. EQUILIBRIUM INTEREST RATE (bisection) ----------
r_low  =  0.01;
r_high =  1/beta - 1 - 0.005;   % just below Bewley bound

fprintf('\nSearching for equilibrium r ...\n');
for iter_eq = 1:30
    r_try = (r_low + r_high) / 2;
    w_try = (1-alpha) * ((alpha/(r_try+delta))^(alpha/(1-alpha)));
    [~, a_pol, excess] = solve_hh(r_try, w_try, beta, sigma, ...
                                   a_grid, z_grid, Pi);
    fprintf('  iter %2d | r = %.5f | excess demand = %+.4f\n', ...
            iter_eq, r_try, excess);
    if excess > 0
        r_high = r_try;
    else
        r_low  = r_try;
    end
    if abs(excess) < 1e-4, break; end
end

r_eq = r_try;
w_eq = (1-alpha) * ((alpha/(r_eq+delta))^(alpha/(1-alpha)));
fprintf('\nEquilibrium:  r = %.5f,  w = %.4f\n', r_eq, w_eq);

%% ---- 5. SOLVE AGAIN AT EQUILIBRIUM ---------------------
[c_pol, a_pol, ~] = solve_hh(r_eq, w_eq, beta, sigma, ...
                               a_grid, z_grid, Pi);

%% ---- 6. PLOT POLICY FUNCTIONS --------------------------
colors = lines(n_z);
lw     = 1.6;

figure('Color','w','Position',[100 100 1100 460]);

% --- consumption policy ---
subplot(1,2,1); hold on; box on; grid on;
for iz = 1:n_z
    plot(a_grid, c_pol(:,iz), 'Color', colors(iz,:), ...
         'LineWidth', lw, 'DisplayName', ...
         sprintf('z = %.2f', z_grid(iz)));
end
xlabel('Current assets  a','FontSize',12);
ylabel('Consumption  c(a,z)','FontSize',12);
title('Consumption policy function','FontWeight','bold','FontSize',13);
legend('Location','northwest','FontSize',9);
xlim([a_min, a_max]);

% --- savings policy ---
subplot(1,2,2); hold on; box on; grid on;
for iz = 1:n_z
    plot(a_grid, a_pol(:,iz), 'Color', colors(iz,:), ...
         'LineWidth', lw, 'DisplayName', ...
         sprintf('z = %.2f', z_grid(iz)));
end
plot([a_min a_max], [a_min a_max], 'k--', 'LineWidth', 1, ...
     'DisplayName', '45° line');
xlabel('Current assets  a','FontSize',12);
ylabel("Next-period assets  a'(a,z)",'FontSize',12);
title('Savings policy function','FontWeight','bold','FontSize',13);
legend('Location','northwest','FontSize',9);
xlim([a_min, a_max]);

sgtitle(sprintf('Aiyagari Model  |  n_a = %d,  n_z = %d,  r^* = %.4f', ...
        n_a, n_z, r_eq), 'FontSize', 14, 'FontWeight','bold');

%% =========================================================
%  LOCAL FUNCTION: solve household problem
%% =========================================================
function [c_pol, a_pol, excess_demand] = ...
         solve_hh(r, w, beta, sigma, a_grid, z_grid, Pi)

n_a = length(a_grid);
n_z = length(z_grid);

% Utility
if sigma == 1
    u  = @(c) log(max(c, 1e-10));
else
    u  = @(c) max(c,1e-10).^(1-sigma) / (1-sigma);
end

% Initialise value function
V  = zeros(n_a, n_z);
V_new = zeros(n_a, n_z);
c_pol = zeros(n_a, n_z);
a_pol = zeros(n_a, n_z);

tol = 1e-6;  maxiter = 2000;

for iter = 1:maxiter

    for iz = 1:n_z
        % Expected continuation value (interpolated)
        EV = V * Pi(iz,:)';          % n_a x 1

        for ia = 1:n_a
            cash = (1+r)*a_grid(ia) + w*z_grid(iz);
            % feasible savings: a' in [a_grid(1), cash-eps]
            a_lo = a_grid(1);
            a_hi = min(cash - 1e-6, a_grid(end));

            if a_hi <= a_lo
                % Constrained: consume everything
                c_opt  = cash - a_lo;
                a_opt  = a_lo;
            else
                % Golden-section search over a'
                [a_opt, ~] = golden_search( ...
                    @(ap) -( u(cash - ap) + ...
                             beta * interp1(a_grid, EV, ap, 'linear','extrap') ), ...
                    a_lo, a_hi, 1e-8);
                c_opt = cash - a_opt;
            end

            V_new(ia,iz) = u(c_opt) + beta * ...
                           interp1(a_grid, EV, a_opt, 'linear','extrap');
            c_pol(ia,iz) = c_opt;
            a_pol(ia,iz) = a_opt;
        end
    end

    diff = max(abs(V_new(:) - V(:)));
    V    = V_new;
    if diff < tol, break; end
end

% Stationary distribution (iterate)
mu = ones(n_a, n_z) / (n_a * n_z);
for t = 1:5000
    mu_new = zeros(n_a, n_z);
    for iz = 1:n_z
        for ia = 1:n_a
            ap = a_pol(ia, iz);
            % distribute to two nearest grid points
            idx  = sum(a_grid <= ap);
            idx  = max(1, min(n_a-1, idx));
            wt   = (ap - a_grid(idx)) / (a_grid(idx+1) - a_grid(idx));
            wt   = max(0, min(1, wt));
            for iz2 = 1:n_z
                mu_new(idx,   iz2) = mu_new(idx,   iz2) + ...
                                     Pi(iz,iz2)*(1-wt)*mu(ia,iz);
                mu_new(idx+1, iz2) = mu_new(idx+1, iz2) + ...
                                     Pi(iz,iz2)*wt    *mu(ia,iz);
            end
        end
    end
    if max(abs(mu_new(:)-mu(:))) < 1e-10, break; end
    mu = mu_new;
end

K_supply  = sum(sum(mu .* a_pol));
K_demand  = (alpha / (r + delta))^(1/(1-alpha));
excess_demand = K_demand - K_supply;
end

%% =========================================================
%  LOCAL FUNCTION: golden-section search (minimise f on [lo,hi])
%% =========================================================
function [x_opt, f_opt] = golden_search(f, lo, hi, tol)
phi = (sqrt(5)-1)/2;
a   = lo;  b = hi;
c   = b - phi*(b-a);
d   = a + phi*(b-a);
fc  = f(c);  fd = f(d);
while (b-a) > tol
    if fc < fd
        b=d; d=c; fd=fc;
        c=b-phi*(b-a); fc=f(c);
    else
        a=c; c=d; fc=fd;
        d=a+phi*(b-a); fd=f(d);
    end
end
x_opt = (a+b)/2;
f_opt = f(x_opt);
end

%% =========================================================
%  LOCAL FUNCTION: Tauchen (1986) discretisation
%  Returns log-income grid z and transition matrix Pi
%% =========================================================
function [z, Pi] = tauchen(n, mu, rho, sigma, m)
if n == 1
    z  = mu;
    Pi = 1;
    return
end
z_max = mu + m * sigma / sqrt(1-rho^2);
z_min = mu - m * sigma / sqrt(1-rho^2);
z     = linspace(z_min, z_max, n)';
d     = z(2) - z(1);
Pi    = zeros(n, n);
for i = 1:n
    Pi(i,1) = normcdf((z(1) - rho*z(i) - mu*(1-rho) + d/2) / sigma);
    Pi(i,n) = 1 - normcdf((z(n) - rho*z(i) - mu*(1-rho) - d/2) / sigma);
    for j = 2:n-1
        lo = (z(j) - rho*z(i) - mu*(1-rho) - d/2) / sigma;
        hi = (z(j) - rho*z(i) - mu*(1-rho) + d/2) / sigma;
        Pi(i,j) = normcdf(hi) - normcdf(lo);
    end
end
end