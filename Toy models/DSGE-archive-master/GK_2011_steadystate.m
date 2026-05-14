addpath('/Applications/Dynare/6.4-arm64/matlab');
dynare gk_ra.mod

%% ---- Plot IRFs to government spending shock (eps_G) ----------------------

vars   = {'Y','c','I','N','r','rk','spread','Q','K','nB','theta','G','B','T'};
labels = {'Output (Y)','Consumption (c)','Investment (I)','Labor (N)', ...
          'Real rate (r)','Return on capital (rk)','Credit spread', ...
          'Tobins Q','Capital (K)','Bank net worth (nB)', ...
          'Leverage (\theta)','Gov. spending (G)','Gov. debt (B)','Taxes (T)'};

T_irf = 100;
t     = 1:T_irf;

figure('Name','IRFs – Government Spending Shock','NumberTitle','off', ...
       'Position',[100 100 1400 900]);

for i = 1:numel(vars)
    field = [vars{i} '_eps_G'];
    if isfield(oo_.irfs, field)
        irf_vec = oo_.irfs.(field);
    else
        irf_vec = zeros(1, T_irf);
    end

    subplot(4, 4, i);
    plot(t, irf_vec * 100, 'b-', 'LineWidth', 1.5);
    hold on;
    yline(0, 'k--', 'LineWidth', 0.8);
    title(labels{i}, 'FontSize', 9);
    xlabel('Quarters', 'FontSize', 8);
    ylabel('% dev. from SS', 'FontSize', 8);
    xlim([1 T_irf]);
    grid on;
end

sgtitle('Impulse Responses to a 1% Government Spending Shock (GK 2011)', ...
        'FontSize', 12, 'FontWeight', 'bold');
