import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

BLUE       = '#002147'
RED        = '#8C1515'
BLUE_MUTED = '#4a6f8a'
RED_MUTED  = '#c0624a'

_COLORS     = [BLUE, RED, BLUE_MUTED, RED_MUTED]
_LINESTYLES = ['-', '--', '-.', ':']
_MARKERS    = ['', '', '', 'o']


def show_irfs(irfs_list, variables, labels=None,
              ylabel='Deviation from SS (pp)', T_plot=100,
              figsize=(18, 5), savepath=None):
    labels = labels or [''] * len(irfs_list)
    n_var  = len(variables)
    fig, axes = plt.subplots(1, n_var, figsize=figsize, sharey=False)
    if n_var == 1:
        axes = [axes]

    for i, (ax, var) in enumerate(zip(axes, variables)):
        for j, (irf, label) in enumerate(zip(irfs_list, labels)):
            data = irf[var][:T_plot] if var in irf else np.zeros(T_plot)
            mkr  = _MARKERS[j % len(_MARKERS)]
            ax.plot(data,
                    color=_COLORS[j % len(_COLORS)],
                    linestyle=_LINESTYLES[j % len(_LINESTYLES)],
                    linewidth=1.8, marker=mkr, markersize=4, markevery=4,
                    label=label)
        ax.axhline(0, color='#888888', linewidth=0.8, linestyle=':')
        ax.set_title(var, fontsize=10, pad=6)
        ax.set_xlabel('Quarter', fontsize=9)
        if i == 0:
            ax.set_ylabel(ylabel, fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=8)
        if any(l for l in labels):
            ax.legend(fontsize=8, frameon=False)

    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=150, bbox_inches='tight')
    plt.close(fig)


def generate_irf_plots(model_results, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    irfs_Z_D     = model_results['irfs_Z_D']
    irfs_def_D   = model_results['irfs_def_D']
    dShock_def_D = model_results['dShock_def_D']
    T            = model_results['T']

    print("Generating IRF plots...")

    # Overview: welfare and spread
    show_irfs(
        [irfs_def_D], ['spread_rb', 'U_D', 'U_F'],
        labels=['Default shock'],
        savepath=output_dir / 'fig_irf_overview_welfare.png'
    )
    print("  Saved fig_irf_overview_welfare.png")

    # Overview: macro variables
    show_irfs(
        [irfs_Z_D, irfs_def_D],
        ['Y_D', 'C_D', 'w_D', 'n_inter_D', 'q_b_D', 'q_b_F'],
        labels=['TFP shock', 'Default shock'],
        savepath=output_dir / 'fig_irf_overview_macro.png'
    )
    print("  Saved fig_irf_overview_macro.png")

    # 1. Output, Consumption & Trade
    show_irfs(
        [irfs_Z_D, irfs_def_D], labels=['TFP shock (D)', 'Default shock (D)'],
        variables=['Y_D', 'Y_F', 'C_D', 'C_F', 'p', 'NX_D'],
        savepath=output_dir / 'fig_irf_goods_trade.png'
    )
    print("  Saved fig_irf_goods_trade.png")

    # 2. Labour, Capital & TFP
    show_irfs(
        [irfs_Z_D, irfs_def_D], labels=['TFP shock (D)', 'Default shock (D)'],
        variables=['N_D', 'N_F', 'K_D', 'K_F', 'I_D', 'I_F', 'Q_D', 'w_D'],
        savepath=output_dir / 'fig_irf_labour_capital.png'
    )
    print("  Saved fig_irf_labour_capital.png")

    # 3. Factor Prices
    show_irfs(
        [irfs_Z_D, irfs_def_D], labels=['TFP shock (D)', 'Default shock (D)'],
        variables=['w_D', 'w_F', 'N_D', 'N_F', 'rk_D', 'rk_F'],
        savepath=output_dir / 'fig_irf_factor_prices.png'
    )
    print("  Saved fig_irf_factor_prices.png")

    # 4. Bond Holdings & External Position
    show_irfs(
        [irfs_Z_D, irfs_def_D], labels=['TFP shock (D)', 'Default shock (D)'],
        variables=['b_D_D', 'b_F_D', 'b_D_F', 'b_F_F', 'nfa_D', 'n_inter_D'],
        savepath=output_dir / 'fig_irf_bonds_nfa.png'
    )
    print("  Saved fig_irf_bonds_nfa.png")

    # 5. Rates & Returns
    show_irfs(
        [irfs_Z_D, irfs_def_D], labels=['TFP shock (D)', 'Default shock (D)'],
        variables=['rb_actual_D', 'rb_actual_F', 'rn_D', 'rn_F', 'rdep_D', 'rdep_F', 'rk_D'],
        savepath=output_dir / 'fig_irf_rates_returns.png'
    )
    print("  Saved fig_irf_rates_returns.png")

    # 6. Fiscal
    show_irfs(
        [irfs_Z_D, irfs_def_D], labels=['TFP shock (D)', 'Default shock (D)'],
        variables=['b_gov_D', 'b_gov_F', 'TAX_D', 'TAX_F', 'def_rate_D'],
        savepath=output_dir / 'fig_irf_fiscal.png'
    )
    print("  Saved fig_irf_fiscal.png")

    # 7. Default decomposition
    irfs_def_D_plot = dict(irfs_def_D)
    irfs_def_D_plot['shock_def_D'] = dShock_def_D
    show_irfs(
        [irfs_def_D_plot], variables=['shock_def_D', 'def_rate_D'],
        labels=['Default shock (D)'], ylabel='Deviation from SS', figsize=(10, 5),
        savepath=output_dir / 'fig_irf_default_decomp.png'
    )
    print("  Saved fig_irf_default_decomp.png")

    # 8. Walras residuals (regression check)
    show_irfs(
        [irfs_Z_D], labels=['TFP shock'],
        variables=['goods_mkt_D', 'deposit_mkt_D', 'rb_D_res', 'rb_F_res', 'b_D_F_res', 'b_F_D_res'],
        ylabel='Residual',
        savepath=output_dir / 'fig_walras_residuals_tfp.png'
    )
    show_irfs(
        [irfs_Z_D], labels=['TFP shock'],
        variables=['global_goods_res', 'goods_mkt_F'],
        ylabel='Walras residual',
        savepath=output_dir / 'fig_walras_untargeted_tfp.png'
    )
    show_irfs(
        [irfs_def_D], labels=['Default shock'],
        variables=['goods_mkt_D', 'deposit_mkt_D', 'rb_D_res', 'rb_F_res', 'b_D_F_res', 'b_F_D_res'],
        ylabel='Residual',
        savepath=output_dir / 'fig_walras_residuals_def.png'
    )
    show_irfs(
        [irfs_def_D], labels=['Default shock'],
        variables=['global_goods_res', 'goods_mkt_F'],
        ylabel='Walras residual',
        savepath=output_dir / 'fig_walras_untargeted_def.png'
    )
    print("  Saved Walras residual figures")
    print("IRF plots done.")
