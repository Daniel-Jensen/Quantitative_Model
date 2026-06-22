"""
TPI figures (7 total) — all saved to output_dir.

Takes the dict returned by tpi.run_tpi() and produces:
  fig_tpi_spread_mitigation.png
  fig_tpi_effectiveness.png
  fig_tpi_welfare_macro.png
  fig_tpi_welfare_bar.png
  fig_tpi_welfare_spread.png
  fig_tpi_bond_price.png
  fig_tpi_utility.png
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

BLUE       = '#002147'
RED        = '#8C1515'
BLUE_MUTED = '#4a6f8a'
RED_MUTED  = '#c0624a'


def _get_var(irf, var, T_plot):
    if var == 'spread_rb' and var not in irf:
        return (irf['rb_actual_D'][:T_plot] - irf['rb_actual_F'][:T_plot]) * 100
    return irf[var][:T_plot] * 100 if var in irf else np.zeros(T_plot)


def generate_tpi_plots(tpi_results, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    irfs_tpi    = tpi_results['irfs_tpi']
    gamma_values = tpi_results['gamma_values']
    gamma_labels = tpi_results['gamma_labels']
    TPI_COLORS  = tpi_results['TPI_COLORS']
    TPI_LSTYLES = tpi_results['TPI_LSTYLES']
    TPI_MARKERS = tpi_results['TPI_MARKERS']
    dW_D        = tpi_results['dW_D']
    dW_F        = tpi_results['dW_F']
    W_D         = tpi_results['W_D']
    W_F         = tpi_results['W_F']
    gammas_fine = tpi_results['gammas_fine']
    peak_arr    = tpi_results['peak_arr']
    cost_arr    = tpi_results['cost_arr']
    frac_closed = tpi_results['frac_closed']
    peak_no_tpi = tpi_results['peak_no_tpi']

    print("Generating TPI plots...")

    # ── Figure 1: Spread Mitigation ───────────────────────────────────────────
    T_plot1  = 20
    fig1_vars   = ['spread_rb', 'def_rate_D', 'q_b_D', 'rb_actual_D', 'n_inter_D', 'cb_buy_D']
    fig1_titles = ['Bond Yield Spread  (rb_D − rb_F)', 'Default Rate  def_rate_D',
                   'Bond Price  q_b_D', 'Bond Return  rb_actual_D',
                   'Bank Net Worth  n_inter_D', 'CB Bond Purchases  cb_buy_D']

    fig1, axes1 = plt.subplots(2, 3, figsize=(18, 9))
    for ax, var, title in zip(axes1.flatten(), fig1_vars, fig1_titles):
        for j, g in enumerate(gamma_values):
            ax.plot(_get_var(irfs_tpi[g], var, T_plot1),
                    color=TPI_COLORS[j], linestyle=TPI_LSTYLES[j], linewidth=1.8,
                    marker=TPI_MARKERS[j], markersize=4, markevery=8, label=gamma_labels[j])
        ax.axhline(0, color='#888888', linewidth=0.8, linestyle=':')
        ax.set_title(title, fontsize=10, pad=6)
        ax.set_xlabel('Quarter', fontsize=9); ax.set_ylabel('pp dev. from SS', fontsize=9)
        ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=8)
    axes1.flatten()[0].legend(fontsize=8, frameon=False, loc='upper right')
    fig1.suptitle('Figure 1: TPI — Spread Mitigation under Default Shock (Country D)', fontsize=12, y=1.01)
    fig1.tight_layout()
    fig1.savefig(output_dir / 'fig_tpi_spread_mitigation.png', dpi=150, bbox_inches='tight')
    plt.close(fig1)
    print("  Saved fig_tpi_spread_mitigation.png")

    # ── Figure 2: Effectiveness ───────────────────────────────────────────────
    fig2, axes2 = plt.subplots(1, 3, figsize=(17, 5))
    axes2[0].plot(gammas_fine, peak_arr * 100, color=BLUE, linewidth=2)
    axes2[0].scatter([g for g in gamma_values],
                     [peak_arr[np.argmin(np.abs(gammas_fine - g))] * 100 for g in gamma_values],
                     color=TPI_COLORS, s=60, zorder=5)
    axes2[0].set_xlabel('Feedback gain  γ', fontsize=9); axes2[0].set_ylabel('Peak spread dev. (pp)', fontsize=9)
    axes2[0].set_title('Peak Spread Remaining', fontsize=10, pad=6)
    axes2[0].axhline(0, color='#888888', linewidth=0.8, linestyle=':')
    axes2[0].spines[['top', 'right']].set_visible(False)

    axes2[1].plot(gammas_fine, frac_closed, color=RED, linewidth=2)
    axes2[1].scatter([g for g in gamma_values],
                     [frac_closed[np.argmin(np.abs(gammas_fine - g))] for g in gamma_values],
                     color=TPI_COLORS, s=60, zorder=5)
    axes2[1].axhline(100, color='#888888', linewidth=0.8, linestyle='--', label='Full closure')
    axes2[1].set_xlabel('Feedback gain  γ', fontsize=9); axes2[1].set_ylabel('%', fontsize=9)
    axes2[1].set_title('Fraction of Peak Spread Closed', fontsize=10, pad=6)
    axes2[1].set_ylim([-5, 115]); axes2[1].legend(fontsize=8, frameon=False)
    axes2[1].spines[['top', 'right']].set_visible(False)

    axes2[2].plot(gammas_fine, cost_arr, color=BLUE_MUTED, linewidth=2)
    axes2[2].scatter([g for g in gamma_values],
                     [cost_arr[np.argmin(np.abs(gammas_fine - g))] for g in gamma_values],
                     color=TPI_COLORS, s=60, zorder=5)
    axes2[2].set_xlabel('Feedback gain  γ', fontsize=9)
    axes2[2].set_ylabel('∑ cb_buy_D × q_b_D  (D-goods·quarters)', fontsize=9)
    axes2[2].set_title('CB Balance-Sheet Cost', fontsize=10, pad=6)
    axes2[2].spines[['top', 'right']].set_visible(False)
    for j, g in enumerate(gamma_values):
        for ax in axes2:
            ax.axvline(g, color=TPI_COLORS[j], alpha=0.25, linewidth=0.8, linestyle=':')
    fig2.suptitle('Figure 2: TPI — Does It Close the Spread? Effectiveness vs Cost', fontsize=12, y=1.01)
    fig2.tight_layout()
    fig2.savefig(output_dir / 'fig_tpi_effectiveness.png', dpi=150, bbox_inches='tight')
    plt.close(fig2)
    print("  Saved fig_tpi_effectiveness.png")

    # ── Figure 3: Welfare & Macro ─────────────────────────────────────────────
    T_plot3 = 60
    fig3_vars   = ['U_D', 'U_F', 'C_D', 'C_F', 'Y_D', 'Y_F', 'TAX_D', 'b_gov_D']
    fig3_titles = ['Welfare  U_D', 'Welfare  U_F', 'Consumption  C_D', 'Consumption  C_F',
                   'Output  Y_D', 'Output  Y_F', 'Lump-sum Tax  TAX_D', 'Govt Debt  b_gov_D']
    fig3, axes3 = plt.subplots(2, 4, figsize=(22, 9))
    for ax, var, title in zip(axes3.flatten(), fig3_vars, fig3_titles):
        for j, g in enumerate(gamma_values):
            data = irfs_tpi[g][var][:T_plot3] * 100 if var in irfs_tpi[g] else np.zeros(T_plot3)
            ax.plot(data, color=TPI_COLORS[j], linestyle=TPI_LSTYLES[j], linewidth=1.8,
                    marker=TPI_MARKERS[j], markersize=4, markevery=8, label=gamma_labels[j])
        ax.axhline(0, color='#888888', linewidth=0.8, linestyle=':')
        ax.set_title(title, fontsize=10, pad=6)
        ax.set_xlabel('Quarter', fontsize=9); ax.set_ylabel('% / pp dev. from SS', fontsize=9)
        ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=8)
    axes3.flatten()[0].legend(fontsize=8, frameon=False, loc='lower right')
    fig3.suptitle('Figure 3: TPI — Welfare & Macro Effects under Default Shock (Country D)', fontsize=12, y=1.01)
    fig3.tight_layout()
    fig3.savefig(output_dir / 'fig_tpi_welfare_macro.png', dpi=150, bbox_inches='tight')
    plt.close(fig3)
    print("  Saved fig_tpi_welfare_macro.png")

    # ── Figure 4: Welfare Bar Chart ───────────────────────────────────────────
    x = np.arange(len(gamma_values)); width = 0.35
    _xlabs = [f'γ={g}' for g in gamma_values]
    fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(14, 6))

    def _annotate_bars(ax, bars, vals):
        for bar, v in zip(bars, vals):
            h = bar.get_height(); sign = 1 if v >= 0 else -1
            ax.text(bar.get_x() + bar.get_width()/2, h + sign * (0.005 * abs(h) + 0.001),
                    f'{v:+.3f}', ha='center', va='bottom' if v >= 0 else 'top', fontsize=7)

    ax4a.bar(x - width/2, W_D, width, color=RED,  label='Country D (Distressed)', alpha=0.85)
    ax4a.bar(x + width/2, W_F, width, color=BLUE, label='Country F (Partner)',    alpha=0.85)
    ax4a.axhline(0, color='#444444', linewidth=0.8)
    ax4a.set_xticks(x); ax4a.set_xticklabels(_xlabs, fontsize=9)
    ax4a.set_title('Discounted Welfare Deviation\n(Σ β^t · U · 100,  t = 0…99)', fontsize=10)
    ax4a.set_ylabel('% of SS quarterly consumption', fontsize=9)
    ax4a.legend(fontsize=9, frameon=False); ax4a.spines[['top', 'right']].set_visible(False)
    _annotate_bars(ax4a, list(ax4a.patches[:len(gamma_values)]), W_D)
    _annotate_bars(ax4a, list(ax4a.patches[len(gamma_values):]), W_F)

    ax4b.bar(x - width/2, dW_D, width, color=RED,  label='Country D (Distressed)', alpha=0.85)
    ax4b.bar(x + width/2, dW_F, width, color=BLUE, label='Country F (Partner)',    alpha=0.85)
    ax4b.axhline(0, color='#444444', linewidth=0.8)
    ax4b.set_xticks(x); ax4b.set_xticklabels(_xlabs, fontsize=9)
    ax4b.set_title('Welfare Gain vs No-TPI (γ = 0)\n(ΔW > 0 = better off with TPI)', fontsize=10)
    ax4b.set_ylabel('Δ% of SS quarterly consumption', fontsize=9)
    ax4b.legend(fontsize=9, frameon=False); ax4b.spines[['top', 'right']].set_visible(False)
    _annotate_bars(ax4b, list(ax4b.patches[:len(gamma_values)]), dW_D)
    _annotate_bars(ax4b, list(ax4b.patches[len(gamma_values):]), dW_F)

    fig4.suptitle('Figure 4: TPI — Discounted Welfare Comparison', fontsize=12, y=1.01)
    fig4.tight_layout()
    fig4.savefig(output_dir / 'fig_tpi_welfare_bar.png', dpi=150, bbox_inches='tight')
    plt.close(fig4)
    print("  Saved fig_tpi_welfare_bar.png")

    # ── Figure 5: Why the Spread Persists ────────────────────────────────────
    T_plot5  = 50
    fig5_sub = [0, 5, 10]
    _c5  = [TPI_COLORS[gamma_values.index(g)] for g in fig5_sub]
    _ls5 = [TPI_LSTYLES[gamma_values.index(g)] for g in fig5_sub]
    _mk5 = [TPI_MARKERS[gamma_values.index(g)] for g in fig5_sub]
    _lb5 = [gamma_labels[gamma_values.index(g)] for g in fig5_sub]

    def _plot_lines(ax, var):
        for g, c, ls, mk, lb in zip(fig5_sub, _c5, _ls5, _mk5, _lb5):
            data = irfs_tpi[g][var][:T_plot5] * 100 if var in irfs_tpi[g] else np.zeros(T_plot5)
            ax.plot(data, color=c, linestyle=ls, linewidth=1.9,
                    marker=mk, markersize=4, markevery=8, label=lb)
        ax.axhline(0, color='#888888', lw=0.8, ls=':')
        ax.set_xlabel('Quarter', fontsize=9)
        ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=8)

    fig5, axes5 = plt.subplots(2, 3, figsize=(18, 10))
    ax5 = axes5.flatten()

    _plot_lines(ax5[0], 'def_rate_D')
    ax5[0].set_title('Default Rate  def_rate_D\n[fundamental spread driver — TPI-invariant]', fontsize=10, pad=6)
    ax5[0].set_ylabel('pp dev. from SS', fontsize=9); ax5[0].legend(fontsize=8, frameon=False)

    _plot_lines(ax5[1], 'spread_rb')
    ax5[1].set_title('Yield Spread  spread_rb = rb_D − rb_F\n[δ_b = 0.10 → insensitive to q_b_D]', fontsize=10, pad=6)
    ax5[1].set_ylabel('pp dev. from SS', fontsize=9)

    for g, c, ls, mk, lb in zip(fig5_sub, _c5, _ls5, _mk5, _lb5):
        tr_sp = (irfs_tpi[g]['rb_actual_D'][:T_plot5] - irfs_tpi[g]['rb_actual_F'][:T_plot5]) * 100
        ax5[2].plot(tr_sp, color=c, linestyle=ls, lw=1.9, marker=mk, markersize=4, markevery=8, label=lb)
    ax5[2].axhline(0, color='#888888', lw=0.8, ls=':')
    ax5[2].set_title('Total-Return Spread  rb_actual_D − rb_actual_F\n[includes capital gain/loss — more TPI-responsive]', fontsize=10, pad=6)
    ax5[2].set_ylabel('pp dev. from SS', fontsize=9); ax5[2].set_xlabel('Quarter', fontsize=9)
    ax5[2].legend(fontsize=8, frameon=False); ax5[2].spines[['top', 'right']].set_visible(False); ax5[2].tick_params(labelsize=8)

    _plot_lines(ax5[3], 'n_inter_D')
    ax5[3].set_title('Bank Net Worth  n_inter_D\n[GK IC relaxes: fewer bonds → lower θ_tgt → n_inter rises]', fontsize=10, pad=6)
    ax5[3].set_ylabel('% dev. from SS', fontsize=9)

    _plot_lines(ax5[4], 'U_D')
    ax5[4].set_title('Domestic Welfare  U_D\n[recapitalisation → output/consumption recovery]', fontsize=10, pad=6)
    ax5[4].set_ylabel('% dev. from SS (÷ C_ss)', fontsize=9)

    x5 = np.arange(len(gamma_values)); w5 = 0.35
    bars_D = ax5[5].bar(x5 - w5/2, dW_D, w5, color=RED,  label='ΔW_D  (Distressed)', alpha=0.85)
    bars_F = ax5[5].bar(x5 + w5/2, dW_F, w5, color=BLUE, label='ΔW_F  (Partner)',    alpha=0.85)
    ax5[5].axhline(0, color='#444', lw=0.8)
    ax5[5].set_xticks(x5); ax5[5].set_xticklabels([f'γ={g}' for g in gamma_values], fontsize=8)
    ax5[5].set_title('Welfare Gain vs No-TPI  (γ = 0)\nΔW = Σ β^t · ΔU · 100,  t = 0…99', fontsize=10, pad=6)
    ax5[5].set_ylabel('Δ% of SS quarterly consumption', fontsize=9)
    ax5[5].legend(fontsize=8, frameon=False); ax5[5].spines[['top', 'right']].set_visible(False); ax5[5].tick_params(labelsize=8)
    for bar, v in zip(list(bars_D) + list(bars_F), list(dW_D) + list(dW_F)):
        ax5[5].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.003 + (0 if v >= 0 else -0.03),
                    f'{v:+.3f}', ha='center', va='bottom', fontsize=6.5)

    fig5.suptitle('Figure 5: TPI — Why the Spread Persists Despite Welfare Gains\n'
                  'Default shock in Country D  |  TPI rule: cb_buy_D = γ × spread_rb  (closed-loop)',
                  fontsize=11, y=1.02)
    fig5.tight_layout()
    fig5.savefig(output_dir / 'fig_tpi_welfare_spread.png', dpi=150, bbox_inches='tight')
    plt.close(fig5)
    print("  Saved fig_tpi_welfare_spread.png")

    # ── Figure 6: Bond Prices ─────────────────────────────────────────────────
    T_plot6 = 60
    fig6, axes6 = plt.subplots(1, 2, figsize=(14, 5))
    for j, g in enumerate(gamma_values):
        kw = dict(color=TPI_COLORS[j], linestyle=TPI_LSTYLES[j], linewidth=1.9,
                  marker=TPI_MARKERS[j], markersize=4, markevery=8, label=gamma_labels[j])
        axes6[0].plot(irfs_tpi[g]['q_b_D'][:T_plot6] * 100, **kw)
        axes6[1].plot(irfs_tpi[g]['q_b_F'][:T_plot6] * 100, **kw)
    for ax, title in zip(axes6, ['Bond Price  q_b_D  (Country D — Distressed)',
                                   'Bond Price  q_b_F  (Country F — Partner)']):
        ax.axhline(0, color='#888888', linewidth=0.8, linestyle=':')
        ax.set_title(title, fontsize=10, pad=6); ax.set_xlabel('Quarter', fontsize=9)
        ax.set_ylabel('% dev. from SS', fontsize=9); ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=8)
    axes6[0].legend(fontsize=8, frameon=False)
    fig6.suptitle('Figure 6: TPI — Bond Prices under Default Shock\n'
                  'CB purchases compress D-bond yields by raising q_b_D', fontsize=11, y=1.02)
    fig6.tight_layout()
    fig6.savefig(output_dir / 'fig_tpi_bond_price.png', dpi=150, bbox_inches='tight')
    plt.close(fig6)
    print("  Saved fig_tpi_bond_price.png")

    # ── Figure 7: Household Utility ───────────────────────────────────────────
    T_plot7 = 60
    fig7, axes7 = plt.subplots(1, 2, figsize=(14, 5))
    for j, g in enumerate(gamma_values):
        lw = 2.5 if g == 0 else 1.9; ls = '--' if g == 0 else TPI_LSTYLES[j]
        kw = dict(color=TPI_COLORS[j], linestyle=ls, linewidth=lw,
                  marker=TPI_MARKERS[j], markersize=4, markevery=8,
                  label=gamma_labels[j], zorder=(10 if g == 0 else 5))
        axes7[0].plot(irfs_tpi[g]['U_D'][:T_plot7] * 100, **kw)
        axes7[1].plot(irfs_tpi[g]['U_F'][:T_plot7] * 100, **kw)
    for ax, title in zip(axes7, ['Utility  U_D  — Country D (Distressed)\n[GHH composite, % dev. from SS]',
                                   'Utility  U_F  — Country F (Partner)\n[GHH composite, % dev. from SS]']):
        ax.axhline(0, color='#888888', linewidth=0.8, linestyle=':')
        ax.set_title(title, fontsize=10, pad=6); ax.set_xlabel('Quarter', fontsize=9)
        ax.set_ylabel('% dev. from SS  (÷ C_ss)', fontsize=9)
        ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=8)
    axes7[0].legend(fontsize=8, frameon=False)
    fig7.suptitle('Figure 7: TPI — Household Utility With vs Without Intervention\n'
                  'Dashed line = no-TPI benchmark (γ = 0)', fontsize=11, y=1.02)
    fig7.tight_layout()
    fig7.savefig(output_dir / 'fig_tpi_utility.png', dpi=150, bbox_inches='tight')
    plt.close(fig7)
    print("  Saved fig_tpi_utility.png")
    print("TPI plots done.")
