# Housekeeping
import dbnomics
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

country_codes = ['ITA', 'DEU', 'ESP', 'GRC', 'PRT']
start_date = pd.to_datetime("2007-01-01")

# --- Fetch data --------------------------------------------------------------

dfs_rates = {
    code: dbnomics.fetch_series(
        'OECD',
        'MEI',
        f'{code}.IRLTLT01.ST.M'
    )
    for code in country_codes
}

# --- Extract and clean -------------------------------------------------------

dfs_rates_extracted = {}

for code, df in dfs_rates.items():
    df = df.copy()
    df['period'] = pd.to_datetime(df['period'])

    df = df.loc[
        df['period'] >= start_date,
        ['period', 'value']
    ].sort_values('period')

    dfs_rates_extracted[code] = df

# Check observation counts
obs_counts = {code: len(df) for code, df in dfs_rates_extracted.items()}
all_equal = len(set(obs_counts.values())) == 1

print(obs_counts, all_equal)

# Coefficient of variation
cvar = {}

for code, df in dfs_rates_extracted.items():
    cvar[code] = df['value'].var() / df['value'].mean()

print(cvar)

# --- Aesthetic settings ------------------------------------------------------

BLUE        = '#002147'
RED         = '#8C1515'
BLUE_MUTED  = '#4a6f8a'
RED_MUTED   = '#c0624a'

_COLORS     = [BLUE, RED, BLUE_MUTED, RED_MUTED]
_LINESTYLES = ['-', '--', '-.', ':']
_MARKERS    = ['', '', '', 'o']

# --- Align data --------------------------------------------------------------

rates_wide = None

for code, df in dfs_rates_extracted.items():
    temp = df.rename(columns={'value': code})

    if rates_wide is None:
        rates_wide = temp
    else:
        rates_wide = rates_wide.merge(temp, on='period', how='outer')

rates_wide = rates_wide.sort_values('period').set_index('period')

# --- Compute spreads relative to Germany ------------------------------------

spread_countries = [code for code in rates_wide.columns if code != 'DEU']

spreads = rates_wide[spread_countries].subtract(
    rates_wide['DEU'],
    axis=0
)

# --- Plot spreads ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(12, 6))

for j, code in enumerate(spread_countries):
    mkr = _MARKERS[j % len(_MARKERS)]

    ax.plot(
        spreads.index,
        spreads[code],
        color=_COLORS[j % len(_COLORS)],
        linestyle=_LINESTYLES[j % len(_LINESTYLES)],
        linewidth=1.8,
        marker=mkr,
        markersize=4,
        markevery=12,
        label=f'{code} - GER'
    )

ax.axhline(0, color='#888888', linewidth=0.8, linestyle=':')

ax.set_title('Long-Term Interest Rate Spreads vs Germany', fontsize=10, pad=6)
ax.set_xlabel('Period', fontsize=9)
ax.set_ylabel('Spread over Germany, pp', fontsize=9)

ax.spines[['top', 'right']].set_visible(False)
ax.tick_params(labelsize=8)
ax.legend(fontsize=8, frameon=False, loc='best')

fig.tight_layout()
plt.show()