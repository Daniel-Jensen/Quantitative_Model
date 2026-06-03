import pandas as pd
import matplotlib.pyplot as plt

start = "1995-01-01"

rates = {
    "ITA": "IRLTLT01ITM156N",
    "DEU": "IRLTLT01DEM156N",
    "ESP": "IRLTLT01ESM156N",
    "GRC": "IRLTLT01GRM156N",
    "PRT": "IRLTLT01PTM156N",
}

recs = {
    "ITA": "ITARECM",
    "ESP": "ESPRECM",
    "GRC": "GRCRECM",
    "PRT": "PRTRECM",
}

cols = {
    "ITA": "#0072B2",
    "ESP": "#E69F00",
    "GRC": "#D55E00",
    "PRT": "#009E73",
}

styles = {
    "ITA": "-",
    "ESP": "-",
    "GRC": "-",
    "PRT": "-",
}

events = {
    "ECB created": "1999-01-01",
    "OMT": "2012-09-06",
    "TPI": "2022-07-21",
}


def fred(series):
    return pd.concat(
        [
            pd.read_csv(
                f"https://fred.stlouisfed.org/graph/fredgraph.csv?id={sid}",
                parse_dates=["observation_date"],
            )
            .rename(columns={"observation_date": "date", sid: c})
            .assign(**{c: lambda x, c=c: pd.to_numeric(x[c], errors="coerce")})
            .loc[lambda x: x["date"] >= start, ["date", c]]
            .set_index("date")
            for c, sid in series.items()
        ],
        axis=1,
    )


plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.dpi": 150,
})


r = fred(rates)
spreads = r.drop(columns="DEU").sub(r["DEU"], axis=0)

all_rec = fred(recs).reindex(spreads.index).ffill().eq(1).all(axis=1)


def shade(ax, x):
    for _, b in x.groupby(x.ne(x.shift()).cumsum()):
        if b.iloc[0]:
            ax.axvspan(
                b.index[0],
                b.index[-1],
                color="#D9D9D9",
                alpha=0.6,
                lw=0,
                zorder=0,
                label="All four in recession",
            )


def plot_spreads(title, ymax=None):
    fig, ax = plt.subplots(figsize=(10.5, 5.6))

    shade(ax, all_rec)

    for c in spreads:
        ax.plot(
            spreads.index,
            spreads[c],
            color=cols[c],
            ls=styles[c],
            lw=2.4,
            label=f"{c} - Germany",
            zorder=2,
        )

    ax.axhline(0, color="#888888", lw=1.2, ls=":", zorder=1)

    if ymax:
        ax.set_ylim(min(0, spreads.min().min()), ymax)

    y0, y1 = ax.get_ylim()

    for label, date in events.items():
        date = pd.to_datetime(date)
        ax.axvline(date, color="#888888", lw=1.1, ls=":", alpha=0.9, zorder=1)
        ax.text(
            date,
            y1 - 0.04 * (y1 - y0),
            f" {label}",
            rotation=90,
            va="top",
            ha="left",
            fontsize=9,
            color="#555555",
        )

    ax.set(
        title=title,
        xlabel="Year",
        ylabel="Spread over Germany, percentage points",
        xlim=(spreads.index.min(), spreads.index.max()),
    )

    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", length=4, width=1)
    ax.grid(axis="y", color="#DDDDDD", lw=0.6, alpha=0.6)
    ax.grid(axis="x", visible=False)

    h, l = ax.get_legend_handles_labels()

    seen = set()
    h2, l2 = [], []

    for handle, label in zip(h, l):
        if label not in seen:
            h2.append(handle)
            l2.append(label)
            seen.add(label)

    if "All four in recession" in l2:
        i = l2.index("All four in recession")
        h2, l2 = h2[:i] + h2[i + 1:] + [h2[i]], l2[:i] + l2[i + 1:] + [l2[i]]

    ax.legend(
        h2,
        l2,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.015, 0.985),
        handlelength=2.8,
        borderaxespad=0,
    )

    fig.tight_layout()
    return fig, ax


plot_spreads("Long-Term Interest Rate Spreads vs Germany")
plot_spreads("Long-Term Interest Rate Spreads vs Germany, Y-Axis Capped at 5", ymax=5)

plt.show()