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


# Reverted back to the original Serif font design
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
        if not b.empty and b.iloc[0]:
            ax.axvspan(
                b.index[0],
                b.index[-1],
                color="#D9D9D9",
                alpha=0.6,
                lw=0,
                zorder=0,
                label="All four in recession",
            )


def plot_spreads(title, start_date, end_date):
    # Added facecolor='none' to make the figure background transparent
    fig, ax = plt.subplots(figsize=(10.5, 6), facecolor="none")
    
    # Make the plot area background transparent as well
    ax.set_facecolor("none")

    s_date = pd.to_datetime(start_date)
    e_date = pd.to_datetime(end_date)

    mask = (spreads.index >= s_date) & (spreads.index <= e_date)
    sub_spreads = spreads.loc[mask]
    sub_rec = all_rec.loc[mask]

    shade(ax, sub_rec)

    for c in sub_spreads:
        ax.plot(
            sub_spreads.index,
            sub_spreads[c],
            color=cols[c],
            ls=styles[c],
            lw=2.4,
            label=f"{c}",
            zorder=2,
        )

    ax.axhline(0, color="#888888", lw=1.2, ls=":", zorder=1)

    y0, y1 = ax.get_ylim()

    for label, date in events.items():
        event_date = pd.to_datetime(date)
        if s_date <= event_date <= e_date:
            ax.axvline(event_date, color="#888888", lw=2.2, ls=":", alpha=0.9, zorder=1)
            
            # Specific rule for TPI: place it strictly above the graph
            if label == "TPI":
                ax.annotate(
                    f"{label}",
                    xy=(event_date, y1),
                    xytext=(0, 6),
                    textcoords="offset points",
                    rotation=90,
                    va="bottom",
                    ha="center",
                    fontsize=10,
                    color="#555555",
                    clip_on=False
                )
            # All other events: place them inside the graph, offset to the right
            else:
                ax.annotate(
                    f" {label}",
                    xy=(event_date, y1 - 0.04 * (y1 - y0)),
                    xytext=(6, 0),
                    textcoords="offset points",
                    rotation=90,
                    va="top",
                    ha="left",
                    fontsize=10,
                    color="#555555",
                )

    ax.set(
        xlabel="Year",
        ylabel="Spread, percentage points", 
        xlim=(sub_spreads.index.min(), sub_spreads.index.max()),
    )
    
    ax.set_title(title, pad=35)

    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", length=4, width=1)
    ax.grid(axis="y", color="#DDDDDD", lw=0.6, alpha=0.6)
    ax.grid(axis="x", visible=False)

    h, l = ax.get_legend_handles_labels()

    seen = set()
    h_countries, l_countries = [], []
    h_rec, l_rec = [], []

    for handle, label in zip(h, l):
        if label not in seen:
            seen.add(label)
            if label == "All four in recession":
                h_rec.append(handle)
                l_rec.append(label)
            else:
                h_countries.append(handle)
                l_countries.append(label)

    # Main Legend (Countries)
    leg1 = ax.legend(
        h_countries,
        l_countries,
        frameon=False,
        loc="upper right",
        handlelength=2.8,
        borderaxespad=1,
    )
    ax.add_artist(leg1) 

    # Secondary Legend (Recession)
    if h_rec:  
        ax.legend(
            h_rec,
            l_rec,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.14),
            handlelength=2.8,
            borderaxespad=0,
        )

    fig.subplots_adjust(top=0.88, bottom=0.2)
    return fig, ax


# Graph 1: 1995 to 2020
plot_spreads("Long-Term Interest Rate Spreads (1995 - 2020)", "1995-01-01", "2020-01-01")

# Graph 2: 2020 to Now
plot_spreads("Long-Term Interest Rate Spreads (2020 - Present)", "2020-01-01", pd.Timestamp.today())

plt.show()