# ONE-LINE SUMMARIES OF SAVED TPI RUNS, FOR THE REPORT TABLES.
#   python3 summarise.py tag1 tag2 ...
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent / "tpi_runs"


def row(tag):
    # THE NUMBERS THE REPORT QUOTES FOR ONE RUN.
    d = json.loads((HERE / f"{tag}.json").read_text())
    r0, r1 = d["rest_off"], d["rest_on"]
    off, on = d["irf_off"], d["irf_on"]
    imp = {k: {r["m"]: r["sov_bp"] for r in d[k]} for k in ("impact_rest", "impact_shock")}
    kern = d.get("kernel", {})
    cum = lambda p, k: sum(p[k])
    trough = lambda p: (min(p["Y"]), p["Y"].index(min(p["Y"])))
    return dict(
        tag=tag, cap=d["cap_solved"], ok=d["solve_ok"],
        rest_sov=(r0["sov_bp"], r1["sov_bp"]), rest_lend=(r0["spread_bp"], r1["spread_bp"]),
        rest_Y=100 * (r1["Y"] / r0["Y"] - 1), rest_I=100 * (r1["I"] / r0["I"] - 1),
        rest_C=100 * (r1["C"] / r0["C"] - 1),
        shock_sov=(off["sov_bp"][0], on["sov_bp"][0]), bought=on["m_cb"][0],
        Y0=(off["Y"][0], on["Y"][0]), trough=(trough(off), trough(on)),
        cumY=(cum(off, "Y"), cum(on, "Y")), cumC=(cum(off, "C"), cum(on, "C")),
        cumI=(cum(off, "I"), cum(on, "I")), n0=(off["n"][0], on["n"][0]),
        d_spread0=(off["d_spread"][0], on["d_spread"][0]), d_rdep0=(off["d_rdep"][0], on["d_rdep"][0]),
        Mpath=on["M_cb"][:8], Pi=on["Pi_cb"][:8],
        impact_rest=imp["impact_rest"], impact_shock=imp["impact_shock"],
        legs={k: (v["exp_loss_bp"], v["risk_prem_bp"], v["liq_prem_bp"], v["Om_ratio"], v["m_cb"])
              for k, v in kern.items()},
        kt=d["kt"])


def main():
    # PRINT EVERY RUN NAMED.
    for tag in sys.argv[1:]:
        try:
            r = row(tag)
        except FileNotFoundError:
            print(f"{tag}: not finished")
            continue
        f = lambda t: f"{t[0]:+.1f} -> {t[1]:+.1f}"
        print(f"\n== {tag}: cap reached {r['cap']}  ok={r['ok']}")
        print(f"   rest: sov {f(r['rest_sov'])}  lending {f(r['rest_lend'])}  dY {r['rest_Y']:+.2f}%  "
              f"dI {r['rest_I']:+.2f}%  dC {r['rest_C']:+.2f}%")
        print(f"   shock: sov on impact {f(r['shock_sov'])}  bought {r['bought']:.1f}% of stock  "
              f"Y0 {r['Y0'][0]:+.3f} / {r['Y0'][1]:+.3f}  trough {r['trough'][0][0]:+.3f}@q{r['trough'][0][1]}"
              f" / {r['trough'][1][0]:+.3f}@q{r['trough'][1][1]}")
        print(f"          cumY {r['cumY'][0]:+.2f} / {r['cumY'][1]:+.2f}  cumC {r['cumC'][0]:+.2f} / "
              f"{r['cumC'][1]:+.2f}  cumI {r['cumI'][0]:+.2f} / {r['cumI'][1]:+.2f}  n0 {r['n0'][0]:+.2f} / "
              f"{r['n0'][1]:+.2f}  d_spread0 {r['d_spread0'][0]:+.0f} / {r['d_spread0'][1]:+.0f}")
        print(f"          book M (% stock) " + " ".join(f"{x:.1f}" for x in r["Mpath"])
              + "   P&L (% qtr GDP) " + " ".join(f"{x:+.2f}" for x in r["Pi"]))
        for k in ("impact_rest", "impact_shock"):
            base = r[k][0.0]
            print(f"   {k}: " + "  ".join(f"{100*m:.0f}%: {v - base:+.1f}bp" for m, v in sorted(r[k].items())
                                         if m > 0))
        for k, v in r["legs"].items():
            print(f"   legs {k:10s} EL {v[0]:6.1f}  RP {v[1]:+6.1f}  LP {v[2]:6.1f}  Om_d/Om_nd {v[3]:.3f}"
                  f"  m_cb {v[4]:.3f}")
        print(f"   KT: binding {r['kt']['n_binding']}/{r['kt']['n']}  corner {r['kt']['n_corner']}  "
              f"min b_DD/B {r['kt']['min_b']:.3f}")


if __name__ == "__main__":
    main()
