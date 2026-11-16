## 5. What to take into the paper, and what to run next

**The faithful TPI, as a specification:**

1. **Funding through a TARGET2 margin.** Use the union deposit market, interpreted as Eurosystem
   intermediation of cross-border funding (§2.1). The TPI's results should not be reported on the national
   market: there the deposit rate insures the banker and the backstop has nothing to act on.
2. **A rollover component of default risk** à la Costain–Nuño–Thomas: π_t = p^d(s_t) − θ_q·M_{t+1}/B_{t+1}, with
   θ_q = 0.0046 mapped from their θ = 2.02×10⁻⁵ per € billion (§2.2). State in the text that p^d(s) is the
   fundamental part and θ_q the backstop-sensitive part.
3. **Pari passu**, as the TPI and OMT texts commit (§2.4).
4. **A fundamentals-based activation rule**: cap = fundamental (expected-loss) spread + 50 bp, the model
   counterpart of "unwarranted, disorderly" (§2.5). Drop the fixed 200 bp cap: below the expected-loss spread
   it is a transfer, not a backstop.
5. **P&L by capital key** as the baseline, PSPP-style national income retention as the sensitivity, and a
   statement that the TPI's actual sharing rule is not public (§2.3).
6. **Reserves as divertible as bonds** (λ_Z = λ), as the regulatory treatment and Diamond–Jiang–Ma suggest,
   with λ_Z = 0 as the Gertler–Karadi bound (§2.6).

**What the faithful TPI does in the model (§3):** it caps the spread at the fundamental level plus 50 bp,
lowers the risk premium (on the quick grid it extracts all of it), and in a realized default becomes a large
cross-border transfer from F to D.

**What it costs:** it needs a very large book, about 80% of the stock on the quick grid, because private
demand is too elastic. On the full grid the headline shock barely reaches the warranted cap, so the book stays
small; the slope, −1.2 bp per 1% of the stock, is still 12–25 times below the evidence.

**Next steps, in order:**

1. **The elasticity gap (§4).** Add CNT's two-type demand: inelastic preferred-habitat holders plus
   mean-variance arbitrageurs, with γ calibrated to the Eser–Schwaab slope, and the 2011–12 Greek spread level
   as a joint target. This is the single most important remaining step toward the real-world magnitudes.
2. **PSPP-style retention.** Narrow the TPI obligation band (`o_band`, currently ±25% of the book's half-width)
   so that the corners with ±20%-of-GDP P&L disappear, then re-run Q3b.
3. **The faithful design at 9 nodes in s.** At 115 points its announcement effect and dynamics are not
   identified (§3.4). Only the baseline (F1c) has been run on the full grid; the other variants have not.
4. **Implement D3** (the explicit Eurosystem facility with a collateral constraint and the ELA penalty) as the
   structural version of the TARGET2 margin.
5. **Welfare.** Weigh the no-default path (D pays the insurance premium) and the default path (F pays the
   claim) with the default probabilities. The two paths in §3.3 are the inputs.

---

## 6. Files

| File | Contents |
|---|---|
| `Claude files/docs/tpi_fix_explained.md` | this report |
| `Claude files/docs/tpi_national_analysis.md` | the previous note: why the shipped TPI fails under national deposits |
| `Claude files/diagnostics/tpi_national/tpi_switches_v3.patch` | every experimental switch as one patch against `domestic-funding-fix` (`0dd0732`): the national-TPI unlock, `lambda_Z_D`, `phi_float`, `theta_cnt`, the warranted cap, and the income/loss key split |
| `Claude files/diagnostics/tpi_national/tpi_national.py` | solve + compare driver (`--theta`, `--warranted`, `--keyD`, `--keyLoss`, `--lamZ`, `--market`, `--s-refine`) |
| `Claude files/diagnostics/tpi_national/default_path.py` | the realized-default paths of §3.3 |
| `Claude files/diagnostics/tpi_national/summarise.py` | the table numbers |
| `results/TPI_NATIONAL/` | every run's JSON and log (git-ignored) |

---

## 7. References

- Altavilla, C., D. Giannone and M. Lenza (2016). The financial and macroeconomic effects of OMT announcements. *International Journal of Central Banking* 12(3).
- Bocola, L. (2016). The pass-through of sovereign risk. *Journal of Political Economy* 124(4).
- Bocola, L. and A. Dovis (2019). Self-fulfilling debt crises: A quantitative analysis. *American Economic Review* 109(12).
- Calvo, G. (1988). Servicing the public debt: The role of expectations. *American Economic Review* 78(4).
- Cole, H. and T. Kehoe (2000). Self-fulfilling debt crises. *Review of Economic Studies* 67(1).
- Corsetti, G. and L. Dedola (2016). The mystery of the printing press. *Journal of the European Economic Association* 14(6).
- Costain, J., G. Nuño and C. Thomas (2025). The term structure of interest rates in a heterogeneous monetary union. *Journal of Finance* 80(4). Working-paper version, November 2021, used for the equations and θ.
- Diamond, W., Z. Jiang and Y. Ma (2024). The reserve supply channel of unconventional monetary policy. *Journal of Financial Economics* 159.
- Eser, F. and B. Schwaab (2016). Evaluating the impact of unconventional monetary policy measures: Empirical evidence from the ECB's Securities Markets Programme. *Journal of Financial Economics*.
- Gabaix, X. and M. Maggiori (2015). International liquidity and exchange rate dynamics. *Quarterly Journal of Economics* 130(3).
- Gertler, M. and P. Karadi (2011). A model of unconventional monetary policy. *Journal of Monetary Economics* 58(1).
- Krishnamurthy, A., S. Nagel and A. Vissing-Jorgensen (2018). ECB policies involving government bond purchases: Impact and channels. *Review of Finance* 22(1).
- Lorenzoni, G. and I. Werning (2019). Slow moving debt crises. *American Economic Review* 109(9).
- Vayanos, D. and J.-L. Vila (2021). A preferred-habitat model of the term structure of interest rates. *Econometrica* 89(1).

**Institutional and data sources consulted:**
- [ECB press release on the TPI, 21 July 2022](https://www.ecb.europa.eu/press/pr/date/2022/html/ecb.pr220721%7e973e6e7273.en.html) · [Bundesbank TPI page](https://bundesbank.de/en/tasks/monetary-policy/outright-transactions/transmission-protection-instrument-tpi--896158)
- [ECB, technical features of OMT, 6 Sep 2012](https://www.ecb.europa.eu/press/pr/date/2012/html/pr120906_1.en.html)
- [ECB FAQ on the PSPP](https://www.ecb.europa.eu/mopo/implement/app/html/ecb.faq_pspp.en.html) · [ECB FAQ on the PEPP](https://www.ecb.europa.eu/mopo/implement/pepp/html/ecb.faq_pepp.en.html) · [Eurosystem monetary income (Springer, 2025)](https://link.springer.com/article/10.1007/s11079-025-09804-w)
- [Eurogroup and the SMP income transfers to Greece: ECB correspondence](https://www.ecb.europa.eu/ecb/access_to_documents/document/correspondence/shared/data/ecb.dr.cor20171010Chountis.fi.pdf) · [Private sector involvement (Greek PSI, SMP exemption)](https://en.wikipedia.org/wiki/Private_sector_involvement) · [Xafa, "Life after debt"](https://mirandaxafa.com/mediaupload/articles/life_after_debt.pdf)
- [Greek ELA and Eurosystem funding: Yale case study](https://elischolar.library.yale.edu/journal-of-financial-crises/vol4/iss2/39/) · [TARGET2 balances 2011: CESifo Forum](https://www.cesifo.org/DocDL/forum-0112-special-5.pdf) · [Greek deposit rates, April 2011 (Bank of Greece data)](https://neoskosmos.com/en/2011/06/15/news/greece/greek-banks-deposits-down-2-4-billion-euros-april/)
- [Costain–Nuño–Thomas, BIS WP 1165](https://www.bis.org/publ/work1165.pdf) · [Diamond–Jiang–Ma, NBER w31693](https://www.nber.org/papers/w31693) · [Altavilla–Giannone–Lenza, ECB WP 1707](https://ECB.eu/pub/pdf/scpwps/ecbwp1707.pdf) · [Krishnamurthy–Nagel–Vissing-Jorgensen, NBER w23985](https://www.nber.org/papers/w23985) · [Eser–Schwaab, EconStor](https://econstor.eu/handle/10419/154075)
