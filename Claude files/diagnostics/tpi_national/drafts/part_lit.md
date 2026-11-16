## 4. How the literature gets central-bank purchases to lower spreads

Five mechanisms appear in the literature. Each maps onto one term of the spread decomposition in §3.1, so
each can be checked against this model.

### M1. Balance-sheet capacity: the central bank is not constrained (Gertler–Karadi)

**Papers.** Gertler & Karadi (2011, *JME*) and Gertler & Kiyotaki (2010, *Handbook of Monetary Economics*);
Gertler & Karadi (2013, *IJCB*) for large-scale asset purchases; Dedola, Karadi & Lombardo (2013, *JME*) for
the two-country case.

**Mechanism.** Private intermediaries face the same incentive constraint as here. The central bank is less
efficient at intermediating (a resource cost τ per unit), but it can raise funds elastically by issuing
riskless government debt or reserves that *households* hold, or that do not count against bank capacity. A
purchase therefore *removes* the asset from the constrained balance sheet:

$$
\lambda\mathcal{A}^{\rm private}_t=\lambda\,(Q_tK_t+q_tb_t-q_tM_t)\quad\Longrightarrow\quad
\frac{\partial\mu_t}{\partial M_t}=-(1-\mu_t)\frac{q_t}{\mathcal{A}_t}<0 .
$$

That lowers the excess return on *every* asset the bank holds: the lending spread and the bond's
liquidity premium. This is what Krishnamurthy, Nagel & Vissing-Jorgensen (2018) call the **market-segmentation**
channel, and it accounts for **50%** of the yield declines they attribute to SMP and OMT for Italy, Spain and
Portugal.

**In this model.** It is switched off because Z enters the divertible base at λ_Z = λ. The concrete
implementation is one parameter: **λ_Z = 0**, reserves that cannot be diverted. They sit at the central bank;
Gertler & Kiyotaki's interbank parameter ω plays the same role for claims that are easier to enforce. It does
not reopen the portfolio-substitution margin that the single-λ rule guards against, because the bank does not
choose Z: the Eurosystem does.

### M2. Limited risk-bearing capacity: risk extraction (Vayanos–Vila; Costain–Nuño–Thomas)

**Papers.** Vayanos & Vila (2021, *Econometrica*) for preferred-habitat arbitrage; Gabaix & Maggiori (2015,
*QJE*) for limited intermediary risk capacity; and, for exactly this question, **Costain, Nuño & Thomas
(2025, *Journal of Finance*), "The Term Structure of Interest Rates in a Heterogeneous Monetary Union."**

**Mechanism.** Risk-averse arbitrageurs price bonds with a market price of risk proportional to the risk they
hold. Central-bank purchases reduce the duration and default risk left in private hands, so they compress the
term premium and the **credit risk premium**. Costain, Nuño & Thomas build a two-country monetary union with one
defaultable issuer, solve it in affine form, and find that both the level and the shifts of the Italian spread
are mainly the credit risk premium. They report that the PEPP announcement lowered Italian yields by about
80 bp, with its flexibility adding about 15 bp. Their extension makes the default probability fall as the
central bank holds more of the debt, which they call **default risk extraction** and justify by rollover crises.

**In this model.** Formula (I) is the GK analogue: the "price of risk" is the k-weighted variance. Under
national deposits Ω^d ≲ Ω^nd, so the price of risk is about zero and there is nothing to extract (§2.4). The
concrete implementations are either (a) a funding closure that keeps Ω^d > Ω^nd (D3, the Eurosystem refinancing
of the companion note, or union parity), or (b) a risk budget on bank sovereign exposure (P2 of the companion
note), which adds a genuine risk-bearing limit.

### M3. Default risk that the central bank can change: multiple equilibria and rollover

**Papers.** Calvo (1988); Cole & Kehoe (2000, *REStud*); Lorenzoni & Werning (2019, *AER*), "Slow Moving Debt
Crises"; Corsetti & Dedola (2016, *JEEA*), "The Mystery of the Printing Press"; Camous & Cooper (2019, *AEJ:
Macro*), "'Whatever It Takes' Is All You Need"; Aguiar, Amador, Farhi & Gopinath (2015, *QJE*); Bacchetta,
Perazzi & van Wincoop (2018, *JIE*). The quantitative counterpart is Bocola & Dovis (2019, *AER*).

**Mechanism.** Default is partly self-fulfilling. High yields raise the debt burden and make default more likely
(slow-moving crises), or investors refuse to roll over (rollover crises). A credible backstop removes the bad
equilibrium:
- Corsetti & Dedola show that above a critical scale of central-bank purchases the government strictly prefers
  fiscal adjustment to default when fundamentals are sound;
- Camous & Cooper show that a state-contingent commitment overcomes the strategic uncertainty, while simple rules
  do not;
- Bacchetta et al. show that monetary policy can preclude a slow-moving crisis only through sustained inflation.

**This is the only mechanism in which an announcement lowers the spread at zero purchases.** The OMT was never
used. Krishnamurthy et al. attribute **37%** of the SMP/OMT yield declines to lower default risk and **13%** to
lower redenomination risk. Bocola & Dovis find that non-fundamental risk played a *limited* role in Italy in
2008–12, so the size of this channel is an empirical question, not a given.

**In this model.** It is absent by construction: π_t = p^d(s_t) is exogenous (Bocola's s-shock), so neither
purchases nor the cap can change the expected-loss term, which is most of the spread. The concrete
implementation is to make the default probability depend on the private float, as in Costain–Nuño–Thomas:

$$
\pi_t=\frac{1}{1+\exp\!\big(-[\,s_t+\phi\,\ln((B_{t+1}-M_{t+1})/B_{t+1})\,]\big)},\qquad \phi\ge0 ,
$$

or on the market-value debt burden, as in Lorenzoni–Werning:

$$
\pi_t=p^d\big(s_t+\phi_b\,(\delta B_{t+1}/q_t-\overline{\delta B/q})\big).
$$

φ = 0 nests the current model. φ would be calibrated to the default-risk share of the OMT effect (KNVJ: 37%), or
bounded by Bocola–Dovis.

### M4. Fiscal backstop: loss-sharing and transfers

**Papers.** Roch & Uhlig (2018, *JIE*), "The dynamics of sovereign debt crises and bailouts" (a bailout guarantee, explicitly motivated by OMT, removes the self-fulfilling run); Corsetti & Dedola
(2016) on central-bank losses and inflation.

**Mechanism.** A third party commits resources in default states, which raises recovery for private holders or
removes the incentive to default.

**In this model.** It is present but mis-targeted. The Eurosystem's book takes the haircut in default and the
capital key sends 92.9% of the loss to F, but private holders' recovery is unchanged, the Eurosystem buys nothing
in default, and π does not respond. So the transfer insures the D banks *on the bonds the Eurosystem has already
bought*: it relieves bank net worth, not the market price. Lowering the spread through this channel would require
raising recovery for the bonds still in private hands. Under national funding the capital key also works
*against* the backstop in no-default states: the book's excess return is a transfer from D to F that D cannot
smooth (§2.5, §5.2 (iii)).

### M5. Redenomination risk

**Paper.** Krishnamurthy, Nagel & Vissing-Jorgensen (2018), who identify redenomination from the gap between
local-law euro bonds and foreign-law or dollar bonds: 13% of the effect on average, and none for Italy.

**In this model.** Absent: there are no currencies and no exit option. It belongs to a nominal extension.

### What the evidence says about the size of a purchase's price impact

| Source | Measured impact |
|---|---|
| Eser & Schwaab (2016, *JFE*), SMP | about **−3 bp** on 5-year yields per purchase of 0.1% of the outstanding stock, i.e. about −30 bp per 1% of the stock, plus large announcement effects |
| Ghysels, Idier, Manganelli & Vergote (2017, *JEEA*) | SMP lowered program-country yields at high frequency (daily regressions understate it, because purchases responded to price falls) |
| Krishnamurthy, Nagel & Vissing-Jorgensen (2018) | SMP + OMT: 2-year yields −400 bp (Italy, Spain), −500 bp (Ireland, Portugal), −1000 bp (Greece); 37% default risk / 13% redenomination / 50% segmentation |
| Costain, Nuño & Thomas (2025) | PEPP announcement about −80 bp on Italian yields, +15 bp from flexibility; mostly credit risk premium |

These set the target for the model. A credible TPI should deliver, at the stress node, a price impact of the
order of **tens of bp per 1% of the stock**, through segmentation and default-risk extraction.
