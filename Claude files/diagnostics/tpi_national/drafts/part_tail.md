## 6. Caveats and what to run next

- **Quick grid.** Every run here uses the 25-point coarse grid. Signs and cross-economy comparisons are
  robust to that, since every economy is solved the same way. Levels are not: the coarse union announcement
  (−49 bp sovereign, −95 bp lending) is several times the full-grid one (−4.3 / −36 bp). The next step is the
  full grid (`S_REFINE = 5`) for `nat`, `nat_phi1_lz0_key1` and `uni_phi1_lz0`.
- **The box.** In the φ economies the Eurosystem's impact purchase (79% of the stock) sits at the edge of
  the TPI book's box (±80%), and every later quarter of those impulse responses is clipped back into the box
  (the run logs print the escapes). Only their impact quarter is clean. Widen `m_band` above 1 before
  reading the dynamics.
- **φ is illustrative.** φ = 1 was chosen to show the channel, not calibrated. §5.4 lists three calibration
  targets.
- **κ_D = 1 is a policy choice, not a free lunch.** It removes the no-default drain to F, but also the
  cross-border loss-sharing in default. The impulse responses condition on no default, so they show only the
  benefit. A welfare comparison needs the default branch.
- **The lending-spread jump in the full-fix economy** (+267 bp on impact on exact reads) comes with a smaller
  fall in the bond price and a rise in bank net worth and investment. The franchise-value reading in §2.5 is an
  interpretation, not a decomposition. Run the output decomposition on that economy before writing it up.
- **Pre-existing open items** that these runs inherit, unrelated to the TPI:
  - impact D consumption rises under national deposits (+0.19% here), the C_D gate;
  - the coarse union TPI impulse response is pathological (impact output −1.8%), as already recorded on
    2026-10-02 to 05.
- **Nothing is committed.** The worktree patches (§1.2) live in the session scratchpad. Porting them means
  adding `lambda_Z_D` and `phi_float` to `calibration/global_projection.py` with their nesting defaults, and
  N-tests as in §5.4.

## 7. References

- Bacchetta, P., E. Perazzi and E. van Wincoop (2018). Self-fulfilling debt crises: What can monetary policy do? *Journal of International Economics* 110.
- Bocola, L. (2016). The pass-through of sovereign risk. *Journal of Political Economy* 124(4).
- Bocola, L. and A. Dovis (2019). Self-fulfilling debt crises: A quantitative analysis. *American Economic Review* 109(12).
- Calvo, G. (1988). Servicing the public debt: The role of expectations. *American Economic Review* 78(4).
- Camous, A. and R. Cooper (2019). "Whatever it takes" is all you need: Monetary policy and debt fragility. *AEJ: Macroeconomics* 11(4).
- Cole, H. and T. Kehoe (2000). Self-fulfilling debt crises. *Review of Economic Studies* 67(1).
- Corsetti, G. and L. Dedola (2016). The mystery of the printing press: Monetary policy and self-fulfilling debt crises. *Journal of the European Economic Association* 14(6).
- Costain, J., G. Nuño and C. Thomas (2025). The term structure of interest rates in a heterogeneous monetary union. *Journal of Finance* 80(4).
- Dedola, L., P. Karadi and G. Lombardo (2013). Global implications of national unconventional policies. *Journal of Monetary Economics* 60(1).
- Eggertsson, G. and M. Woodford (2003). The zero bound on interest rates and optimal monetary policy. *Brookings Papers on Economic Activity* 2003(1).
- Eser, F. and B. Schwaab (2016). Evaluating the impact of unconventional monetary policy measures: Empirical evidence from the ECB's Securities Markets Programme. *Journal of Financial Economics*.
- Gabaix, X. and M. Maggiori (2015). International liquidity and exchange rate dynamics. *Quarterly Journal of Economics* 130(3).
- Gertler, M. and P. Karadi (2011). A model of unconventional monetary policy. *Journal of Monetary Economics* 58(1).
- Gertler, M. and P. Karadi (2013). QE 1 vs. 2 vs. 3...: A framework for analyzing large-scale asset purchases as a monetary policy tool. *International Journal of Central Banking*.
- Gertler, M. and N. Kiyotaki (2010). Financial intermediation and credit policy in business cycle analysis. *Handbook of Monetary Economics* 3A.
- Ghysels, E., J. Idier, S. Manganelli and O. Vergote (2017). A high-frequency assessment of the ECB Securities Markets Programme. *Journal of the European Economic Association*.
- Krishnamurthy, A., S. Nagel and A. Vissing-Jorgensen (2018). ECB policies involving government bond purchases: Impact and channels. *Review of Finance* 22(1).
- Lorenzoni, G. and I. Werning (2019). Slow moving debt crises. *American Economic Review* 109(9).
- Roch, F. and H. Uhlig (2018). The dynamics of sovereign debt crises and bailouts. *Journal of International Economics* 114.
- Vayanos, D. and J.-L. Vila (2021). A preferred-habitat model of the term structure of interest rates. *Econometrica* 89(1).
- Wallace, N. (1981). A Modigliani–Miller theorem for open-market operations. *American Economic Review* 71(3).
- Aguiar, M., M. Amador, E. Farhi and G. Gopinath (2015). Coordination and crisis in monetary unions. *Quarterly Journal of Economics* 130(4).

**Web sources consulted** (publication details and headline findings):
[Costain–Nuño–Thomas, Banco de España summary](https://www.bde.es/wbe/en/areas-actuacion/analisis-e-investigacion/investigacion-destacada/research-features/term-structure-interest-rates-heterogeneous-monetary-union.html) ·
[Camous–Cooper, AEJ: Macro](https://www.aeaweb.org/doi/10.1257/mac.20170167) ·
[Krishnamurthy–Nagel–Vissing-Jorgensen, NBER w23985](https://www.nber.org/papers/w23985) ·
[NBER Digest summary of KNVJ](https://www.nber.org/digest/feb18/quantifying-impact-ecb-policies-during-debt-crisis) ·
[Bocola–Dovis, AER](https://www.aeaweb.org/doi/10.1257/aer.20161471) ·
[Bacchetta–Perazzi–van Wincoop, NBER w21158](https://www.nber.org/papers/w21158) ·
[Eser–Schwaab, EconStor](https://econstor.eu/handle/10419/154075) ·
[Corsetti–Dedola, IDEAS](https://ideas.repec.org/a/bla/jeurec/v14y2016i6p1329-1371.html) ·
[Lorenzoni–Werning, AER](https://www.aeaweb.org/doi/10.1257/aer.20141766) ·
[Aguiar–Amador–Farhi–Gopinath, NBER w20277](https://www.nber.org/papers/20277) ·
[Roch–Uhlig, IDEAS](https://ideas.repec.org:443/p/imf/imfwpa/2016-136.html) ·
[Bocola (2016) and the LTRO, Minneapolis Fed](https://minneapolisfed.org/article/2015/risky-business) ·
[Gertler–Karadi (2011), IDEAS](https://ideas.repec.org/a/eee/moneco/v58y2011i1p35-38.html) ·
[DNB WP 778, fragmentation](https://www.dnb.nl/en/publications/research-publications/working-paper-2023/778-identifying-financial-fragmentation-do-sovereign-spreads-in-the-emu-reflect-differences-in-fundamentals/)
