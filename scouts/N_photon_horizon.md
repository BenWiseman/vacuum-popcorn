# N. UHE photon attenuation, 1e19 to 1e20 eV (scout, 2026-10-01)

**Result.** The claim holds except for one gap, Galactic radio emission (last section). No published central estimate I found gives a pair-production interaction length below 1.6 Mpc from 1e19 to 1e20 eV, CMB plus radio. A 100 kpc path then absorbs at most 6.1% of photons, an 8.5 kpc path 0.5% (my arithmetic). `N_photon_check.py` reproduces four published CMB lengths to 1.5%.

**1. Interaction length at 1e19, 5e19, 1e20 eV, in Mpc** (distance for 1/e survival). "Plot" means measured off the figure.
- CMB only: 2.5, 10, 19 (Heiter et al. 1710.11406 Fig. 5 and Protheroe and Biermann 1996, PB96, Fig. 6, plots; my calculation agrees).
- CMB + radio, Auger 1612.04155 Fig. 1: 2.4, 6.4, 7.2; Heiter et al. Fig. 5: 2.4, 6.7, 7.7 (plots). Radio: PB96, case not stated.
- Same, Bhattacharjee and Sigl astro-ph/9811011 Fig. 11, dashed (plot): ~2, 3.2, 3.4. Radio: Clark et al. 1970, observed.
- PB96 astro-ph/9605119 Fig. 6, radio only (plot): 19, 7.8, 6.5 without galaxy evolution; 9.5, 3.5, 2.3 with it. With CMB added: 2.2, 4.5, 4.9; 2.0, 2.6, 2.1.
- Niţu et al. 2004.13596 Fig. 15, radio only (plot): 7.6, 2.5, 1.7; their Table 2 fit, my calculation: 7.7, 2.5, 1.8. With CMB added: 1.9, 2.0, 1.6; lower edge of their uncertainty band at 1e20 eV, 1.3.

Energy-loss lengths (cascade included) run longer. Risse and Homola (astro-ph/0702632): "Typical energy loss lengths assumed for UHE photons range between 7–15 Mpc at 10^19 eV and 5–30 Mpc at 10^20 eV." Auger (2509.05113): "approximately 10 Mpc at 10^19 eV". Gelmini, Kalashev, Semikoz (astro-ph/0506128), photons near 5e19 eV: "cannot reach Earth from beyond 10 to 40 Mpc".

**2. Shortest length**
- PB96: "interactions in the cosmic microwave background radiation give a mean interaction length of less than 10 kpc at 10^6 GeV" (1e15 eV).
- Vernetto and Lipari (1608.01587): "most effective for absorption at gamma ray energy Eγ ≃ 2.2 PeV, where the absorption length is 7.07 kpc"; "galactic sources beyond the Galactic center are very strongly attenuated".
- Auger (2509.05113): "about 30 kpc at 10^15 eV"; Risse and Homola Fig. 1 (plot) has the same minimum, 30 kpc, as an energy-loss length.
- Auger 1612.04155: "the attenuation length varies between 90 kpc at 10^17.3 eV and 900 kpc at 10^18.5 eV". Galactic attenuation matters from about 1e14 to 1e18 eV, not above.

**3. Statements on the Galactic halo**
- Kalashev and Kuznetsov (1606.07354, Sec. 2): "for photons with E ≳ 10^18 eV the attenuation length in interstellar medium exceeds the size of our Galaxy halo. This implies that for higher energy photon we can neglect the absorption and cascaded radiation." (their cascade check: CMB and EBL.)
- Bhattacharjee and Sigl (Sec. 6.13.1; Eq. 97 scales R_H to 100 kpc): "R_H is too small compared to EHE photon interaction length of ≳ 10 Mpc (see Fig. 11)". No energy named; their Fig. 11 interaction length (plot) is 3.4 Mpc at 1e20 eV, 12 Mpc at 1e21 eV.

**Gaps**
- NOT FOUND: any attenuation length in TA photon papers 1811.03920 and 2512.01638.
- NOT FOUND: a published length with the Galaxy's own radio emission at 1e19 to 1e20 eV. Every curve above uses the CMB, infrared and extragalactic radio only. Next: compute with a Galactic radio sky map.
- NOT FOUND: the 1997 PB96 erratum; PB96 values here use the 1996 arXiv text.
- My calculation, unpublished: an ARCADE-2 type power law (1.26 K at 1 GHz, index -2.6, values from 1110.5257, cut at 1 MHz) gives 1.2 Mpc at 1e20 eV.
