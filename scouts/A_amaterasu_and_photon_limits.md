# Scout A: Amaterasu and UHE photon limits (2026-09-30)

**Summary.** TA prints E = 244 ± 29 (stat.) +51/-76 (syst.) EeV at RA 255.9°, Dec 16.1° (arXiv:2311.14231), excludes a photon primary at 99.986% CL, and finds no candidate source. Later papers restate the energy as 212 EeV (iron), 154 ± 18 EeV (Auger scale) or about 170 to 175 EeV (Yakutsk-side calibration); I found no TA revision. Two papers propose SHDM decay: 2406.03174 finds severe tension with TA and Auger photon limits, and 2504.15272 reaches Amaterasu only at 2σ. The strongest printed photon-based lifetime limit is τ ≳ 3 to 4 × 10^30 s (b b̄, 10^12.2 to 10^13.8 GeV; 2302.02993, 95% CL). Auger SD photon flux limits above 10, 20, 40 EeV are 2.11, 0.312, 0.172 × 10^-3 km^-2 sr^-1 yr^-1. Neutrino-based limits at these masses are much weaker (10^25 to 10^27 s); an IceCube or KM3NeT limit at 10^11 to 10^13 GeV is NOT FOUND.

## 1. The TA event (2311.14231, Science 382, 903)

- Energy: "244 ± 29 (stat.) +51 -76 (syst.) EeV in the detector frame" (stacked in print; the abstract gives the same in exa-electron volts, "~40 joules"). "The original TA SD reconstructed energy of 309±37(stat.) EeV". Systematics named: "-10% in the direction of lower energies" for unknown primary; migration "-3%"; SD "systematic uncertainty of 21%".
- Direction (Table 1): 27 May 2021 10:35:56 UTC; RA 255.9±0.6°, Dec 16.1±0.5°; zenith 38.6±0.4°; azimuth 206.8±0.6°; S800 530±57 m^-2. The arXiv text (one version) has no Galactic coordinates; the Science page returned 403, so its supplement is unread. Printed elsewhere: 2312.13273 "(ℓ, b) = (36.2, 30.9)◦ in Galactic coordinates"; 2501.16158 "(l, b) = (36.2◦, 30.9◦)". My conversion of the printed RA/Dec gives (36.15°, 30.94°).
- Composition: "excludes a photon as the primary particle at the 99.986% confidence level, instead favoring a primary proton. However, the classifier cannot distinguish between protons and heavier nuclei because the fluorescence detectors were not operating during this event due to bright moon light."
- Source: "The arrival direction of this event is consistent with the location of the Local Void ... Even taking into account the range of possible GMF deflections and primary mass, we do not identify any candidate sources for this event." Intro: BSM production "is constrained by upper bounds on the flux of ultra-high-energy (UHE) photons (5, 6)" (TA Astropart. Phys. 110, 8; Auger ApJ 933, 125).
- Energy re-estimates (third party; none is a TA re-analysis):
  - 2312.13273 (Unger, Farrar): "E_nom = (2.12±0.25) × 10^20 eV, including the quoted corrections for resolution effects (-3%) and heavy primaries (-10%)"; "E_low = (1.64±0.19) × 10^20 eV" with the 20% scale uncertainty.
  - 2501.16158 (Korochkin, Semikoz, Tinyakov): "rescaled to E = 212 EeV to account for the systematic bias of the energy reconstruction for the heavy primaries".
  - 2502.15876 (Unger, Farrar): "ΔE_TA/E_TA = -0.09 - 0.2 (lg(E_TA/eV) - 19)" (their ref: EPJ Web Conf. 283 (2023) 2003); "The correction amounts to -36.7% and the energy of the Amaterasu particle at the Auger energy scale is 154±18 (stat.) EeV". I checked the arithmetic. The newer working-group proceeding 2509.05530 prints "an energy shift of 16-20% per decade" above 10^19 eV and does not mention Amaterasu.
  - 2404.16948 (Glushkov et al., CORSIKA + QGSJet-II.04 protons): TA's nomogram "10^20.36 ≈ 2.29×10^20 eV"; their calibration "10^20.22 ≈ 1.75 × 10^20 eV"; vertical conversion "10^20.23 ≈ 1.70 × 10^20 eV".
  - 2405.12004 (Farrar, PRL 134, 081003): "Amaterasu at 210-250 EeV".
  - TA Collaboration re-analysis with another hadronic model: NOT FOUND.

## 2. Exotic-origin proposals

- 2406.03174 (Sarmah et al., PRD 111, 083048), SHDM decay in the Galactic halo. Mass: "allowed range of mχ ... between ∼ (10^11.8 - 10^13.6) GeV"; intro: "mass larger than ∼ 5 × 10^11 GeV". Lifetime (b b̄, fit to TA spectrum plus Amaterasu, CL not stated): "τχ = 5.5 × 10^28 s and 7.5 × 10^29 s for mχ = 10^12 GeV and mχ = 10^14 GeV". Photon limits: inconsistent, "in severe tension with the SHDM parameter space required to explain the TA Amaterasu event". Check: it calls the Amaterasu direction "outside the FOV of PAO" while printing PAO's declination range as -π/2 to +π/4; Dec 16.1° lies inside it, and inside Auger's 30° to 60° zenith photon range (latitude 35.2° S reaches Dec +24.8°, my arithmetic).
- 2504.15272 (Murase, Narita, Yin, JCAP 10 (2025) 109), scalar inflaton DM. Mass: "we predict that the DM mass is heavier than 10^13 GeV". Lifetime: no number printed; Fig. 3 (right top) marks Amaterasu with a green star at m_φ = 10^13 GeV on decay-rate upper limits from photon fluxes that "just sastify [sic] the photon limits". Photon limits: partial. "can be explained within the 2σ level for the ϕ → gg or ϕ → Hqq" (3σ for ϕ → H̄ll), "only explain it within 2σ level due to the severe multi-messenger γ-ray bound"; a dark-sector extension is offered "to further alleviate γ-ray bounds".
- 2408.07172 (Fargion et al., Universe 10, 323): Z-burst, "ZeV (Zeta, 10^21 eV)" neutrinos on relic neutrinos; "neutrino masses needed for this Z boson model could range from 0.1 to 0.4 eV"; "extreme tens of ZeV neutrino energy". No lifetime; UHE photon limits not discussed.
- Monopoles (no photon-limit discussion): 2312.08115 (mass "probably between 4 to 10 TeV"); 2403.12322 ("M ≤ 10^8 GeV", "M ≥ 10^4 GeV"); 2604.21099 (same bounds).

## 3. Lifetime limits

Photon-based:
- 2302.02993 (Das, Murase, Fujii, PRD 107, 103013; Auger limits; 95% CL): "For the bb channel, τχ ≳ 3 × 10^30 s at 10^12.2 GeV ≲ mχ ≲ 10^13.8 GeV"; "τχ ≳ 4 × 10^30 s at 10^13 GeV for the bb channel"; "τχ ≳ 10^30 s for 10^11 GeV < mχ < 10^15 GeV" (all SM two-body channels).
- 2406.03174 (TA limit, b b̄): "maximum of ∼ 3 × 10^30 s at around mχ ∼ 10^12 GeV"; 8 × 10^29 s at 10^14 GeV.
- 2508.08779 (Aloisio, Ambrosone, Evoli; νν̄ and b b̄; gamma, IceCube, Auger, KM3-230213A; 95% CL): "τχ ≳ 5 · 10^29 - 10^30 s for 10^7 GeV ≲ mχ ≲ 10^17 GeV".
- 2512.01638 (TA): "τX = 7.1 × 10^22 yr for the X → q q̄ model and τX = 1.7 × 10^22 yr for the X → ν ν̄ model", M_X printed "10^12 eV" (sic), "normalized to the Pierre Auger hybrid limits".
- 2601.11703 (Chianese et al., 10^7 to 10^15 GeV): figures only; abstract: "marginally less stringent than earlier evaluations". Printed values: NOT FOUND.

Neutrino-based:
- 2210.01303 (νν̄ channel): "τ > 10^27 s for mχ ∼ 10^11 GeV"; above 10^7.5 GeV KASCADE-Grande photon limits are "outperforming existing neutrino telescopes"; above 10^10 GeV "Auger-SD supersedes all other experiments".
- 2406.03174 (TA neutrino limit, b b̄): 7 × 10^25 s at 10^12 GeV; 4.8 × 10^25 s at 10^14 GeV.
- 2205.12950 (IceCube HESE): 1 PeV only, "larger than 10^28 s" (b b̄: "10^27 s").
- KM3NeT 2606.09986 prints best-fit lifetimes "10^26-10^27 s" above about 100 PeV, not a limit at 10^11 to 10^13 GeV.

Lifetime at which SHDM reaches the observed flux above 10^20 eV, against photon limits:
- 2302.02993: "The cosmic-ray flux constrains τχ to ≳ 4 × 10^29 s at 10^13 GeV for the qq decay channel"; the photon limit above is "an order of magnitude longer". SHDM cannot supply the full flux.
- 2406.03174: flux fit 5.5 × 10^28 s (10^12 GeV), 7.5 × 10^29 s (10^14 GeV); photon limits 3 × 10^30 s, 8 × 10^29 s.
- 1504.01319 (Aloisio, Matarrese, Olinto): "τX = 2.2 × 10^22 yr" at M_X = 4.5 × 10^13 GeV, from a photon limit "at the level of 2%". That is 6.9 × 10^29 s (my conversion), below the 3 × 10^30 s of 2302.02993 at that mass. Auger 2209.05926: "excludes that SHDM particles could explain the bulk of UHECRs, they can still contribute in a subdominant way".

## 4. Auger photon limits (95% CL, E^-2 photon spectrum)

- 2209.05926 (JCAP 05 (2023) 021, data to June 2020), Table 1, units 10^-3 km^-2 sr^-1 yr^-1, above 10, 20, 40 EeV: 2.11, 0.312, 0.172 (abstract: 3.12 × 10^-4, 1.72 × 10^-4). Index 1: 1.48, 0.273, 0.162. Index 3: 2.86, 0.332, 0.166. Fractions: "1.6%, 1.2% and 3.2% above 10 EeV, 20 EeV and 40 EeV". Above 40 EeV: NOT FOUND (paper: "will be discussed in a future work"). ICRC2025 2509.05113 repeats the three flux values.
- 2406.07439 (PRD 110, 062005), above 10 EeV: 0.0021 km^-2 sr^-1 yr^-1; fraction 0.77%.
- 0712.1147 (older SD): 3.8, 2.5, 2.2 × 10^-3 above 10, 20, 40 EeV; fractions 2.0%, 5.1%, 31%.
- TA for higher energies, 2512.01638 (14 yr, 95% CL, Table 2), above 10^19, 10^19.5, 10^20 eV: 2.3 × 10^-3, 8.2 × 10^-4, 3.0 × 10^-4 km^-2 sr^-1 yr^-1; fractions 6.8 × 10^-3, 2.1 × 10^-2, 0.21.

Method: quotes read from PDF text or arXiv HTML; the ≳ signs in 2302.02993 and 2406.03174 were confirmed in the HTML.
