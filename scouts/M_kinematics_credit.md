# Scout M: kinematics credit (1 Oct 2026)

**Summary.** (a) is in print in general form (Bödeker and Moore 2009); the at-rest value is a substitution into it and is in no paper I read. (b) is in print at rest (Garcia Garcia, Koszegi, Petrossian-Byrne, Dec 2022). The threshold in (c) is in print since 1992. The product relation and the one-bounce ceiling are NOT FOUND. Check script: `M_kin_check.py` in this folder (explicit boosts, all pass, one control fails on purpose). Correction to (c): the product relation needs the transverse mass.

## 1. Fact (a)

- Bödeker and Moore (B&M) 0903.4099, Eq. (3.8), plasma frame: "Since the wall passes at the speed of light and purely in the +z direction, it leaves each particle's p⊥ and E − pz unchanged. However E² − p² must change by Δm² ... (E, pz, p⊥) → (E + [m²(h2) − m²(h1)]/[2(E − pz)], pz + [m²(h2) − m²(h1)]/[2(E − pz)], p⊥)". Footnote 9: "E − pz conservation is the same as energy conservation in the wall rest frame." Setting E = m_in, pz = 0 gives (m_out² + m_in²)/(2 m_in); checked numerically (50.5 for m_in = 1, m_out = 10). B&M do not write the at-rest case.
- Wall frame, B&M 0903.4099 Eq. (3.6): "pz,in − pz,out = −(m²(h2) − m²(h1))/(2E) + O(m²/E³, p⊥²/E³)"; the minus sign comes from their axes (plasma moves along −z). B&M 1703.08215 Eq. (1) has the letter's sign: "Δp1→1 = pz,in − pz,out ≃ (m²a,h − m²a,s)/(2E)". Azatov and Vanvlasselaer 2010.02590 Eq. (8): "(pz^s − pz^h) ≃ (m2² − m1²)/(2p0)".
- Scaling form only: Baldes et al. 2403.05615 §2.3.2, "the change in momentum in the z-direction is given pX ≃ Δm²a/Tn" (bath frame).
- The at-rest formula itself: NOT FOUND.

## 2. Facts (b) and (c)

**(b)**
- Garcia Garcia, Koszegi, Petrossian-Byrne 2212.10572 Eq. (55): "γdr = (1+|veq|²)/(1−|veq|²) ≃ 2γeq² ≫ γeq", for dark photons at rest far from the wall.
- Kawana, Lu, Xie 2206.09923 Eq. (2.8), general case: "Ẽ = [(1 + vw²)E + 2vw pz]/(1 − vw²)"; pz = 0 gives γ²(1+β²)E.
- Garcia Garcia and Petrossian-Byrne 2407.09603 Eqs. (75), (76): "ωr = (1+v²)γ²m ≃ 2γ²m", "γshell ≡ ωr/m ≃ 2γ² ≫ γ".
- Baker et al. 2105.07481, non-relativistic: "On each reflection, the particle's energy increases by dE = 2vw E."

**(c) threshold**
- Dine, Leigh, Huet, Linde, Linde hep-ph/9203203 (1992), App. A: "sufficient momentum (p²x > (m²1 − m²o)) to go through the slice, or are reflected."
- Arnold hep-ph/9302258 (1993), Eq. (1.2) and text: "The momentum transferred to the wall will be 2p1z if the particle is reflected and p1z − p2z if it is transmitted".
- Lewicki, Vaskonen, Veermäe 2205.05667 Eq. (4); Lewicki et al. 2305.07702 Eq. (2): "F_i(u) = 2, u² < m_j² − m_i²".
- Baker, Kopp, Long 1912.02830, massless-in form, citing B&M: "will only have sufficient energy to enter the bubble if γw(pz + vw|p|) > m_χ^in". Their "in" is inside the bubble: reversed labels.

**(c) product relation and ceiling: NOT FOUND.** Closest, qualitative only: Baker et al. 2110.00005, "the overdensity does not extend far beyond pr ≈ m∞χ since particles with such large momenta are able to pass through the wall and escape the bubble".
My check, 5,789 random reflections: E·E' = γ²(m_in² + p⊥²) + u² holds to 1e-56. The letter's form without p⊥ is exact only at p⊥ = 0 and is off by up to 100% otherwise. The ceiling holds: for m_in = 1, m_out = 10, E' → 199 = (2m_out² − m_in²)/m_in as γ → m_out/m_in from below; just above that γ the particle is transmitted.

## 3. IDs

All nine IDs exist and match, with one wrong pairing. 2010.02590 is Azatov, Vanvlasselaer, "Bubble wall velocity: heavy physics effects" (JCAP 01 (2021) 058, Oct 2020). "Dark matter production from relativistic bubble walls" is 2101.05721, Azatov, Vanvlasselaer, Yin (JHEP 03 (2021) 288). Filtered DM: arXiv Dec 2019, PRL 2020, so both years are right. 2306.15555 is "Bubbletrons", Baldes, Dichtl, Gouttenoire, Sala (arXiv Jun 2023, PRL 2025). 2412.18752 is Shakya, "Cosmic Colliders", a Dec 2024 proceeding. 2407.09603 (Garcia Garcia, Petrossian-Byrne, 2024) and 2305.07702 (Lewicki et al., 2023) match; Eq. (2) is as you say.

## (d) credit

- Bubbletrons 2306.15555: "In this letter we point out that cosmological particle accelerators may also have existed"; Eq. (4) s_coll ≃ 4γ²coll E²V.
- Shakya 2412.18752: "Collisions of vacuum bubbles in the early Universe can act as cosmic-scale high-energy colliders with energy reach close to the Planck scale." Wall-wall collision (p_max ≈ M_P, Eq. 7), a different mechanism; credits Watkins and Widrow 1992, Falkowski and No 2012.
- Baker, Kopp, Long: "most of them are reflected and quickly annihilate away".

## 4. Present-day sources: searches

Ran over 180 INSPIRE queries (Fermi ball, Q-ball, soliton, false vacuum remnant or pocket, bubble wall, vacuum bubble, domain wall, texture, KM3-230213A, Amaterasu, against cosmic ray, UHECR, neutrino, burst, GZK), about 20 arXiv API queries, 15 web searches; Google Scholar gave a captcha.
NOT FOUND: a proposal that collapsing false-vacuum pockets, Fermi balls, Q-balls or solitons give present-day UHECR, photons or neutrinos by a wall mechanism. "Fermi ball" AND "cosmic ray": 0 hits. Nearest:
- Sakharov, Konoplich, Gogberashvili 2506.23387: merging black holes nucleate true-vacuum bubbles; "Collisions between these bubbles produce microscopic black holes that rapidly evaporate via Hawking radiation, emitting intense, short-lived bursts of neutrinos with energies exceeding 100 PeV." Scout I did not list it.
- Brandenberger, Cyr, Jiao 2005.11099: textures collapse "at close to the speed of light"; the process "will continue to the present time"; output includes "high energy photons and ultra-high energy neutrinos". A global defect, no vacuum pocket.
- Sengupta, Stojkovic, Dai 2501.15848: "bubbles of true vacuum might already exist in our universe"; "a macroscopically large burst of high energy neutrinos and photons from Higgs decays". Expanding wall, no EeV figure.
- Yin 2505.15764: walls "effectively making the walls a source of cosmic rays. For a typical injected momentum pX ∼ mϕ γw", flux ≃ 10 km⁻² yr⁻¹ (Eq. 23).
- Q-balls: Addazi et al. 2510.00254 (gamma rays beyond PeV); Kasuya et al. 2403.01675 (decay now, MeV).

## 5. Extra items

1. Kusenko, Loveridge, Fong, Q-balls as UHECR sources: NOT FOUND. 32 INSPIRE records by Kusenko mention Q-balls, none on UHECR, none with Fong. Real: Kusenko, Loveridge, Shaposhnikov hep-ph/0405044, astro-ph/0507225.
2. Batell, Caren, Fong, "Domain wall annihilation as a dark matter source": NOT FOUND (no title or author-pair match).
3. Batell, Verhaaren on collapsing bubbles or walls: NOT FOUND. Their joint papers (1904.10468, 2004.10761) are on twin Higgs. The only INSPIRE title with "catastrophic particle production" is Yoshimura hep-th/9506176.
4, 5. Long and Schleich, "Catastrophic particle production by a thin wall"; B. Link, "Acceleration of particles by moving cosmological domain walls": NOT FOUND.
6. Tye, Wexler: NOT FOUND; hep-ph/0011033 is Hosek on superfluid QCD. Real: Davoudiasl, Hewett, Rizzo hep-ph/0010066, "Gravi burst: Super GZK cosmic rays from localized gravity".

**Simultaneous showers.** A direct UHE (above 1e19 eV) limit on several showers within microseconds to seconds from one direction: NOT FOUND for Auger and TA. Auger multiplet searches (1111.2472, 2004.10591) look for direction and 1/E alignment; the source "should be steady", so they do not test simultaneity. What exists:
- CHICOS, Carlson et al. astro-ph/0411212, E > 1e14 eV: "coincidence times ranging from 1 µsec up to 1 second. The results are consistent with the absence of excess coincidences except for a 2.9σ excess observed for coincidence times less than 10 µsec. We report upper limits for the coincidence probability as a function of coincidence time."
- Albin, Whiteson 2102.03466: the Gerasimova-Zatsepin effect "was considered by The Large Area Air Shower (LAAS) observatory [51] but has not been observed"; of mechanisms including top-down "collapse processes", "To date, no evidence of these mechanisms have been experimentally confirmed".
- CREDO, Homola et al. 2010.08351: ensembles "not yet probed"; limits only promised.
- TA, Abbasi et al., Phys. Lett. A 381, 2565 (2017): "several short-time bursts of air shower like events", correlated with lightning.

## Table

| Fact | Earliest or clearest source | Letter should say |
|---|---|---|
| (a) E = (m_out² + m_in²)/(2 m_in) | B&M 0903.4099 Eq. (3.8), fn 9 | "following Bödeker and Moore" for E − pz; at-rest value as our substitution |
| (b) E' = γ²(1+β²)m ≈ 2γ²m | Garcia Garcia et al. 2212.10572 Eq. (55); general Kawana et al. 2206.09923 Eq. (2.8) | "following Garcia Garcia et al." |
| (c) reflect if u² < m_out² − m_in² | Dine et al. 1992 App. A; Arnold 1993 Eq. (1.2); Lewicki et al. 2305.07702 Eq. (2) | "following Dine et al. and Arnold" |
| (c) E·E' = γ²(m_in² + p⊥²) + u²; ceiling (2m_out² − m_in²)/m_in | none found; script verifies | "new as stated", a short consequence of the sources above |
| (d) walls as accelerators or colliders | Baldes et al. 2306.15555; Shakya 2412.18752 (wall-wall); Baker, Kopp, Long 1912.02830 | "following" all three; present-day use new, nearest 2506.23387, 2505.15764 |
