# Data

Two tables of published cosmic-ray events, used by `calc/p12_halo_directions.R`. Every row names its
source in the last column, and numbers are copied as printed.

- `uhecr_above_100EeV.csv`: events at or above 100 EeV.
- `uhecr_57_to_100EeV_control.csv`: events between 57 and 100 EeV, the control sample.

## Sources and terms

- **Pierre Auger** rows are from the Pierre Auger Open Data (catalogue of the highest-energy events,
  Abdul Halim et al. 2023, ApJS 264, 50, arXiv:2211.16020; opendata.auger.org). The Open Data are
  released under CC BY-SA 4.0, and these rows remain under that licence. Cite: Pierre Auger
  Collaboration, Auger Open Data, DOI 10.5281/zenodo.4487612.
- **Telescope Array** rows are from Abbasi et al. 2014, ApJ 790, L21 (arXiv:1404.5890), Table 1, and
  the 244 EeV event from Science 382, 903 (arXiv:2311.14231), Table 1.
- **AGASA** rows are from the AGASA results page (Table 1, modified 24 June 2003), Takeda et al. 1999
  (ApJ 522, 225) and Hayashida et al. 2000 (astro-ph/0008102).
- **Fly's Eye**: Bird et al. 1995, ApJ 441, 144. **Yakutsk**: Glushkov et al. 2025, arXiv:2506.00517.

Energy scales differ between experiments, so the rows should not be pooled by energy without
rescaling.
