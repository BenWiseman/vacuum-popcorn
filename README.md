# Vacuum popcorn

Code, data and manuscript for the letter

**Vacuum popcorn: ultra-high-energy particles from collapsing false-vacuum pockets**
B. H. Wiseman (ORCID 0009-0002-1023-9026)

Preprint: https://doi.org/10.5281/zenodo.23096148
Archive of this repository (v1, as submitted): https://doi.org/10.5281/zenodo.23099969

A pocket of false vacuum that traps particles which are light inside it and heavy outside throws
them out when it collapses. This repository holds every calculation behind that claim. Each number
in the letter carries a note naming the script that produces it.

## Layout

    letter/   the manuscript (vacuum_popcorn.md), its figure, its reference records and its checks
    letter/preprint/   a plain typeset copy of the manuscript (PDF, LaTeX, reference list)
    calc/     the calculations (R, with one Python reference check) and the output of each
    data/     published cosmic-ray events above 57 EeV, with the source of every row
    scouts/   notes on the literature, each fact with its source

## Running it

R (base packages only) and Python 3 (standard library only).

    cd calc
    ./run_all.sh              # every script, about an hour and a half on 14 cores
    LOADED=0 ./run_all.sh     # the same without the hour-long loaded-collapse scan (p10)

Each script writes its output to the `.txt` file of the same name, and those files are in the
repository as run. `run_all.sh` names any script that fails.

| scripts | what they do |
|---|---|
| `pk_lib.R`, `p1_*` | one particle meeting one moving wall, exactly; checked against 100-digit arithmetic |
| `pt_lib.R`, `pt3_lib.R`, `p2`, `p4*` | following a particle through a collapse; checked against brute-force time stepping |
| `p3*`, `p6`, `p7` | the escape spectra and the mechanism |
| `p8`, `p15` | the two figures: the spectra, and the sketch of the sequence |
| `p5*` | what Earth would see: the photon cap, simultaneous showers |
| `pl_lib.R`, `p9`, `p10`, `p11` | a collapse whose particles slow the wall |
| `p12`, `p13`, `p14` | which published events could be pops: arrival directions, photon searches, photons per proton |

The checks plant errors on purpose and require them to be caught, so a check that passes has been
shown able to fail.

## Manuscript

`letter/vacuum_popcorn.md` is the source. From `letter/`:

    python3 cite_check.py vacuum_popcorn.md      # citations against the reference list and INSPIRE records
    python3 style_check.py vacuum_popcorn.md     # the letter's style sheet (STYLE.md)

`letter/preprint/vacuum_popcorn.pdf` is the same text typeset; `vacuum_popcorn.tex` there compiles on
its own with pdflatex.

## Use of AI tools

The letter's acknowledgements say how language models were used in writing the code, searching the
literature and drafting and checking the text.

## Licence

Code (`.R`, `.py`, `.sh`): MIT. Manuscript, figures, literature notes and documentation: CC BY 4.0.
The Pierre Auger rows in `data/` are from the Pierre Auger Open Data and stay under CC BY-SA 4.0;
`data/README.md` gives the source of every table. See `LICENSE`.

## Citing

`CITATION.cff` has the details. Until the letter is published, cite it as B. H. Wiseman (2026),
"Vacuum popcorn: ultra-high-energy particles from collapsing false-vacuum pockets", with this
repository's address.

## Contact

benjamin.h.wiseman@gmail.com
