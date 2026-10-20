# Profex.str to TOPAS.str Converter

A Python script for the batch conversion of Profex/BGMN .str crystal structure files to TOPAS format .str files.

Code written by prompting Microsoft Copilot. Code review and bug fixes (v2 updates, October 2026) aided by Claude (Anthropic). All crystallographic requirements and conversion rules were defined, tested and validated by the author.

[![DOI](https://zenodo.org/badge/1187204487.svg)](https://doi.org/10.5281/zenodo.19134548)

# Purpose

Batch conversion of Profex .str files to a format readable by TOPAS, e.g. the structures supplied with a Profex install. TOPAS users without ICDD PDF-5+ can then use starting structures comparable to those available to Profex users, where it matters (e.g. QPA).

Every Profex .str produces one uniquely named TOPAS .str.

# Design rules

* **The source must prove itself.** The Profex .str is assumed to possibly be wrong. Every site is checked against SPACEGRP.DAT from the Profex install, which is assumed to be correct.
* **Invent nothing.** Free coordinates are never guessed. A site that cannot be resolved is skipped and flagged, never written with made-up values.
* **Preserve everything.** The original Profex line is kept as a comment beside every TOPAS line it produced, so the user can always check the conversion.
* **Name parameters meaningfully.** Parameter names are derived from the phase name (e.g. a_Calcite, scale_Calcite, csl_Calcite), so no two phases in a refinement collide.
* **Leave Hermann–Mauguin notation to the user.** The Profex HM symbol is kept verbatim as a comment; the user confirms a form TOPAS accepts.
* **Report every outcome.** Each output folder gets a CSV report listing every file with its status and the reason for any problem.

# How to use

*Ensure Profex is installed.* The script needs SPACEGRP.DAT from the Profex install.

1. Download and install Python 3.
2. Download `ProfexStr_to_TOPAS_v2.py`.
3. Open Command Prompt in the folder containing the .py file.

Convert one folder:

```
py ProfexStr_to_TOPAS_v2.py ^
  --in "C:\Program Files\Profex5\Structures\SUB-FOLDER" ^
  --out "C:\STRING2\output\SUB-FOLDER"
```

Convert every Profex structure folder in one go (paste into Command Prompt as is; in a .bat file use `%%F` and `%%~nxF`):

```
for /d %F in ("C:\Program Files\Profex5\Structures\*") do py ProfexStr_to_TOPAS_v2.py --in "%F" --out "C:\STRING2\output\%~nxF"
```

* Use `python` instead of `py` if that is how Python is installed on your PC.
* Write the output somewhere you own (not inside Program Files).
* Use a new, empty output folder for each run. Existing files are not overwritten; new copies are saved as `name_1.str`.
* SPACEGRP.DAT is looked for at `C:\Program Files\Profex5\BGMNwin\SPACEGRP.DAT`. Use `--dat "path\to\SPACEGRP.DAT"` if yours is elsewhere.

# THE USER SHOULD VERIFY THE OUTPUT IS CORRECT AND MATCHES THE INTENT OF THE PROFEX STR. THERE MAY BE FRINGE CASES I HAVE MISSED.

Most issues during testing came from formatting differences or typos in the Profex .str files (for example `SpacegropuNo` instead of `SpacegroupNo`, handled by fuzzy matching). The original Profex data is preserved as a comment for ALL atomic sites for this reason.

# What the converter does

**Space group and setting**
* Uses the Setting given in the Profex .str; otherwise a unique HM match, a rhombohedral-axes check, or Setting 1 with a NOTE.
* If any site fails against that entry, **every** SPACEGRP.DAT entry for the same space group number is tested. Only if exactly one entry matches ALL sites (and suits the given cell) is it used, with a WARNING naming both entries. If it contradicts the Setting/OriginChoice stated in the Profex .str, the warning says so. If several entries match, nothing is changed and the file is flagged as ambiguous.

**Atom sites**
* Special-position coordinates are checked against all equivalent positions in SPACEGRP.DAT, including centring translations. Coordinates fixed by symmetry are filled in from SPACEGRP.DAT as exact fractions (1/3, not 0.3333).
* General positions must have x, y and z given, otherwise the site is skipped.
* Coordinates are output modulo 1.
* Beq = 100 × TDS. If no TDS is given, beq 1 is inserted (the TOPAS GUI default when a structure without Beq is imported; note that in launch mode TOPAS uses Beq 0 if none is given).

**Unit cell**
* Profex min/max limits are carried over; if none are given, ±1% is applied.
* Cell parameters fixed or linked by symmetry (e.g. b = a, γ = 120°) are written as constraints, not as refinable values. If the Profex file states such a value (e.g. `GAMMA=120` in a trigonal file), it is kept as a comment only.

**Occupancies**
* Mixed sites are expanded into one TOPAS site per element with shared coordinates.
* Refinable occupancies (`PARAM=p=...`) become named `prm` parameters with the Profex limits, or 0–1 if no limits are given. Each site gets its own name, so sites are never linked by accident.
* Occupancy values defined without `PARAM=` (on the atom line or in the header) are fixed in BGMN, so they are written as fixed parameters (`prm !name value min 0 max 1`) with a comment: they can be refined by removing the `!`, but that is probably not the intent of the source .str.
* Values defined as formulas in the header (e.g. `pAL=(1-pMG-pFE)`) become linked TOPAS parameters (`prm pAL_phase = (1-pMG_phase-pFE_phase);`).

**Also inserted:** `scale`, `MVW()` and `CS_L(200 min 32 max 5000)` with derived parameter names.

# CSV report

Each output folder gets `_canary_report_<date_time>.csv` with one row per input file:

| Column | Meaning |
|---|---|
| status | `OK`, `SKIPPED_SITES` (some sites not written), `UNDEFINED_SYMBOLS` (an occupancy symbol could not be resolved) or `ERROR` (file not converted; see `<name>_ERROR.txt`) |
| setting_note | How the space-group entry was settled when sites failed: resolved by site match, ambiguous, or no entry matches |
| skip_reasons | Why sites were skipped, with counts |
| error | The reason a file could not be converted |

Other columns record the space group, setting, origin, HM symbol and lattice used.

# Results on the structures supplied with Profex

1,113 structure files in 11 folders (Profex 5):

* **1,108 OK**, including CoFe2O4, Cu2O and ZnAl2O4, whose origin was found by testing every SPACEGRP.DAT setting (all sites fit origin 2; for CoFe2O4 and Cu2O this contradicts the Profex .str).
* **5 with skipped sites:**
  * Arsenolite.str: As given on Wyckoff e (x, x, x) as x=0.3529, y=0, z=0. No equivalent position matches.
  * CRYOLITE.STR: Wyckoff a and b are both given as 0, 0, 0, which is impossible for P12_1/c1.
  * CSH-0625.str (CaSO4·0.625H2O): one Ca given on Wyckoff b as (0.2858, 0.2858, 1/6), which matches no position for b. Stoichiometry suggests the coordinates, not the Wyckoff letter, are in error.
  * nontronite15a.str and smectitedi2wfix1.str: turbostratic/layer-disorder models with coordinates given as formulas. Out of scope for v2.
* **0 not converted.**

# Limitations

* Validated on the 1,113 structures supplied with Profex. Other .str sources (hand-written, CIF-derived) have not been tested yet.
* Only `SpacegroupNo=` is read; a file giving only a Hermann–Mauguin symbol will fail with ERROR.
* Coordinates or TDS given as formulas are not evaluated.
* The TOPAS space_group is written as a number (Setting 1, cubic/tetragonal/orthorhombic) or as the SPACEGRP.DAT HM string otherwise. Always confirm it is in a form TOPAS accepts and matches the intended origin.

# Planned for v3

* Turbostratic/layer-disorder models (BGMN supercell approach → TOPAS equivalent)
* .cif input
* Reverse conversion, TOPAS .str → Profex .str
* Hermann–Mauguin translation to TOPAS-ready symbols
