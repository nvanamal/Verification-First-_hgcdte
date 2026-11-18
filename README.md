# Verification-First HgCdTe Screening — Paper 2

Code and derived simulation data supporting **Verification-First Surrogate-Assisted Screening of Square-Footprint ITO Nanostructures for HgCdTe Infrared Photodetectors** (V7 computational record).

**Author:** Dr. Nagendrababu Vanamala · [ORCID](https://orcid.org/0000-0002-0909-0343)

This repository contains the analysis and figure-generation workflow, not a claim of journal acceptance or an independently reproducible COMSOL material model.

## Study scope

- 300 square-footprint geometries × 20 wavelengths = 6,000 Main300 responses over 2.5–5.0 µm.
- Documented geometry-level allocation: 210 training, 45 validation and 45 retrospective held-out audit geometries. All wavelengths from each geometry remain together.
- Separate HgCdTe and ITO absorbed-power predictions; full-penalty objective `J = mean(P_HgCdTe) - mean(P_ITO)`.
- Independent COMSOL verification of 12 previously frozen guided candidates and 12 random controls.
- Candidate-specific BO06/M201 mesh checks, BO06 spectral-grid check, objective-weight sensitivity and native-ITO-range sensitivity.

The guided/random comparison characterizes these two frozen sets. The complete historical screening pool cannot be regenerated, and no general repeated-run optimization superiority is claimed. Original screening predictions remain distinct from later audit refits.

## Reproduce the analysis

Use Python 3.12 with the versions in `requirements.txt`. From the repository root:

```bash
python -m pip install -r requirements.txt
python analyze.py
python plot_results.py
python style_spectra.py
python plot_paired_spectra.py
python plot_compact_supp.py
python refinement_analysis.py
python plot_refinement.py
python reviewer_analyses.py
python native_range_analysis.py
python plot_workflow.py
```

`analyze.py` refits models and learning curves and may take several minutes. It writes derived tables to `data/`, model artifacts locally and plots to `figures/`. The later scripts regenerate spectra, refinement results and sensitivity analyses. Scripts resolve paths from their own location. This reproduces analysis of the archived outputs; it does not run COMSOL or reconstruct the unavailable historical screening pool.

For sensitivity checks without retraining, the supplied derived CSVs allow `python native_range_analysis.py`, `python reviewer_analyses.py` and `python refinement_analysis.py` directly.

## Contents

| Path | Contents |
|---|---|
| `inputs/` | Three original corrected response/analysis workbooks; no optical-source tables |
| `data/` | Derived responses, exact partitions, model metrics, learning subsets, candidate records, sensitivity and configuration tables |
| `support/` | Documented reference partition recovery |
| `refinement_exports/` | Nine BO06/M201 port and integrated-power CSV exports |
| `provenance/` | Author-provided production, BO06/M201 and material configuration transcriptions |
| `Material_Provenance.json` | Unit handling, material domains, ITO hash/range and disclosure restrictions |
| `environment.txt` | Recorded analysis environment |
| `PUBLIC_FILE_SHA256SUMS.txt` | Checksums of the published files |

Historical filenames containing “Blind” are retained for data compatibility; final publication terminology is retrospective held-out audit.

## Restricted materials and interpretation

**HgCdTe refractive-index/extinction-coefficient source data are restricted and are not included.** The measurements originate from the U.S. Army Research Laboratory framework; access or redistribution requires authorization from the data owner. Files named `*_hgcdte.csv` contain simulated integrated absorbed powers, not HgCdTe n/k values. Derived optical responses are supplied for analysis.

The raw `ito.txt` file and COMSOL `.mph` models are excluded. The archived ITO input uses explicit micrometer conversion, linear interpolation and constant extension beyond its native 2.623526293–4.712033503 µm interval. Its original bibliographic and film-specific provenance remains unresolved; background ITO literature is not asserted as the proven source. The four-wavelength n/k check is **ITO only** and is an offline reconstruction, not a direct COMSOL numerical export.

Historical nominal 20%-Cd, approximately 5 µm HgCdTe film details provide provenance and are not substituted for present-model geometry/composition. BO06 and M201 share the documented physical configuration but use candidate-specific refined meshes, with 267,780 and 303,002 solved degrees of freedom.

The engineering objective is not measured responsivity or detectivity. Candidate rankings and group advantage depend on the ITO penalty weight. No new simulation is performed by these Python scripts.

## Citation and release status

Use `CITATION.cff` for this software record. The manuscript remains a research manuscript; no journal publication DOI or archived-release DOI is assigned here. No software/data license grant is specified in this repository; third-party and restricted-source rights remain with their owners.
