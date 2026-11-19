# How to publish this code and get a DOI

You need your own GitHub and Zenodo accounts; nothing here has been uploaded anywhere yet.
Two routes — pick one.

---

## Route A (recommended): GitHub repository + Zenodo DOI

1. Create a **public** repository on GitHub, e.g. `pulpitis-compartment-immunity`.
   Do **not** initialise it with a README (this folder already has one).
2. Upload the contents of this folder (drag-and-drop in the GitHub web interface is enough, or
   `git init && git add . && git commit -m "v1.0.0" && git remote add origin <url> && git push`).
3. Go to Zenodo → Settings → GitHub → flip the toggle for that repository to ON.
4. On GitHub, create a Release tagged `v1.0.0`. Zenodo then archives the snapshot automatically
   and mints a DOI (a citable DOI, not the "concept" DOI).
5. Copy the DOI badge, paste it into `README.md` §7, and use it in the manuscript.

Resulting manuscript text for the Code availability section:

> All analysis code and the numerical result tables are available at
> https://github.com/<user>/pulpitis-compartment-immunity and archived at Zenodo
> (https://doi.org/10.5281/zenodo.XXXXXXX).

## Route B (simplest): Zenodo upload only

1. Zip this folder.
2. Zenodo → New upload → Upload the zip → fill title/creators/description → Publish.
3. Zenodo mints a DOI immediately; use the same manuscript sentence, replacing the GitHub URL
   with "available at Zenodo".

Either route satisfies Scientific Reports' data/code availability policy. Route A is preferable
because reviewers can browse individual scripts.

---

## Two consistency points after the DOI exists

1. The manuscript currently reads "on reasonable request"; replace the whole sentence with the
   text above (one place, Code availability section).
2. `README.md` §7 contains `DOI: <to be inserted after archiving>`. Replace it with the real DOI
   so the repository and the paper point at each other.

---

## What is in this folder

- `R/bulk/` bulk-transcriptome scripts
- `R/singlecell/` GSE185222 scripts
- `R/emtab/` E-MTAB-14896 scripts, including the leave-one-patient-out and permutation tests added
  in the final revision (`56_endo_ratio_lopo.R`, `57_LR_patient_perm.R`)
- `R/manuscript/` figure assembly and manuscript generation
- `results/` 59 numerical tables behind the figures
- `README.md`, `LICENSE`

Raw GEO/ArrayExpress matrices and the intermediate expression RDS files are deliberately excluded
(regenerable from the accessions and far too large for a repository).
