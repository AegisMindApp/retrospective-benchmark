# Provenance of the pre-registration

The paper's tamper-evidence argument rests on the pre-registration being frozen **before**
any post-2024 result was scored. That claim needs to be checkable, and this file states
exactly how far it can be checked from this repository — and where it cannot.

## What is checkable here

`PREREGISTRATION.md` contains the frozen design, the primary endpoint, the decoy protocol,
the 0.60 gate threshold, the commitment to publish regardless of outcome, and all twenty
amendments with their dates and rationale.

## What is not checkable here, and why

The original freeze is git commit **`92c6f8bf` (2026-07-14)** in a **private** monorepo
that also holds unrelated work, credentials and patent drafts, and cannot be opened.

**This repository is an extract. Its commit SHAs are new and do not include `92c6f8bf`.**
A reader therefore cannot verify the original timestamp from this repository alone.

We state this rather than leave it implicit, because a pre-registration whose timestamp
cannot be independently verified is exactly the kind of unfalsifiable claim the paper
criticises. Presenting an extract as though it carried the original chain of custody would
reproduce the problem being reported.

## Independent timestamp

To provide a verifiable date that does not depend on the private repository, this extract
is deposited to Zenodo, which issues its own timestamp and DOI:

- **DOI `10.5281/zenodo.21824633` was reserved but never minted and has never resolved.**
  The draft deposit was deleted on 7 Aug 2026 (Amendment 26). It is named here only so that
  anyone who saw it referenced in an earlier revision of these files knows there is nothing
  behind it. There is no archived copy of this study.

  That is the *version* DOI for this deposit. Once published, Zenodo also mints a
  **concept DOI** that always resolves to the latest version; cite the concept DOI in
  preference, so that any later correction or retraction is visible to a reader following
  the citation rather than hidden behind a pinned version.

That timestamp establishes the date of **publication of the extract**, not the date of the
original freeze. The two are different claims and are not conflated here.

The original commit and its full history can be shown to editors or reviewers on request.

## Amendments and self-reporting

Amendments 19 and 20 record our own procedural failures — a silent ligand-preparation
defect that made every ligand rigid, and calibration controls initially placed where they
could not serve their purpose. They are included because a pre-registration that logs only
protocol changes, and not the mistakes that forced them, is worth less than one that logs
both.

## One deliberate difference from the private repository

`PREREGISTRATION.md` refers throughout to a module `run_vina_docking.py` living in a
sibling directory (`amr_glass/docking/`). In this extract the protocol-neutral parts of
that module are provided as **`vina_tools.py`** instead, and the imports are repointed.

The original module is a target-specific script that also contains candidate compound
structures and patent-evidence generation for unrelated work; it is not published here.
The docking functions themselves — receptor preparation, ligand preparation, the Vina
invocation — are byte-identical to the versions used in the study.

**`PREREGISTRATION.md` has not been edited to match.** It records what was true when each
amendment was written, and silently rewriting a pre-registration to tidy up a later
refactor is precisely the practice this study reports on. The discrepancy is documented
here instead.
