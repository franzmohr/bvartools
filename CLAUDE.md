# bvartools

R package assisting in the set-up of algorithms for Bayesian inference of vector autoregressive and error correction models.

## Branching workflow

- `main` **is the development version**. It carries a development version number
  between releases, such as `1.2.3.9000`.
- New work goes on a short-lived branch created from `main` and is merged back
  into `main` **through a pull request**. Nothing is pushed to `main` directly,
  including small changes: the repository requires a pull request, and a direct
  push only succeeds by bypassing that rule.
- A release is a **tag plus a GitHub release**, not a branch and not a state of
  `main`. The released version is whatever the tag points at, so nothing is held
  back from `main` on the grounds that it is unreleased. A release candidate is
  the same thing, tagged `vX.Y.Z-rcN` and marked as a pre-release on GitHub.

### Cutting a release

On a branch off `main`, in one pull request:

1. `DESCRIPTION`: set `Version` to the release version, without the
   development suffix.
2. `CITATION.cff`: set `version` to the same, and set `date-released` to the
   day the release is published. The file names the release it is heading for
   and is never set to a release candidate, even while a candidate is the
   newest thing archived on Zenodo.
3. `.zenodo.json`: set `version` to the same. Zenodo reads this file in
   preference to `CITATION.cff`, so it is the one that decides how the archived
   record is labelled, and the `release-version` workflow does **not** check it.
4. `NEWS.md`: the section for the release carries its version as its heading.
5. `cran-comments.md`: describe the submission that is actually being made,
   and confirm the check results against a fresh run rather than an older one.

Merge that pull request, then tag the merge commit and create the GitHub release
from the tag. The `release-version` workflow checks the tag against `DESCRIPTION` and `CITATION.cff`,
ignoring any `-rcN` suffix, so a mismatch fails the tag.

Immediately afterwards, open a second pull request setting `DESCRIPTION` back to
a development version, `X.Y.Z.9000`.
