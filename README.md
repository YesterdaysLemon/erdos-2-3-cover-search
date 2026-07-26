# Erdős \(2^k3^\ell m+1\) cover search

> **Research status: unresolved.** This repository does not contain a value
> of \(m\), a proof that one exists, or a proof that none exists.

The problem studied here is to determine whether there is an integer
\(m\geq 1\), with \(\gcd(m,6)=1\), for which

\[
2^k3^\ell m+1
\]

is composite for every pair \(k,\ell\geq 0\).

The computational approach searches for finite systems of prime divisors.
For a suitable prime \(p\), the condition
\(p\mid 2^k3^\ell m+1\) selects a periodic affine fibre in the
\((k,\ell)\)-lattice. A finite fibre cover of the complete lattice would
produce congruences for \(m\), which could then be combined by the Chinese
remainder theorem.

## Current certified checkpoint

As of 2026-07-26:

- no candidate \(m\) has been found;
- no global impossibility theorem has been proved;
- 4,609 of 4,637 ranked finite divisor-period families have independently
  replayed no-cover certificates;
- 28 ranked finite families remain in the intersection of the aggregate
  block-star and separately verified single-anchor frontiers.
- for the closest remaining family, at period `2533395664800`, 22
  independently replayed conditional-overlap bounds reduce the exact
  period-level upper bound to
  `730799865304831267447/729328817548572040800`, or
  `1.0020169883883154`. Because this is still greater than one, that family
  remains unresolved.

The fraction \(4{,}609/4{,}637\) measures only this deliberately chosen
finite menu. It is **not** a probability, a completeness claim, or evidence
that the original infinite problem is 99.4% solved.

The exact checkpoint and the prominent correction to an earlier invalid
parallel-class argument are recorded in [RESEARCH_LOG.md](RESEARCH_LOG.md).

## Repository contents

- `*.py`: search programs, exact certificate generators, independent
  verifiers, and regression tests.
- `order_pool_1050000_component_core_corrected_max32_stable.json`: stable
  finite candidate pool used by the current frontier.
- `order_pool_1050000_max32_period_divisor_family_ranking_all.json`: the
  4,637 ranked divisor-period families.
- `order_pool_1050000_max32_all4637_pairanchor_star_scan_v35.json`: aggregate
  finite scan.
- `order_pool_1050000_max32_pairanchor_survivors_singleanchor_star_scan_v26.json`
  and its verification: the independent single-anchor frontier.
- the `order_pool_1050000` exact rational certificates and replay reports
  needed by the current checkpoint and its dependency closure.
- [PUBLICATION_MANIFEST.md](PUBLICATION_MANIFEST.md): what the public snapshot
  includes and deliberately excludes.

The files remain flat because certificates refer to one another by basename.

## Reproducing checks

Python 3.11 or later is recommended.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m unittest discover -v
```

An independently replayed conditional-fibre certificate can be checked with:

```powershell
python verify_projected_pair_conditional_fibre_overlap.py `
  order_pool_1050000_component_core_corrected_max32_stable.json `
  order_pool_1050000_conditional_fibre31_vs_full_pair29_239_53_67_23_71_47block_certificate.json `
  --output replay.json
```

Successful replay prints `verified=True`. This proves only the finite
intersection claim encoded by that certificate.

## Evidence policy

Discovery programs do not certify their own output. Important results are
stored as exact integer or rational data and replayed by a separately written
verifier using a different enumeration orientation or contraction. Finite
checks are always labeled finite; absence of a found cover is never promoted
to a global theorem.
