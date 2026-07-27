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
- 4,617 of 4,637 ranked finite divisor-period families have independently
  replayed no-cover certificates;
- 20 ranked finite families remain after intersecting the aggregate
  block-star and separately verified single-anchor frontiers and applying
  23 exact period certificates;
- the former closest family, at period `2533395664800`, is now eliminated.
  A 16-anchor block and 40 independently replayed conditional-overlap bounds
  give the exact period-level upper bound
  `40107466081993334654251/40113084965171462244000`, or
  `0.9998599239329758`, which is strictly less than one.
- period `101264763600` is also eliminated. Promoting `p=71`, `p=31`, and
  `p=191` into a 17-anchor block gives the independently replayed upper bound
  `2102349443216882789/2103497917325925120`, or
  `0.9994540169973154`.
- period `13127595717600` is eliminated by another 17-anchor block, built by
  promoting `p=31`, `p=71`, and `p=191`, plus one conditional `p=37` edge.
  Its independently replayed upper bound is
  `3065470386035078889247/3068621248660075836000`, or
  `0.9989731992417205`.
- period `4247163320400` is eliminated by a 16-anchor block promoting
  `p=31`, `p=601`, `p=191`, and `p=59`, plus 18 conditional edges. Its
  independently replayed upper bound is
  `17632397122197491960353/17632840928789005920000`, or
  `0.9999748306813799`.
- period `216497080800` is eliminated by a 17-anchor block promoting
  `p=71`, `p=191`, and `p=601`, plus 14 conditional edges. Its independently
  replayed upper bound is
  `1301657254484087906011/1301971833192270240000`, or
  `0.999758382861931`.
- period `330442912800` is eliminated by a 16-anchor block promoting
  `p=71`, `p=601`, `p=31`, and `p=103`, plus 21 conditional edges. Its
  independently replayed upper bound is
  `2599707959748681404489/2599821871898771520000`, or
  `0.9999561846327545`.
- period `1659810952800` is eliminated by a 19-anchor block promoting
  `p=191`, `p=601`, `p=31`, `p=311`, `p=599`, and `p=647`, plus 40
  conditional edges. Its independently replayed upper bound is
  `35318499010624116759961/35320158961985177280000`, or
  `0.9999530027211133`.
- period `776363187600` is eliminated by an 18-anchor block promoting
  `p=31`, `p=191`, `p=59`, `p=599`, and `p=601`, plus 38 conditional
  edges. Its independently replayed upper bound is
  `602703301209279912260669/602739799060432793760000`, or
  `0.999939446754288`. The exact no-cover margin is approximately
  `0.000060553245712`.

The fraction \(4{,}617/4{,}637\) measures only this deliberately chosen
finite menu. It is **not** a probability, a completeness claim, or evidence
that the original infinite problem is 99.4% solved.

The exact checkpoint and the prominent correction to an earlier invalid
parallel-class argument are recorded in [RESEARCH_LOG.md](RESEARCH_LOG.md).
The 20 unresolved periods are indexed in
[CURRENT_FINITE_FRONTIER.md](CURRENT_FINITE_FRONTIER.md).

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
- [CURRENT_FINITE_FRONTIER.md](CURRENT_FINITE_FRONTIER.md): a compact table
  of the 20 ranked finite families not yet eliminated.
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

The assembled current period checkpoint can be replayed without manually
listing its 40 dependency reports:

```powershell
python replay_ranked_period_conditional_star.py `
  order_pool_1050000_period2533395664800_conditional_star_v3_certificate.json `
  --output replay-period2533395664800.json
```

Successful replay prints `verified=True`. These checks prove only the finite
claims encoded by their certificates.

## Evidence policy

Discovery programs do not certify their own output. Important results are
stored as exact integer or rational data and replayed by a separately written
verifier using a different enumeration orientation or contraction. Finite
checks are always labeled finite; absence of a found cover is never promoted
to a global theorem.
