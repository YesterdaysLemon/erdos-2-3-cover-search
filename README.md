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
- 4,620 of 4,637 ranked finite divisor-period families have independently
  replayed no-cover certificates;
- 17 ranked finite families remain after intersecting the aggregate
  block-star and separately verified single-anchor frontiers and applying
  26 exact period certificates;
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
- period `3139207671600` is eliminated by a 24-anchor block, 165 promoted
  projected-conditional edges, and 1,007 explicit pairwise star-subset
  witnesses. The complete independent replay gives exact upper bound
  `11731653453882582685001/11731696196943206484000`, or
  `0.9999963566171587`. The exact no-cover margin is approximately
  `0.00000364338284134`.
- period `14440355289360` is eliminated by a 19-anchor block, 136
  independently replayed conditional edges, and 1,105 explicit pairwise
  star-subset witnesses. Its complete independent replay gives exact upper
  bound
  `12234784315852076932153/12235161701751295507200`, or
  `0.9999691556263482`. The exact no-cover margin is approximately
  `0.0000308443736518`.
- period `931635825120` is eliminated by a 22-anchor block, 107
  independently replayed conditional edges, and 923 explicit pairwise
  star-subset witnesses. Its complete independent replay gives exact upper
  bound
  `6026778698707389062689/6027397990604327937600`, or
  `0.999897253856821`. The exact no-cover margin is approximately
  `0.000102746143179`.

The fraction \(4{,}620/4{,}637\) measures only this deliberately chosen
finite menu. It is **not** a probability, a completeness claim, or evidence
that the original infinite problem is nearly solved.

## Exploratory direct-cover checkpoint

A separate exact counterexample-guided search now attacks a 14,629-fibre
family directly.  Its component cap is 128 and its raw density sum is
approximately `2.092255659167`.  A cover by this finite family would yield
the requested integer `m`; failure of the search would not prove that no
such integer exists.

The learner evaluates candidate columns on demand and stores enormous exact
checker coordinates as residues modulo the relevant prime-power components.
It also reuses exact cover counts between rounds.  An optional validated
NumPy checkpoint now persists those coordinate residues and counts across
process restarts; candidate, phase, and ordered point-prefix fingerprints
prevent stale counts from being trusted.  Appended lesson batches reuse the
cached prefix.  This reduced the learner from multi-gigabyte dense matrices
to less than 1 GB of private memory during the main experiment and removes
most repeated continuation startup work.

Each exact hole is expanded into the union of four structural tiles: one
900-point Cartesian tile over the lowest binary, ternary, and quinary digits,
plus separate 49-, 121-, and 169-point tiles for components 7, 11, and 13.
Because all four include the original point, their union has 1,236 points
per hole.  Exact witness-diversity constraints make each 100-hole batch
range across additional component residues.

One five-round pass began with 344,000 accumulated lessons.  Its completed
repairs used `0, 98, 64, 62, 107` phase changes at
`344000, 467600, 591200, 714800, 838400` points, respectively.  After every
repair, the separate exact SAT checker returned 100 genuine uncovered
points.

A continuation added separate lowest-digit tiles for components 17, 19, and
23, raising the union to 2,412 points per exact hole.  Repairs at
`962000, 1203200, 1444400` points took `45, 138, 149` phase changes, and the
exact checker again returned 100 holes after each repair.  The final
expansion left 1,685,600 local lessons.

The next repair replaced the singleton 7- and 13-component slices by
Cartesian `(7,11)` and `(11,13)` digit tiles.  It eliminated 146,427 retained
misses with 61 phase changes, but the exact checker returned 10 fresh holes.
Their 28,332-point tile unions leave a 1,968,920-point continuation
checkpoint locally.  The low repair count still shows substantial phase
flexibility, so this branch was paused in favour of the structurally
different perfect-power construction.
Consequently this experiment has found neither a cover nor a no-cover proof.
The method, the earlier max-32 experiment, and the measured null results are
recorded in
[RESEARCH_LOG.md](RESEARCH_LOG.md).

The exact checkpoint and the prominent correction to an earlier invalid
parallel-class argument are recorded in [RESEARCH_LOG.md](RESEARCH_LOG.md).
The 17 unresolved periods are indexed in
[CURRENT_FINITE_FRONTIER.md](CURRENT_FINITE_FRONTIER.md).

## Exploratory perfect-power checkpoint

The strongest structurally different construction writes
`m = M^1616615`.  Algebraic factorization then handles every exponent pair
whose two coordinates are divisible by one of `7,11,13,17,19`.  The current
conditioned pool has 12,577 additional prime fibres and reciprocal-density
sum `7.737722692865`.

An ordinary phase covered 77,988 retained exact lessons.  Four rounds that
varied every single binary--ternary digit pair grew the checkpoint to 357,488
points; each new batch needed only one phase change and each exact checker
still found 25 holes.  Stronger four-digit blocks, varying two binary and two
ternary digits at once, then reached 830,162 points.  Again, each retained
batch needed one phase change and all three exact checks returned five holes.

A separate triple-coverage phase then required every retained lesson to lie
on at least three selected fibres.  Starting from 87,107 lessons, three
four-digit-block rounds grew that checkpoint to 587,225 points.  The new
batches began with `128638` and `131585` points below triple coverage and
were repaired with only three phase changes each; every exact checker still
returned five holes.  The search has therefore found no conditioned-cell
cover.  More point accumulation by itself is no longer the active tactic:
the next experiments impose structural conditions on complete CRT cells.

The first complete-cell condition freezes the 17-adic coordinate and requires
the compatible fibres in every remaining CRT cell to have total line density
at least one.  This is a necessary condition for a cover, not a sufficient
one.  Exact holes supplied 5, 10, 15, and then 115 cumulative density cuts.
All four masters were satisfiable in `1.887`, `3.296`, `5.416`, and `115.325`
seconds, respectively, and each round needed one new full-weight phase move.
The resulting full checker still returned its complete requested batches of
5, 5, 100, and 100 holes, all independently replayed with 17-adic density
below the required threshold.

The last 100-hole batch capped target reuse on the four earlier mover rows,
so it could not be repaired by reusing any of them.  This removed one known
degeneracy but did not remove the holes.  The combined 215-cut master was
also satisfiable, after `863.551` seconds, by one additional change:
the unrestricted modulus-2 row for `p=41` flipped from target 1 to target 0.
Exact integer replay gave minimum scaled density 25 on the 215 cuts, above
the required 17.

The checker then supplied one genuinely uncovered point at a time under a
target cap on every accumulated small mover.  Exhaustive one-row scans
repaired the growing exact core with

```
cuts   moved row   target move   exact replay range
315    p=17, h=4   3 -> 1        23..98
316    p=19, h=3   1 -> 0        23..115
317    p=37, h=3   1 -> 2        30..115
318    p=41, h=2   0 -> 1        25..98
319    p=73, h=3   2 -> 0        35..115
```

Each replay has zero density violations, but every returned phase still has
an exact full-domain hole.  In particular, `p=41` flipped back once the
other movers protected the older half-plane.  This remains an oscillating
necessary-condition search, not a cover.

The binary oscillation was then quotiented out explicitly.  With every other
phase frozen, the checker generated 100 coordinate-diverse holes for each of
the eight target pairs of `(p=41,p=17)`.  Their union with the earlier core
contains 1,121 distinct exact coarse cells.  The exhaustive one-change scan
tested 1,198,544 legal retargets in `130.255` seconds and found none.

An exact two-change scanner replaces the million-variable MILP for the next
radius.  It anchors one failing cell, intersects the possible second-row
target classes over all remaining deficits, and exact-integer replays every
surviving pair.  The first 1,121-cell run found

```
p=41: 0 -> 1
p=107: 31 -> 0
```

in `17.568` seconds, with density range `23..99`.  The full checker still
returned 100 holes.  Ninety-nine had density 12; the remaining hole had
density 18 but its 18 active rows represented only 17 distinct affine lines,
covering 177 of the 289 points in the residual `F_17^2` plane.

The scanner therefore now exact-replays complete residual-plane line unions,
not only their total density.  Eleven further checker/repair batches grew the
joint core to 2,221 density cells and 1,100 residual planes.  Exact
radius-two repairs successively used

```
cells   changed pair
1221    p=41, p=293
1321    p=41, p=439
1421    p=41, p=179
1521    p=41, p=263
1621    p=19, p=97
1721    p=19, p=193
1821    p=19, p=557
1921    p=37, p=97
2021    p=97, p=193
2121    p=97, p=2131
2221    p=277, p=59
```

Every listed pair exactly satisfies the accumulated density and finite-plane
constraints, but every corresponding full-domain checker again produced its
complete 100-hole batch.  The latest batch consists entirely of density
violations, with scaled density `8..11`.  Thus the method has exhausted whole
one- and two-row repair families in succession, but it has not found a cover
or proved that no larger coordinated repair exists.

The checker and master now support
algebraic sublattices, restricted perfect-power targets, unioned cut files,
hard phase-change budgets, and exact integer replay of every accepted MILP
answer.  They also support phase overrides for adversarial branch queries,
exhaustive exact one- and two-change scans, and first-level affine-plane
bitset replay.  No conditioned-cell cover, global `m`, or global
impossibility proof has been obtained.

## Exact five-anchor quotient obstruction

The most conspicuous phase oscillation was also moved into a different
finite space.  The rows

```
p = 41, 17, 19, 37, 73
h =  2,  4,  3,  3,  3
```

have exactly `2*4*3*3*3 = 216` joint target assignments.  Freeze all other
12,572 row phases at the seed-13 checkpoint and regard each exact uncovered
exponent pair as a set of anchor assignments that it defeats.  Of 2,722
candidate points, 2,402 remain outside every frozen nonanchor row and every
algebraic origin sublattice.  Their incidence sets defeat all 216 anchor
assignments.  The resulting self-contained certificate records one witness
per branch, using only 28 distinct witnesses, and has SHA-256
`faba66c4a77cc74cfed7bcfdf75947b0b99967a6176b128e42eaacf4764c6c84`.

The independent verifier does not trust the search files.  It replays every
recorded branch witness against all 12,577 embedded affine rows, checks all
five algebraic sublattices, reconstructs the complete legal target product,
and returns:

```
verified = true
legal_anchor_branches = 216
verified_branch_witnesses = 216
distinct_recorded_witnesses = 28
```

The proof object and replay are:

```
power1616615_anchor41_17_19_37_73_frozen_quotient_certificate.json
power1616615_anchor41_17_19_37_73_frozen_quotient_verification.json
```

This is stronger than a small-radius statement for those five rows, but its
scope is still local: the other phases are frozen.  Freeing 11 additional
rows suggested by the radius-two repairs gives a 16-anchor exact SAT master
with 708 target options.  It solved each of the first three accumulated
masters in about one millisecond, while the first two resulting phases each
returned 100 new exact full-domain holes.

The quotient was then split hierarchically.  Fix the low-anchor branch
`(0,1,0,0,0)` in the displayed order and free only these 11 repair rows:

```
p = 97, 107, 179, 193, 263, 277, 293, 439, 557, 2131, 59.
```

Their legal target space has exactly
`693636364523341088` assignments.  Coordinate-diverse counterexamples left
the master SAT through 8,300 point clauses, but they were highly redundant
as phase clauses.  Asking the exact checker for 100 distinct 11-anchor target
fingerprints made the 8,400-point master UNSAT.  Independent Z3 core
extraction reduced this to 258 points.  PySAT and the separately encoded Z3
integer master both return UNSAT on that compact core; Z3 uses 269 assertions.

The self-contained proof object and replay are:

```
power1616615_anchor11_lowbranch_0_1_0_0_0_unsat_certificate.json
power1616615_anchor11_lowbranch_0_1_0_0_0_unsat_verification.json
```

The certificate SHA-256 is
`7330119c171612c05fd634b8b51686e3d3156ad259f1c08dcafb7b2e29501f3f`.
Thus this one low-anchor branch remains impossible even after all 11 repair
rows are freed.  The other 215 low-anchor branches and all rows outside this
16-coordinate quotient are not eliminated by this result.

The hierarchical audit has now accumulated 4,044 distinct joint
low/repair fingerprints.  An exact parity relation reduces the 216 nominal
low branches to 162 genuinely different residual systems:

```
p=41 target = l (mod 2)
p=17 target = 2k + 3l (mod 4), hence target = l (mod 2).
```

After fixing a `p=41` phase, every surviving point has the opposite
`l`-parity.  The two `p=17` phases with the fixed `p=41` parity can therefore
never cover a survivor and are indistinguishable.  The reusable
`audit_two_level_anchor_quotient.py` checks this coefficient relation and
audits all branches with one assumption-based SAT model.  On the current
finite fingerprint corpus, 66/216 raw branches, or 43/162 quotient classes,
are UNSAT.  Except for the independently certified branch above, these are
exploratory finite-corpus results, not new proof certificates.

A construction-oriented continuation freed three additional repair rows,
`p=4297,653,577`, and switched the SAT master to minimum-Hamming-change
MaxSAT.  The 14-anchor master covers 440 accumulated exact lessons.  Its
last three 10-point repairs required only 1, 1, and 2 phase changes, instead
of the 12--14 changes made by unconstrained solves.  Every resulting exact
checker nevertheless returned all ten requested holes.  This controlled
column-generation path is still active but has not produced a conditioned
cover, an integer `m`, or a global obstruction.

The continuation now also has a sparse full-pool repair path.  Stable MaxSAT
column generation promoted six further rows, giving a 20-anchor phase that
covers 526 accumulated exact points.  Checker-side target caps produced
finite packs of four mutually target-distinct holes at cap one and eight at
cap two after excluding the tiny `h < 7` anchors.  These are packing
statements, not bounds on the total number of holes.

Direct MILP repair of the growing pack repeatedly found four-row repairs
with `p=41,p=17` frozen, but became solver-hard at 29 supplied points.
The replacement `finite_sample_mask_repair.py` projects the current base
misses to bit masks, enumerates a tiny mask cover, assigns distinct source
rows, and exact-replays every supplied point.  On the 29-point instance,
57,454 deficit-hitting moves collapsed to 198 masks; after six mask models
the first exact finite repair was found in about 0.003 seconds.  Subsequent
single-hole rounds grew the branch corpus to 38 points.  Every repaired
phase checked so far still has a genuine full-domain hole.

The mask engine is intentionally one-sided.  A returned phase is an exact
finite repair, but exhaustion does not exclude repairs that use a row with
zero gain on the base-missed points solely to compensate for a loss
elsewhere.  `finite_sample_sat_repair.py` and
`finite_sample_z3_repair.py` retain complete bounded-repair encodings for
that purpose; their current 29-point runs timed out rather than proving
UNSAT.  A resumable `sparse_anchor_quotient_sweep.py` enumerates all 162
low-anchor quotient representatives.  The first cross-branch sample also
showed why each branch needs its own exact adversary: the phase-one corpus is
already covered by the opposite `p=41` anchor in phase-zero branches.

An experimental 62-row extension with prime-power components above 16,384
was admitted by widening finite target matrices to 64 bits.  None of those
rows covered any of the latest three adversarial holes at its base phase,
and no one retargeting covered more than one.  A full Z3 bit-vector check
reached a five-minute cap without a model or proof.  The extension is
therefore not part of the active exact-checker loop.

## Exact small affine-subpool obstruction

As a deliberately different strategy, the 14 derived rows with modulus
`h <= 12` were treated as a complete affine-cover problem on their
period-5,544 torus.  Exact CEGIS reached an UNSAT master on 13,908 explicit
witness points.  A separately written verifier then enumerated all 746,496
legal phase assignments using exact bitset unions and found that none covers
those points.  It also reconstructs every row from its prime, multiplicative
orders, perfect-power restriction, and period-60 coordinate change, then
checks that every witness avoids the five declared algebraic sublattices.

The self-contained proof object and replay are:

```
power1616615_hle12_structured_noncover_certificate.json
power1616615_hle12_structured_noncover_verification.json
```

This proves a finite no-cover statement for exactly the embedded 14 rows.
It rules out that attractive low-modulus affine template, but says nothing
by itself about the other 12,563 rows, the roughly 1,711 neighboring
period-60 cells needed for a global construction, or the original infinite
problem.

## Bounded homogeneous obstructions

The complete raw pool through subgroup index `1050000` contains 129,497 leaf
lattices.  The explicit exponent pair

```
(k,l) = (552897490806962158, 6049004616530593493)
```

lies in none of them.  The discovery check uses the stored signature
equations, while an independent replay recomputes
`2^k * 3^l mod p` for every source prime and returns `verified=True`.
Consequently, even the union of every homogeneous fibre in this bounded pool
is not a cover.

A separate refinement-specific certificate groups all 429,197 possible
parent reductions and finds no complete
set of prime-\(q\) siblings, and an independently oriented replay returns
`verified=True`.  Every nontrivial finite refinement tree has a deepest
internal node whose children are all leaves, so no homogeneous refinement of
the trivial cover can use only this finite pool.

These are bounded finite obstructions.  The explicit witness rules out every
homogeneous subfamily of this pool, including non-refinement families, but
does not exclude homogeneous fibres beyond the bound, general affine covers,
or fibres of higher index.

## Repository contents

- `*.py`: search programs, exact certificate generators, independent
  verifiers, and regression tests.
- `local_phase_cegis.py` and `exact_uncovered.py`: direct finite-cover CEGIS,
  including streaming low-memory repair, unions of exact component-digit
  tiles (including multiple digits of one prime component), and adversarial
  witness diversity.
- `finite_sample_mask_repair.py`: fast gain-mask construction search with
  exact full-corpus replay of every emitted phase.
- `finite_sample_sat_repair.py` and `finite_sample_z3_repair.py`: complete
  bounded-change finite repair encodings in SAT and pseudo-Boolean spaces.
- `sparse_anchor_quotient_sweep.py`: resumable driver over the 162 exact
  low-anchor quotient representatives.
- `certify_homogeneous_refinement_obstruction.py` and its independent
  verifier: the deepest-node sibling test for bounded homogeneous
  refinement trees.
- `certify_bounded_homogeneous_noncover.py` and its modular-exponentiation
  verifier: an explicit point outside the union of all bounded homogeneous
  fibres.
- `certify_small_derived_pool_noncover.py` and
  `verify_small_derived_pool_noncover.py`: a self-contained 14-row affine
  no-cover certificate and exhaustive phase replay.
- `search_projected_conditional_designs.py`: checkpointable projected
  conditional-overlap discovery, including a proof-structural shared
  residual pair.
- `search_mitm_pairwise_conditional_bounds.py`: fast full-anchor subset
  selection whose chosen subset is reevaluated with exact rational
  arithmetic.
- `search_sixterm_conditional_bounds.py`: discovery-only sixth-order
  Bonferroni screening.
- `promote_projected_conditional_designs.py`: resumable regeneration and
  independent replay of selected discovery designs.
- `certify_pairwise_star_witnesses.py` and
  `verify_pairwise_star_witnesses.py`: proof-producing replay of explicit
  star subsets without an exhaustive `2^20` optimization.
- `order_pool_1050000_component_core_corrected_max32_stable.json`: stable
  finite candidate pool used by the current frontier.
- `order_pool_1050000_component_core_corrected_max128.json`: 14,629-row
  candidate pool used by the current exploratory direct-cover search.
- `order_pool_1050000_period3139207671600_conditional_fibre6553_paired_`
  `autodesign_{certificate,verification}.json`: permanent proof object and
  independent replay for the strongest paired conditional edge in the
  now-eliminated family.
- `order_pool_1050000_period3139207671600_conditional_star_`
  `{certificate,verification}.json`: assembled period proof and independent
  replay; its promoted conditional manifest and pairwise-star witness
  certificate identify the complete dependency closure.
- `order_pool_1050000_period14440355289360_conditional_star_`
  `{certificate,verification}.json`: assembled period proof and independent
  replay for the newly eliminated 1,124-row family; its block chain,
  promoted conditional manifest, and 1,105 pairwise-star witnesses are
  included.
- `order_pool_1050000_period931635825120_conditional_star_`
  `{certificate,verification}.json`: assembled period proof and independent
  replay for the eliminated 945-row family; its 22-anchor block chain,
  promoted conditional manifest, and 923 pairwise-star witnesses are
  included.
- `order_pool_1050000_max32_period_divisor_family_ranking_all.json`: the
  4,637 ranked divisor-period families.
- `order_pool_1050000_max32_all4637_pairanchor_star_scan_v35.json`: aggregate
  finite scan.
- `order_pool_1050000_max32_pairanchor_survivors_singleanchor_star_scan_v26.json`
  and its verification: the independent single-anchor frontier.
- the `order_pool_1050000` exact rational certificates and replay reports
  needed by the current checkpoint and its dependency closure.
- [CURRENT_FINITE_FRONTIER.md](CURRENT_FINITE_FRONTIER.md): a compact table
  of the 17 ranked finite families not yet eliminated.
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

The bounded homogeneous non-cover witness can be replayed directly from the
complete source pool:

```powershell
python verify_bounded_homogeneous_noncover.py `
  order_pool_1050000_bounded_homogeneous_noncover_certificate.json `
  --pool order_pool_1050000.json `
  --output replay-bounded-homogeneous.json
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

The small affine-subpool obstruction can be replayed without its discovery
pool:

```powershell
python verify_small_derived_pool_noncover.py `
  power1616615_hle12_structured_noncover_certificate.json `
  --output replay-small-derived-pool.json
```

## Evidence policy

Discovery programs do not certify their own output. Important results are
stored as exact integer or rational data and replayed by a separately written
verifier using a different enumeration orientation or contraction. Finite
checks are always labeled finite; absence of a found cover is never promoted
to a global theorem.
