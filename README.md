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

As of 2026-07-28:

- no candidate \(m\) has been found;
- no global impossibility theorem has been proved;
- one class-27 neighborhood has a packaged radius-five Hamming-ball
  obstruction candidate: an exhaustive 34-node, 130-leaf phase tree over
  12,579 rows records `repair_exists=false`, but its independent full replay
  is still pending and the claim is not yet promoted;
- a determinant-one shear along exponent direction `(3,1)` produced a
  99-row, period-15,120 layered candidate whose weakest exact column capacity
  is about `1.10120`, but new exact certificates rule out both the chosen
  placement and every other active-class placement of the same 99-row
  family;
- a wider shear `(1,-3)` produced an independently replayed 81-row,
  period-110,880 placement with weakest capacity about `1.09784` and no
  forced coprime-pair violation; an exact modulo-30 weighted projection then
  ruled out that placement's weakest column in all 30 anchor-phase branches;
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

The mask engine now closes the secondary-loss gap explicitly.  After
selecting one distinct source row for every gain mask, it recursively tries
all remaining zero-gain or duplicate-mask moves that could cover a newly
created deficit.  A 40-point subcore in low-anchor branch `(1,1,0,1,0)` is
therefore a complete radius-four obstruction.  It has 25 base misses,
76,239 deficit-hitting moves, 660 gain masks, and relaxed mask-cover number
three.  An independent scalar verifier enumerates 416 full mask skeletons,
rejects 415 by row ownership and the last by exact replay, and returns
`verified=true` with no repair.

The proof object and replay are:

```
power1616615_lowbranch_1_1_0_1_0_radius4_obstruction_certificate.json
power1616615_lowbranch_1_1_0_1_0_radius4_obstruction_verification.json
```

The certificate SHA-256 is
`9531f57cef120ed007eb5816ebaaa3cdbd432d47adcab1eda4acc97fbecf9c45`.
This rules out only four-row repairs around the declared base phase with all
five low anchors fixed.  Radius five, distant full-pool assignments, the
other 161 quotient classes, and the original problem remain open.

`finite_sample_sat_repair.py` and `finite_sample_z3_repair.py` retain
independent complete bounded-repair encodings; their current larger runs
timed out rather than contributing another result.  A resumable
`sparse_anchor_quotient_sweep.py` enumerates all 162 low-anchor quotient
representatives.  The first cross-branch sample also showed why each branch
needs its own exact adversary: the phase-one corpus is already covered by
the opposite `p=41` anchor in phase-zero branches.

An experimental 62-row extension with prime-power components above 16,384
was admitted by widening finite target matrices to 64 bits.  None of those
rows covered any of the latest three adversarial holes at its base phase,
and no one retargeting covered more than one.  A full Z3 bit-vector check
reached a five-minute cap without a model or proof.  The extension is
therefore not part of the active exact-checker loop.

That generic negative result did not hold uniformly across the quotient.
`select_targeted_signature_rows.py` now scores an extension row only when
several branch-specific exact holes share one legal target, and it can reject
large prime-power components before they enter the checker.  In quotient
class 27, the row

```
p=63700993, h=884736, largest component=32768
```

hits two of six accumulated branch holes at target 763570.  In class 54,
`p=19131877`, with `h=177147`, similarly hits two of six at target 84275.
Starting those rows at the selected targets reduced both finite masters from
four ordinary phase changes to three.  The resulting radius-three searches
remain feasible on 54 and 44 accumulated points, respectively.

The large component is kept out of the expensive exact SAT encoding.
`verify_layered_exact_misses.py` instead takes explicit holes from the stable
12,577-row checker and scalar-replays every augmented row.  In the latest
audited batches, four of ten base holes in class 27 and all ten in class 54
remained exact augmented-pool holes.  This proves that the checked phases are
not covers while preserving a cheap base checker.  It does not prove that
either targeted row belongs to a global construction, and the newest finite
repairs have not produced an integer `m`.

A recursive scan of the enlarged verified hole corpus found a stronger row,

```
p=11337409, h=157464, largest component=19683.
```

Its selected target hit 101/108 training holes in class 27 and 70/117 in
class 54.  Large independent batches showed that this was residue-family
leverage rather than a cover: one later class-27 phase intercepted 0/100
fresh base holes, and a changed-row-diverse check retained another exact
augmented miss.  Even so, adding this one row reduced the class-27 finite
repair to two ordinary phase changes and it exactly covers 252 accumulated
points.

Class 54 moved in the opposite direction.  On its 151-point corpus the
complete search found no repair within two changes of the augmented base
phase, even when either targeted row is allowed to move.  The independent
scalar verifier needs no mask skeleton enumeration because the relaxed
two-mask cover is already impossible.  The proof object and replay are:

```
power1616615_lowbranch_0_3_0_0_0_targeted_radius2_obstruction_certificate.json
power1616615_lowbranch_0_3_0_0_0_targeted_radius2_obstruction_verification.json
```

The certificate SHA-256 is
`23f14c18a283c4b03471673d8503a712006a57e9171053f2d83b5d28ee4fa17d`.
This is another finite Hamming-ball theorem, not a full quotient-class
obstruction.  The class-27 two-change phase still has an exact full-domain
hole, and neither branch yields an integer `m`.

The selector now accepts a separately declared validation corpus.  A target
is chosen from training points only and must hit a requested number of held-
out points at that same target before promotion.  This rejected a tempting
third row, `p=3439853569`, which hit 52/101 training points but 0/2
changed-row-diverse validation holes.

The audit also exposed a checker artifact: one 100-hole batch occupied a
single residue pair modulo almost every tested small prime.  Three subsequent
checker runs therefore forced different coordinate cells modulo 127, 113,
and 109.  Both promoted rows intercepted 0/10 holes in each cross-section.
After all 30 exact augmented misses were added, class 27 nevertheless retained
a two-change finite repair covering 281 accumulated points, now using
`p=233` and `p=823`.  This is robust finite repair flexibility; the promoted
rows have not passed out-of-sample validation and the phase is not accepted
as a construction.

`exact_common_phase_misses.py` strengthens the adversary in a different
space.  It unions every distinct fibre selected by several saved phase maps
and asks for one exact point outside that union.  Because recent repairs
differ on only two rows, four historical phases enlarged the stable checker
from 12,577 to only 12,585 distinct fibres.  Every returned point is scalar-
replayed and is simultaneously missed by all supplied phases.

This produced one exact augmented point common to two early radius-two
repairs.  Later 100-point pair and triple batches eliminated several repair
responses at once.  The master nevertheless remains feasible at radius two
on 582 accumulated points, most recently by changing `p=853` and `p=12841`.
The promoted rows sometimes intercept almost every common base-pool hole;
for the latest pair they caught 98/100, leaving two exact augmented common
misses.  Thus common-phase cuts improve Benders efficiency but have not
produced a cover or a full radius-two obstruction.

Coordinate-diverse common cuts ultimately closed class 27 at radius two.
The 741-point discovery corpus was already UNSAT in its relaxed mask layer.
`minimize_relaxed_radius_obstruction.py` guarded every base-miss clause by a
SAT assumption, extracted an UNSAT core, greedily deleted redundant
assumptions, and complete-replayed the result.  Only nine points are needed:

```
base misses                 9
legal deficit-hitting moves 33349
distinct gain masks         110
relaxed two-mask cover      none
```

The proof object and independent scalar replay are:

```
power1616615_lowbranch_0_1_0_0_0_targeted_radius2_obstruction_certificate.json
power1616615_lowbranch_0_1_0_0_0_targeted_radius2_obstruction_verification.json
```

Certificate SHA-256:
`ce270923ce476d05dcacd6a2e1f8bd94c7e3ad68d8bc49cc082f22b094dc6cc8`.
Only the five low anchors are fixed; either promoted row may be among the
two changed rows.  This rules out a radius-two neighborhood of one augmented
base phase.  Radius three, distant assignments, the rest of quotient class
27, the other 161 classes, and the original problem remain open.

The relaxed obstruction has also been compressed into a small hypergraph.
All 33,349 legal one-row moves induce 110 distinct gain masks on the nine
points, but only 14 masks are inclusion-maximal.  Exhaustively unioning the
105 unordered pairs of maximal masks gives maximum cardinality eight, so
every pair omits at least one of the nine witnesses.  This is an exact
finite combinatorial proof of the relaxed radius-two obstruction, independent
of row ownership and secondary losses.  A fractional point-weight LP lands
on equality and is therefore not used as the proof.

The hypergraph proof object and separately implemented scalar reconstruction
are:

```
power1616615_lowbranch_0_1_0_0_0_targeted_radius2_hypergraph_dual_certificate.json
power1616615_lowbranch_0_1_0_0_0_targeted_radius2_hypergraph_dual_verification.json
```

The hypergraph certificate SHA-256 is
`d281ad38fdb8777f00bdce4bed78ef9d1ac582905ba18ca49870c016f4254aa4`.
The verifier reconstructs all 110 masks directly from scalar affine
congruences, reduces them again to the same 14 maximal masks, and confirms
maximum pair union eight.

Radius three remains feasible on the tracked adversarial corpus.  The first
three-change response covered all 741 points by retargeting
`p=15121,223,569`.  A direct exact checker over the augmented component
domain found a genuine full-domain hole.  After that cut, the response
swapped `p=223` for `p=1777`; one exact point was then found outside both
responses simultaneously.  On the resulting 743-point corpus, the current
finite response retargets `p=15121,709,569` and has exact sampled minimum
coverage one.  It has not passed a full-domain cover check and is not a
candidate integer `m`.  This oscillation motivates common-phase cuts and a
quotient-wide symbolic engine rather than unbounded single-hole sampling.

A single augmented-domain query against those three responses then returned
ten coordinate-diverse exact holes common to all three.  Adding the complete
batch gives a 753-point corpus.  Radius three is still feasible there.
`enumerate_relaxed_mask_covers.py` enumerated 3,851 distinct relaxed
three-mask covers in 60 seconds before its declared time limit and found no
mask common even to that incomplete sample.  Thus the apparent two-move
backbone was solver-order bias, not a valid structural conclusion.

The mask repair engine now accepts previous phases and a required Hamming
distance.  On the same 753 points it produced six exact finite repairs at
pairwise distance six, the maximum possible for radius-three phases.  Their
18 changed rows are pairwise disjoint.  This creates deliberately separated
responses for common-phase exact cuts and avoids spending each expensive
checker call on a cosmetic one-row oscillation.  It also shows that this
finite radius-three cover space remains broad; it is not evidence of a
full-domain cover.

The exact augmented-domain checker then unioned five of those separated
phases.  Their 15 alternative fibres enlarged the 12,579-row base assignment
to 12,594 distinct selected fibres.  One checker call returned ten
coordinate-diverse points missed by all five phases, and scalar replay
confirmed every point.  Adding the batch gives 763 accumulated points.
Radius three remains feasible: the current finite response retargets
`p=15121,269,167`, and an independent scalar audit reports sampled minimum
coverage one.  Thus diversified common cuts are more efficient, but they
had not yet closed radius three at that checkpoint.

Subsequent branch exhaustion closes radius three for this finite class-27
Hamming ball.  On 773 points, a complete relaxed enumeration found 3,921
three-mask skeletons.  Expanding every mask to its legal row owners gave
14,302 raw products, 9,625 distinct-row trials, and 9,468 exact finite
repairs; 157 owner assignments failed full-corpus replay.  Of the exact
repairs, 9,462 use `p=97: 1 -> 2`.  The remaining six all use
`p=109: 2 -> 6`, `p=4177: 28 -> 9`, and one of six explicitly replayed
third moves.

After two exact 100-hole batches were merged, all six exceptional repairs
failed the 973-point corpus.  Fixing the dominant `p=97` move reduces the
remaining problem to a two-change residual search.  Its complete relaxed
enumeration has 26 skeletons and exactly 89 legal owner assignments, all of
which pass finite replay.  Those 89 repairs use 50 distinct alternative
fibres.  One exact augmented-domain query against their complete fibre union
returned 100 exponent pairs missed by every repair, with independent scalar
replay.

Adding that batch gives 1,073 points and an exact radius-three UNSAT result:

```
initial base misses          869
legal deficit-hitting moves  210047
distinct gain masks          15835
relaxed three-mask skeletons 14
distinct-row matching fails  14
exact repairs                0
```

The separately implemented scalar verifier reconstructs the 12,579 rows,
all 1,073 points, all 210,047 legal moves, and all 14 terminal skeletons.  It
returns `verified=true` and `repair_exists=false`:

```
power1616615_lowbranch_0_1_0_0_0_targeted_radius3_1073_sparse_obstruction_certificate.json
power1616615_lowbranch_0_1_0_0_0_targeted_radius3_1073_sparse_obstruction_verification.json
```

Certificate SHA-256:
`38fa3e2774e84f0c621c612100c7fdaa4a6c713b144dc7f5c3204e31837de26d`.
This is an exact finite Hamming-ball obstruction with the five low anchors
fixed.

Radius four has since been closed on a 1,373-point extension.  A monolithic
five-minute mask search examined only eight terminal skeletons and a separate
600-second Z3 encoding returned `TIME_LIMIT`, so neither contributed a
negative result.  The successful proof instead partitions the Hamming ball
by the complete legal phase sets of three high-leverage rows:

```
p=97,  h=4
  base target 1 -> partition p=109, h=9
    base target 2 -> partition p=193, h=8
```

The other three `p=97` targets leave four terminal relaxed skeletons in
total, all rejected by distinct-row matching.  Eight alternate `p=109`
targets and seven alternate `p=193` targets have 8--10 and 7--9-point
relaxed UNSAT cores respectively.  With all three rows at their base targets,
an 11-point relaxed core closes the final branch.  Thus the complete decision
tree has 19 leaves.  The combined certificate embeds only 3,377 point
appearances across those leaves rather than repeating the full corpus.

The independently implemented scalar verifier authenticates the pool and
base phase, checks that every legal target occurs exactly once in the tree,
derives each leaf phase and remaining Hamming budget from its path, and then
reconstructs all sparse-radius searches.  It returns:

```
verified=true
repair_exists=false
partition_count=3
leaf_count=19
total_full_skeleton_count=4
elapsed_seconds=437.196
```

The public proof objects and the reproducible discovery corpus are:

```
power1616615_lowbranch_0_1_0_0_0_targeted_radius4_1373_partitioned_obstruction_certificate.json
power1616615_lowbranch_0_1_0_0_0_targeted_radius4_1373_partitioned_obstruction_verification.json
power1616615_lowbranch_0_1_0_0_0_targeted_radius4_1373_points.json
```

Certificate SHA-256:
`010807fa45cf948fe1f79b49713543347266986ca42589fc5d8201afe8126ed7`.
This is still only a finite Hamming-ball obstruction.  It does not rule out
radius five, distant assignments, other quotient branches, fibres outside
the bounded pool, or the original infinite problem.

Radius five initially remained feasible on the same corpus.  Quotienting
first by `p=97:1 -> 2` found a five-change response in 10.6 seconds, with the
changes

```
p=97:1->2, p=109:2->6, p=193:1->5, p=433:12->0, p=577:7->1.
```

A 200,000-draw reproducible audit found 19 uncovered points among 191,831
algebraically eligible draws.  More importantly, the complete-domain exact
checker returned ten genuine holes, and a separate scalar affine replay
gave coverage zero on all ten.  Those counterexamples started a 1,392-point
adversarial continuation.

An exhaustive radius-five candidate certificate has now been packaged on
that continuation.  The monolithic search was replaced by a complete phase
tree with 34 partition nodes using 15 distinct primes.  Its 130 leaves
comprise 80 residual radius-four
obstructions, 49 residual radius-three obstructions, and one residual
radius-five obstruction.  Compact leaf cores embed 16,603 point appearances
instead of repeating the full 1,392-point corpus at every leaf.

The deepest base branch ends in a 31-point core.  An exact integer-weight
dual has total weight 410 and maximum one-mask weight 82, so any relaxed
five-mask cover must be an equality case.  Equality analysis isolates one
five-mask skeleton; the complete distinct-row replay rejects it because its
owners cannot be assigned to five different rows.  This illustrates why the
dual equality alone was not reported as a proof.  The ordinary sparse-radius
certificate supplies the actual leaf proof.

The independently implemented scalar tree verifier authenticates the source
pool and base phase, reconstructs every legal target branch, derives each
leaf phase and remaining budget from its path, and replays every embedded
sparse obstruction.  The first full replay was stopped after 5.2 hours to
preserve the host's 15% free-memory floor, so no `verified=true` report
exists yet.  The candidate certificate itself records:

```
repair_exists=false
row_count=12579
max_changes=5
partition_count=34
leaf_count=130
total_leaf_point_count=16603
total_full_skeleton_count=28
```

The public artifacts are:

```
power1616615_lowbranch_0_1_0_0_0_targeted_radius5_1392_partitioned_obstruction_certificate.json
power1616615_lowbranch_0_1_0_0_0_targeted_radius5_1392_points.json
```

Certificate SHA-256:
`b578027be29cabec886031c0a297e4e5150c182f7c23cac712f264ad563d5feb`.
If the independent replay succeeds, this will close only the declared
radius-five Hamming ball with the five low anchors fixed.  Until then it is
public work in progress, not a certified result.  Either way, it is neither
a value of \(m\) nor a global nonexistence proof.

## Unimodular layered construction pivot and obstruction

A separate constructive route changes lattice coordinates before selecting
phases.  The determinant-one basis

```
(k, l) = x * (3, 1) + y * (-1, 0)
```

means `k=3x-y`, `l=x`; advancing `x` multiplies by
\(2^3\cdot3=24\).  Restricting the bounded source pool to rows whose
base-24 order divides 15,120 leaves 99 primes.  A max-min binary allocation
chooses the active `x` class of each row.  Exact rational replay gives:

```
raw reciprocal density = 1.1861303234043865
weakest column capacity = 1.1011988792930392
minimum multiplicity    = 18
forced coprime-pair violations = 0
```

The verifier independently reconstructs the basis transform from all 129,497
source rows, confirms exactly which 99 rows qualify, authenticates every
target restriction, replays all 15,120 column capacities, and repeats the
complete pair screen.  It returns `verified=true`.  These are necessary
preprocessing facts, not a cover: the remaining residual phase of each prime
is globally coupled across every column in which that row is active.

The reproducible candidate and its independent replay are:

```
order_pool_1050000_direction3_1_period15120_layered_pool.json
order_pool_1050000_direction3_1_period15120_layered_pool_verification.json
```

Candidate SHA-256:
`13cc7468117b88256a9e4158200da56602e529cc708993478eadff8abf19302c`.

The capacity and pair screens were not sufficient.  A stronger
covering-system argument now rules out the placement before phase synthesis.
For any prime \(q\), if fewer than \(q\) residue-class moduli are divisible
by \(q\), those classes are redundant in a cover of the integers.  To see
this, hold all prime-to-\(q\) coordinates fixed and vary through \(q\)
translates: every \(q\)-free class has constant truth value, while each
\(q\)-divisible class can hit at most one translate.

In weakest column 454, 78 residual classes are active but only 6 of their
moduli are divisible by 11.  Removing those six classes leaves 72 classes
with exact reciprocal density

```
1331626919328932373170013200623287663407217762378152211484151
----------------------------------------------------------------
1332735606562889452651421880641636334058989019743327654591360
```

or `0.9991681116430765`, which is strictly below one.  Those 72 classes
cannot cover even by the union bound, so the independently phased column
relaxation is impossible.  A fortiori, no globally coupled choice of the 99
row phases covers this layered placement.

The proof object and separately implemented exact-rational replay are:

```
order_pool_1050000_direction3_1_period15120_layered_noncover_certificate.json
order_pool_1050000_direction3_1_period15120_layered_noncover_verification.json
```

Their SHA-256 hashes are respectively
`5d88f3db107a0411bf3433bdf41b9e3b2c8c739056d502ef27c787505b405e1c`
and
`2778bfc99cc7158d1ec48bbf2406c8d1223bcb8f95ba294f2783f86faebdc6e6`.

The same lemma also gives a placement-independent reduction.  Repeatedly
remove every row containing the smallest residual prime \(q\) that occurs in
fewer than \(q\) current residual moduli.  If any active-class placement of
the 99 rows covered every column, every one of these removals would preserve
the cover column by column.  Twenty-eight exact pruning rounds remove 48
rows, leaving 51 rows with raw reciprocal density
`1.0627209309649785`.

Of those survivors, 27 are active in every column and 24 choose one residue
class modulo their active-coordinate modulus.  A reduced multi-valued
decision diagram uses exact integer weights with common scale 6,350,400 and
exhausts all choices on a 68-column core modulo 2,520.  The independent
replay constructs 1,152,581 MDD nodes including its two terminals and reaches
terminal zero.  Therefore no active-class allocation of this entire 99-row
family can even satisfy the necessary per-column capacity inequalities.

The stronger proof object and independent replay are:

```
order_pool_1050000_direction3_1_period15120_layered_family_noncover_certificate.json
order_pool_1050000_direction3_1_period15120_layered_family_noncover_verification.json
```

Their SHA-256 hashes are respectively
`b83eb87058d58f403b7d60b8d3896bed78146f47a4330e70c775737728b35df6`
and
`6f21d8bf8032a1f6dbc6929162b8540c08ce79c7f2267fa1a621192d35fb0ed6`.
The point-only CEGIS checkpoints are therefore superseded.  Future layered
families will receive prime-deficit pruning and an exact capacity-placement
test before any BDD over residual phases, SAT, or pointwise search.  Other
shears, periods, and larger source pools remain open; neither obstruction is
a global nonexistence proof.

## Wider shear and exact modulo-30 projection obstruction

A wider scan applied the proof-safe prime-deficit reduction before allocating
active layer classes.  The first candidate to survive both exact capacity
and the complete unavoidable-pair screen uses

```
(k, l) = x * (1, -3) + y * (0, 1),
```

so `k=x` and `l=-3x+y`.  At layer period 110,880, 140 source rows qualify.
Thirty-three prime-deficit rounds remove 59 rows.  A max-min allocation of
the 81 survivors has capacity pattern period 9,240 and exact weakest capacity

```
140596447 / 128066400 = 1.0978402375642635,
```

attained in 2,640 columns of the declared full period.  The independent
verifier reconstructs the shear and source selection, replays every pruning
round and target restriction, checks all exact column capacities, and finds
zero forced coprime-pair violations.

That promising necessary-condition survivor still does not cover.  Column
24 contains 56 residual classes.  Project their transverse coordinate
modulo 30.  The `p=7` and `p=11` classes occupy complete cells modulo 6 and
5; translation symmetry fixes both targets at zero.  The complete-cell
`p=31` class then has 30 possible targets modulo 30.  For every target, an
exact nonnegative integer weighting of the remaining projection cells bounds
the maximum total contribution of all 53 tail rows strictly below the
required weight.  Across the four distinct exact ratios, the weakest gap
still has tail/required ratio

```
42580991 / 51226560 = 0.8312287805388455 < 1.
```

The separately implemented verifier reconstructs all 56 rows and replays
all 30 weighted branches using rational arithmetic.  It returns
`verified=true` and `proved_no_declared_layered_cover=true`.  Thus the
modulo-30 projection replaces an already slow point-by-point column CEGIS
with a small finite proof.

The pool, its necessary-condition replay, and the projection obstruction are:

```
order_pool_1050000_direction1_neg3_period110880_pruned_layered_pool.json
order_pool_1050000_direction1_neg3_period110880_pruned_layered_pool_verification.json
order_pool_1050000_direction1_neg3_period110880_pruned_layered_projection_noncover_certificate.json
order_pool_1050000_direction1_neg3_period110880_pruned_layered_projection_noncover_verification.json
```

Their SHA-256 hashes in the same order are:

```
88a3b076961bdb7f9ac281afc21adb43d6193778b8c12f717d280c4b929c27cb
e60856156a60e911dcd85e85252582b2791e83531876455bc5909832ae197e4a
8d94523c81a9c98d0bca87182d108dbec2c9bf46667edc8135315a8b3dc0f6f9
d4f24ef8d6d0b6b5f6c27436da29f948c138ec983264ea1540c16b00b7d55a53
```

This obstruction rules out only the chosen 81-row placement.  It motivates
the next allocator: enforce projected-cell feasibility while choosing the
active classes, then pass surviving placements to phase-aware column BDDs or
CRT-component transfer matrices.

An experimental quantified-SMT encoding moved phase targets outside a
universal pair of integer exponent coordinates.  It correctly handles tiny
parity examples, but Z3 returned `UNKNOWN` after 120 seconds on the already
certified 14-row, period-5,544 no-cover instance.  That admission test is
recorded in `power1616615_hle12_quantified_z3_result.json`; the encoding is
not being scaled to the 12,579-row branch.

The sparse-radius proof dependencies were also promoted from ignored
`_tmp` files to permanent tracked pool and base-phase artifacts.  The
standalone and partitioned sparse-radius verifiers read only files present
in a fresh clone.

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
  exact full-corpus replay, complete bounded secondary-loss search, and
  optional Hamming separation from previous repair phases.
- `certify_sparse_radius_obstruction.py` and its independent scalar verifier:
  self-contained finite Hamming-ball obstruction certificates.
- `finite_sample_sat_repair.py` and `finite_sample_z3_repair.py`: complete
  bounded-change finite repair encodings in SAT and pseudo-Boolean spaces.
- `sparse_anchor_quotient_sweep.py`: resumable driver over the 162 exact
  low-anchor quotient representatives.
- `select_targeted_signature_rows.py`: verification-cost- and
  validation-aware promotion of extension rows whose one training-selected
  target also hits separately declared held-out points.
- `verify_layered_exact_misses.py`: independent scalar replay of base-pool
  holes against a small targeted extension.
- `exact_common_phase_misses.py`: exact counterexamples shared by several
  saved phase assignments, with distinct-fibre union, optional coordinate
  diversity, and scalar replay.
- `minimize_relaxed_radius_obstruction.py`: assumption-core extraction and
  complete replay for obstructions already UNSAT in the relaxed mask layer.
- `analyze_relaxed_mask_dual.py` and `verify_relaxed_mask_dual.py`: exact
  compression of a relaxed obstruction to inclusion-maximal gain masks and
  independent scalar replay of the resulting hypergraph certificate,
  including exact handling of tight dual-equality cases.
- `build_axis_layered_pool.py` and `verify_axis_layered_pool.py`:
  max-min layered-pool discovery in ordinary or unimodularly sheared
  coordinates, followed by independent exact source-selection, capacity,
  and unavoidable-pair replay.
- `certify_axis_layered_column_noncover.py` and
  `verify_axis_layered_column_noncover.py`: a complete column scan for the
  prime-deficit covering-system obstruction and a separately implemented
  exact-rational replay.
- `certify_axis_layered_family_noncover.py` and
  `verify_axis_layered_family_noncover.py`: iterative proof-safe
  prime-deficit pruning followed by an exact reduced MDD over every remaining
  active-class placement.
- `materialize_axis_layered_column.py`: exact one-dimensional extraction of
  a selected layered column for structural or adversarial checks.
- `certify_axis_layered_projection_noncover.py` and
  `verify_axis_layered_projection_noncover.py`: branch on translated
  full-cell anchors and replay exact weighted projection duals.
- `enumerate_relaxed_mask_covers.py`: complete-or-explicitly-limited
  enumeration of exact-radius relaxed mask covers and observed backbones.
- `quantified_phase_cover_z3.py`: experimental quantified-SMT encoding; its
  public 14-row timeout is an inconclusive scaling result, not a certificate.
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
- `order_pool_1050000_direction3_1_period15120_layered_pool{,_verification}.json`:
  the 99-row base-24 layered candidate and its independent necessary-condition
  replay.
- `order_pool_1050000_direction3_1_period15120_layered_noncover_`
  `{certificate,verification}.json`: the exact column-454 prime-deficit
  obstruction proving that this particular layered placement cannot cover.
- `order_pool_1050000_direction3_1_period15120_layered_family_noncover_`
  `{certificate,verification}.json`: the stronger 68-column exact MDD
  obstruction for every active-class placement of the declared 99-row
  family.
- `order_pool_1050000_direction1_neg3_period110880_pruned_layered_`
  `pool{,_verification}.json`: the independently replayed 81-row
  necessary-condition survivor in the wider shear.
- `order_pool_1050000_direction1_neg3_period110880_pruned_layered_`
  `projection_noncover_{certificate,verification}.json`: the exact
  30-branch modulo-30 obstruction for column 24 of that placement.
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
