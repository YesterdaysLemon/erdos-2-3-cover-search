# Erdos Problem 203 research log

Status: **open / not solved**. Nothing in this directory is a proof unless an
exact all-residue verifier returns UNSAT for the uncovered-point formula.
Random checks, finite samples, density bounds, and CEGIS progress are only
discovery evidence.

## Critical component-bound correction (2026-07-26)

The original parallel-class union bound in `component_core.py` was applied to
the raw number `s` of same-direction rows even when `s >= q`.  Its overlap
inequality assumes a proper parallel class `s < q`.  Once `q` rows of one
direction are available, they may choose all `q` distinct offsets and cover
`F_q^2` by themselves; extra rows merely duplicate offsets.  Applying the
formula beyond its domain could therefore report a false local obstruction.
The nominally independent component-core verifier repeated the same formula.

Both implementations now return feasible immediately when a direction has at
least `q` rows.  A permanent regression checks the formerly failing
`q=3`, capacities `(6,2,2,1)` example against both implementations and exact
plane enumeration.  The corrected verifier rejects the old
`2^43*3^5` four-row, density-`7/16` artifact at its first recorded peel.

Every historical component-core certificate containing a
`parallel_class_union_bound` record with `max_parallel >= prime`, and every
claim depending on such a certificate, is **withdrawn**.  A complete scan of
3,877 JSON files found 106 affected artifacts among 238 component-core
artifacts.  In particular, the prior empty
bounded ordinary cores, the 801/868-row equivalence reductions, and the
claimed `2^43*3^5` obstruction are not proofs.

Rebuilding the current artifacts gives:

* the `2^43*3^5` family retains all 75 rows, density `2.153234145944`;
* the 898-row parity-child branch still independently reduces to the same
  872 rows, because its actual eliminations remain valid under the corrected
  bound;
* the last valid 12,639-row parent pre-core now retains all 12,639 rows;
* filtering that parent to components at most 16,384 retains 12,577 rows,
  density `7.737722692865`.

The corrected 12,577-row construction pool is a strict superset of the former
9,396-row pool: it restores 3,181 legitimate fibres.  Constructive CEGIS was
restarted on this larger pool while preserving the accumulated exact
witnesses and phase checkpoints.  No result below that relies on an affected
component-core reduction should be treated as current unless it has been
explicitly rebuilt under this correction.

The correction also reopens the direct non-perfect-power route.  Rebuilding
the complete independently verified `h <= 1,050,000` pool reduces 129,497
rows to a nonempty 26,032-row core of density `2.270053134203`.  The
independent verifier passes 19,402 eliminations (19,399 tangent bounds, two
valid proper-parallel-class bounds, and one independently enumerated
inessential row).  Exact whole-lattice CEGIS is now active on this core.  A
cover here would directly construct the required `m`, without any
conditioned-cell lifting.

For faster exact iteration, restricting that corrected global core to
prime-power components at most 256 retains 23,729 rows and density
`2.248860507386`.  The independent component verifier removes no further
row.  This is still a genuine complete-period global construction family,
not a local or sampled problem.  A fast exact CEGIS branch is running
200-coordinate-diverse counterexamples per round on this subpool while a
separate 26,032-row branch requests batches of 2,000.

## Corrected normalized-leaf and conditional-fibre frontier (2026-07-26)

The finite divisor-family sweep now has independently verified no-cover
certificates for **4,609 of 4,637** ranked families.  The remaining 28 are
the intersection of the aggregate block-star survivors with the separately
verified single-anchor survivors.  This is a finite checkpoint only; it does
not prove that no unrestricted affine prime-fibre cover exists.

The main normalized-leaf correction is that the `p=71`, `h=35` row does not
remain active on every period-720 base cell after jointly normalizing its
target with `p=5` and `p=7`.  Its projection modulo five is a fixed
target-zero activity mask.  On the residual modulus seven, compatible and
incompatible target configurations were evaluated separately; compatibility
is pointwise worst for the union.  The certificate generators use algebraic
residual weights, while the verifiers directly enumerate the residual
equations and transpose the base grid.  This is distinct from, and replaces,
the earlier invalid unconditional-leaf experiment.

The corrected normalized-leaf blocks and exact period replays certify:

* period `137005268400`, upper `0.9959650988716595`;
* period `8454045600`, upper `0.9971118301032595`;
* periods `1266697832400` and `1520037398880`, uppers
  `0.9978572271640542` and `0.9979678602785422`;
* period `131037706800`, upper `0.998769347828364`;
* period `1255683068640`, upper
  `104982827676739454263/105963707585293477920`;
* periods `53542288800` and `81722440800`, uppers
  `4952507641673617/4959992888928000` and
  `221111775048971251/221437577264904000`.

A second exact refinement replaces selected generic Bonferroni star edges by
the minimum intersection of an outside fibre with a subunion of verified
block anchors.  Two proof modes are implemented:

1. a finite kernel-image subgroup with a complement transform over every
   adversarial anchor-target tuple; and
2. a projected base tensor whose residual matrices exactly count independent
   congruence lines, including a two-line shared residual pair.

The generators parameterize the outside fibre or traverse its kernel image.
Their independent verifiers instead scan the complete base plane or enumerate
disjoint cyclic cosets, and replay the target transforms through a separate
matrix-contraction implementation.  The period-level verifier treats the
exact block as the centre of a star, checks that every conditional anchor is
a recorded block anchor, and independently recomputes every remaining edge.

This conditional-fibre method additionally certifies:

* period `2092805114400`, upper
  `3531818089585662619499/3532123586176449264000`;
* period `465817912560`, upper
  `1682832143791918121/1682925587213270400`;
* period `2831442213600`, upper
  `2154205577113076005151/2155125002407545168000`;
* period `829905476400`, upper
  `17986428066422989931/17989894208141856000`;
* period `995886571680`, upper
  `6377772380734115701/6378235219250294400`;
* period `9626903526240`, upper
  `9872363416620982493917/9873990145272975321600`.
* period `36100888223400`, upper
  `65212664589123122303/65213366504514228000`;
* period `369318549600`, upper
  `10237100158509227/10237510194912000`;
* period `11107965607200`, upper
  `65279951250600888761279/65284074037427547456000`;
* period `292907815200`, upper
  `24898388077728761/24899507554521600`;
* period `122583661200`, upper
  `2656923511077515879/2657250927178848000`;
* period `274010536800`, upper
  `2501623385216309010677/2502725809247471808000`;
* period `447069823200`, upper
  `94974513775765669/95011278826464000`;
* period `262075413600`, upper
  `4343824929882296897/4344308818065216000`;
* period `80313433200`, upper
  `1740461905358613053/1740957504013728000`.

The last two families require a stronger conditional replay.  When an outside
modulus divides the period-720 component of a projected-pair/start-leaf block,
`certify_projected_pair_conditional_fibre_overlap.py` restricts the exact
base-plane partition to that fibre and minimizes over every remaining anchor
target.  On the shared residual triangle, incompatible targets are pointwise
minimal for the union.  Factorized extra rows are retained as independent
residual components; when the start leaf cannot be normalized, omitting it
gives a still-valid anchor subunion.  The independent verifier transposes the
base grid, directly enumerates the shared residual equations, and contracts
the target tensor in the opposite orientation.

For period `122583661200`, 17 independently replayed conditional edges reduce
the old upper `1.0075654877188904` to approximately `0.99987679`.  For period
`274010536800`, 19 projected-subunion edges reduce the old upper
`1.0081171072672201` to approximately `0.99955951`.  Both final period
certificates have `verified=True`.  Period `447069823200` reuses 15 of the
same independently checked conditional fibres, adds seven new fibres plus the
previously unused `p=59` replay, and reduces its old upper
`1.0082490596819373` to approximately `0.99961305`.  Period
`262075413600` combines 29 independently replayed full-block and projected
subunion edges, reducing its old upper `1.0085328464605865` to approximately
`0.99988862`.  Period `80313433200` combines 30 replayed conditional edges,
including a new `p=191` edge, and reduces its old upper
`1.0090495177160266` to approximately `0.99971533`.

All listed upper bounds are strictly below one and their period certificates
have `verified=True`.  The current source compiles and all 30 regression tests
pass.  The whole-lattice exact CEGIS branch remains active; it has not produced
a cover or a global UNSAT result.

## Exact forced-overlap frontier (2026-07-26)

A new phase-independent obstruction packages unavoidable overlaps among a
small set of low-modulus fibres.  If a disjoint block of rows has maximum
possible union density `U`, its loss relative to the sum of its individual
densities can be subtracted before applying the union bound to every other
row.  The calculation is exact for every target assignment, not sampled.

The independently replayed base blocks are:

* `(p=5,7,11)`: maximum union `53/120`, loss `3/40`;
* `(p=5,7,11,13)`: maximum union `353/720`, loss `79/720`;
* `(p=5,7,11,13,17)`: maximum union `1529/2880`, loss `379/2880`;
* `(p=5,7,11,13,17,19)`: maximum union `983/1728`, loss `1289/8640`.

The first two targets, `p=5` and `p=7`, are fixed to zero only after an
exhaustive check that their target map is jointly surjective; translation
then preserves the union size and removes a factor of 24 from each target
enumeration.  The certificate generator uses bitset unions.  Its independent
verifier instead constructs all intersection tables and applies exact
inclusion-exclusion.

A modulus coprime to the base block can be adjoined without another grid
search.  If the base maximum is `U`, one independent modulus-`h` line changes
it exactly to `U + (1-U)/h`.  This gives verified extensions by the
`p=23`, `h=11` row.  A separate nested-component certificate handles
`p=101`, `h=100`: relative to the period-720 base it projects to a
modulus-20 line and covers exactly one fifth of the top-digit lifts over each
compatible base cell.  Four further disjoint, jointly-surjective row pairs
close the remaining period-277,200 gap.  Each pair map is exhaustively
enumerated again by the final verifier.

The same projected-row theorem was independently replayed for `p=97`,
`h=48` in the period-332,640 family.  After fixing a target-selection bug in
the certificate generator, both implementations recover maximum block union
`20167/34560`, hence loss `5369/34560`.  Adding the coprime `p=23`, `h=11`
row gives extended block union `23623/38016` and a full-pool upper bound below
one.

A further phase-independent inequality replaces the rows outside an exact
block by the vertices of a forest.  For any forest `F` on sets `A_v`,

```
|union A_v| <= sum_v |A_v| - sum_(uv in F) |A_u intersect A_v|.
```

This is pointwise: if a point lies in `r` sets, their induced subgraph has at
most `r-1` edges.  The exact anchor union can itself be one forest vertex.
An edge from an outside row to that block is witnessed by its
phase-independent intersection with one anchor row, which is a subset of its
intersection with the whole anchor union.  Separate verifiers reconstruct
every selected pair map by exhaustive enumeration and independently reject
cycles.

For larger families, a two-anchor version strengthens each block-to-row
edge.  For an outside row `C` and anchors `A_i,A_j`, the second Bonferroni
bound gives

```
|C intersect (A_i union A_j)|
  >= |C intersect A_i| + |C intersect A_j|
     - |C intersect A_i intersect A_j|.
```

The pair densities are fixed when their two-target maps are surjective.  The
maximum nonempty triple-fibre density is the cokernel index of the
three-target homomorphism divided by the product of its three moduli.  The
generator obtains that index from maximal presentation minors.  A separately
implemented verifier checks the lattice formula, uses direct image
enumeration on all manageable grids, and compares against deterministic
small-grid unit tests.

The following corrected divisor-period families are now independently
certified not to cover:

* period `5040`, 31 rows: upper union density `823/840`;
* period `7920`, 28 rows: corrected component cascade reaches the empty core;
* period `15120`, 41 rows: upper `6043/6048`;
* period `25200`, 47 rows: upper `8321/8400`;
* period `55440`, 55 rows: upper `21803/22176`;
* period `83160`, 57 rows: upper `15473/16632`;
* period `110880`, 60 rows: upper `27659/27720`;
* period `138600`, 62 rows: upper `43559/46200`;
* period `166320`, 72 rows: upper `7321/7392`;
* period `277200`, 81 rows: upper `2216537/2217600`;
* period `332640`, 80 rows: upper `265585/266112`;
* period `554400`, 90 rows: upper `1330313/1330560`;
* period `655200`, 93 rows: upper `275172031/275184000`;
* period `720720`, 99 rows: upper `17279789/17297280`;
* period `831600`, 106 rows: upper `6648713/6652800`;
* period `942480`, 108 rows: upper `22497661/22619520`;
* period `982800`, 107 rows: upper `3140957/3144960`;
* period `3603600`, 148 rows: upper `86443361/86486400`;
* period `82882800`, 246 rows: upper
  `10940115133/10940529600`.

The original block-catalog scan eliminates 192 of the 200 smallest periods
of reciprocal density at least one.  Eight separately replayed projected
block and overlap-forest certificates eliminate the former exceptions, so
all 200 are now resolved.

The ranking was then extended from the first 200 to all 4,637 divisor-period
families of reciprocal density at least one in the corrected max-component-32
pool.  Two further independently replayed catalog blocks handle important
missing-component cases: `(p=5,7,13,17,19,23,47)` for families without the
modulus-5 component, and a coprime `p=239`, `h=119` extension of the
`(p=5,7,11,13,17,19,97,23,47)` block.  An independent aggregate replay
verifies:

* 4,403 block plus pair-anchor-star obstructions;
* 183 additional single-anchor-star obstructions.

Thus 4,586 of the 4,637 ranked finite families are now certified not to
cover, with 51 still unresolved by these inequalities.  Applying the
strongest verified projected/factorized block and pair-anchor star to the
entire 1,577-row max-component-32 core gives upper `1.088154030194`, so this
method does not rule out that full finite family.

Three further exact extension principles account for the newer catalog
entries.  First, an order-`23` row and an order-`23*13` row have an exact
two-line union on their shared component when their determinant is a unit.
Second, the residual path

```
7 -- 17 -- 19 -- 29
```

from the `p=29,239,647,1103` rows has product intersection densities because
each adjacent two-target map is jointly surjective.  Third, the `p=53` row
adds an independent order-`13` residual above a separately enumerated
modulus-`4` base projection.  Bitset generators exhaust the normalized base
targets.  The `p=311`, `h=5*31` row supplies a second independent projected
residual; its exhaustive calculation separates all eight base-projection
activity categories.  Independent inclusion-exclusion verifiers reconstruct
the maxima.  Small-grid and four-component path regressions agree with the
closed formulas.

A third independent projected row, `p=59`, `h=2*29`, is handled by an exact
160-cell projection histogram over all 34,560 normalized base assignments.
The generator and transposed-grid verifier both check 5,529,600 target
combinations and obtain maximum block union
`68101091/111260240`.  After the coprime `p=23` and `p=47` extensions, an
independent 796-row star replay gives upper
`63885794597480972802877/64180935944274339590400 =
0.9954014172206895` for period `760018699440`, eliminating the former
closest finite survivor.

For the next survivor, the projected `p=139`, `h=6*23` row and the
base-independent `p=47`, `h=23` row form a jointly-surjective two-line pair
on their shared order-23 component.  Combining this pair with the
`7--17--19` chain and the `p=53,311` projections gives 480 projection cells.
Both implementations check all 16,588,800 normalized target combinations
and recover maximum `9581660713/15268823700`.  After adjoining the coprime
`p=23` row, the independently replayed 815-row period
`497943285840` certificate has upper
`5161457148184908294739/5165573248190332177200 =
0.9992031668495136`.

The same histogram engine combines the `7--17--19` chain with all three
independent projected residuals `p=53,311,59`.  Its two implementations agree
on maximum `1574272253/2566932680` after 5,529,600 target combinations.
Following the coprime `p=23,47` extensions, period `4813451763120` has the
independently replayed upper
`578945742365550948639547/579232946897075914803360 =
0.9995041640275065`.

The four-residual endpoint path

```
p29 -- p239 -- p647 -- p1103 -- p59
```

has conditional projected events at both ends.  Unit determinants on the
shared components `7,17,19,29` make every enabled event subset independent.
The generator and transposed verifier check 1,105,920 target combinations
and agree on maximum `39058710917/63875221592`.  With the coprime
`p=23,47` factors, period `310545275040` has verified upper
`740282634304323847021/741833598843609592320 =
0.9979092824297747`.

A related block omits the order-13 `p=53` projection, retaining the
`7--17--19` chain and the independent projected `p=311,59` residuals.  The
generator and transposed verifier check all 1,382,400 target combinations
and agree on maximum `16075768847/26656608600`.  After the coprime
`p=23,47` factors, period `5553982803600` has independently replayed upper
`255730731253272758898047/257056041522370968108000 =
0.9948442749633531`.

Keeping compatible sub-blocks in the catalog is essential when a family
lacks one of those coprime factors.  The endpoint path with only the
order-23 `p=47` extension gives independently replayed upper
`8567188230163970041/8588483885158632000 =
0.9975204407111409` for period `423470829600`.  The shorter
`p=29,239` residual pair with `p=53,311,59` and only the order-11 `p=23`
extension gives independently replayed upper
`363813167679822985463/363975062065828012800 =
0.9995552047299995` for period `165221456400`.

The `p=139,47` shared order-23 pair can also coexist with the independent
`p=59` order-29 residual.  The generator and transposed verifier enumerate
all 33,177,600 normalized targets and agree on maximum
`93593928823/147598629100`.  The resulting period-`6563797858800`
star certificate has independently replayed upper
`41573082968544303719633/41626514329649724384000 =
0.9987164103945316`.

On the order-11 component, the projected `p=67`, `h=6*11` row and the
base-independent `p=23`, `h=11` row form another jointly-surjective pair.
Adding this pair to the full `p=29--239--647--1103--59` endpoint path and
the independent `p=53` residual requires 192 projection cells.  Generator
and transposed verifier agree across 6,635,520 target combinations on
maximum `229688420981/351313718756`.  Period `67509842400` then has
independently replayed upper
`5842167542599736773/5843049770349792000 =
0.9998490124533027`.

The projected `p=53`, `h=4*13` and `p=79`, `h=6*13` rows likewise form a
jointly-surjective pair on their shared order-13 component.  Combining that
pair with the order-11 `p=67,23` pair and the `p=29,239,647` chain requires
576 projection cells.  Generator and transposed verifier agree across
19,906,560 targets on maximum `2052849507/3143992852`.  Period
`6983776800` has independently replayed upper
`2415239476483751/2419302499614000 =
0.9983205807744598`.

Removing the unsupported `p=647` tail leaves the `p=29,239` residual pair
and the same order-13 and order-11 projected pairs.  Its two exhaustive
implementations again agree across 19,906,560 targets, now on maximum
`88831121/136272136`.  Period `10659448800` has independently replayed
upper `65241301247409269/65301062482656000 =
0.999084835177642`.

Adding the independent `p=311`, `h=5*31` projection to the full
`p=29,239,647` chain and both shared pairs gives 2,880 projection cells.
The generator and independently transposed verifier agree across
99,532,800 targets on maximum `414653009/632881678`.  The separate 741-row
period replay then gives period `108248540400` the strict upper
`390416007535904678527/390591549957681072000 =
0.9995505729148636`.

The same block can be extended by the independent `p=59`, `h=2*29`
projection without constructing a 5,760-square dense score matrix.  The
exact scorer separates the plain `p=311,59` residuals from the 576-state
shared-pair block and combines them by
`1-(1-U_base)(1-U_extra)`.  A dense 480-state regression reproduces an
older verified certificate's maximum, loss, maximizing base targets, and
maximizing projections exactly.  The generator and independently
transposed verifier then check all 199,065,600 targets and agree on maximum
`36408872615/55060705986`.  A separate 842-row replay gives period
`1046402557200` upper
`109390335866483365592779/109495831171469927184000 =
0.9990365358766823`.

Extending the path through the projected `p=1103` endpoint while retaining
both shared pairs and the factorized `p=311` residual produces a 16-row
block.  The generator and independently transposed verifier again check all
199,065,600 normalized target combinations and agree on maximum
`912724574738/1379020409013`, with forced-overlap loss
`972133091188361/3640613879794320`.  All 27 regression tests pass after
replay.  Adding the block as the 86th distinct verified catalog entry does
not resolve another ranked family, but it improves the closest remaining
survivor, period `1255683068640`, to exact upper
`647379757524147129373/647265997072728633600 =
1.0001757553338704`.  The bound remains above one and is therefore not an
obstruction.

This is a rigorous finite-family advance only.  It neither rules out larger
periods nor proves that a hypothetical `m` must arise from a finite prime
cover.

## Problem

Find an integer `m >= 1`, coprime to 6, such that

`2^k * 3^l * m + 1`

is composite for every `k,l >= 0`, or prove no such `m` exists. The maintained
problem page still marks it open: https://www.erdosproblems.com/203

## Exact affine-line reduction

For a prime `p >= 5`, put

* `r2 = ord_p(2)`,
* `r3 = ord_p(3)`, and
* `h = lcm(r2,r3) = |<2,3>|`.

Choose a generator `g` of `<2,3>` and write `2 = g^a`, `3 = g^b`. A fixed
nonzero `m (mod p)` selects one fibre

`p | 2^k 3^l m + 1  <=>  a*k + b*l = c (mod h)`.

Consequently, one affine fibre for each of finitely many distinct primes
covers `Z^2` if and only if their CRT residues produce a solution. The residue
is `m = -g^(-c) (mod p)`. Also impose `m = 1 (mod 6)` and choose the positive
CRT representative larger than every selected prime, so every divisor is
proper.

## Correctness repair

The original multiplicative-order routine could skip a remaining prime factor
after shrinking its trial-division bound. For example it reported
`ord_34511(3)=493`; the correct value is `17`. The routine now factors the
fixed group order `p-1` first and only then reduces the candidate order. A
500-prime randomized comparison with `sympy.n_order` found zero mismatches.

Every exact loader recomputes signatures from `p`, so stale JSON metadata
cannot enter a proof. Older candidate artifacts are exploratory only. The
corrected smooth-29 density is `1.567801157677`.

## Exact checker

`exact_uncovered.py` decomposes every line modulus into maximal prime-power
components. CRT makes the two residues on distinct components independent.
It SAT-checks whether a pair exists outside every selected affine line and can
therefore return either a genuine missed pair or an exact UNSAT cover result.
The implementation uses linked one-hot selectors for lower powers and only
`O(q)` matching pairs for a modulus `q`, since at least one line coefficient
is a unit on each prime-power component.

`exact_uncovered_z3.py` is an independent direct-modular checker. It is slower
but is reserved as a second verifier for any final certificate.

## Rigorous finite-family obstructions

Capacity bounds using normalized `p=5` and `p=7` anchor fibres rigorously rule
out complete-period covers for

`720720, 1441440, 2882880, 5765760, 11531520, 23063040`.

Some are certified by exact rational LP duals. For `L=720720`, for example,
the exact upper bound is `26597/732160`, below the required `1/24` capacity in
every uncovered anchor cell. A three-anchor (`p=5,7,11`) relaxation also rules
out the first 75 smooth-order candidates, while the first 100 and the larger
unrestricted pools remain inconclusive.

These are impossibility results only for the declared finite families, not for
the original problem.

## Corrected prime pools and active exact searches

Scanning by subgroup order through `h <= 300000` found 40,796 corrected prime
signatures, no unresolved cofactors, total reciprocal density
`2.751770465527`, and largest prime `9277494431`.

Active construction methods are:

* `exact_cegis.py`: SAT master plus exact uncovered-point checker. An ordinary
  4,979-prime run uses the `h <= 30000` pool (density `2.436995880416`).
* The same exact loop with `m=M^3`, using power-compatible residues and the
  algebraic `X^3+1` cover of pairs with `3|k` and `3|l`.
* `exact_greedy.py`: a memory-light exact loop that greedily spends an unused
  prime fibre on batches of genuine checker witnesses. Ordinary and `M^3`
  variants use the full 100,000-order pool.

The most promising current construction is a disjoint induced parity split.
The three parent lines `k=0`, `l=0`, and `k+l=0 (mod 2)` cover `Z^2`. An
even-order prime whose signature reduces to the corresponding parity direction
induces a line of half the original index on that parent. Odd-order fibres
induce full-index lines on any parent and can be allocated disjointly.

At order 100,000, an alternating allocation of the 2,811 odd-order fibres
between the first two parents gives derived densities

* `1.522946681957` for the `k`-even slice,
* `1.630117437051` for the `l`-even slice, and
* `1.589907891058` for the `k+l`-even slice.

All three survive exact two- and three-anchor capacity tests. Their respective
first three anchor moduli are `(5,6,11)`, `(3,8,15)`, and `(2,9,14)`; each
triple map is jointly surjective, so all three targets normalize to zero. The
exact three-anchor upper capacities are respectively about `0.00332995245`
versus required `1/330`, `0.00327406825` versus `1/360`, and `0.00454402429`
versus `1/252`. These are feasibility screens, not covers. Normalized exact
CEGIS and greedy searches are active on all three disjoint pools.

The same disjoint allocation through order 300,000 raises the derived
densities to `1.614399391921`, `1.721828690888`, and `1.666576385716`.
The memory-light searches have been resumed on these strict supersets; the
100,000-order CEGIS masters remain active because loading all 40,796 original
fibres into their bit-equality encoding at once would be unnecessarily large.

`incremental_greedy.py` now retains one exact SAT complement checker while
fibres are added.  Refining a prime-power residue component links the old and
new one-hot levels in both directions, so every reported model is a genuine
uncovered CRT class.  If the incremental instance becomes UNSAT, the script
rebuilds the independent primary checker before emitting a cover.  Regression
tests include the three parity lines and a later 2-adic refinement.

Batch model enumeration in `exact_uncovered.py` can additionally block the
target fingerprint on a declared jointly-surjective normalization family.
This prevents CEGIS batches from differing only in irrelevant high-order CRT
coordinates.  On the `k`-slice anchors `(p,h)=(11,5),(13,6),(23,11)`, a test
batch returned 100 witnesses with 100 distinct nonzero anchor fingerprints.
The default checker path is unchanged when no diversity primes are supplied.

`local_phase_cegis.py` replaces the binary target master with a dense
multivalued target matrix and min-conflicts repair.  It can retain all random
tests, not only misses, and still accepts a construction only after the exact
checker returns UNSAT.  This reduced the memory for a 5,000-fibre / 8,000-point
master from about 30 GB to about 0.2 GB.  `phase_hillclimb.py` performs
population-level coordinate descent and supports frozen bands plus split
validation.  A corrected freeze regression is essential: an early version
mistakenly zeroed frozen phases; all results from that run were discarded.

On the direct first-5,000 pool (density `2.437694261408`), sample CEGIS can
absorb tens of thousands of counterexamples but plateaus near a `6%` fresh
miss rate.  Optimizing the first 1,000 broad fibres on one million fixed points
and merging them back lowers the full-pool miss rate to about `5.7%` on
100,000--200,000 fresh points.  Higher-index band optimization easily
overfits, even with split validation.  A validated residual assignment over
all 40,796 fibres appeared to leave only 59/100,000 training points, but an
independent audit missed 4,294/100,000; it is explicitly not a certificate.
The practical conclusion is that extra density cannot simply memorize the
roughly random residual.  A successful phase design must align the residual
with later fibres or with algebraically covered sublattices.

No `cover_*exact.json` file should be trusted unless both exact checkers return
UNSAT and the CRT construction is independently verified.

## Maximal-component core reduction

`component_core.py` now applies a rigorous covering-system reduction before
synthesis.  At a largest surviving `q`-adic exponent, freeze every lower
digit and every other CRT component.  The fibres at that exponent become
affine lines in `F_q^2`.  If those lines cannot cover the top-digit plane,
then all lower-exponent fibres already cover every frozen context, so the
maximal group is redundant and may be deleted.  Iterating preserves the
existence or nonexistence of a cover within the declared finite pool.

The cheap screen uses a parallel-class union bound, Blokhuis's theorem that a
nontrivial blocking set in `PG(2,p)` has at least `3(p+1)/2` points, and the
Blokhuis--Brouwer theorem that every essential point of a blocking set `B` lies
on at least `2p+1-|B|` tangents.  The latter gives a much stronger
direction-capacity obstruction.  Dualize an affine line cover and adjoin the
point `P` dual to the line at infinity.  If fewer than all `p+1` directions are
present, `P` is essential.  Reduce to a minimal blocking subset containing
`P`.  If that subset has size `b`, then at most `b-p` direction rays through
`P` can contain its other `b-1` points.  Therefore the sum of the largest
`b-p` available direction capacities must be at least `b-1` for some
`3(p+1)/2 <= b <= N+1`.  Failure for every such `b` is an exact obstruction.

The tangent result is Result 3.11 in Bartoli et al., *On the metric dimension
of affine planes, biaffine planes and generalized quadrangles*, Australasian
Journal of Combinatorics 72(2) (2018), 226--248:
https://ajc.maths.uq.edu.au/pdf/72/ajc_v72_p226.pdf .  Its short proof reduces
to Jamison's `2p-1` lower bound for an affine blocking set.  Independent SAT
checks of all triggered profiles in `PG(2,3)` and randomized triggered
profiles in `PG(2,5)` and `PG(2,7)` found no counterexample.

The two trivial projective-line exceptions are exactly `p` parallel affine
lines and `p+1` concurrent affine lines using every direction.  Tiny top
planes are also enumerated to remove a row that can never be essential in any
cover; such removals are performed one at a time before recomputing.

For the order-300,000 pool the tangent-capacity cascade reduces 40,796 fibres
of density `2.751770465527` to an equivalent 970-fibre core of density
`1.146387846807` without using the optional finite-plane enumeration or
essential-row routines.  The optional exact routines remove one further row,
leaving 969 fibres of density `1.146382539402`.

The three-anchor exact rational LP on the tangent-only core uses the fibres
for `p=5,7,23`, whose jointly realized anchor space has 264 cells.  For every
target of the third anchor it certifies residual capacity at most
`0.002618516158422`, strictly below the required `1/264 =
0.003787878787879`.  Thus **no cover exists using any of the 40,796 fibres
with subgroup order at most 300,000**.  This is a complete finite-pool
impossibility result, not a proof about primes or subgroup orders beyond the
scan.

An independently structured fourth-anchor LP adds the `p=47`, order-23 fibre.
Because 23 is coprime to the first three anchor moduli, a candidate fibre is
uniform across the new coordinate unless its order is divisible by 23 and its
23-direction is parallel to the anchor; in that exceptional case its target
selects one of the 23 coordinate values.  A memory-efficient dual formulation
with marginal variables exactly reproduces the direct formulation.  On the
complete `h <= 400000` pool, tangent pruning leaves 3,212 fibres of density
`1.445831666587`, but the four-anchor exact ratio is only
`0.965372476397517 < 1`; hence that entire pool is impossible too.

The complete `h <= 500000` scan has 65,363 fibres of raw density
`2.814603178392`.  Tangent pruning leaves 6,024 fibres of density
`1.606144872761`.  This pool survives the separate fourth-anchor tests at
orders 23, 43, and 83, with exact normalized upper ratios about `1.1373`,
`1.1578`, and `1.1686`.  Five-anchor tests at order pairs `(23,43)` and
`(43,83)` are also inconclusive, with exact ratios about `1.1140` and
`1.1454`; the `(23,83)` calculation is still running.  This is feasibility
evidence, not a cover.

The complete `h <= 700000` scan has 89,114 fibres of raw density
`2.854600039781`, with every 25,000-order interval fully factored and no
unresolved cofactors.  Adding the higher-order fibres triggers a stronger
tangent-pruning cascade: only 2,515 fibres of density `1.187887889881`
remain.  On this equivalent core, normalize the `p=5` and `p=7` targets and
condition on each of the 12 possible targets of the surviving `p=13`,
order-12 fibre.  For all 12 targets the exact three-anchor dual upper bound is
at most `0.005300965084432`, strictly below the required `1/144 =
0.006944444444444`.  Therefore **no cover exists using any prime fibre with
subgroup order at most 700,000**.  As before, this is a bounded finite-pool
theorem, not a solution of the original problem.

The complete `h <= 900000` scan has 112,417 fibres of raw density
`2.883901207176`, again with no unresolved cofactors.  The same
tangent-capacity cascade leaves an equivalent 5,284-fibre core of density
`1.315612167778`.  With the `p=5`, `p=7`, and `p=13` fibres as anchors, all
12 possible targets of the third anchor have exact rational dual capacity at
most `0.006303707640108`, below `1/144 = 0.006944444444444`.  Consequently
**no finite prime-fibre cover exists using any fibre with subgroup order at
most 900,000**.  This strengthens the bounded theorem only; an arbitrary
solution need not have a finite covering set.

The complete `h <= 1000000` scan adds 11,449 independently rechecked
signatures, for 123,866 fibres of raw density `2.895965407621`; all four
intervals have zero unresolved cofactors.  Unexpectedly, the enlarged pool
triggers a complete tangent-capacity cascade: 20 rounds reduce the equivalent
ordinary core to **zero rows**.  Thus no finite prime-fibre cover exists using
any subgroup order at most one million, without needing an anchor LP.
An independently coded replay rebuilt every maximal prime-power group from
the merged pool and verified all 18,696 recorded eliminations (18,662
Blokhuis--Brouwer tangent bounds and 34 parallel-class bounds), ending with
the same empty core.

Applying the tangent-capacity cascade to the alternating allocation of all
odd-order fibres between the `k` and `l` branches reduces each of the `k`,
`l`, and `k+l` finite branch pools to the empty set.  Hence that entire
disjoint parity construction is impossible through subgroup order 300,000.
Earlier compact local repair on the weakest pre-tangent core fit 11,000 of
11,200 accumulated points, while exact binary CEGIS reached 8,200 genuine
counterexamples before its master became solver-hard.  Those searches are now
superseded by the rigorous tangent-capacity obstruction.

## Perfect-power route

If `m=M^D`, then whenever an odd prime `q|D` divides both `k` and `l`, the
expression is `X^q+1` and factors.  If `4|D`, Sophie Germain's identity also
factors every cell with `k=2 (mod 4)` and `l=0 (mod 4)`.  For the remaining pairs, the
selected residue `-g^(-c)` must be a `D`-th power modulo `p`.

Writing a primitive root as `r`, `g=r^((p-1)/h)`, this is exactly

`gcd(D,p-1) | (p-1)/2 - ((p-1)/h)*c`.

The CEGIS search encodes this congruence with a small binary remainder
automaton; it does not enumerate or relax invalid targets. Translations by
`D(s,t)` preserve the algebraic sublattices and the power condition, so the
`p=5,p=7` fibres can still be normalized to zero for odd `D`.  For even `D`,
the permitted target tuples form an affine translation orbit rather than a
subgroup; the exact LP now normalizes to a canonical tuple in that orbit and
keeps the absolute mod-4 Sophie cells fixed.

Component pruning for a perfect-power search skips prime components dividing
`D`, because the algebraically covered set changes across those top-digit
planes.  This yields safe supersets rather than an unsafe reuse of the ordinary
core.  Power-aware exact-dual LPs track both anchor cells and the excluded
`qZ^2` algebraic cells while enforcing the permitted target congruence for
every prime.  Through `h <= 300000` they certify:

* `D=3`: exact upper `0.000551845948526 < 1/792`;
* `D=15`: exact upper `0.000402950419033 < 1/1800`.

Thus neither the cube nor fifteenth-power construction works in that finite
pool.  A coarser exact measure LP for `D=105` is inconclusive (ratio about
`1.3250`), and its local phase assignment missed `13.370%` of 100,000 fresh
pairs, close to the random-density baseline.  It remains exploratory.

The power-aware LPs now also handle even `D`.  If `4|D`, they include the
Sophie Germain cells exactly and normalize the `p=5,p=7` targets to the
canonical affine orbit (for example `(2,1)` at `D=4` and `(2,3)` at `D=12`).
A CRT-factorized counter computes local prime-power cell counts and multiplies
them, avoiding an explicit `D^2` grid.  It exactly reproduces the older grid
LPs at `D=60` and `D=420`.

On safe tangent cores through `h <= 300000`, exact two-anchor normalized
ratios rule out `D=4` (`0.383154`), `D=12` (`0.338738`), `D=60`
(`0.794600`), `D=420` (`0.986595`), and `D=4620` (`0.978621`).  At
`D=60060` the two-anchor ratio rises to `1.086921`, but adding the fixed
power-compatible `p=13`, order-12 anchor target gives the stronger exact
ratio `0.930089 < 1`.  Extending the odd primorial to 17 and 19 lowers that
three-anchor ratio further, to `0.924051` and `0.921891`.  Thus every tested
mixed fourth-power construction is still impossible in the 300,000-order
pool.  The ratios above are normalized residual capacities; values below one
are rationally post-checked impossibility certificates.

Among smaller exact power cores, `D=660`, `780`, `924`, and `1092` are also
excluded.  The strongest surviving design is currently `D=1540 =
4*5*7*11`: fifth-, seventh-, and eleventh-power identities cover
`5Z^2`, `7Z^2`, and `11Z^2`, while Sophie Germain covers
`k=2 (mod 4), l=0 (mod 4)`.  Through order 300,000 its safe core has 5,938
rows, but a fourth anchor at `p=11` gives exact ratio
`0.988099496982 < 1`, so that pool is impossible.  Through order 700,000 the
safe core has 12,750 rows of density `2.020896782249` (11,871 rows and
density `1.749269848739` after enforcing the 1540th-power residue condition).
Three- and four-anchor tests are inconclusive.  The combined five-anchor
test at `p=5,7,13,11,17` also survives, but narrowly: its exact rational
normalized capacity is `1.023803931727 > 1`.  This is not a cover, but it is
the first rigorously screened power space here with enough residual capacity
to warrant exact CEGIS synthesis.  The synthesis fixes the canonical anchor
tuple `(2,1,2,5,0)`, retains all counterexamples, and writes no trusted cover
unless an exact checker returns UNSAT.

The nearby mixed powers were re-screened on their own safe order-900,000
cores.  `D=420` is already impossible by the two-anchor exact ratio
`0.963468672475`.  `D=4620` and `D=60060` survive two anchors with ratios
`1.077733857392` and `1.070785438261`, but the power-compatible `p=13`
anchor drops them to `0.851770130858` and `0.845434195573`; both are
therefore impossible through order 900,000.  This leaves `D=1540` as the
only currently surviving exact power construction.

The complete order-one-million pool gives additional safe-core exclusions.
Two anchors already rule out `D=660` (`0.851094264995`), `D=780`
(`0.742919632398`), `D=924` (`0.624739498679`), `D=1092`
(`0.540740590431`), and `D=4004` (`0.780390934544`).  `D=2860` barely
survives two anchors at `1.001057639131`, but the `p=13` anchor lowers its
exact ratio to `0.945974693271`.  `D=1820=4*5*7*13` survives three anchors
at `1.026905567803`, then a fourth `p=11` anchor rules it out at
`0.916926729671`.  These are exact rational dual certificates for the full
declared pools, not heuristic synthesis failures.

Combining the odd factors from the two strongest branches gives
`D=20020=4*5*7*11*13`.  Its safe one-million core has 13,486 rows of density
`1.936069112328`, of which 12,606 rows are power-compatible.  Its exact
two-anchor ratio is `1.210262656613`, so this larger algebraic family remains
promising and has been promoted to higher-anchor screening.  It retains
substantial slack at three anchors (`1.150345639257`) and remains just above
the obstruction threshold at four anchors (`1.051359922007`).

The one-million `D=1540` core does not survive the corresponding five-anchor
test: its exact ratio is `0.966161171397 < 1`.  This retrospectively rules out
the entire order-700,000 synthesis pool as well, so its unfinished Z3 and
local-repair assignments are not construction candidates.

The same fifth anchor also closes the apparently stronger `D=20020` branch:
its ratio falls from `1.051359922007` at four anchors to
`0.959469739883 < 1`.  All of its finite-period masters are therefore
superseded by an exact obstruction on the full one-million core.

Adding the next odd algebraic exponent instead gives
`D=340340=4*5*7*11*13*17`.  Its safe core has 16,982 rows of density
`2.057319340200`, including 15,800 power-compatible rows.  Exact ratios remain
comfortably feasible at two, three, and four anchors:
`1.305813505690`, `1.250403928928`, and `1.148768763031`.  Its five-anchor
ratio is still `1.054023398569 > 1`, so this is the first algebraic family
here to survive the complete five-anchor, order-one-million screen.

Within `D=20020`, restricting all fibre orders to divisors of `180180` leaves
60 compatible rows of density `0.954911754912`; the exact binary master is
already UNSAT on 90 genuine residual points.  The divisor family for period
`720720` has 96 compatible rows and density `1.127279664780`, but is likewise
master-UNSAT on 173 residual points.  The larger period `10810800` (190 rows,
density `1.204661172161`) is UNSAT on 259 points, and period `302702400`
(315 rows, density `1.243387333566`) becomes UNSAT after adding its first
100 counterexamples, at 526 points.  The 90-point `180180` obstruction was
replayed with both binary target-index and independently one-hot target
encodings.  These are finite-family obstructions, not evidence against
arbitrary smooth periods.

For `D=340340`, divisor-period CEGIS is substantially more resilient.  A
smooth period with 1,099 compatible rows and density `1.365099659518`
absorbed more than two thousand residual points.  Adding a `29` component
raises this to 1,434 rows and density `1.407936504042`; adding `23` and then
`37` gives 1,858 rows/density `1.444017517569` and 2,263
rows/density `1.470079894059`.  These are active synthesis families, not
certificates: every assignment still has fresh genuine counterexamples.
Extending the period through ten additional CRT prime components gives 5,392
rows and density `1.590528116969`.  Point-wise repair can memorize thousands
of witnesses, but a retained 100,000-point audit still missed
17,407/86,521 non-algebraic pairs.  Global best-target coordinate descent
halved the training misses from 16,871 to 8,289; on a fresh audit it missed
16,383/86,609 (`18.92%`), so most of the apparent gain was overfit.  An
independent exact whole-lattice pair-intersection objective reduced the excess
pair overlap among the 2,000 densest rows from `0.267071707810` to
`0.237024316941`, but its fresh miss rate was still `19.33%`.  These negative
audits are retained explicitly; neither phase assignment is a cover.

The exact checker can now quotient its counterexamples by a declared tuple of
small fibre maps.  With `p=5,7,13,11,17`, an exact pass through the `D=340340`
assignment found all `1,351` residual anchor fingerprints, one genuine CRT
witness per cell.  Coordinate repair covered that complete cross-section in
33 phase moves.  A quota extension, independently smoke-tested on a complete
four-cell toy enumeration, retains a chosen number of distinct full-CRT
representatives from every anchor fingerprint before blocking the cell.  The
current searches fix the five even-power anchor phases (without loss, since
their valid target tuples form one translation orbit) and request two
representatives per residual cell.  Checker UNSAT remains the only acceptance
condition.

An odd-power branch avoids the fourth-power target restrictions.  At order
300,000, `D=15015=3*5*7*11*13` had exact two-anchor ratio
`1.432252431505` and three-anchor float ratio `1.328633654642`; all tested
translation-orbit representatives gave the same value.  The authoritative
one-million tangent core is stricter: it has 13,486 rows of density
`1.936070802969`, with exact two-anchor ratio `1.335002608778` and
three-anchor float ratio `1.220438619942`.  A smooth 55-digit period whose
prime-power components are all at most 256 retains 8,156 of those rows at
density `1.830494678664`.  Its exact synthesis fixes `p=5,p=7` at zero by
translation and uses the two-witness-per-five-anchor-cell checker.  Four- and
five-anchor screening is still in progress, so this is a construction branch,
not a certificate.

Omitting the cubic factor gives substantially more target freedom at the
order-6 and order-12 anchors.  For
`D=85085=5*7*11*13*17`, the exact one-million safe core has 16,980 rows of
density `2.057311716542`; adding 19 gives
`D=1616615` with 20,275 rows of density `2.139794360484`.  Conservative
two-anchor supersets give exact rational ratios `1.608176507534` and
`1.685705541648`.  At three anchors the two target orbits for `D=85085`
have float ratios `1.567366788709` and `1.565385410181`; the first
`D=1616615` orbit is higher still at `1.648037361388`.  These values are
prioritization evidence until the higher-anchor screens finish.

The period `lcm(1,...,256)` retains 18,214 rows of the `D=1616615` core at
density `2.119103858785`, with no prime-power checker component above 256.
The concurrent assignment `c=0` for every row is not a cover: the primary
exact checker produced a genuine CRT counterexample (and an independent
20,000-point audit missed about 30.6 percent of non-algebraic points).
A general-phase exact CEGIS run on the same period is active, with `p=5,p=7`
normalized to zero.  Its first assignment repaired all 2,319 misses on a
retained 20,000-point sample in 590 coordinate moves, but then missed
2,446/18,411 (`13.2855%`) on an independently generated fresh sample.  This
is useful progress, not a certificate: the exact whole-torus checker is now
supplying genuine counterexamples.  A separately reproduced 32-bit LCG
initialization has zero power-compatibility errors across all 18,214 fibres.
It initially appeared to miss 5,154/46,073 (`11.1866%`) points when the audit
coordinates came from the same LCG family, but a later Python-64 audit found
11,711/92,123 (`12.7124%`).  The optimistic LCG figure is therefore treated
as generator correlation, not as a robust benchmark.

An attempted clause-reduction in the primary checker briefly returned a false
UNSAT result on that LCG branch.  It was rejected immediately: an independent
100,000-draw audit found 9,304 misses among 91,988 non-algebraic points.  The
defect was an incorrectly oriented implication when `a` is a unit but `b` is
not on a prime-power component.  The artifact is explicitly quarantined as
`cover_power1616615_lcm256_lcg32retained64_rejected_encoding_bug.json`.
The corrected encoding now finds genuine witnesses again, agrees with Z3 on
the targeted regression, and preserves the exact diversity-quota count.  It
still removes redundant local predicates and pair variables, reducing a
typical check from 77,702 component instances to 35,741 distinct equalities.
No result from the rejected encoding is used as mathematical evidence.

The corrected checker now enumerates an entire five-anchor cross-section in
roughly one to two CPU minutes: the first LCG assignment had 1,351 residual
fingerprints, and successive one-representative rounds were absorbed in
14--24 phase moves.  A two-representative round returned exactly 2,702
genuine CRT witnesses.  The untouched older checker independently returned
1,361 genuine witnesses on a different assignment before it was retired for
speed, providing a second SAT-side implementation check.

A separate random-retention experiment deliberately kept every tested point,
not merely misses.  Its first repair moved the fresh miss rate from
2,000/19,956 (`10.02%`) to 2,000/21,019 (`9.52%`), but after growing the
retained matrix to 63,326 points the next audit regressed to 2,000/19,952
(`10.02%`).  The branch was stopped as sample memorization; this apparent
short-lived improvement is not treated as progress toward an exact cover.

The four-anchor fractional LP completed at `1.589795452560`.  Eight rounded
mixtures had robust 100,000-point miss rates between `10.0192%` and
`10.3329%`, improving on the raw random seed but not nearly enough.  Starting
from the best rounding, a power-compatible coordinate sweep accepted only
moves non-worsening on each half of a 100,000-point sample.  On a separate
one-million-draw audit, it reduced misses from 93,198 to 91,399 among 921,700
eligible points (`10.1115%` to `9.91635%`).  This is a statistically real
heuristic improvement, still far from an exact cover.

Repeated fresh-sample streaming epochs continue to improve that same fixed
million-draw benchmark without retaining a giant target matrix.  The best
completed checkpoint, epoch 8 on a 400,000-point split sample, leaves 87,301
misses among 921,700 eligible points (`9.47174%`), a 6.33 percent reduction
in hole count relative to the rounded fractional seed.  It stops at 120
training holes because every available move would worsen one of the two
200,000-point halves.  A soft-margin `rho=0.1` branch is not materially
better: from the epoch-5 base it gives `9.63567%`, statistically tied with
the next zero-hole epoch.

Pure exact-witness repair and measure optimization pull in opposite
directions.  Fitting 5,404 exact CRT witnesses from the epoch-6 checkpoint
with 63 moves raises a common robust audit from `9.70346%` to `11.87136%`.
Combining those witnesses with 50,000 random validation points and requiring
every move to be non-worsening on each side is much safer: it covers all
55,404 retained points while changing a separate 100,000-point audit only
from `9.67598%` to `9.68032%`.  Nevertheless, an exact replay still finds
holes in 1,375 five-anchor fingerprints, so this cycle has shifted exact
holes rather than eliminated them.

Increasing the retained prime-power component bound gives a small but
reproducible improvement.  On the same million-draw audit, two fresh
component-512 hill-climb epochs reach `9.34002%` missed, and two
component-1024 epochs reach `9.21981%`, the best mature fixed-phase result so
far.  A component-2048 epoch is statistically tied rather than clearly
better.

The five power-compatible anchor targets have 9,216 possibilities, split
into exactly 16 translation orbits of size 576.  Exact enumeration of the
residual anchor cells ranges from 1,351 to 1,450 cells.  The best orbit has
representative `(0,0,1,0,3)` and 1,351 cells; the previously optimized orbit
has representative `(0,0,0,0,3)` and 1,375.  A first frozen-anchor hill
epoch in the best orbit reaches `9.35066%`, so the lower exact cell count has
not yet overcome its one-epoch optimization deficit.

For a fixed phase assignment, summing every fibre's exact conditional
density in each residual anchor cell gives a rigorous necessary condition:
any cell below 1 is impossible to cover under those phases, regardless of
subsequent within-cell search.  The component-256 epoch-9 phase has minimum
capacity `0.889002517902` with 8 of 1,375 cells below 1.  The first
component-1024 phase has minimum `0.912532681229` with 17 cells below 1.
This motivated a deterministic exact-capacity coordinate optimizer rather
than further fitting of finite witness samples.

An exact weak-cell pass selected the 100 component-1024 cells whose original
capacity was below 1.2.  After 1,027 phase moves, every selected cell has
capacity at least `1.111800657920`; the subsequent full 1,375-cell
recomputation has the same minimum and no cell below 1.  An independent
million-draw audit misses `9.36518%`.  The latter is slightly worse than the
mature measure-only phase, so the result removes a rigorous fixed-phase
obstruction but is not evidence of a cover.

Exact SAT still finds one hole in every one of the 1,375 five-anchor cells.
The streaming phase optimizer now supports mixed 64-bit and arbitrary-size
exact CRT witnesses.  It simultaneously fitted four exact witnesses per cell
(5,500), 3,796 additional six-anchor fingerprints, and 50,000 random
validation points while keeping the missed random measure near 9.4 percent.
Thus the accumulated finite constraints remain feasible, but a fresh exact
counterexample pass is still required after every repair.

The accumulated exact-witness repair has since fitted 17,943 distinct
adversarial exponent pairs plus the fixed 50,000-point validation set.  The
latest hard repair started with exactly 1,712 misses and eliminated all of
them in 450 phase moves; an independent 20,000-point audit still missed
`9.57%`.  A separate six-anchor
capacity profile using anchors `5,7,13,11,17,19` has minimum normalized
coverage `1.039079272030` over all 3,797 required cells, so it gives no local
capacity obstruction.  On the algebraic-aware component-1024 core, the final
full recheck of the accumulated weak-cell repair has minimum
`1.091216469709` and no cell below 1.

An exact decomposition with anchors `5,7,13,11,31` provides a small but real
improvement.  There are 3,600 anchor cells and 1,780 cells outside the five
anchor fibres; algebraic identities remove 36, leaving 1,744 required cells.
The current phase fully covers 32 of those cells and the exact checker returns
1,712 witnesses, one in every remaining cell.  This is progress in phase
coupling, not an exact cover.  After repairing against those 1,712 witnesses,
an exact recheck returned a coordinate-disjoint second witness in precisely
the same 1,712 fingerprints: all 32 closed cells persisted, but no additional
cell closed.  A different `5,7,13,19,31` basis likewise returned a witness in
all 1,086 required cells.  Combining both new layers gives 20,740 distinct
exact adversarial points for the next repair.

The component cascade now treats algebraic primes exactly instead of skipping
them wholesale.  At exponent `q^1`, the identity `X^q+1` covers only the
origin of the top `F_q^2` plane, so an exact SAT relaxation tests whether the
available maximal fibres can cover the other `q^2-1` points.  At higher
q-adic exponent there are lower-digit states on which the identity covers
nothing, so the ordinary full-plane obstruction remains valid.  The
punctured-plane SAT encoding matches exhaustive enumeration on 60 small
finite-field cases.

Applying this proof-preserving algebraic-aware peel to the main
`D=1616615` million-order core reduces it from 20,275 rows and density
`2.139794360484` to 15,571 rows and density `2.009693161638`.  The
component-256/512/1024 subcores have respectively 14,437/14,745/15,246 rows
and densities `1.998348171298`, `2.001018353743`, and `2.006178369093`.
The component-1024 core can still fit the cumulative exact tests: a first
repair removed 3,363 misses and covered all 17,943 exact plus 50,000
validation points; a second repair covered the expanded 20,740 exact plus
50,000 validation set.  The latter required 464 moves from 2,585 initial
misses and had a fresh 20,000-point miss rate of `10.81%`.  This preserves
finite-sample feasibility after every proof-safe row elimination, but it is
not an exact cover.

A stronger within-cell test enumerated exactly four new CRT representatives
in each of the 1,712 persistent `5,7,13,11,31` fingerprints.  After
deduplication the strict-core set contained 29,299 exact witnesses plus
50,000 validation points.  A 1,471-move repair covered all 8,559 newly
introduced coordinates, but the exact recheck still returned the identical
1,712 open-fingerprint set and no new closed cell.  A bounded soft-margin
pass (`rho=0.25`) then made 1,168 moves while preserving every hard point and
improved a 20,000-point audit from 2,176 to 2,071 misses.  Nevertheless, its
exact recheck again left precisely the same 1,712 cells open.  Thus point
multiplicity and sampled coverage margin improve generalization without yet
changing exact cell closure.

The complete subgroup-order scan has been extended through `h=1,025,000`.
The new interval contains 2,817 rows of density `0.002782430787` and no
unresolved cofactors.  An independent verifier rechecked every prime, both
multiplicative orders, the stored subgroup order, coefficient orders,
signature surjectivity, and the logarithmic relation; all rows passed
(`sha256=8ab369011999f4233332c7e9f977ab75640aaec847bf0791d79ce4a20cb90cab`).
The merged pool has 126,683 unique rows and density `2.898747838408`.

The optimized scan has since been extended through `h=1,050,000`.  Four
independently verified shards add 2,814 fibres of density
`0.002712370303`, with no unresolved cofactors; the merged pool now has
129,497 fibres and density `2.901460208711`.  Re-running the ordinary
finite-cover component cascade eliminates all 129,497 rows in 20 rounds.
An independent replay reconstructed every live maximal prime-power group and
validated all 19,480 recorded eliminations (19,445
Blokhuis--Brouwer tangent-capacity bounds and 35 parallel-class bounds),
ending with zero survivors.  Thus no ordinary finite prime-fibre cover using
subgroup orders at most 1,050,000 exists.  This remains a bounded theorem and
does not exclude an infinite or higher-order construction.

For `D=1616615`, the new algebraic-aware two-stage core contains 16,342 rows
of density `2.015442942737`; its component-1024 portion contains 16,003
usable rows.  Projecting the strongest saved `p=31`-orbit phase onto this
core preserves target 25 on that anchor.  Exact SAT again returns precisely
1,712 open `5,7,13,11,31` fingerprints.  Those fingerprints are identical
to the million-order set, while all 1,712 full CRT coordinates are disjoint.
An independent direct replay confirms that every returned coordinate avoids
all 16,003 fibres and all six algebraic sublattices.  The additional rows
therefore create a new CRT component but close no coarse cell under the
projected phase.

Conditioning the projected phase on its 1,744 required five-anchor cells
gives a mean union-bound capacity of `1.48684992962946`, but 96 cells have
capacity below one and hence cannot close under that assignment.  A targeted
coordinate pass moved 1,918 phases and raised all 132 cells below `1.10` past
that threshold.  On the full profile this reduced the sub-unit cells from 96
to 14, yet exact SAT still returned the identical 1,712 open fingerprints at
coordinate-disjoint witnesses.  A second pass repaired all 42 then-current
cells below `1.05`, but the full profile created 16 new sub-unit cells
elsewhere (minimum `0.772161616940027`) and exact SAT again returned 1,712
open fingerprints.  Thus weak-cell capacity optimization is oscillatory and
does not change exact closure.

The first persistent cell, `(k,l)=(0,1) mod 60`, has also been transformed
into canonical local coordinates.  The transformation preserves every
perfect-power target congruence and moves the residual algebraic
`7,11,13,17,19` exclusions to their origin sublattices.  An independent
replay checked 24,069 allowed-target samples across all 12,068 transformed
rows and verified the line equivalence and all five algebraic translations.
The component cascade leaves an equivalent 1,716-row cell core of density
`2.436435680431`; a separate replay validates all 38 eliminations over ten
rounds (23 tangent-capacity and 15 parallel-class bounds).  Exact cell
synthesis is active on this much smaller, rigorously equivalent problem.

That component-1024 cell branch is now rigorously obstructed.  Conditioning
again on the admissible nested cell `(178152,229352) mod 323323`, only 411 of
the 1,716 rows can meet the cell under any allowed target.  Their exact
conditional union-bound capacity is
`532580301465724074037255 / 680183364370876925855557 =
0.782995188302384... < 1`; an independent rational replay agrees.  This
rules out the component-1024 conditioned pool, not the unrestricted original
problem.

Restoring every prime-power component through subgroup order 1,050,000
changes the conditioned problem substantially.  The first ordinary cascade
leaves 10,480 rows.  A new proof-safe lower-digit peel uses the fixed derived
target modulo an algebraic prime `q`: at a maximal `q^e`, freeze the lowest
digit of `(k,l)` and retain only rows available in that state.  For the
`19^2`, `13^3`, and `11^3` groups, every required non-origin state has fewer
than `q` available top-plane lines, so none can cover `F_q^2`.  The peel
removes 263 rows; an independent verifier replays all 648 states
(`PASS input=10480 removed=263 survivors=10217 states=648`).  The subsequent
ordinary bounds remove another 956 rows, independently replayed as one
tangent and two parallel-class certificates.  The resulting equivalent
full-component cell core has 9,261 rows, density `7.478463834429`, and maximum
prime-power component `2^14=16384`.

The projected full-component phase initially misses only 56 of 95,974
admissible random points.  An algebraic-preserving translation simultaneously
normalizes nine composite high-order anchors to zero; an independent line
replay verifies the translation on all target restrictions.  A fixed
100,000-point hill pass then closed its two retained holes in two moves and
missed zero of a separate 100,000-point audit.  The exact modular checker is
tractable at component 16,384 and returned a genuine CRT hole, so this was not
yet a cover.  Exact CEGIS is now active: each exact batch supplies diverse
full-period witnesses, and the first batches were repaired in only a handful
of phase moves.  No exact cell cover has yet been certified.

A rational-dual one-axis capacity relaxation was also checked on both
coordinate axes.  Relative to the 24 persistent residue classes modulo 60,
the floating optima and rational-dual ratios agree near
`8.357204935629255` on the `k` axis and `9.088795517975106` on the `l`
axis.  These values are far above one, so this relaxation supplies no
axis-capacity obstruction; any successful obstruction must use phase
coupling or genuinely two-dimensional overlap.

An SPF/grouped implementation of the component cascade exactly reproduces
the established one-million safe core, including all 20,275 ordered rows and
round totals.  The extended safe core has 20,552 rows/density
`2.140067983153`.  The algebraic-aware second stage also exactly reproduces
the old 15,571-row baseline; on the extension it leaves 15,763 rows/density
`2.009882884952`.  Its component-1024 subcore has 15,432 rows/density
`2.006362161949`, a net gain of 186 usable high-order fibres.  Projecting the
best saved phase onto those rows changed a matched 200,000-point audit by only
three misses, and exact SAT still returned the identical 1,712 open
`5,7,13,11,31` fingerprints.  Thus this interval adds real but very small
construction capacity and does not solve or rule out the branch.

The three-anchor LP on the tightened full core is still inconclusive at
`1.524639275751`, and the four-anchor ratio is
`1.457526667290`.  Applying the same peel to
`D=4849845=3*1616615` produces exactly the same final row set; its extra
cube identity affects the residual domain, not this final candidate core.
For that exponent the three- and four-anchor ratios are much tighter,
`1.278158161841` and `1.169616085023`, but still above the rigorous
impossibility threshold.  The completed alternate five-anchor LP using
`5,7,11,13,23` in the all-zero orbit has ratio `1.091100032570 > 1`
over 3,520 required cells, so that stronger relaxation is close but still
does not exclude the orbit.  Eight random component-1024 phases missed
`17.54%`--`18.02%` of non-algebraic audit points.  Projecting the best
`D=1616615` margin phase onto the stricter cube-compatible targets improved
this to `15.84%`; a hard repair reached `15.29%`.

The repaired cube branch has 1,672 required `5,7,13,11,31` cells after its
extra algebraic exclusions.  Exact SAT leaves 1,648 open and proves 24
closed.  An alternate `5,7,13,19,31` basis leaves a witness in all 1,104
required cells.  After jointly repairing both exact layers, the original
decomposition again had precisely the same 1,648 open fingerprints, with a
coordinate-disjoint witness in each.  Its weaker audit and smaller closed-cell
count rank this branch below `D=1616615`; it remains exploratory, not ruled
out globally.

The fractional solution of the cube branch's `5,7,11,13,23` LP was rounded
with twelve independent seeds and compared on the same 82,032 admissible
audit points.  Seed 409 was best at 12,161 misses (`14.8247025551004%`).
This average improvement did not translate into exact cell closure: the
corrected component-1024 SAT check returned a witness in every one of the
3,520 required five-anchor fingerprints.  An earlier invocation that passed
the unfiltered core to the checker was rejected by its component-size guard
and produced no mathematical result.

On the stronger `D=1616615` phase, a deeper exact enumeration has now
produced sixteen distinct CRT representatives in each of the same 1,712
persistent `5,7,13,11,31` fingerprints, for 27,392 new adversarial points.
The enumeration is exact and exhaustive at the fingerprint level.  A joint
repair covered those points together with all 79,299 prior exact/validation
points in 1,621 moves.  The ensuing exact recheck nevertheless returned one
witness in precisely the same 1,712 fingerprints, and all 1,712 coordinates
were disjoint from the 27,392 training coordinates.  Increasing sampled depth
from four to sixteen representatives per cell therefore moved the holes
without closing a cell; further point-only repair is now deprioritized in
favour of cell-level relaxations and new fibres.

The long-running five-anchor fractional LPs have completed.  On the full
million-order core the optimum normalized residual capacity is
`1.542122263315`; on the component-256 truncation it is
`1.519265352386`.  Both are comfortably above 1, so neither is an
impossibility certificate.  Eight independent roundings of the full
solution miss between `9.92413%` and `10.31979%` on a common 100,000-draw
audit before hill optimization.

Replacing the exponent factor 19 by 23 or 29 collapses the authoritative
million-order core to densities `2.057324701474` and `2.057352848798`.
Adding 23 instead gives `D=37182145` and only seven additional full-pool
rows; its component-256 subpool is exactly the same 18,214 rows and density
as `D=1616615`.  On a common robust audit seed its LCG assignment is slightly
worse (`12.7799%` versus `12.7124%`), so the extra algebraic 23-cell does not
offset the stronger phase restriction in this baseline.

Removing the factor 5 does not improve the authoritative core despite giving
some fibres more phase freedom.  For `D=323323=7*11*13*17*19`, the exact
one-million cascade leaves 10,174 rows of density `1.593556273048`, and its
exact two-anchor capacity is `1.204004569975`.  Removing 7 as well gives
`D=46189=11*13*17*19`: only 5,299 rows of density `1.359070850659` survive,
and the exact two-anchor ratio is `0.980211663129 < 1`.  Thus that entire
finite prime-fibre pool is impossible already at two anchors.  The factor 5
in the current `D=1616615` branch is structurally useful because it prevents
this tangent-capacity collapse.

For comparison, the completed full five-anchor float screen for the strongest
even exponent `D=6466460` is `1.108816858158 > 1`, improving materially on
the `D=340340` value `1.054023398569`.

Multiplying the strongest even exponent by 3 is not helpful.  On the safe
order-300,000 cores, `D=1021020` and `D=19399380` survive two anchors at
`1.081718065158` and `1.079689729670`, but the `p=13` anchor gives exact
upper ratios `0.924050883777` and `0.921891471308`.  Both finite pools are
therefore impossible.

`exact_cegis.py` now accepts even powers, removes primes having no compatible
power target, skips both odd-power algebraic cells and the Sophie Germain
cells, and passes the same exclusions to both exact checkers.  Out-of-range
binary target codes are replaced only by a valid power-compatible target.
This was compile-checked and exercised on a small `D=4` round before use on a
real synthesis search.

## Certified low-cell locking

For the `D=1616615` branch, an exact optimization of the 21 rows whose
moduli divide 432 leaves exactly 622 of the 186,624
`(k mod 16,l mod 16,k mod 27,l mod 27)` cells open.  This is the current
low-adic base assignment, recorded in `tmp/low432_multistart_best.json`.

Primary PySAT unsatisfiable cores initially closed 135 of those cells using
the phases of only 24 distinct rows.  Independent Z3 replay proved every one
of the 135 core claims UNSAT.  The core bundle is
`power1616615_order1050000_cell60_0_1_low432_checkpoint1_closed_cell_cores_pysat.json`;
its independent audit is
`verify_low432_checkpoint1_closed_cell_cores_z3.log`.

One remaining cell,
`(k,l)=(0,5) mod 16` and `(k,l)=(2,0) mod 27`, was then conditioned exactly.
The conditioned pool has 9,237 rows, density `19.659935002332`, and largest
prime-power component 1,024.  Exact CEGIS found a cover, from which the
primary checker extracted a 22-row core.  The child-coordinate core was
independently proved UNSAT by the structurally separate Z3 bit-vector checker.

The 22 child phases were lifted algebraically back to the parent pool.
Both the primary PySAT checker and independent Z3 bit-vector checker then
proved the lifted 22-row core UNSAT directly under the original fixed
coordinate conditions.  The parent proof object and verifier outputs are:

* `power1616615_order1050000_cell60_0_1_focus_a_parent_core.json`
* `power1616615_order1050000_cell60_0_1_focus_a_parent_core_pysat_verification.json`
* `power1616615_order1050000_cell60_0_1_focus_a_parent_core_z3bv_verification.json`

All 24 phases supporting the earlier 135 certificates are unchanged.  The
new core overlaps seven of those rows, so the residual-cell certificates use
39 distinct prime phases.  Fourteen of these are also among the 21 low-base
rows.  Preserving both the original low-base cover and every residual-cell
certificate therefore freezes 46 distinct prime phases in total.  A full
exact recount of the 9,261-row assignment returns exactly 486 open low cells,
down from 487, in
`power1616615_order1050000_cell60_0_1_low432_checkpointA_open_cells_exact.json`.
Thus 136 of the base's 622 residual cells are now rigorously closed.

This is a certified local advance, not a solution of Erdos problem 203.
There is still no globally covering phase assignment and no constructed
integer `m`; the remaining 486 low cells must be closed (or the construction
abandoned) before CRT assembly and a global exact verification can begin.

## Other explored directions

* A common-projective-form reduction would turn the problem into a classical
  one-dimensional covering system, but scans found maximum compatible density
  below `0.493`, far short of 1.
* The public Vela repository contains an exact but partial 20-prime witness
  with `m=8168305011630835886634520238999`; it leaves about 25.3 percent of
  sampled exponent pairs uncovered and is explicitly not a solution.
* Cochrane--Myerson's composite one-dimensional cover yields a homogeneous
  cover of `Z^2`, but its small exact directions do not occur in the corrected
  prime-signature pool. Higher-order fibres can refine several of its early
  classes, so recursive refinement remains a structural lead.

## Main files

* `search_cover.py`: prime orders and finite-torus SAT search.
* `exact_uncovered.py`: primary exact all-residue checker.
* `exact_uncovered_z3.py`: independent exact modular checker.
* `exact_cegis.py`: exact counterexample-guided SAT synthesis.
* `exact_greedy.py`: memory-light exact greedy synthesis.
* `generate_order_pool.py`: order-first candidate discovery.
* `anchor_capacity*.py`, `triple_anchor_capacity_lp.py`: rigorous obstruction
  tools.
* `candidate_*.json`, `order_pool_*.json`: search pools or exploratory phases,
  not certificates by themselves.

## Reduced-block finite synthesis update (2026-07-24)

The active `q=19`, direction `(1,1)`, mod-16 cell `(0,4)` still has two
unresolved mod-9 blocks.  Proof-replayed component peeling reduces `(2,2)`
from 3,468 to 298 equivalent rows and the correctly 36-lock-conditioned
`(7,6)` block from 3,463 to 289 equivalent rows.  The obsolete 25-lock,
299-row artifacts are not evidence.

Weighted min-conflicts with clause breakout and basin restarts materially
improved the finite stress tests:

* on the 298-row `(2,2)` block, the best phase now misses 24 of 921,432
  retained adversarial/random points (previous best 42);
* on the 289-row `(7,6)` block, the saved phase misses only one of 299,501
  retained adversarial points.

These are finite sample results, not covers of the exponent lattice.  Exact
critical-core masters on the 289-row block were SAT: the 26,703-point
one-hole-plus-singletons core solved in 2.510 seconds but its phase missed
2,460 points in the full finite universe; after adding those misses, the
29,163-point master solved in 190.327 seconds and its phase missed 1,417 full
finite points.  Thus neither core is a finite obstruction, and the decoded
phases overfit.  Fast WalkSAT CEGIS rounds on the same accumulating core also
oscillated rather than closing the full universe.

`finite_component_cegis.py` now uses two proof-safe exact strengthenings.
Every direction group is nonempty, and in fact may be assumed to select its
full number of distinct flats: whenever two interchangeable rows duplicate a
flat while an unused flat exists, moving one duplicate to that flat preserves
the old coverage and only adds new coverage.  The master therefore enforces
exactly `min(number of rows, number of flats)` selections per group.
Independently on each CRT component, two simple groups with unit-determinant
directions are translated so that a selected flat in each has target zero.
The strengthened exact `23 x 29` master remains undecided; no SAT or UNSAT
claim is made from its runtime.

The new finite-search utilities are `sample_phase_walksat.py` (weighted
breakout, basin restarts, and rare global moves),
`finite_sample_phase_misses.py`, `finite_sample_phase_core.py`,
`finite_sample_coordinate_cegis.py`, `finite_sample_soft_repair.py`,
`finite_sample_milp_repair.py`, and `merge_point_sets.py`.  Every positive
finite result above was replayed by rebuilding the complete target matrix and
counting uncovered points.  There is still no integer `m`, no exact local
cover, and no global impossibility proof.

## Exact-coordinate residual synthesis (2026-07-24)

The one-miss finite phase for the 289-row `(7,6)` block is not an infinite
cover.  The exact CRT checker produced fresh genuine uncovered exponent pairs
on every repair round.  Ten 100-witness rounds, then batches of 1,000, 5,000,
50,000, and 100,000, were all satisfiable as finite point sets; the phase
search could cover even the accumulated 83,703-point core quickly, while the
exact checker still found more holes.  These results rule out treating the
one-miss stress score as near-proof, but do not rule out another phase
assignment.

An exact coordinate decomposition gives substantially stronger structure.
For the repaired phase, PySAT proved 426 of the 529 mod-23 cells closed, using
a union of 31 rows, and the independent Z3 bit-vector checker replayed all 426
cell certificates.  Only 103 mod-23 cells remained open.  Conditioning the
first open cell gave a 258-row child in which 712 of 841 mod-29 cells were
dually certified closed.  Deeper exact grids similarly exposed sparse
residuals rather than relying on random lattice samples.

`close_coordinate_residual.py` now supports recursive CRT-coordinate
descent.  On one branch it descended through the full component schedule
`23,29,31,37,23,29,31,32`.  The deepest mod-32 grid had 156 open cells;
three verified full-cell retargets reduced this to 68, then 40, then zero.
The exact cover minimized to four rows with primes
`59393,380929,151553,471041`.  Lifting that core closed the surrounding
mod-31, mod-29, and mod-23 grids completely, with both PySAT and Z3 checks at
every completed level.

That greedy branch is not a cover of the parent block.  A second higher
coordinate branch exhausted the remaining mutable 2-adic fibres and left 44
cells open.  Reordering the descent to put the mod-32 coordinate first
eventually produced a leaf with only 22 compatible mod-31 lines, scaled
density `22/31 < 1`; that leaf is impossible under the phases locked along
that particular path.  Neither failure is a global obstruction because a
different coupled choice of higher-component phases can avoid the locked
branch.

`component_density_cegis.py` implements that stronger coupling condition.
For a chosen residual prime power it asks, exactly, whether some assignment of
all other CRT coordinates leaves compatible line density below one.  A Z3
bit-vector pseudo-Boolean checker finds violating cells, and a conservative
PySAT master retargets only unprotected phases while enforcing every learned
density cut.  The original repaired phase already passes the exact residual
density condition for components 23, 29, and 31; component 37 has violations
and remains under CEGIS.  Density feasibility is necessary but not sufficient
for an affine cover, so a successful density run will still require the full
exact uncovered-point checker.

The primary checker was also corrected so a fixed coordinate modulus is added
to the ambient CRT domain even when a supplied proof core uses only a proper
divisor.  This is required, for example, when a modulus-4 row certifies a
single cell of a mod-32 grid.  The independent Z3 checker already had this
semantics.

There is still no exact cover of the 289-row block, no lifted global phase
assignment, and no integer `m`.

## Coupled affine-plane closures and a locked-branch obstruction (2026-07-24)

The attempted `q=31` density repair inside the first open mod-37 branch was
not merely difficult: the 66 protected mod-37 phases make that branch
impossible.  In the conditioned 139-row pool the modulus counts are

```
23:19, 29:18, 31:20, 37:75, 23*37:1, 29*37:4, 31*37:2.
```

The 20 modulus-31 rows and two `31*37` rows contribute scaled residual
density at most 22, below the threshold 31, unless a full-density row is
active.  The 66 fixed mod-37 lines leave 175 cells of `F_37^2` open.  Each
of the nine remaining mod-37 fibres can hit at most 8 or 9 of those cells,
and the sum of the nine individual maxima is 77, so at least 98 mod-37
cells remain for every extension.  Independently, 18 mod-29 lines leave at
least `29*(29-18)=319` cells.  Nineteen mod-23 lines leave at least 92
cells, more than the single `23*37` row can project onto.  Fixing a mod-23
hole outside that projection forces the four `29*37` rows to cover at least
`319*98=31,262` pairs, while their union has capacity at most
`4*29*37=4,292`.  This is a contradiction.

The replayable certificate and independent verifier are:

* `reduced289_recursive32_guard31_v1_locked66_q31_product_obstruction.json`
* `reduced289_recursive32_guard31_v1_locked66_q31_product_obstruction_verification.json`
* `component_product_obstruction.py`
* `verify_component_product_obstruction.py`

This obstruction applies only while those 66 phases are protected.  It is
not an obstruction to the 289-row block or to the Erdos problem.

Backtracking and retargeting whole affine-plane layers jointly was much more
effective than greedily closing one coordinate cell at a time.  All 75
mod-37 fibres in the same conditioned pool were jointly retargeted to cover
all 1,369 cells of `F_37^2`.  Core minimization retained all 75 rows, and
both the primary PySAT checker and independent Z3 bit-vector checker proved
the 75-row core globally covering.

One level higher, joint retargeting of all 63 mod-31 fibres covered 960 of
the 961 cells of `F_31^2`, leaving only `(7,1)`.  Conditioning that one cell
produced 75 mod-37 fibres, which were jointly retargeted to cover
`F_37^2`.  After exact lifting, the union consists of 138 disjoint rows:
63 mod-31 rows plus the 75-row lifted exceptional-cell core.  Primary PySAT
and independent Z3 checks both prove that this 138-row union globally
covers the entire 202-row conditioned parent branch.  The main proof
objects are:

* `reduced289_r31_r37_twolevel_parent_cover.json`
* `reduced289_r31_r37_twolevel_parent_core138.json`
* `reduced289_r31_r37_twolevel_cell7_1_verification.json`

The 138-row core was then lifted into the surrounding mod-29 cell `(1,28)`;
both exact checkers verified that lift.  At the mod-29 level, joint search
over the 51 mod-29 fibres currently leaves 21 affine-plane cells before
deeper refinement.  This is exact local progress, but the surrounding
mod-29 grid, the earlier recursive levels, the other low-adic blocks, and
the remaining global cells are still unresolved.  There is still no
integer `m` and no global impossibility proof.

## Phase-independent three-component obstruction and repair pivot (2026-07-25)

The formerly search-only 289-row `(7,6)` block now has a rigorous
phase-independent obstruction.  A 233-row subproblem splits exactly into 189
rows on the squarefree coarse torus
`F_29^2 x F_31^2 x F_37^2`, 22 rows with one extra 29-adic digit, and 22 rows
with one extra 31-adic digit.  The two refining families cannot fill a missed
coarse cell: their affine-line hole lower bounds are respectively 203 and
279.  Thus the 233-row problem covers if and only if its 189-row coarse
subproblem covers.

The 189 coarse rows have family counts

```
X:51, Y:46, Z:49, XY:17, XZ:13, YZ:12, XYZ:1.
```

An Alon--Furedi grid bound gives at least six holes after the 51 `X`-only
lines.  If the full three-component family covered, every such hole would
have incidence at least 11 among the 13 `XZ` projections or at least 15 among
the 17 `XY` projections.  Exact direction multiplicities and pair-incidence
counts bound the first high-incidence set by one point and the second by two
points, giving at most three `X` holes, a contradiction.  The independently
recomputed certificate is:

* `reduced289_q23density_squarefree189_obstruction.json`
* `reduced289_q23density_squarefree189_obstruction_verification.json`
* `certify_three_component_cover_obstruction.py`
* `verify_three_component_cover_obstruction.py`

The obstruction composes upward without phase locks through the 233-row and
253-row pools.  A finite-affine-plane line-cover bound then lifts it through a
287-row parent.  Finally, the two remaining period-16 and period-32 rows can
occupy at most `64+32=96` of the 1,024 mod-32 cells, so at least 928 cells
inherit the impossible 287-row geometry.  The resulting all-phase
obstruction for the 289 present rows and its independent replay are:

* `reduced289_locked36_reduced289_obstruction.json`
* `reduced289_locked36_reduced289_obstruction_verification.json`
* `lift_period_cell_obstruction.py`
* `verify_period_cell_obstruction.py`

The preceding 24-step component peeling was replayed independently and
re-materialized a byte-identical 289-row residual from the 3,463-row
conditioned pool.  This proves that the current locked `(7,6)` branch cannot
cover; it does not prove that the surrounding mod-16 parent is impossible.

Auditing that boundary exposed a constructive alternative.  Direct peeling
of the surrounding mod-16 pool leaves 314 rows, including thirteen rows whose
moduli divide 9.  With the actual 46-prime global certificate set frozen,
seven of those rows have fixed projected targets and four unprotected rows
may be retargeted.  PySAT found a complete mod-9 cover, minimized to a
nine-row core.  A structurally separate brute-force verifier checked all 81
cells and exhausted 1,813 assignments for every proper subset of the four
mutable rows, proving that all four are necessary relative to the fixed
targets:

* `checkpointA_min17_q19closure_d02_1_1_mod16_0_4_mod9_global46safe_cover.json`
* `checkpointA_min17_q19closure_d02_1_1_mod16_0_4_mod9_global46safe_cover_verification.json`
* `solve_period_divisor_cover.py`
* `verify_period_divisor_cover.py`

The four changed primes are
`8709121,139969,530713,1492993`.  Their targets were lifted through seven
conditioning levels into the global 9,261-row phase file.  Independent
forward reprojection confirms that every other phase is unchanged and none
of the 46 protected primes is touched:

* `power1616615_order1050000_cell60_0_1_low432_checkpointA_mod9repair_lift_audit.json`
* `power1616615_order1050000_cell60_0_1_low432_checkpointA_mod9repair_lift_verification.json`

The authoritative full exact recount completed over the same component
domain as the prior checkpoint.  It still finds exactly 486 distinct
mod-16/mod-27 open cells: no previously closed cell reopened, but no coarse
cell closed.  The repair therefore closes a genuine deep branch without yet
changing the coarse ledger:

* `power1616615_order1050000_cell60_0_1_low432_checkpointA_mod9repair_exact_misses.json`

The reason is now exact.  At the surrounding `q=7` level, the prescribed
target-zero modulus-7 fibres cover 43 of 49 cells and leave precisely
`(0,1),...,(0,6)`.  Every other modulus-7 row has
`target_modulus=7`, so it is not actually retargetable.  Both PySAT and the
corrected HiGHS MILP prove the six-cell parent cover infeasible.  During this
audit `solve_component_grid_cover_milp.py` was fixed to enforce each row's
`target_residue/target_modulus`; the earlier invalid MILP phase/result files
were deleted, and the corrected infeasibility record is
`checkpointA_min17_screen_r001_12_1_2_0_q7_global46safe_milp_corrected_result.json`.

There is still no integer `m` and no global impossibility proof.  The exact
next construction boundary is to close the five sibling `q=7` cells without
changing the protected 46-prime certificate set, or to prove that their
required deeper repairs cannot coexist.

## Phase-invariant 486-cell barrier

Three further exact repair lifts tested whether the remaining coarse ledger
was an artifact of the particular saved phase.  A `p=760321` mod-9 repair, a
`p=8317` repair of the conditioned `q=11` child `(5,10)`, and a `p=64513`
repair covering an entire mod-16 child across eight `q=19` siblings were each
lifted back through the full conditioning chain.  Every lift passed the
independent forward-reprojection verifier and preserved all protected
prime phases.  The authoritative global recount after each repair still
returned exactly the same 486 open mod-16/mod-27 cells.  In every comparison
the exact witness coordinates moved; no old witness survived unchanged.

The phase search was then widened substantially.  Whole-screen exact CEGIS
retained 35,001 distinct counterexamples, including coordinate-diverse
batches and high-order phase anchors.  A double-coverage variant retained
40,001 lessons, and a focused `q=7` child search retained 27,192 lessons.
Every feasible repair acquired a fresh batch of 500 exact misses; none
closed the selected coarse cell.  These are synthesis failures, not
impossibility certificates.

Finally, the double-coverage whole-screen assignment was lifted in bulk.
All 9,215 unlocked leaf rows were considered, 4,533 root phases changed,
none of the protected primes changed, and every one of the 9,261 root target
restrictions was checked independently.  The full exact recount again found
486 misses in precisely the existing coarse cells:

* `power1616615_order1050000_cell60_0_1_low432_checkpointA_bulk_margin2_lift_audit.json`
* `power1616615_order1050000_cell60_0_1_low432_checkpointA_bulk_margin2_lift_verification.json`
* `power1616615_order1050000_cell60_0_1_low432_checkpointA_bulk_margin2_exact_misses.json`

This large, globally valid phase displacement makes further phase-only
search in the same 9,261-row pool a low-value direction.  It does not prove
that those cells are impossible, and it says nothing directly about
higher-order fibres.  The next construction experiment is therefore to
extend the complete subgroup-order cutoff, rerun every proof-preserving
core reduction, and test a genuinely strengthened row set.  There remains
no integer `m` and no global impossibility proof.

## Higher-order fibres and unlocked low-adic synthesis (2026-07-25)

The complete subgroup-order scan was extended from `h=1,050,000` through
`h=1,100,000` in four exact contiguous shards.  They contain respectively
1,395, 1,480, 1,424, and 1,453 rows, all with zero unresolved cofactors.
Independent recomputation of every prime, order, and signature passed with
SHA-256 values

```
1a427a698026191e28e8c929fb3d23f72e427f17826e0df138a077a97c53843a
eb3a5c5c3d1270928f58d766e973fc3ce1bca6a62a8363b4554655f0588be909
a13cf05a647095d839e23a637a5fe11081122f9de5dcb9e946a822812ba66391
472eb601297cb5e3a66a4402a67ea6f490b469e9348c604f2888f6d0b1f5e2f2
```

The merged complete pool has 135,249 rows.  A separate construction-only
scan on exponents divisible by 25,920 through 20,000,000 found 1,278
independently verified rows.  After conditioning on the period-60 cell, 411
remain compatible and 278 are new relative to the full conditioned pool.

Four denser selected-exponent scans, using every exponent divisible by 60
between 1,050,001 and 3,000,000, found 17,999 rows.  Each shard has zero
unresolved cofactors and independently passed row replay:

```
949284834d89ec22c0710b9a29e8720143c838b2a5afc56529c2c7aca37fcb44
568b726403c98c62933350ae35169133b03561e4a4011869eec1a2a04eff3ae2
75a57d3b6dc04ca16064955efc7e133fffae68564543b328e5ab782236c4f60c
f41adfc035f15bee020076b996b89faacd25d03dc0e17e7780f21ae715aba80f
```

Conditioning and the exact 20,000-component guard retain 3,132 rows, of
which 2,392 are new in the construction pool.  Two farther 25,920-step
shards covering 20,000,001 through 60,000,000 emitted 2,319 and 3,146 rows.
Their emitted rows independently pass all primality, order, and signature
checks, but their source scans retain respectively 38 and 145 unresolved
cofactors.  They are therefore valid construction rows but not complete
interval evidence.  This provenance is preserved explicitly in the merged
pool.  After conditioning, they contribute another 982 unique rows, giving
18,982 rows in
`power1616615_cell60_0_1_augmented_step60_3000000_p25920_60000000_max20000.json`.

Exact margin-two, margin-three, and margin-four assignments were built on
supports of 9,503, 9,697, and 9,846 rows.  They covered respectively 2,916,
3,402, and 3,888 accumulated exact witnesses with the requested finite
multiplicity.  Every authoritative all-residue recount still returned 486
misses in the same 486 mod-16/mod-27 cells, while sharing no exact witness
with the preceding recount.  These are genuine moving holes, not a
certificate that the cells are impossible.

The invariance had an important algorithmic cause: every greedy extension
preserved all phases already present in the preceding support, including the
low 2-adic and 3-adic phases.  Re-running the complete proof-preserving
component cascade on the 18,982-row pool changes the synthesis problem
dramatically:

* the algebraic-aware component peel leaves 12,883 rows;
* the derived lower-digit peel leaves 12,570 rows;
* the final ordinary component cascade leaves only 801 rows, of density
  `3.618338566880`.

Independent replays verify all three stages.  The first replay checks 1,376
tangent-capacity eliminations plus two independently enumerated
never-essential rows.  The lower-digit replay checks 696 finite states and
313 removals.  The final replay checks 29 tangent and eight parallel-class
eliminations and reconstructs the same 801 survivors.  Thus the 801-row pool
is cover-existence-equivalent to the declared 18,982-row construction pool;
this does not make the construction pool complete over all possible primes.

An unlocked exact CEGIS pass is now active on the 801-row core.  Under the
projected phase, 3,482 of the first 3,888 lessons were missed.  A single
legal move changed the unique modulus-2 row, `p=41`, from target 0 to target
1 and covered all of them.  The next exact batch moved into previously
untouched even-`l` mod-16 cells, disjoint from the old 486-cell set.  Further
rounds are jointly retargeting the remaining low rows rather than protecting
the old certificates; a second branch enforces double coverage of every
retained lesson.

This is the first phase search in the strengthened family that genuinely
escapes the old low-adic lock.  It has not yet produced even a complete
conditioned-cell cover.  Such a local cover would still require exact
lifting and the surrounding period-60/global branches before CRT assembly.
There remains no integer `m` and no global impossibility proof.

## Recursive parity reduction and small-period exclusion (2026-07-25)

The unlocked 801-row problem has an exact recursive low-adic structure.  The
unique `p=41` row is

```
h=2, a=0, b=1.
```

Choosing its phase covers one parity of `l`, and the complementary parity is
an affine rectangular lattice.  Conditioning on `l=0 mod 2` or `l=1 mod 2`,
excluding `p=41`, and independently replaying every row produces two
800-row branches.  After algebraic-preserving canonical translation, their
ordered `(h,p,a,b,target_residue,target_modulus)` data are exactly identical,
and every child target modulus is one.  An independent symbolic lift check
verifies arbitrary child phases, so a cover of this canonical 800-row branch
plus the complementary `p=41` phase is equivalent to a cover of the 801-row
parent.

The canonical branch contains another exact splitter,

```
p=17, h=2, a=1, b=1.
```

Writing `k=x`, `l=r-x+2y` conditions on one parity of `k+l`.  Both residual
parities independently replay to canonically identical 799-row branches of
density `3.241161768428`, again with unrestricted child targets.  The
general affine-lattice constructor, independent replay, phase/lesson
projector, and parity checks are:

* `condition_derived_lattice.py`;
* `verify_conditioned_lattice.py`;
* `project_lattice_state.py`;
* `verify_lattice_projection.py`;
* `verify_parity_split_equivalence.py`.

On the synthesis side, changing one prime-power CRT digit at a time proved
far more stable than replacing a whole phase.  On 6,804 exact lessons, one
component sweep made 189 moves and achieved double coverage of every point.
After merging independent streams, two sweeps double-covered 24,786 points.
A 200,000-point random training tranche plus 28,244 exact lessons was reduced
from 2,676 holes to one by component descent; a sparse MILP changed one phase
and closed the last point.  Independent recount covered all 228,244 points.
On a disjoint 500,000-draw audit, the admissible random hole rate fell from
`0.0133006373222` to `0.00184245738796`.  A second random tranche and
full-sample point repair covered all 428,244 retained points, but a new
500,000-draw audit still found rate `0.00176214390009`.  These are heuristic
construction improvements, not cover certificates.  Exact 10,000-witness
batches still reach their cap.

The two parity reductions also expose complete small-period subproblems.
Every divisor-only period through 216 tested so far is exact UNSAT.  Verified
component-core cascades additionally reduce the period
`276,348,372,2088,2484,4176,4968` subsets to zero, giving phase-independent
bounded impossibility certificates for those row families.

The first nontrivial stable family is period 432.  It has 19 rows of density
`2.037037037037` on the complete `432 x 432` torus.  Translation symmetry
sets the independent `p=19,p=37` order-3 phases to zero; the `p=73` phase has
two orbits (concurrent target 0, or triangular targets 1/2), and a remaining
translation sets the order-4 `p=97` phase to zero.  Complete SAT returns
UNSAT in both representatives.  The concurrent result was independently
reproduced by Glucose, and the triangle result by MiniCard.  The symmetry
and result replay is
`power1616615_augmented799_period432_symmetry_unsat_audit.json`.

The period-648 subset first reduces from 19 rows to an equivalent 17-row
core of density `2.027777777778`.  The same two symmetry representatives are
both UNSAT under CaDiCaL and independently under MiniCard; the replay is
`power1616615_augmented799_period648_symmetry_unsat_audit.json`.
Period 864 and unrestricted 799-row exact searches remain active.

These results rule out several concrete finite-period constructions and
sharpen the first persistent period-60 cell.  They do not cover the cell,
do not settle the surrounding 1,711 persistent coarse cells, and do not
produce an integer `m` or a global impossibility proof.

## Period-864 closure and deeper pure-(2,3) search (2026-07-25)

The period-864 problem is now exactly closed.  After fixing the order-3
anchors `p=19,p=37,p=73` and the order-4 anchor `p=97`, translations
preserving those anchors independently normalize the two order-9 anchors
`p=109,p=15121` to their residues modulo 3.  Exhaustive replay checks all
8,748 assignments on these six anchors and reduces them to 18 canonical
cases: two `p=73` orbits times nine order-9 residue pairs.  CaDiCaL returned
UNSAT on the full `864 x 864` plane in every case.  The symmetry/result audit
is:

* `power1616615_augmented799_period864_symmetry_unsat_audit.json`.

MiniCard independently reproduced all 18 UNSAT results, including both the
concurrent and triangle orbits.  The same exhaustive normalization replay
passed against that second result set:

* `power1616615_augmented799_period864_symmetry_unsat_mc_audit.json`.

Thus the complete period-864 finite theorem is duplicated across CaDiCaL and
MiniCard, with independent encodings of each row's exactly-one target choice.

Independently replayed component reductions transfer the already certified
period-432, period-648, and period-864 obstructions to larger divisor-period
families having exactly the same stable core:

* period 432 also excludes `1296,15984`;
* period 648 also excludes `7992,14904,18792`;
* period 864 also excludes `1728,2592,59616,75168`.

The transfer records are
`power1616615_augmented799_period432_identical_core_transfers.json`,
`power1616615_augmented799_period648_identical_core_transfers.json`, and
`power1616615_augmented799_period864_identical_core_transfers.json`.  Each
source-to-core reduction also has a separate
`period*_component_core_verification.json` replay.

The strongest live periodic construction family uses every one of the 53
surviving rows whose modulus has no prime factor other than 2 or 3.  Their
least common period is `53,747,712`; the stable core has reciprocal density
`2.156166201084`.  Both `p=73` symmetry orbits are under simultaneous local
component-repair CEGIS and exact SAT-master CEGIS.  The exact checker is
diversifying counterexamples by the maximal coordinate components
`8192,6561`.  A separate 31-row period-15,552 family of density
`2.137860082305` remains active as a smaller exact target.

`exact_cegis.py` now exposes coordinate-diverse counterexample batches and
emits a JSON certificate when its finite master becomes UNSAT.
`lift_lattice_phases.py` and `verify_lattice_phase_lift.py` implement and
independently replay the exact 799-to-800 phase lift; the smoke lift with
the complementary `p=17` phase passed.

These are exact finite-family obstructions and active construction searches,
not a solution of Erdős problem 203.  There is still no integer `m` and no
global impossibility proof.

## Higher pure-(2,3) fibres and structural-density synthesis (2026-07-25)

Three selected-exponent gcd scans were completed and independently replayed:

* `107,495,424 = 2^14 3^8` contributed the new primes `71,663,617`
  and `11,609,505,793`;
* `161,243,136 = 2^13 3^9` contributed `120,932,353` and
  `483,729,409`;
* `214,990,848 = 2^15 3^8` contributed no row beyond the first scan.

All reported cofactors were resolved.  After the period-60 power cell and
the two independently replayed parity transformations, the first pair
survives the exact component cascade and expands the former 53-row
pure-(2,3) core to 55 rows.  The second pair is algebraically valid but is
removed as inessential by the same cascade.  The 55-row stable core has
period `53,747,712` and density `2.156166870880`.

A new construction scan, `scan_pure23_multipliers.py`, avoids factoring the
enormous common divisor directly.  It exhaustively tests the declared finite
box `p=t*h+1`, where `h` ranges over all 143 pure-(2,3) values from
`1,100,001` through `500,000,000` and `1 <= t <= 1000`.  It found 20 exact
rows.  `verify_order_range.py` independently checked primality, both
multiplicative orders, the affine logarithm relation, and the absence of
duplicates against the 1.05-million pool:

```
PASS rows=20 sha256=2a362bf00989baab1eca1eb629abc6e009845730ceec8befc0dedefe5757f476
```

After power and parity conditioning, 17 of those rows are new relative to
the 801-row branch.  Seven survive the unrestricted component cascade.  In
the old period, the genuinely new row

```
p=6879707137, h=995328
```

survives while the two rows from the `161,243,136` scan peel away.  This
gives a verified 56-row core of density `2.156167875574`.  Expanding only
the binary period component to `2^16` retains all seven stable new rows and
initially gave a 62-row core of period `429,981,696`.

A second complete declared multiplier box covers all 89 pure-(2,3) values
from `500,000,001` through `10,000,000,000`, again with
`1 <= t <= 1000`.  It found nine more exact rows, independently replayed
with hash
`3cb709650644324ab60e778646a3bc8356aaf589a74135b5e69682a214a0b840`.
After the same power/parity chain,

```
p=3439853569, h=35831808
```

fits the existing period and expands that core to 63 rows, density
`2.156178792318`, without increasing its largest checker component.  A
larger 67-row core of period `6,879,707,136` is also stable, but its
2-primary checker component is `2^20`; the 63-row period-`429,981,696`
family is therefore the stronger practical exact-synthesis target.
Independent component-core replay passes for the 56-, 63-, and 67-row
artifacts.

The current 55-row local phase passes the exact necessary residual-density
condition on the 3-adic component but fails it on the 2-adic component: an
exact Z3 witness has scaled residual density `5797 < 8192`.  This proves
that phase assignment cannot be a cover.  `component_density_cegis.py` now
has an exact non-expanded Z3 weighted master for synthesizing against these
structural violations; its SAT outcomes remain only necessary-condition
hints, while an eventual master-UNSAT result would be a finite-family
obstruction.

Local exact-counterexample searches have been migrated to the stronger
63-row family.  None has yet returned checker-UNSAT.  The results above are
verified prime discovery and finite-family strengthening, not a cover of the
conditioned cell, not a resolution of the surrounding 1,711 cells, and not
an integer `m` or a global impossibility proof.

## Sparse density masters and higher pure-(2,3) boxes (2026-07-26)

The next two multiplier boxes are now closed on their stated finite domains.
For the 29 pure-(2,3) values from `10,000,000,001` through
`25,000,000,000`, the range `1001 <= t <= 10000` tested 256,500 odd
candidates and found no additional joint-order row.  Its independent replay
hash is:

```
ab14a6cb6cfcab5d2a3dd6f100c7fa31fd30bba031c39ba9c8d6e5587e7798e9
```

For all 129 pure-(2,3) values from `25,000,000,001` through
`1,000,000,000,000`, the range `1 <= t <= 1000` tested 126,938 odd
candidates and found 17 exact rows.  Independent primality, order, and
affine-log replay passed with hash:

```
a9a9e937b278a13aaac342d26459146b5c0e6ba6d385ecd891c44e25e2b7cf8d
```

All 17 rows survive the power-cell and two parity transformations, with
separate replays passing at every stage.  Nine survive the unrestricted
component cascade.  Together with the previous 79-row pure-(2,3) family they
give a verified 88-row stable core of period
`190,210,142,896,128 = 2^30 3^11` and density `2.156195732810`.
This density increase is real but very small; these high-order fibres are
structurally useful rather than a substantial reciprocal-density gain.

The preceding 79-row core has period
`11,888,133,931,008 = 2^26 3^11`, density `2.156195729291`, and a
separately replayed component-core audit.  A finite seed phase fails the exact
2-primary residual-density condition with scaled density
`30,623,234 < 67,108,864`; after adding the nine stable rows, the inherited
88-row seed likewise fails with
`472,514,624 < 1,073,741,824`.  These are exact counterexamples to those
phase assignments, not phase-independent obstructions.

`component_density_cegis.py` now also has an exact sparse-target
SciPy/HiGHS MILP master.  For each learned coarse-cell cut it materializes
only targets that occur in at least one cut, plus one sentinel representing
all still-unrepresented targets; those targets have identical zero
coefficients on the accumulated finite cut set.  On the 79-row core, the
first 20-cell master used 908 binary variables and solved optimally in
`0.316` seconds.  On the 88-row core it used 1,002 variables and solved in
`0.153` seconds.  Exact Z3 checkers remain responsible for finding every
new violated cell.

As an independent bounded benchmark, the complete period-432 concurrent
2-primary density master contains all 729 cells of its `27 x 27` coarse
plane.  The weighted Z3 formulation returned UNSAT in `152.2` seconds, and
the sparse SciPy/HiGHS MILP formulation independently returned infeasible in
`109.9` seconds.  This is a stronger necessary-condition obstruction for
that already-closed finite symmetry branch; it does not extend by itself to
the triangular branch or to a larger period.

The selected-exponent scan at `483,729,408` also completed with no unresolved
cofactor.  Its four exact rows independently replay with hash
`8f7a5751a7b66269edc791c1eb9c8bf97445b0cbbc023115b2289de5d8ea78fb`,
but all four duplicate earlier discoveries and therefore do not enlarge the
pool.

The 79- and 88-row residual-density CEGIS searches, both symmetry orbits of
the 63-row exact/local search, and the direct period-15,552 density master
remain active.  None has produced an exact cover or an all-phase
obstruction.  There is still no integer `m` and no global impossibility
proof.

## Period-15,552 density obstruction and full-core restart (2026-07-26)

The residual-density master has now closed the complete 31-row
period-15,552 family.  For the 2-primary component, every phase assignment
must give scaled line density at least `32` in each of the `243^2 = 59,049`
cells of the complementary 3-primary plane.  This is a necessary condition
for a cover, with no union-bound relaxation in the cell enumeration itself.

The six small anchors have moduli `3,3,3,4,9,9`.  Exhaustive translation and
negation replay checks all 8,748 assignments of their targets and reduces
them to 18 canonical cases:

* concurrent or triangular `p=73` orbit;
* each of the nine residue pairs for `p=109,p=15121` modulo 3.

An independent symmetry audit verifies every normalization and the fixed
targets recorded in every case artifact:

```
PASS assignments=8748 cases=18 engine=ortools-cp-sat
```

OR-Tools CP-SAT returned `INFEASIBLE` on the complete weighted all-cell
master in all 18 cases.  The concurrent cases took 56.7--129.6 seconds and
the triangular cases 211.0--288.4 seconds under parallel load.  Therefore no
phase assignment of the supplied 31 rows can meet even this necessary
2-primary density condition, and hence no affine cover exists in the entire
period-15,552 family.  The audit is:

* `power1616615_augmented799_period15552_density2_cpsat_symmetry_unsat_audit.json`.

This strengthens and supersedes the stalled local/exact construction searches
on that period.  HiGHS and Z3 replays of a canonical case are still running;
until one finishes, the new obstruction has one exact master engine plus the
separate symmetry replay, rather than two independent UNSAT engines.

The CP-SAT implementation was first regressed on the unrestricted
period-432 master.  It reproduced the established UNSAT result in `0.510`
seconds with 249 target variables and all 729 cells, matching the earlier
Z3 and HiGHS formulations.

The strongest live construction search has also been restored.  The
proof-preserving full-component conditioned core has 9,261 fibres, reciprocal
density `7.478463834429`, and largest exact-checker component `16,384`.
Two saved phases have been resumed against exact counterexamples.  The
ordinary branch absorbs batches of 2,000 coarse-coordinate-diverse holes.
The stronger triple-margin branch first had only about 1,100 hole cells at
the `(16,27)` fingerprint level; it is now learning batches of 5,000 at the
finer `(256,243)` level.  These are active exact-checker CEGIS searches.
Every emitted point is genuine, but no checker has returned UNSAT and no
cell-cover artifact exists.

Finally, the formerly missing multiplier box on all 89 pure-(2,3) orders
from `500,000,001` through `10,000,000,000` and
`1001 <= t <= 10000` is complete.  It tested 792,000 odd candidates and
found no row; independent replay passed with hash
`0804981c69833a1bf9d31ddf2166a83281cd0b156be98e3ac193adaabfa2f7ec`.

The period-15,552 result is a rigorous finite-family theorem.  It does not
cover the first persistent cell, settle the 1,711 surrounding cells, produce
an integer `m`, or prove global impossibility.

## Maximal binary-period obstruction and independent Boolean-Z3 replay (2026-07-26)

The complete divisor family available in the verified 834-row pure-(2,3)
core with period

```
2^30 * 3^5 = 260,919,263,232
```

contains 60 rows of reciprocal density `2.153234145754`.  The
proof-preserving component cascade removes 56 rows in seven rounds: eight
parallel-class capacity records and one Blokhuis--Brouwer tangent-capacity
record.  The four survivors have moduli `4,8,32,32` and exact total density

```
1/4 + 1/8 + 1/32 + 1/32 = 7/16 < 1.
```

The independent component verifier now optionally combines replay of every
elimination with an exact rational survivor-density check.  On this family it
reports:

```
PASS input=60 rounds=7 records=9 survivors=4 tangent=1 parallel=8
     inessential=0 density=7/16 proved_no_cover=True
```

The combined certificate is
`power1616615_augmented846_pure23_multiplier_t1000_period260919263232_component_core_density_obstruction_verification.json`.
Thus no phase assignment of any of these 60 rows covers the exponent lattice.
This single theorem subsumes all divisor subfamilies with binary exponent at
most 30 and ternary exponent at most 5.

The obstruction ends sharply at the next ternary layer.  For periods
`2^30*3^b`, independently replayed stable cores have:

* `b=6`: 65 rows, density `2.154649425911`;
* `b=7`: 73 input rows reduce to the same 65-row `b=6` core;
* `b=8`: 78 input rows reduce to 76 rows, density `2.156177184894`;
* `b=9` and `b=10`: 81 rows, density `2.156186513985`.

The complete all-cell 2-primary density master for the first surviving
65-row `b=6` family has `729^2=531,441` coarse cells.  Its canonical
concurrent symmetry case is active under OR-Tools CP-SAT.  SAT would be only
a necessary-density phase hint; INFEASIBLE would be a finite-family
obstruction after all symmetry cases are checked.

A separate direct solver,
`component_density_direct_z3bool.py`, represents every row target by Boolean
variables and uses Z3 native exactly-one and weighted pseudo-Boolean
constraints.  It independently reproduced the unrestricted period-432
density obstruction in `4.603` seconds.  On the canonical concurrent
period-15,552 case it returned `INFEASIBLE` in `244.338` seconds, independently
confirming the CP-SAT encoding on that case.  The remaining symmetry cases
are being replayed in controlled parallel waves.

The two full-component construction searches also remain active.  The
ordinary 9,261-row branch has repaired 27,988 exact witnesses; the
triple-coverage branch has repaired 37,107.  Each exact checker continues to
find a full fresh batch, so neither phase is a conditioned-cell cover.

These are rigorous finite-family advances and live construction searches.
There is still no integer `m` and no global impossibility proof.

## Smooth-log extension through order 10^15 and strengthened cell core (2026-07-26)

The pure-(2,3) multiplier scanner previously used a generic
baby-step/giant-step logarithm.  At orders near `10^15` this requires roughly
`sqrt(h)` memory and reached 3.3 GB before the incomplete run was stopped.
The scanner now uses Pohlig--Hellman digit lifting specialized to
`h=2^a*3^b`.  It computes each logarithm with `a+b` modular exponentiations
and a final exact CRT check.  Replaying all 48 rows from the four earlier
`t<=1000` boxes reproduced every stored logarithm exactly:

```
POHLIG_HELLMAN_REPLAY_PASS rows=48
```

With that correction, the complete declared box

```
10^12 < h <= 10^15,  1 <= t <= 1000,  p=t*h+1
```

contains 290 pure-(2,3) orders and 286,819 tested odd candidates.  It produced
37 exact fibres.  Independent primality, multiplicative-order, signature, and
affine-log replay passed with hash:

```
20f8381ac182128b849efef6a8c3d1b23e474288809fcdc68bcaf13bf08561b0
```

All 37 fibres survive the `D=1,616,615` power cell and both parity
transformations, with a separate verifier passing at each stage.  Thirty of
the new fibres survive the unrestricted component cascade; their presence
also makes two older rows essential again.  The final parity-child core grows
from 834 to 866 rows, density `3.241192138354`.  An entirely separate route
through the parent cell core, followed by the two verified parity
transformations, produces exactly the same ordered 866-row data.

The strongest pure-(2,3) divisor core now has 120 rows, period

```
2^43 * 3^15 = 126,214,320,739,011,526,656
```

and density `2.156196571010`.  Its independent component replay removes no
row.  The maximum current family with ternary exponent at most five has 74
rows and period

```
2^43 * 3^5 = 2,137,450,604,396,544.
```

The component cascade removes 70 rows in eleven rounds and again leaves the
same four-row `7/16` core.  The combined verifier reports
`proved_no_cover=True`, extending the earlier 60-row obstruction through all
newly available binary layers.

At ternary exponent six, the maximal current family has 80 rows and period
`2^43*3^6`.  Its density exceeds the former 65-row family by only
`1.92e-10`.  The former 65-row canonical all-cell density master returned
`INFEASIBLE` over all `729^2` cells in `673.064` seconds.  The corresponding
80-row superset master is active; remaining 65-row cases were stopped as
dominated.

The complete full-cell proof chain was also rerun from the pre-core rather
than appending to a reduced artifact:

* 19,065 conditioned rows reduce to a 12,952-row algebraic component core;
* the independently replayed lower-digit peel removes 313 rows, leaving
  12,639;
* the ordinary cascade leaves an equivalent 868-row parent core, density
  `3.618364190774`;
* the two parity transformations give the independently identical 866-row
  child core.

Filtering the 868-row parent to prime-power components at most 16,384 and
replaying the cascade leaves 804 rows.  Unioning those with the older
9,261-row construction pool adds 135 genuinely absent usable fibres, giving
9,396 rows and density `7.480677913166`.  Exact CEGIS has restarted on this
stronger pool while preserving 37,988 ordinary witnesses and 62,107
triple-coverage witnesses.

A live status refresh of the Erdős Problems database still lists Problem 203
as open, with no partial or complete solution claimed in its comments.
No external claim is being treated as a solution.

There is still no integer `m`, no complete conditioned-cell cover, and no
global impossibility proof.

## 64-bit order extension and high-multiplier row (2026-07-26)

The deterministic pure-(2,3) scan was extended through every declared order
in

```
10^15 < h <= 18*10^15,  1 <= t <= 1000,  p=t*h+1.
```

It checked 138,445 odd candidates over 140 orders and found 14 exact fibres.
The independent replay hash is

```
6c7f521ae8039e2976a3ae4d6b13d8be188fe898881dce31afc379e8e6421d76
```

All fourteen rows pass the power-cell, rectangle, and lattice verifiers; five
survive the unrestricted child component cascade.

A separate complete box

```
10^12 < h <= 1.8*10^15,  1001 <= t <= 10000,  p=t*h+1
```

checked 2,825,985 candidates over 317 orders and found one exact fibre:

```
p = 22,799,473,113,563,137
h = 2,783,138,807,808
t = 8192.
```

Its independent replay hash is

```
4e2299e6970f248df49393579e98ca296651727c82a08f1f7f1b2caeb97742bc
```

The fibre passes all three transformation verifiers and survives the child
component cascade with child order `2^33*3^3 = 231,928,233,984`.

After incorporating this row, the maximal `2^43*3^15` divisor core has 126
stable rows.  The `2^43*3^5` family has 75 input rows; an independent replay
of eleven proof-safe cascade rounds leaves four rows with exact reciprocal
density `7/16`, again certifying that this entire finite family cannot cover.
The `2^43*3^6` family now has 82 stable rows.  The superseded 81-row exact
master was stopped and the corrected 82-row master started.

On the constructive side, the 9,396-row exact CEGIS pool has reached 41,988
ordinary witnesses and 67,107 triple-coverage witnesses.  Both candidates
repair the current witness sets, but their exact checkers continue to produce
fresh misses; neither is a cover.

These computations strengthen both the finite obstruction and the positive
search.  They still do not supply an integer `m` or a global impossibility
proof.

## Closest-family conditional-overlap checkpoint (2026-07-26)

For the surviving period

```
2533395664800
```

the corrected 14-anchor block was combined with 22 conditional-fibre
intersection certificates.  Each conditional certificate has a separate
verification artifact, and the assembled period certificate was replayed by
`verify_ranked_period_conditional_star.py`.

The exact current upper bound is

```
730799865304831267447 / 729328817548572040800
= 1.0020169883883154
```

so the remaining gap to the no-cover cutoff is exactly

```
1471047756259226647 / 729328817548572040800
= 0.00201698838831534.
```

The newest `p=103` conditional certificate contributes an independently
verified improvement of

```
152102 / 1882330695
= 0.0000808051424779003.
```

This is the strongest certified bound currently assembled for that finite
family, but it does not eliminate the family because the upper bound remains
greater than one.  It also does not imply a corresponding percentage of
progress on the unrestricted Erdős problem.

## Period-2,533,395,664,800 elimination (2026-07-26)

The preceding checkpoint has now been superseded by a strict no-cover
certificate.  The corrected 14-anchor block was first extended by the
`p=31` and `p=37` fibres.  The resulting 16-anchor exact block has forced
overlap loss

```
146904738209 / 456011219664
= 0.32215158722902243.
```

The conditional projected-pair generator and its independent verifier now
support appended base-period anchors.  They quotient adversarial base targets
by the exact translation stabilizer; the verifier reconstructs the
stabilizer independently by residue lifting and checks that its orbits are
disjoint and exhaustive.  On the `p=31` regression this reduces 7,372,800
raw target tuples to 460,800 orbit representatives without changing the exact
answer.

Forty independently replayed conditional-fibre edges were assembled around
the enlarged block.  The decisive new `p=599`, `h=299` certificate uses a
compact period-17,940 decomposition with 1,076,400 outside base points and
63,590,400 target states.  Its exact forced intersection is

```
370879 / 185472690
= 0.001999642103643399,
```

improving the previous generic edge by

```
18630527429 / 81575340406560
= 0.00022838430506263383.
```

The complete period certificate and the independently recomputed replay now
give

```
40107466081993334654251 / 40113084965171462244000
= 0.9998599239329758 < 1.
```

The no-cover margin is approximately `0.00014007606702417014`.  Therefore no
choice of phases from this 933-row divisor-period family can cover the full
exponent lattice.  The ranked finite frontier is reduced from 28 to 27
families: 4,610 of 4,637 are now certified no-cover.  All four changed or new
certificate programs compile, all 30 regression tests pass, and the complete
period replay reports `verified=True`.

This is a rigorous finite-family elimination only.  It does not produce an
integer `m`, rule out all possible finite prime-fibre systems, or prove a
global impossibility theorem for the original Erdős problem.

## Period-101,264,763,600 elimination (2026-07-26)

The next ranked survivor began with aggregate block-star upper bound
`1.0096749113138785`.  Reusing every compatible verified conditional edge
lowered it only to `1.0061146977205686`.

A stronger transformation promotes an outside fibre into the exact block.
For this purpose, `conditional_block_extension_v2` permits the intersection
witness to use any recorded subunion of the existing block: an intersection
lower bound against a subunion is automatically a lower bound against the
whole block.  The independent verifier checks that every conditional anchor
is an actual base-block anchor, that all recorded affine rows match the stable
pool, and that the supplied conditional report independently verified the
same exact fraction.  Legacy complete-block extension certificates remain
replayable under the original schema.

Starting from the 14-anchor projected endpoint-path block, the independently
verified promotion chain adds:

```
p=71  with intersection 659/41580
p=31  with intersection 4/297
p=191 with intersection 659/112860.
```

Promoting `p=71` reduced the assembled upper bound to
`1.0011385289539028`.  Promoting `p=31` next gave
`1.0010877783324745`.  The final `p=191` promotion, together with 13
independently verified conditional-fibre edges around the 17-anchor block,
gives

```
2102349443216882789 / 2103497917325925120
= 0.9994540169973154 < 1.
```

The exact no-cover margin is

```
1148474109042331 / 2103497917325925120
= 0.0005459830026845619.
```

The complete 759-row period certificate independently replays with
`verified=True`.  The ranked finite frontier is therefore reduced again:
4,611 of 4,637 families are certified no-cover and 26 remain.

This remains a finite divisor-family theorem, not a value of `m` and not a
global impossibility proof.
