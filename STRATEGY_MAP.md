# Strategy map

> **Status: unresolved.** None of the routes below has produced the requested
> integer \(m\) or a proof that no such integer exists.

For every usable prime \(p\), divisibility by \(p\) selects one affine fibre

\[
  a_p k+b_p\ell\equiv c_p\pmod {h_p}.
\]

A finite collection of fibres that covers \(\mathbb Z^2\), with no prime
reused, gives compatible congruences for \(m\) and therefore a construction
by the Chinese remainder theorem. This is the proof boundary for every
constructive route below.

## Ranked routes

The numbering is execution order, not an estimated probability of resolving
the original problem. The pseudo-Boolean route is closest to a finite theorem;
the arithmetic-matched reverse design is the preferred new route for actually
finding \(m\).

### 1. Proof-producing pseudo-Boolean closure

This is the near-term finite-family route. The current 81-row layered family
has one active-class variable per legal layer residue. Exact integer
inequalities encode:

- one active class per row;
- capacity at every distinct layer-column pattern;
- exact modulo-30 dual obstructions, lifted monotonically to supersets; and
- unavoidable coprime-residual overlap cuts.

The generated OPB instance contains only necessary conditions. A SAT answer
is merely a placement that survives those conditions. An UNSAT answer becomes
a theorem about this one finite family only after a proof log and the
independent exact reconstruction both verify.

This route is likely to settle the current family efficiently, but it cannot
by itself settle the original problem.

### 2. Symbolic CRT game instead of point CEGIS

The exact checker already decomposes \(k\) and \(\ell\) into independent
prime-power CRT components. Include the phase choices as existential
variables and the exponent-component digits as universal variables:

\[
  \exists(c_i)\ \forall(k_q,\ell_q)\quad
  \bigvee_i [a_i k+b_i\ell=c_i\pmod {h_i}].
\]

This can be attacked as QBF, or as a BDD/MDD and transfer-matrix game over
the component factor graph. The crucial change from the current direct
search is the learned object: a symbolic component cube or an exact
separator replaces a single uncovered exponent pair.

For the smallest new survivor,
\(G=(\mathbb Z/5040\mathbb Z)^2\), the universal point has only 28 binary
CRT-coordinate bits:

\[
5040=2^4\cdot3^2\cdot5\cdot7.
\]

Each row phase is likewise split into independent one-hot prime-power
residues. Precomputed lookup gates evaluate
\(a_i k+b_i\ell\bmod q^e\) on each local component; a row matches when all
its component gates match, and the cover gate is the disjunction of the 31
rows. This instance has 83 phase components, 531 one-hot phase variables,
and only 2,037 pairwise at-most-one clauses before circuit sharing. Invalid
binary encodings of the 3-, 5-, and 7-adic digits satisfy the universal
implication automatically. This represents 25,401,600 quotient cells with
28 universal bits rather than 25 million point clauses.

`build_crt_cover_qcir.py` implements this deterministic lookup-gate
translation from an authenticated arithmetic inventory. Its manifest keeps
all solver claims `UNRUN`/false. An independent reconstruction verifier and
a proof-capable QBF run are still required. The emitter uses the cleansed
integer-identifier form of the
[QCIR-G14 specification](https://www.qbflib.org/qcir.pdf).

The conservative proof route is QCIR-to-PCNF conversion, then
[DepQBF](https://lonsing.github.io/depqbf/) with an ASCII or binary QRP trace,
followed by
[QRPcheck/QBFcert](https://fmv.jku.at/qbfcert). The conversion, QCIR
reconstruction, and QRP proof all remain separate gates; an unlogged solver
answer is telemetry only.

In a Boolean circuit encoding, phase bits belong before the universal CRT
digits, while match/Tseitin auxiliaries that depend on those digits must
belong after them. Getting this dependency order wrong would create a
different and unsound synthesis problem.

The admission benchmark is deliberately strict:

1. reproduce the known 14-row no-cover result;
2. emit a checkable UNSAT proof or an independently replayable symbolic
   certificate;
3. decide the first survivor not already removed by exact arithmetic;
4. demonstrate that learned cubes remove substantially more phase space
   than point clauses before scaling to the 12,579-row branch.

The initially selected 31-row `5040 x 5040` target is not a genuine QBF
benchmark after all. For a primitive row on a square quotient, descent is
equivalent to `h | 5040` and is invariant under a unimodular shear. It is
therefore the older 31-row divisor-period family of density `143/140`.
The modulus-4 `p=5` and modulus-6 `p=7` maps are jointly surjective despite
their noncoprime moduli, forcing overlap `1/24` and the exact union bound

\[
143/140-1/24=823/840<1.
\]

`verify_square_quotient_period_bridge.py` authenticates the source-row
identity, transports the two rows through the shear, and independently
enumerates their 144-cell target map. The real replay returns
`verified=true`. The independently reconstructed maximum affine forest uses
30 forced edges and improves the same finite-family upper bound to
`21229/25200 < 1`. The QBF target is therefore rigorously retired.

### 3. Reverse-designed finite-group cover with arithmetic matching

Work first in a small sheared finite abelian group, initially
\(\mathbb Z/120\mathbb Z\times\mathbb Z/60\mathbb Z\), and synthesize a
cover by affine line slots. This reverses the current workflow: geometry
proposes the covering system, then arithmetic supplies its primes.

An unconstrained abstract cover is not enough. For each line type
\((h,a,b)\), its multiplicity is capped by the distinct source primes that
actually realize that signature. A bipartite matching must assign every
selected slot to a distinct prime before the design is promoted. Rows that
do not descend to the chosen quotient may appear only through an explicitly
labelled relaxation.

The quotient is then lifted one prime-power component at a time. At each
lift, the design must cover every child cell and preserve the prime matching.
Only a full lift that passes both exact lattice checkers may proceed to CRT
construction of \(m\).

This is the highest-upside route for finding \(m\), but it may fail quickly
because useful abstract line multiplicities are arithmetically scarce.

`inventory_reverse_group_lines.py` implements the exact preflight. It applies
the declared unimodular shear, retains only predicates that genuinely descend
to the finite group, canonicalizes projective directions, and reports the
actual distinct-prime capacity and reciprocal density of every line type. It
also detects fibres that are CRT rectangles. If \(h=uv\), \((u,v)=1\), and
the coefficient pair isolates \(k\bmod u\) and \(\ell\bmod v\), then crossing
two ordinary one-dimensional covering systems becomes a finite
distinct-prime matching problem.

The same inventory constructs a maximum-weight forest of phase-independent
intersections. Two affine rows are joined whenever the Smith-normal-form
index of their joint map is one; this includes every coprime-modulus pair
and also independent rows with noncoprime moduli. Hunter's inequality gives
the exact bound

\[
 \left|\bigcup_i L_i\right|/|G|
 \leq \sum_i 1/h_i-\sum_{\{i,j\}\in F}1/(h_i h_j).
\]

This can reject a quotient family without choosing phases or enumerating
cells. `verify_reverse_group_density_obstruction.py` reconstructs that
particular finite-family claim directly from the authenticated source pool.

There is also a proof-level probabilistic preflight. For \(N\) quotient cells
and unrestricted distinct-prime rows of indices \(h_i\), independent uniform
phases leave

\[
  N\prod_i(1-1/h_i)
\]

cells uncovered in expectation. If this exact quantity is below one, some
phase assignment is a complete cover. Conditional expectation can then make
the construction deterministic. A passing result still requires independent
inventory replay and an explicit phase/CRT artifact before an integer \(m\)
is claimed.

The simple product is used only for unrestricted, independently phased,
distinct-prime rows whose \(h_i\) targets partition every quotient point.
Restricted phase domains require the cell-dependent expression
\(\sum_g\prod_i(1-p_{i,g})\); the current preflight conservatively ignores
restricted rows rather than applying the uniform formula to them.

`construct_reverse_group_cover.py` implements that deterministic phase
choice. It also emits exploratory greedy assignments when the sufficient
first-moment inequality does not pass, but only a zero-miss assignment that
survives the independent full-domain and arithmetic checks can be promoted.

#### First reverse-design pilot

A 64-configuration scan crossed six small coordinate bases/periods with
eight transverse quotient sizes. Sixty-three families had raw reciprocal
density below one and therefore could not cover before phase search. The
only survivor was the `(3,1)` basis with
\(G=\mathbb Z/15120\times\mathbb Z/2520\): it had 28 descending rows and
raw density

\[
337/336.
\]

Its modulus-4 and modulus-35 rows have coprime indices, so every pair of
their phases intersects in density \(1/140\). Even subtracting only this
one forced overlap gives

\[
337/336-1/140=239/240<1.
\]

Using all nine coprime neighbors of the modulus-35 row as a Hunter star
strengthens the same bound to \(59/60\). In that pilot the star already is
the complete positive-overlap forest.

Thus the entire 38,102,400-cell pilot family is ruled out by exact quotient
arithmetic, with no graph placement or cell enumeration. This is a useful
eliminator, not a result about all possible quotient designs or the original
infinite problem.

`scan_reverse_group_quotients.py` makes the search reproducible and keeps it
in this arithmetic space. The first broader exact preflight over 16 small
bases and the declared smooth-period grid found only seven distinct modulus
multisets surviving the older Hunter-star test. The strongest had 33 rows,
raw density \(2609/2520\), and star-adjusted upper bound \(341/336\); its
random phase expectation was still enormous. The scanner now uses the
strictly stronger maximum joint-surjectivity forest. Its full rescan remains
deferred at this stopping checkpoint; the square `5040 x 5040` candidate is
already eliminated, so the next scan should focus QBF work only on surviving
nonsquare families.

### 4. Library of one-dimensional covering templates

Choose a layer period \(L\). For each residue of \(\ell\bmod L\), the active
rows induce an ordinary covering system in \(k\). Instead of allocating rows
by max-min density, enumerate or certify a library of complete one-dimensional
covers and solve a matching problem that assigns compatible templates to all
layers without reusing primes.

This is mathematically a structured subcase of the layered search, but its
state space is a library of whole covers rather than individual row phases.
It is promising if a small number of reusable template types exists; it
should be abandoned if the library has high entropy and no repeated cores.

### 5. Fourier or exact LP hierarchy on quotient groups

An affine-fibre indicator has sparse Fourier support on a finite abelian
quotient. Exact dual weights combining several CRT projections can certify
that a finite family cannot cover even when every single projection looks
feasible. This generalizes the successful modulo-30 weighted obstruction.

This is primarily an obstruction engine. It can prune finite families and
produce rational certificates, but the primal side must still return an
actual phase assignment before it contributes to a construction.

### 6. Perfect-power algebraic coverage

Setting \(m=M^D\) makes some exponent sublattices automatically composite by
identities such as \(X^q+1\), leaving a smaller affine-cover remainder. This
route already has extensive exact experiments in the repository. It remains
valid, but further point-by-point continuation is low priority; it becomes
attractive again only if the remaining cells can be closed symbolically or
matched to a small reverse-designed template.

There is also a useful completeness warning here. Capelli's reducibility
criterion for \(X^n-a\), specialized to a positive binomial \(X^n+A\) over
the rationals, leaves exactly two universal mechanisms: an odd-power
factorization or the exceptional \(A=4B^4\) case behind Sophie Germain's
identity. Both are already encoded by the perfect-power checker. Thus a new
algebraic route must introduce genuinely non-binomial structure rather than
searching for another identity of the same form.

Reference:
[Koley--Reddy, *Irreducibility of \(x^n-a\)*](https://arxiv.org/abs/2006.03787).

### 7. Irredundant coset-cover theorem screens

The affine fibres are cosets of subgroup kernels in a finite abelian group,
so classical coset-cover theorems supply cheap necessary conditions before
any solver runs. Lettl and Sun prove that if a \(k\)-coset, multiplicity-one
cover contains an irredundant coset of index \(h\), then

\[
h\leq 2^{k-1},
\qquad
k\geq 1+\sum_{p^\alpha\parallel h}\alpha(p-1).
\]

A row is automatically irredundant in any hypothetical cover when deleting
it leaves raw density below one. This makes the theorem an exact algebraic
screen for sparse, high-index quotient families. It does not eliminate the
current 31--34-row survivors: their largest relevant Mycielski value is only
18 and their indices are at most 3,780. It therefore remains a supporting
filter, not a main route.

Primary references:
[Lettl--Sun, *On covers of abelian groups by cosets*](https://arxiv.org/abs/math/0411144)
and
[Szegedy, *Coverings of abelian groups and vector spaces*](https://arxiv.org/abs/math/0411244).

## Current decision policy

1. Finish and independently reconstruct the OPB master while respecting the
   host-memory floor.
2. If the current family is proof-certified UNSAT, preserve the proof and
   move construction effort to arithmetic-matched reverse design.
3. If it is SAT, pass the survivor to exact projected and full-column
   checkers; never call a necessary-condition survivor a cover.
4. Prototype the symbolic CRT game on the known 14-row benchmark before
   spending resources on a large QBF or decision diagram.
5. Do not resume unbounded point-graph accumulation unless a measured
   symbolic compression ratio justifies it.

## What would count as progress

- **Construction:** a distinct-prime affine cover accepted by both exact
  full-domain checkers, followed by a verified CRT value of \(m\).
- **Finite obstruction:** a self-contained exact certificate for one
  explicitly declared row family.
- **Global nonexistence:** a theorem covering all possible primes and
  signatures, not merely a bounded source pool.

Finite obstruction counts, solver telemetry, and sampled holes must not be
reported as probabilities that the original problem has or lacks a solution.
