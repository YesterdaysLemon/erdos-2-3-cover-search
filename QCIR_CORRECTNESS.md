# Correctness boundary for the symbolic CRT cover encoding

> **Status: encoding theorem only.** This note proves what a generated QCIR
> instance means. It does not report a solver result, a covering phase
> assignment, an integer \(m\), or a global theorem about Erdős Problem 203.

## Declared finite family

Let

\[
G=\mathbb Z/U\mathbb Z\times\mathbb Z/V\mathbb Z.
\]

For every inventory row \(i\), let

\[
\phi_i(x,y)=a_i x+b_i y\pmod {h_i}.
\]

The inventory admits the row only after checking

\[
h_i\mid a_iU,\qquad h_i\mid b_iV,
\qquad \gcd(a_i,b_i,h_i)=1.
\]

The first two conditions make \(\phi_i\) a well-defined homomorphism on
\(G\). The last makes it surjective. The current QCIR prototype additionally
requires every target modulo \(h_i\) to be arithmetically legal.

A phase assignment is a tuple \(c_i\in\mathbb Z/h_i\mathbb Z\). It covers
the declared finite group exactly when

\[
\forall (x,y)\in G\quad
\bigvee_i \left[\phi_i(x,y)=c_i\right].
\tag{1}
\]

## CRT decomposition of coordinates

Write

\[
U=\prod_q q^{u_q},\qquad V=\prod_q q^{v_q}.
\]

The Chinese remainder theorem gives

\[
G\cong
\prod_q
\left(
  \mathbb Z/q^{u_q}\mathbb Z
  \times
  \mathbb Z/q^{v_q}\mathbb Z
\right).
\tag{2}
\]

The exporter uses binary variables for each local coordinate in (2). A
domain such as \(\mathbb Z/9\mathbb Z\) needs four bits and therefore has
seven invalid bit strings. For each prime \(q\), the circuit's local-valid
gate is the disjunction of the equality gates for the genuinely legal
coordinate pairs. The global-valid gate is their conjunction.

Consequently:

- every element of \(G\) has exactly one valid binary encoding;
- every valid encoding denotes exactly one element of \(G\); and
- an invalid encoding makes the implication
  `valid -> covered` true without imposing a spurious cover condition.

## CRT decomposition of phases

For

\[
h_i=\prod_{q\mid h_i}q^{e_{i,q}},
\]

CRT gives a bijection

\[
c_i\pmod {h_i}
\longleftrightarrow
\left(c_{i,q}\pmod {q^{e_{i,q}}}\right)_{q\mid h_i}.
\tag{3}
\]

The exporter creates one existential one-hot block for every component on
the right of (3). An at-least-one gate and all pairwise at-most-one gates
make each block choose exactly one residue. Because all targets are
unrestricted, every tuple of component choices corresponds to exactly one
legal \(c_i\bmod h_i\), and vice versa.

## Lookup and match gates

Fix a row \(i\) and a component \(q^{e_{i,q}}\). For every valid local
coordinate pair, the exporter evaluates

\[
a_i x_q+b_i y_q\pmod {q^{e_{i,q}}}
\]

with integer arithmetic while constructing the file. It groups the local
point-equality gates by this residue. The bucket for residue \(r\) is
therefore true exactly when the local affine form equals \(r\).

The component-match gate is the disjunction, over \(r\), of

```
phase_component_is_r AND local_affine_bucket_is_r.
```

It is true exactly when the selected phase component equals the local affine
value. The row-match gate is the conjunction over all prime-power components
of \(h_i\). By (3), it is true exactly when

\[
\phi_i(x,y)=c_i\pmod {h_i}.
\tag{4}
\]

The coverage gate is the disjunction of all row-match gates. By (4), it is
true exactly when the represented group element is covered by at least one
declared row.

## Quantifier theorem

The generated formula has prefix

\[
\exists\{\text{phase bits}\}
\forall\{\text{coordinate bits}\}
\]

and matrix

\[
\text{phase-valid}\ \wedge\
\left(\text{coordinate-valid}\Longrightarrow\text{covered}\right).
\tag{5}
\]

### Theorem

Subject to the authenticated inventory checks above, the generated QCIR
formula is true if and only if the declared affine rows admit a phase
assignment covering \(G\).

### Proof

Suppose the QCIR formula is true. Its existential assignment satisfies every
one-hot constraint, so (3) yields one legal target \(c_i\) for every row.
Every element of \(G\) has a valid binary encoding. Applying the universal
part of (5) to that encoding makes `covered` true, and (4) proves (1).

Conversely, suppose targets \(c_i\) satisfy (1). Encode their prime-power
components using the existential one-hot blocks. Phase-valid is true by
construction. For any universal bit assignment, an invalid coordinate makes
the implication in (5) true. A valid assignment represents an element of
\(G\), which satisfies (1), so one row-match gate and hence `covered` is
true. Thus (5) holds for every universal assignment. \(\square\)

## From quotient to the complete exponent lattice

Let

\[
\pi:\mathbb Z^2\longrightarrow G
\]

be coordinate reduction modulo \(U,V\). The descent conditions imply that
every affine predicate factors through \(\pi\). Therefore the declared rows
cover \(G\) if and only if their periodic lifts cover every
\((k,\ell)\in\mathbb Z^2\), and hence every \((k,\ell)\in\mathbb Z_{\ge0}^2\).

This equivalence concerns only the declared finite row family. A false QCIR
formula rules out that family, not every possible prime divisor. A true
formula supplies phase data only; it still requires:

1. extraction of the outer existential phase assignment;
2. conversion from canonical projective targets back to each source row;
3. distinct-prime and arithmetic target replay;
4. Chinese-remainder construction of \(m\);
5. independent full-lattice replay; and
6. direct verification that every listed prime divides the intended terms.

Only after those steps could a true instance contribute a candidate \(m\).

## Proof-artifact gate

The accepted negative-result chain is:

1. authenticate and independently reconstruct the arithmetic inventory;
2. independently reconstruct the cleansed QCIR file with
   `verify_crt_cover_qcir.py`;
3. convert the circuit to a quantified CNF with recorded hashes;
4. solve with proof logging;
5. verify the QBF proof independently; and
6. promote only the explicitly declared finite-family noncover statement.

An unlogged `UNSAT`, timeout, sampled search, or `UNKNOWN` result is not a
proof.
