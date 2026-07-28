# Correctness boundary for the affine Hunter forest

> **Status: finite-family obstruction only.** This note proves the overlap
> bound emitted by `inventory_reverse_group_lines.py`. It does not prove a
> global statement about Erdos Problem 203 and does not construct an integer
> \(m\).

## Joint affine maps

Let

\[
G=\mathbb Z/U\mathbb Z\times\mathbb Z/V\mathbb Z
\]

and let two descending inventory rows define surjective homomorphisms

\[
\phi_i(x,y)=a_i x+b_i y\pmod {h_i},\qquad
\phi_j(x,y)=a_j x+b_j y\pmod {h_j}.
\]

Their joint map is

\[
\Phi=(\phi_i,\phi_j):
G\longrightarrow
\mathbb Z/h_i\mathbb Z\times\mathbb Z/h_j\mathbb Z.
\tag{1}
\]

Before quotienting the domain, the cokernel of the corresponding map from
\(\mathbb Z^2\) has presentation matrix

\[
M=
\begin{pmatrix}
a_i & b_i & h_i & 0\\
a_j & b_j & 0 & h_j
\end{pmatrix}.
\tag{2}
\]

The Smith normal form theorem says that the order of this cokernel is the
greatest common divisor of the \(2\)-by-\(2\) minors of \(M\). Thus the
joint-image index is

\[
\begin{split}
d_{ij}=\gcd(&a_i b_j-b_i a_j,\ h_i a_j,\ h_i b_j,\\
            &h_j a_i,\ h_j b_i,\ h_i h_j).
\end{split}
\tag{3}
\]

The descent conditions make both rows factor through \(G\), so (1) is
surjective exactly when \(d_{ij}=1\). Coprime row moduli imply this condition,
but (3) also detects independent maps with noncoprime moduli.

If \(d_{ij}=1\), every target pair
\((c_i,c_j)\) has exactly

\[
\frac{|G|}{h_i h_j}
\]

preimages. Therefore every legal choice of the two row phases forces an
intersection of density \(1/(h_i h_j)\). This remains true when the phase
domains are restricted: surjectivity supplies the intersection for every
target pair, not merely for unrestricted or random targets.

## Forest inequality

Let \(A_i\) be the fibre selected by row \(i\). For any forest \(F\) on the
rows, Hunter's inequality gives

\[
\frac{|\bigcup_i A_i|}{|G|}
\leq
\sum_i \frac1{h_i}
-
\sum_{\{i,j\}\in F}
\frac{|A_i\cap A_j|}{|G|}.
\tag{4}
\]

One direct proof roots each tree and adds vertices after their parents.
When \(A_i\) is added, its new contribution is at most
\(|A_i|-|A_i\cap A_{\operatorname{parent}(i)}|\). Summing over every tree
proves (4); summing the component bounds proves the forest form.

Build a graph whose vertices are source-prime row slots. Join \(i,j\) exactly
when (3) equals one, and give that edge weight \(1/(h_i h_j)\). All weights
are nonnegative, so a maximum-weight spanning forest gives the strongest
bound obtainable from these certified pairwise intersections:

\[
B=
\sum_i\frac1{h_i}
-
\max_F\sum_{\{i,j\}\in F}\frac1{h_i h_j}.
\tag{5}
\]

The inventory producer computes the forest with Prim's algorithm. The
independent verifier:

1. reconstructs every descending row from the authenticated source;
2. recomputes (3) directly;
3. checks that every emitted edge is forced and that the emitted edges are
   acyclic;
4. recomputes the optimum with Kruskal's algorithm; and
5. checks every exact fraction and finite-group fibre count.

If the independently verified value \(B\) is strictly below one, the
declared rows cannot cover the declared group for any phase assignment.
This eliminates only that finite arithmetic family. It says nothing about
rows omitted from the source pool or about all possible primes.
