# Public snapshot manifest

This repository is a curated snapshot of the active workspace, updated on
2026-07-27.

## Included

- all Python source and regression tests present at publication time;
- the complete research log, including superseded approaches and corrections;
- the streaming low-memory repair engine, validated restart cache, exact
  component-digit tile-union generator, adversarial witness-diversity
  checker, and their regression tests;
- the stable 1,050,000-order candidate pool used by the current finite sweep;
- the 14,629-row max-128 pool used by the current exploratory direct-cover
  search;
- the exact bounded homogeneous-refinement sibling obstruction over all
  129,497 raw rows, together with its independently oriented replay;
- the explicit exponent pair outside all 129,497 bounded homogeneous fibres,
  together with an independent modular-exponentiation replay over every
  source prime;
- the self-contained 14-row, period-5,544 affine-subpool no-cover
  certificate, its 13,908 witness points, and the independent exhaustive
  replay of all 746,496 legal phase assignments;
- the self-contained five-anchor frozen-remainder quotient certificate and
  independent replay of all 216 legal anchor phase assignments against all
  12,577 embedded rows;
- the self-contained 258-point, 11-anchor UNSAT certificate for low branch
  `(0,1,0,0,0)`, eliminating 693,636,364,523,341,088 joint assignments, and
  its independent Z3 integer replay;
- the nine-point class-27 radius-two obstruction, its 14-mask maximal
  hypergraph compression, and an independent scalar reconstruction of all
  110 legal gain masks;
- the reproducible class-27 radius-three 741/742/743-point finite-repair
  progression and exact augmented-domain counterexamples that refute the
  first two responses;
- the 753/763-point continuation, incomplete 3,851-model relaxed-cover-space
  audit, six pairwise distance-six exact finite repair phases, and ten exact
  augmented-domain points common to five separated responses;
- the complete 1,073-point class-27 radius-three Hamming-ball obstruction,
  covering all 210,047 legal moves and 14 terminal gain-mask skeletons,
  together with an independent scalar reconstruction returning
  `repair_exists=false`;
- the complete 1,373-point class-27 radius-four Hamming-ball obstruction,
  represented by an exhaustive 19-leaf phase tree over `p=97,109,193`,
  its compact embedded leaf cores, and an independent scalar replay
  returning `repair_exists=false`;
- the partitioned-radius composer and verifier, their regression test, and
  the random-audit witness-export path used for fast constructive
  falsification at radius five;
- the experimental quantified-Z3 source and its explicit `UNKNOWN` timeout
  on the known 14-row no-cover benchmark;
- the 4,637-family ranking and the two scans defining the current frontier;
- a compact index of the 17 currently unresolved ranked finite families;
- the exact period-`2533395664800` conditional-star checkpoint and its
  independently replayed dependency certificates;
- the exact period-`776363187600` conditional-star certificate, all five
  verified block extensions, all conditional dependencies, and the complete
  independent period replay;
- the exact period-`330442912800` conditional-star certificate, all four
  verified block extensions, all conditional dependencies, and the complete
  independent period replay;
- the exact period-`1659810952800` conditional-star certificate, all six
  verified block extensions, all conditional dependencies, and the complete
  independent period replay;
- the complete exact period-`3139207671600` no-cover certificate and replay,
  including the verified 24-anchor extension chain, all 165 promoted
  conditional certificates and replays, and 1,007 independently checked
  explicit pairwise star-subset witnesses;
- the complete exact period-`14440355289360` no-cover certificate and replay,
  including the verified 19-anchor extension chain, 129 promoted projected
  conditional certificates and replays, 136 merged conditional edges, and
  1,105 independently checked explicit pairwise star-subset witnesses;
- the complete exact period-`931635825120` no-cover certificate and replay,
  including the verified 22-anchor extension chain, 17 promoted projected
  conditional certificates and replays, 107 merged conditional edges, and
  923 independently checked explicit pairwise star-subset witnesses;
- the exact `order_pool_1050000` certificate and independent-verification
  dependency closure for the current checkpoint.

## Excluded

- transient standard-output and standard-error logs;
- PID files, bytecode, caches, and temporary files;
- local-search checkpoints and intermediate solver phase assignments;
- the max-32 and max-128 direct-cover CEGIS lesson sets and current phase
  assignments, including the 1,968,920-point continuation checkpoint; their
  parameters and measured repair curves are recorded in `RESEARCH_LOG.md`;
- the local `D=1616615` perfect-power phase assignments, stream caches, and
  lesson sets, including the 830,162-point ordinary and 587,225-point
  triple-coverage multi-digit continuations, the transient 17-adic phase
  assignments, the 2,221-cell/1,100-plane exact radius-two core, the
  16-anchor CEGIS phase assignments, and their newest 100-hole checker
  batches; their exact checker outcomes and structural measurements are
  recorded in `README.md` and `RESEARCH_LOG.md`;
- superseded discovery-only period-`3139207671600` design JSON, CEGIS witnesses,
  phases, and checker logs; their exact aggregate checkpoint and caveats are
  recorded in `RESEARCH_LOG.md`;
- superseded discovery-only period-`14440355289360` projected-design JSON and
  checker logs; the promoted proof artifacts and exact aggregate result are
  included;
- superseded discovery-only period-`931635825120` projected-design JSON and
  checker logs; the promoted proof artifacts and exact aggregate result are
  included;
- large sampled or training point sets;
- obsolete exploratory pools and redundant generated scans.

The excluded local material is approximately 3.7 GB and is not needed to
replay the public exact certificates. Its omission is a publication and
usability decision, not evidence that the omitted experiments succeeded.

No claim in this repository resolves the original infinite problem. The
authoritative status is the beginning of `README.md` and the latest checkpoint
in `RESEARCH_LOG.md`.
