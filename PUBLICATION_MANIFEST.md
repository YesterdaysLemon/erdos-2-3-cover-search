# Public snapshot manifest

This repository is a curated snapshot of the active workspace, published on
2026-07-26.

## Included

- all Python source and regression tests present at publication time;
- the complete research log, including superseded approaches and corrections;
- the stable 1,050,000-order candidate pool used by the current finite sweep;
- the 4,637-family ranking and the two scans defining the current frontier;
- a compact index of the 18 currently unresolved ranked finite families;
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
- the exact `order_pool_1050000` certificate and independent-verification
  dependency closure for the current checkpoint.

## Excluded

- transient standard-output and standard-error logs;
- PID files, bytecode, caches, and temporary files;
- local-search checkpoints and intermediate solver phase assignments;
- superseded discovery-only period-`3139207671600` design JSON, CEGIS witnesses,
  phases, and checker logs; their exact aggregate checkpoint and caveats are
  recorded in `RESEARCH_LOG.md`;
- superseded discovery-only period-`14440355289360` projected-design JSON and
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
