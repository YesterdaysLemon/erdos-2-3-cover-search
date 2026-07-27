# Public snapshot manifest

This repository is a curated snapshot of the active workspace, published on
2026-07-26.

## Included

- all Python source and regression tests present at publication time;
- the complete research log, including superseded approaches and corrections;
- the stable 1,050,000-order candidate pool used by the current finite sweep;
- the 4,637-family ranking and the two scans defining the current frontier;
- a compact index of the 23 currently unresolved ranked finite families;
- the exact period-`2533395664800` conditional-star checkpoint and its
  independently replayed dependency certificates;
- the exact `order_pool_1050000` certificate and independent-verification
  dependency closure for the current checkpoint.

## Excluded

- transient standard-output and standard-error logs;
- PID files, bytecode, caches, and temporary files;
- local-search checkpoints and intermediate solver phase assignments;
- large sampled or training point sets;
- obsolete exploratory pools and redundant generated scans.

The excluded local material is approximately 3.7 GB and is not needed to
replay the public exact certificates. Its omission is a publication and
usability decision, not evidence that the omitted experiments succeeded.

No claim in this repository resolves the original infinite problem. The
authoritative status is the beginning of `README.md` and the latest checkpoint
in `RESEARCH_LOG.md`.
