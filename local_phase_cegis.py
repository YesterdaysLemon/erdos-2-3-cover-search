#!/usr/bin/env python3
"""Compact multivalued CEGIS master for affine-line covers.

Each prime owns one integer target c modulo h.  Unlike the binary SAT master,
this keeps sampled target values in a dense uint32 matrix and repairs uncovered
sample points by min-conflicts moves.  It is heuristic only on the synthesis
side: a construction is emitted solely after exact_uncovered proves UNSAT.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import os
import random
import sys
import time
from pathlib import Path

import exact_greedy
import exact_uncovered

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def build_targets(points, candidates, np):
    started = time.monotonic()
    # Repair is column-oriented: each move inspects all retained points for
    # one candidate.  Fortran order keeps those candidate columns contiguous.
    targets = np.empty(
        (len(points), len(candidates)), dtype=np.uint32, order="F"
    )
    uint64_flags = [
        all(0 <= coordinate < (1 << 64) for coordinate in point)
        for point in points
    ]
    uint64_count = sum(uint64_flags)
    uint64_prefix = (
        uint64_count == len(points)
        or uint64_flags == [True] * uint64_count + [False] * (
            len(points) - uint64_count
        )
    )
    if uint64_count:
        uint64_points = (
            points[:uint64_count]
            if uint64_prefix
            else [point for point, flag in zip(points, uint64_flags) if flag]
        )
        ks = np.fromiter(
            (point[0] for point in uint64_points),
            dtype=np.uint64,
            count=uint64_count,
        )
        ls = np.fromiter(
            (point[1] for point in uint64_points),
            dtype=np.uint64,
            count=uint64_count,
        )
        uint64_selector = (
            slice(0, uint64_count)
            if uint64_prefix
            else np.flatnonzero(uint64_flags)
        )
    big_points = (
        points[uint64_count:]
        if uint64_prefix
        else [point for point, flag in zip(points, uint64_flags) if not flag]
    )
    big_selector = (
        slice(uint64_count, len(points))
        if uint64_prefix
        else np.flatnonzero(np.logical_not(uint64_flags))
    )
    for column, (h, _p, a, b, _ord2, _ord3) in enumerate(candidates):
        if uint64_count:
            targets[uint64_selector, column] = (
                (a * (ks % h) + b * (ls % h)) % h
            ).astype(np.uint32)
        if big_points:
            targets[big_selector, column] = np.fromiter(
                ((a * (k % h) + b * (l % h)) % h for k, l in big_points),
                dtype=np.uint32,
                count=len(big_points),
            )
    return targets, time.monotonic() - started


def streaming_coordinate_arrays(points, candidates, np):
    """Return coordinate arrays safe for on-demand target evaluation."""
    max_h = max((int(candidate[0]) for candidate in candidates), default=1)
    uint64_limit = (1 << 64) - 1
    safe_products = 2 * max(0, max_h - 1) ** 2 <= uint64_limit
    safe_coordinates = all(
        0 <= int(coordinate) <= uint64_limit
        for point in points
        for coordinate in point
    )
    dtype = np.uint64 if safe_products and safe_coordinates else object
    ks = np.asarray([int(point[0]) for point in points], dtype=dtype)
    ls = np.asarray([int(point[1]) for point in points], dtype=dtype)
    return ks, ls


def streaming_target_values(ks, ls, candidate, np):
    """Evaluate one candidate column without retaining the full matrix."""
    h, _p, a, b, _ord2, _ord3 = candidate
    h = int(h)
    a = int(a) % h
    b = int(b) % h
    return (
        (a * (ks % h) + b * (ls % h)) % h
    ).astype(np.uint32, copy=False)


def streaming_coordinate_state(
    points,
    candidates,
    np,
    force_components=False,
):
    """Prepare either raw uint64 coordinates or compact CRT residues.

    Exact checker witnesses can be far larger than 64 bits.  Retaining them
    as NumPy object arrays would make every candidate evaluation execute
    Python big-integer arithmetic point by point.  Instead, store residues
    modulo the prime-power factors that actually occur in candidate moduli.
    """
    ks, ls = streaming_coordinate_arrays(points, candidates, np)
    if ks.dtype != object and not force_components:
        return {
            "mode": "raw",
            "ks": ks,
            "ls": ls,
        }

    component_moduli = sorted(
        {
            prime**exponent
            for candidate in candidates
            for prime, exponent in exact_uncovered.factor(
                int(candidate[0])
            ).items()
        }
    )
    component_index = {
        modulus: index
        for index, modulus in enumerate(component_moduli)
    }
    max_residue = max(component_moduli, default=1) - 1
    residue_dtype = np.min_scalar_type(max_residue)
    k_components = np.empty(
        (len(points), len(component_moduli)),
        dtype=residue_dtype,
        order="F",
    )
    l_components = np.empty_like(k_components, order="F")
    for column, modulus in enumerate(component_moduli):
        k_components[:, column] = np.fromiter(
            (int(point[0]) % modulus for point in points),
            dtype=residue_dtype,
            count=len(points),
        )
        l_components[:, column] = np.fromiter(
            (int(point[1]) % modulus for point in points),
            dtype=residue_dtype,
            count=len(points),
        )

    spec_by_h = {}
    for candidate in candidates:
        h = int(candidate[0])
        if h in spec_by_h:
            continue
        spec = []
        for prime, exponent in exact_uncovered.factor(h).items():
            modulus = prime**exponent
            cofactor = h // modulus
            coefficient = (
                cofactor * pow(cofactor, -1, modulus)
            ) % h
            spec.append(
                (component_index[modulus], modulus, coefficient)
            )
        spec_by_h[h] = tuple(spec)
    return {
        "mode": "components",
        "k_components": k_components,
        "l_components": l_components,
        "spec_by_h": spec_by_h,
    }


def append_streaming_coordinate_state(state, addition, np):
    """Append compatible state rows without recomputing old residues.

    A raw uint64 checkpoint can become unsuitable when an exact checker first
    returns a larger coordinate.  That one-time mode transition needs a full
    rebuild.  Once both states use the same representation, their row arrays
    can be concatenated directly.
    """
    if state["mode"] != addition["mode"]:
        return None
    if state["mode"] == "raw":
        return {
            "mode": "raw",
            "ks": np.concatenate((state["ks"], addition["ks"])),
            "ls": np.concatenate((state["ls"], addition["ls"])),
        }
    if state["spec_by_h"] != addition["spec_by_h"]:
        raise ValueError("incompatible streaming component specifications")
    return {
        "mode": "components",
        "k_components": np.concatenate(
            (
                state["k_components"],
                addition["k_components"],
            ),
            axis=0,
        ),
        "l_components": np.concatenate(
            (
                state["l_components"],
                addition["l_components"],
            ),
            axis=0,
        ),
        "spec_by_h": state["spec_by_h"],
    }


def streaming_candidate_digest(candidates) -> str:
    """Return a stable fingerprint of the ordered streaming columns."""
    payload = [
        [int(value) for value in candidate]
        for candidate in candidates
    ]
    return hashlib.sha256(
        json.dumps(payload, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def streaming_assignment_digest(assignment) -> str:
    """Return a platform-independent fingerprint of a phase assignment."""
    digest = hashlib.sha256(b"erdos203-stream-assignment-v1\0")
    for value in assignment:
        digest.update(int(value).to_bytes(8, "little", signed=False))
    return digest.hexdigest()


def streaming_points_digest(points) -> str:
    """Fingerprint an ordered point prefix without a giant JSON copy."""
    digest = hashlib.sha256(b"erdos203-stream-points-v1\0")
    for point in points:
        for raw_value in point:
            value = int(raw_value)
            magnitude = abs(value)
            encoded = magnitude.to_bytes(
                max(1, (magnitude.bit_length() + 7) // 8),
                "little",
                signed=False,
            )
            digest.update(b"\x01" if value < 0 else b"\x00")
            digest.update(len(encoded).to_bytes(8, "little"))
            digest.update(encoded)
    return digest.hexdigest()


def save_streaming_cache(
    path,
    state,
    cover,
    candidates,
    assignment,
    points,
    np,
) -> None:
    """Atomically save reusable streaming coordinates and cover counts."""
    metadata = {
        "version": 1,
        "candidate_digest": streaming_candidate_digest(candidates),
        "assignment_digest": streaming_assignment_digest(assignment),
        "point_count": len(points),
        "point_digest": streaming_points_digest(points),
        "mode": state["mode"],
    }
    arrays = {
        "metadata": np.frombuffer(
            json.dumps(metadata, separators=(",", ":")).encode("utf-8"),
            dtype=np.uint8,
        ),
        "cover": np.asarray(cover),
    }
    if state["mode"] == "raw":
        arrays["ks"] = state["ks"]
        arrays["ls"] = state["ls"]
    else:
        arrays["k_components"] = state["k_components"]
        arrays["l_components"] = state["l_components"]
        specs = {
            str(h): [
                [int(value) for value in component]
                for component in spec
            ]
            for h, spec in state["spec_by_h"].items()
        }
        arrays["spec_by_h"] = np.frombuffer(
            json.dumps(specs, separators=(",", ":")).encode("utf-8"),
            dtype=np.uint8,
        )
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as stream:
        np.savez(stream, **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def load_streaming_cache(
    path,
    candidates,
    assignment,
    points,
    np,
):
    """Load a validated point-prefix cache.

    Coordinates remain reusable when phases change, but cover counts do not.
    The returned prefix length lets the caller append any newer lessons.
    """
    if not path.exists():
        return None, None, 0, "missing"
    try:
        with np.load(path, allow_pickle=False) as cache:
            metadata = json.loads(
                cache["metadata"].tobytes().decode("utf-8")
            )
            if int(metadata.get("version", -1)) != 1:
                return None, None, 0, "version-mismatch"
            if metadata.get(
                "candidate_digest"
            ) != streaming_candidate_digest(candidates):
                return None, None, 0, "candidate-mismatch"
            point_count = int(metadata["point_count"])
            if not 0 <= point_count <= len(points):
                return None, None, 0, "point-count-mismatch"
            if metadata.get("point_digest") != streaming_points_digest(
                points[:point_count]
            ):
                return None, None, 0, "point-prefix-mismatch"
            mode = metadata.get("mode")
            if mode == "raw":
                state = {
                    "mode": "raw",
                    "ks": cache["ks"].copy(),
                    "ls": cache["ls"].copy(),
                }
                row_count = len(state["ks"])
                if len(state["ls"]) != row_count:
                    return None, None, 0, "coordinate-shape-mismatch"
            elif mode == "components":
                raw_specs = json.loads(
                    cache["spec_by_h"].tobytes().decode("utf-8")
                )
                state = {
                    "mode": "components",
                    "k_components": cache["k_components"].copy(),
                    "l_components": cache["l_components"].copy(),
                    "spec_by_h": {
                        int(h): tuple(
                            tuple(int(value) for value in component)
                            for component in spec
                        )
                        for h, spec in raw_specs.items()
                    },
                }
                row_count = state["k_components"].shape[0]
                if (
                    state["l_components"].shape
                    != state["k_components"].shape
                ):
                    return None, None, 0, "coordinate-shape-mismatch"
            else:
                return None, None, 0, "mode-mismatch"
            if row_count != point_count:
                return None, None, 0, "coordinate-count-mismatch"
            cover = cache["cover"].copy()
            if len(cover) != point_count:
                return None, None, 0, "cover-count-mismatch"
            if metadata.get(
                "assignment_digest"
            ) != streaming_assignment_digest(assignment):
                cover = None
                status = "coordinates-only"
            else:
                status = "exact"
            return state, cover, point_count, status
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        json.JSONDecodeError,
    ):
        return None, None, 0, "invalid"


def streaming_state_target_values(
    state,
    candidate,
    np,
    indices=None,
):
    """Evaluate one target column from a compact streaming state."""
    if state["mode"] == "raw":
        ks = state["ks"] if indices is None else state["ks"][indices]
        ls = state["ls"] if indices is None else state["ls"][indices]
        return streaming_target_values(ks, ls, candidate, np)

    h, _p, a, b, _ord2, _ord3 = candidate
    h = int(h)
    a = int(a) % h
    b = int(b) % h
    k_components = state["k_components"]
    l_components = state["l_components"]
    length = (
        len(k_components)
        if indices is None
        else len(indices)
    )
    result = np.zeros(length, dtype=np.uint64)
    for column, modulus, coefficient in state["spec_by_h"][h]:
        kres = (
            k_components[:, column]
            if indices is None
            else k_components[indices, column]
        )
        lres = (
            l_components[:, column]
            if indices is None
            else l_components[indices, column]
        )
        local_target = (
            a * kres.astype(np.uint64, copy=False)
            + b * lres.astype(np.uint64, copy=False)
        ) % modulus
        result += local_target * coefficient
    return (result % h).astype(np.uint32, copy=False)


def streaming_cover_counts(
    coordinate_state,
    candidates,
    assignment,
    np,
    indices=None,
):
    """Count selected-fibre coverage on all or selected retained points."""
    if coordinate_state["mode"] == "raw":
        total = len(coordinate_state["ks"])
    else:
        total = len(coordinate_state["k_components"])
    length = total if indices is None else len(indices)
    cover = np.zeros(length, dtype=np.int32)
    for candidate_index, candidate in enumerate(candidates):
        values = streaming_state_target_values(
            coordinate_state,
            candidate,
            np,
            indices,
        )
        cover += (
            values == assignment[candidate_index]
        ).astype(np.int32)
    return cover


def repair_streaming(
    points,
    candidates,
    assignment,
    mutable,
    rng,
    np,
    max_steps,
    sample_size,
    valid_moduli=None,
    valid_residues=None,
    candidate_moduli=None,
    required_coverage=1,
    coordinate_descent=False,
    coordinate_state=None,
    initial_cover=None,
    return_state=False,
):
    """Point-wise min-conflicts with candidate columns computed on demand.

    The dense learner retains one uint32 target for every point/candidate
    pair.  This variant keeps only the coordinates and current cover counts.
    Candidate columns are evaluated over the uncovered and fragile points
    while scoring, then once over all points for the selected move.
    """
    if coordinate_state is None:
        coordinate_state = streaming_coordinate_state(
            points,
            candidates,
            np,
        )
    if initial_cover is None:
        cover = streaming_cover_counts(
            coordinate_state,
            candidates,
            assignment,
            np,
        )
    else:
        if len(initial_cover) != len(points):
            raise ValueError("initial streaming cover has the wrong length")
        cover = initial_cover.copy()

    def outcome(result):
        if return_state:
            return (*result, coordinate_state, cover)
        return result

    start_uncovered = int(np.count_nonzero(cover < required_coverage))
    best_uncovered = start_uncovered
    best_deficit = int(
        np.maximum(required_coverage - cover, 0).sum()
    )
    best_assignment = assignment.copy()
    stagnant = 0
    for step in range(max_steps + 1):
        uncovered = np.flatnonzero(cover < required_coverage)
        if not len(uncovered):
            return outcome((True, step, start_uncovered, 0))
        if step == max_steps:
            break

        point = int(uncovered[rng.randrange(len(uncovered))])
        selected = rng.sample(mutable, min(sample_size, len(mutable)))
        fragile = np.flatnonzero(cover <= required_coverage)
        moves = []
        best_score = None
        for candidate_index in selected:
            candidate = candidates[candidate_index]
            h, _p, a, b, _ord2, _ord3 = candidate
            h = int(h)
            old_value = int(assignment[candidate_index])
            uncovered_values = streaming_state_target_values(
                coordinate_state,
                candidate,
                np,
                uncovered,
            )
            if coordinate_descent:
                observed_targets, observed_counts = np.unique(
                    uncovered_values,
                    return_counts=True,
                )
                usable = observed_targets != old_value
                if valid_moduli is not None:
                    usable &= (
                        observed_targets
                        % int(valid_moduli[candidate_index])
                        == int(valid_residues[candidate_index])
                    )
                if not np.any(usable):
                    continue
                observed_targets = observed_targets[usable]
                observed_counts = observed_counts[usable]
                top_count = int(observed_counts.max())
                top_targets = observed_targets[
                    observed_counts == top_count
                ]
                new_value = int(
                    top_targets[rng.randrange(len(top_targets))]
                )
                gains = top_count
            else:
                new_value = (
                    int(
                        streaming_state_target_values(
                            coordinate_state,
                            candidate,
                            np,
                            np.asarray([point], dtype=np.int64),
                        )[0]
                    )
                )
                if new_value == old_value:
                    continue
                if (
                    valid_moduli is not None
                    and new_value % int(valid_moduli[candidate_index])
                    != int(valid_residues[candidate_index])
                ):
                    continue
                gains = int(
                    np.count_nonzero(uncovered_values == new_value)
                )
            fragile_values = streaming_state_target_values(
                coordinate_state,
                candidate,
                np,
                fragile,
            )
            losses = int(np.count_nonzero(fragile_values == old_value))
            score = gains - losses
            if best_score is None or score > best_score:
                best_score = score
                moves = [(candidate_index, new_value)]
            elif score == best_score:
                moves.append((candidate_index, new_value))
        if not moves:
            continue

        candidate_index, new_value = moves[rng.randrange(len(moves))]
        old_value = int(assignment[candidate_index])
        values = streaming_state_target_values(
            coordinate_state,
            candidates[candidate_index],
            np,
        )
        cover -= (values == old_value).astype(np.int32)
        cover += (values == new_value).astype(np.int32)
        assignment[candidate_index] = new_value

        current = int(np.count_nonzero(cover < required_coverage))
        current_deficit = int(
            np.maximum(required_coverage - cover, 0).sum()
        )
        if current_deficit < best_deficit:
            best_deficit = current_deficit
            best_uncovered = current
            best_assignment[:] = assignment
            stagnant = 0
        else:
            stagnant += 1
        if (step + 1) % 250 == 0:
            print(
                f"stream_step={step + 1} "
                f"mode={'coordinate' if coordinate_descent else 'point'} "
                f"misses={current} best={best_uncovered} "
                f"deficit={current_deficit} best_deficit={best_deficit}",
                flush=True,
            )
        if stagnant >= 2500:
            # Preserve the dense learner's cycle-breaking neighborhood while
            # still materializing only one candidate column at a time.
            for candidate_index in rng.sample(
                mutable, min(25, len(mutable))
            ):
                old_value = int(assignment[candidate_index])
                if valid_moduli is None:
                    new_value = rng.randrange(
                        int(candidate_moduli[candidate_index])
                    )
                else:
                    modulus = int(valid_moduli[candidate_index])
                    residue = int(valid_residues[candidate_index])
                    choices = (
                        int(candidate_moduli[candidate_index])
                        - 1
                        - residue
                    ) // modulus + 1
                    new_value = residue + modulus * rng.randrange(choices)
                if new_value == old_value:
                    continue
                values = streaming_state_target_values(
                    coordinate_state,
                    candidates[candidate_index],
                    np,
                )
                cover -= (values == old_value).astype(np.int32)
                cover += (values == new_value).astype(np.int32)
                assignment[candidate_index] = new_value
            stagnant = 0
    assignment[:] = best_assignment
    return outcome(
        (False, max_steps, start_uncovered, best_uncovered)
    )


def merge_congruences(residue_a, modulus_a, residue_b, modulus_b):
    """Return the intersection of two congruences, or None if it is empty."""
    common = math.gcd(modulus_a, modulus_b)
    if (residue_b - residue_a) % common:
        return None
    left = modulus_a // common
    right = modulus_b // common
    if right == 1:
        multiplier = 0
    else:
        multiplier = (
            ((residue_b - residue_a) // common) * pow(left, -1, right)
        ) % right
    modulus = math.lcm(modulus_a, modulus_b)
    return (residue_a + modulus_a * multiplier) % modulus, modulus


def expand_component_digit_tiles(
    points, candidate_moduli, component_digits
):
    """Expand exact holes across selected independent component digits.

    For a prime q whose largest exponent in the common candidate period is
    q^e, adding a multiple of (period / q^e) * q^j preserves every other CRT
    component and every lower q-adic digit while cycling digit j.  Applying
    this independently to both coordinates and to several distinct
    prime-adic digits returns the complete Cartesian digit tile around each
    point.  Multiple distinct digits of the same prime may be varied
    together.
    """
    component_digits = tuple(component_digits)
    if len(set(component_digits)) != len(component_digits):
        raise ValueError("tile component digits must be distinct")
    components = tuple(sorted({prime for prime, _digit in component_digits}))
    period = math.lcm(1, *(int(modulus) for modulus in candidate_moduli))
    factorization = exact_uncovered.factor(period)
    missing = [prime for prime in components if prime not in factorization]
    if missing:
        raise ValueError(
            f"tile components do not divide the candidate period: {missing}"
        )

    coordinate_offsets = [0]
    for prime, digit in component_digits:
        exponent = factorization[prime]
        if not 0 <= digit < exponent:
            raise ValueError(
                f"digit {digit} is outside component {prime}^{exponent}"
            )
        prime_power = prime**exponent
        step = (period // prime_power) * prime**digit
        coordinate_offsets = [
            (offset + digit * step) % period
            for offset in coordinate_offsets
            for digit in range(prime)
        ]

    expanded = []
    seen = set()
    for raw_k, raw_l in points:
        k = int(raw_k)
        l = int(raw_l)
        for k_offset in coordinate_offsets:
            for l_offset in coordinate_offsets:
                point = ((k + k_offset) % period, (l + l_offset) % period)
                if point not in seen:
                    seen.add(point)
                    expanded.append(point)
    return expanded


def expand_top_component_tiles(points, candidate_moduli, components):
    """Expand holes over the top digit of each declared component."""
    period = math.lcm(1, *(int(modulus) for modulus in candidate_moduli))
    factorization = exact_uncovered.factor(period)
    missing = [prime for prime in components if prime not in factorization]
    if missing:
        raise ValueError(
            f"top-tile components do not divide the candidate period: {missing}"
        )
    return expand_component_digit_tiles(
        points,
        candidate_moduli,
        [
            (prime, factorization[prime] - 1)
            for prime in components
        ],
    )


def expand_component_digit_tile_groups(
    points,
    candidate_moduli,
    component_digit_groups,
):
    """Return the union of several independent component-digit tiles."""
    expanded = []
    seen = set()
    for component_digits in component_digit_groups:
        for point in expand_component_digit_tiles(
            points,
            candidate_moduli,
            component_digits,
        ):
            if point not in seen:
                seen.add(point)
                expanded.append(point)
    return expanded


def repair(
    targets,
    assignment,
    mutable,
    rng,
    np,
    max_steps,
    sample_size,
    valid_moduli=None,
    valid_residues=None,
    candidate_moduli=None,
    required_coverage=1,
):
    cover = np.count_nonzero(targets == assignment, axis=1).astype(np.int32)
    start_uncovered = int(np.count_nonzero(cover < required_coverage))
    best_uncovered = start_uncovered
    best_deficit = int(
        np.maximum(required_coverage - cover, 0).sum()
    )
    best_assignment = assignment.copy()
    stagnant = 0
    for step in range(max_steps + 1):
        uncovered = np.flatnonzero(cover < required_coverage)
        if not len(uncovered):
            return True, step, start_uncovered, 0
        if step == max_steps:
            break
        point = int(uncovered[rng.randrange(len(uncovered))])
        size = min(sample_size, len(mutable))
        selected = rng.sample(mutable, size)
        selected_array = np.asarray(selected, dtype=np.int32)
        new_values = targets[point, selected_array]
        changed = new_values != assignment[selected_array]
        if valid_moduli is not None:
            changed &= (
                new_values % valid_moduli[selected_array]
                == valid_residues[selected_array]
            )
        if not np.any(changed):
            continue
        selected_array = selected_array[changed]
        new_values = new_values[changed]
        slab = targets[:, selected_array]
        old_matches = slab == assignment[selected_array]
        new_matches = slab == new_values
        losses = np.count_nonzero(
            old_matches & (cover[:, None] <= required_coverage),
            axis=0,
        )
        gains = np.count_nonzero(
            new_matches & (cover[:, None] < required_coverage),
            axis=0,
        )
        scores = gains - losses
        best_score = int(scores.max())
        choices = np.flatnonzero(scores == best_score)
        choice = int(choices[rng.randrange(len(choices))])
        candidate = int(selected_array[choice])
        old_value = int(assignment[candidate])
        new_value = int(new_values[choice])
        old_mask = targets[:, candidate] == old_value
        new_mask = targets[:, candidate] == new_value
        cover -= old_mask.astype(np.int32)
        cover += new_mask.astype(np.int32)
        assignment[candidate] = new_value

        current = int(np.count_nonzero(cover < required_coverage))
        current_deficit = int(
            np.maximum(required_coverage - cover, 0).sum()
        )
        if current_deficit < best_deficit:
            best_deficit = current_deficit
            best_uncovered = current
            best_assignment[:] = assignment
            stagnant = 0
        else:
            stagnant += 1
        if stagnant >= 2500:
            # A small exact-domain-neutral perturbation breaks repair cycles.
            for candidate in rng.sample(mutable, min(25, len(mutable))):
                old_value = int(assignment[candidate])
                if valid_moduli is None:
                    new_value = rng.randrange(int(candidate_moduli[candidate]))
                else:
                    modulus = int(valid_moduli[candidate])
                    residue = int(valid_residues[candidate])
                    choices = (int(candidate_moduli[candidate]) - 1 - residue) // modulus + 1
                    new_value = residue + modulus * rng.randrange(choices)
                if new_value == old_value:
                    continue
                old_mask = targets[:, candidate] == old_value
                new_mask = targets[:, candidate] == new_value
                cover -= old_mask.astype(np.int32)
                cover += new_mask.astype(np.int32)
                assignment[candidate] = new_value
            stagnant = 0
    assignment[:] = best_assignment
    return False, max_steps, start_uncovered, best_uncovered


def repair_coordinate(
    targets,
    assignment,
    mutable,
    rng,
    np,
    max_steps,
    sample_size,
    valid_moduli,
    valid_residues,
    candidate_moduli,
    required_coverage=1,
):
    """Coordinate descent using the best valid target on all missed points."""
    cover = np.count_nonzero(targets == assignment, axis=1).astype(np.int32)
    start_uncovered = int(np.count_nonzero(cover < required_coverage))
    best_uncovered = start_uncovered
    best_deficit = int(
        np.maximum(required_coverage - cover, 0).sum()
    )
    best_assignment = assignment.copy()
    for step in range(max_steps + 1):
        uncovered_mask = cover < required_coverage
        if not np.any(uncovered_mask):
            return True, step, start_uncovered, 0
        if step == max_steps:
            break
        selected = rng.sample(mutable, min(sample_size, len(mutable)))
        best_move = None
        best_score = None
        for candidate in selected:
            values = targets[:, candidate]
            h = int(candidate_moduli[candidate])
            modulus = int(valid_moduli[candidate])
            residue = int(valid_residues[candidate])
            observed = values[uncovered_mask]
            valid_observed = observed[
                observed % modulus == residue
            ]
            if not len(valid_observed):
                continue
            # Large subgroup orders made the old minlength=h bincount allocate
            # mostly empty arrays for every coordinate.  Count only targets
            # that actually occur on the current missed points unless the
            # domain is already compact.
            if h <= max(256, 4 * len(valid_observed)):
                counts = np.bincount(valid_observed, minlength=h)
                valid_counts = counts[residue::modulus]
                target_index = int(np.argmax(valid_counts))
                new_value = residue + modulus * target_index
                gains = int(valid_counts[target_index])
            else:
                observed_targets, observed_counts = np.unique(
                    valid_observed, return_counts=True
                )
                target_index = int(np.argmax(observed_counts))
                new_value = int(observed_targets[target_index])
                gains = int(observed_counts[target_index])
            old_value = int(assignment[candidate])
            if new_value == old_value:
                continue
            losses = int(
                np.count_nonzero(
                    (cover <= required_coverage) & (values == old_value)
                )
            )
            score = gains - losses
            if best_score is None or score > best_score:
                best_score = score
                best_move = (candidate, old_value, new_value, values)
        if best_move is None:
            continue
        candidate, old_value, new_value, values = best_move
        cover -= (values == old_value).astype(np.int32)
        cover += (values == new_value).astype(np.int32)
        assignment[candidate] = new_value
        current = int(np.count_nonzero(cover < required_coverage))
        current_deficit = int(
            np.maximum(required_coverage - cover, 0).sum()
        )
        if current_deficit < best_deficit:
            best_deficit = current_deficit
            best_uncovered = current
            best_assignment[:] = assignment
    assignment[:] = best_assignment
    return False, max_steps, start_uncovered, best_uncovered


def repair_component(
    targets,
    assignment,
    mutable,
    rng,
    np,
    max_steps,
    sample_size,
    valid_moduli,
    valid_residues,
    candidate_moduli,
    required_coverage=1,
):
    """Greedily retarget one prime-power CRT digit, then fall back.

    A component move preserves every coprime digit of a phase.  This gives a
    much smoother neighborhood than replacing the entire phase by one point's
    target, especially when a candidate modulus has several components.
    """
    cover = np.count_nonzero(targets == assignment, axis=1).astype(np.int32)
    start_uncovered = int(np.count_nonzero(cover < required_coverage))
    if not start_uncovered:
        return True, 0, 0, 0
    components = {
        candidate: tuple(
            prime**exponent
            for prime, exponent in exact_uncovered.factor(
                int(candidate_moduli[candidate])
            ).items()
        )
        for candidate in mutable
    }
    moves = 0
    while moves < max_steps:
        selected = rng.sample(
            mutable, min(sample_size, len(mutable))
        )
        rng.shuffle(selected)
        sweep_moves = 0
        for candidate in selected:
            h = int(candidate_moduli[candidate])
            values = targets[:, candidate]
            powers = list(components[candidate])
            rng.shuffle(powers)
            for component in powers:
                if moves >= max_steps:
                    break
                other = h // component
                old_value = int(assignment[candidate])
                old_hits = values == old_value
                base_cover = cover - old_hits
                needed = base_cover < required_coverage
                other_residue = old_value % other if other > 1 else 0
                compatible = needed
                if other > 1:
                    compatible = compatible & (
                        values % other == other_residue
                    )
                component_constraint = math.gcd(
                    component, int(valid_moduli[candidate])
                )
                legal_residue = (
                    int(valid_residues[candidate])
                    % component_constraint
                )
                component_values = values % component
                if component_constraint > 1:
                    compatible = compatible & (
                        component_values % component_constraint
                        == legal_residue
                    )
                observed = component_values[compatible]
                if not len(observed):
                    continue
                choices, counts = np.unique(
                    observed, return_counts=True
                )
                best_score = int(counts.max())
                old_component = old_value % component
                old_score = int(
                    np.count_nonzero(
                        compatible
                        & (component_values == old_component)
                    )
                )
                if best_score <= old_score:
                    continue
                best_components = choices[counts == best_score]
                new_component = int(
                    best_components[
                        rng.randrange(len(best_components))
                    ]
                )
                if other == 1:
                    new_value = new_component
                elif component == 1:
                    new_value = other_residue
                else:
                    multiplier = (
                        (new_component - other_residue)
                        * pow(other, -1, component)
                    ) % component
                    new_value = (
                        other_residue + other * multiplier
                    ) % h
                if (
                    new_value % int(valid_moduli[candidate])
                    != int(valid_residues[candidate])
                ):
                    raise AssertionError(
                        "component move violates target restriction"
                    )
                new_hits = values == new_value
                cover = base_cover + new_hits
                assignment[candidate] = new_value
                moves += 1
                sweep_moves += 1
                if not np.any(cover < required_coverage):
                    return True, moves, start_uncovered, 0
        if not sweep_moves:
            break

    # Whole-phase min-conflicts remains useful for escaping a component-local
    # optimum.  Component moves were strictly deficit-improving, so the
    # current assignment is a safe warm start for this fallback.
    remaining = max(0, max_steps - moves)
    solved, fallback_steps, _before, after = repair(
        targets,
        assignment,
        mutable,
        rng,
        np,
        remaining,
        sample_size,
        valid_moduli,
        valid_residues,
        candidate_moduli,
        required_coverage,
    )
    return solved, moves + fallback_steps, start_uncovered, after


def random_misses(
    candidates,
    assignment,
    coordinate_bits,
    rng,
    checks,
    batch,
    algebraic_primes=(),
    sophie_germain=False,
):
    misses = []
    tested = []
    for _ in range(checks):
        # A huge exact LCM is expensive to construct for large pools.  Fresh
        # fixed-width integers are ample for heuristic discovery; the exact
        # checker remains the only acceptance test.
        k = rng.getrandbits(coordinate_bits)
        l = rng.getrandbits(coordinate_bits)
        if any(k % prime == 0 and l % prime == 0 for prime in algebraic_primes):
            continue
        if sophie_germain and k % 4 == 2 and l % 4 == 0:
            continue
        tested.append((k, l))
        if all(
            (a * (k % h) + b * (l % h) - int(assignment[index])) % h
            for index, (h, _p, a, b, _ord2, _ord3) in enumerate(candidates)
        ):
            misses.append((k, l))
            if len(misses) >= batch:
                break
    return misses, tested


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_pool", type=Path)
    parser.add_argument("--derived-pool", action="store_true")
    parser.add_argument(
        "--derived-targets",
        action="store_true",
        help=(
            "read target_residue and target_modulus from each derived-pool "
            "row instead of recomputing a perfect-power condition"
        ),
    )
    parser.add_argument(
        "--period",
        type=int,
        default=0,
        help="if nonzero, retain only candidates whose h divides this period",
    )
    parser.add_argument(
        "--min-parent-common",
        type=int,
        default=0,
        help=(
            "for a conditioned derived pool, retain only rows whose recorded "
            "parent_common is at least this value"
        ),
    )
    parser.add_argument("--candidate-count", type=int)
    parser.add_argument(
        "--fixed-prefix",
        type=int,
        default=0,
        help=(
            "keep this many leading compatible candidate phases fixed at "
            "their saved values"
        ),
    )
    parser.add_argument(
        "--fixed-primes",
        default="",
        help="comma-separated primes to keep fixed at their saved phases",
    )
    parser.add_argument(
        "--drop-fixed-covered-points",
        action="store_true",
        help=(
            "omit training points already covered by a fixed candidate; "
            "useful for residual high-band repair"
        ),
    )
    parser.add_argument("--normalize-primes", default="")
    parser.add_argument(
        "--diversity-primes",
        default="",
        help=(
            "comma-separated small fibres whose exact target fingerprint is "
            "blocked between checker witnesses, without fixing their phases"
        ),
    )
    parser.add_argument(
        "--diversity-quota",
        type=int,
        default=1,
        help=(
            "retain this many exact CRT representatives per diversity "
            "fingerprint before blocking the whole fingerprint"
        ),
    )
    parser.add_argument(
        "--diversity-target-cap",
        type=int,
        default=0,
        help=(
            "after this many exact witnesses share one target value on a "
            "declared diversity prime, forbid that row-target value for the "
            "rest of the current checker batch; zero disables the cap"
        ),
    )
    parser.add_argument(
        "--diversity-coordinate-moduli",
        default="",
        help=(
            "comma-separated prime-power moduli whose (k,l) residue cells "
            "are diversified by the exact checker"
        ),
    )
    parser.add_argument(
        "--diversity-coordinate-schedule",
        default="",
        help=(
            "semicolon-separated coordinate-modulus stages, for example "
            "'7;7,11;7,13;none'; reuses the retained target matrix"
        ),
    )
    parser.add_argument(
        "--diversity-stage-rounds",
        type=int,
        default=0,
        help="number of exact rounds per scheduled diversity stage",
    )
    parser.add_argument(
        "--diversity-schedule-cycle",
        action="store_true",
        help="repeat the coordinate-diversity schedule after its last stage",
    )
    parser.add_argument(
        "--freeze-target-coordinate-moduli",
        default="",
        help=(
            "comma-separated pairwise-compatible coordinate moduli; preserve "
            "each saved target modulo gcd(h,lcm(moduli)) while repairing the "
            "remaining CRT components"
        ),
    )
    parser.add_argument(
        "--checker-coordinate-cell",
        default="",
        help=(
            "comma-separated modulus:k:l restrictions for a focused exact "
            "cell search, for example 16:0:1,27:8:0; random checks must be 0"
        ),
    )
    parser.add_argument("--points-file", type=Path, required=True)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--initial-phase-file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--checker-checkpoint-file",
        type=Path,
        help=(
            "overwrite a JSON record of the latest random and exact checker "
            "witnesses; with --rounds 1 the phase file is the tested phase"
        ),
    )
    parser.add_argument("--rounds", type=int, default=1000)
    parser.add_argument(
        "--seed",
        type=int,
        help="override the deterministic repair/random-audit seed",
    )
    parser.add_argument(
        "--randomize-mutable-initial",
        action="store_true",
        help=(
            "replace saved mutable phases by seed-deterministic random valid "
            "targets before the first repair; fixed phases are unchanged"
        ),
    )
    parser.add_argument("--max-steps", type=int, default=50000)
    parser.add_argument("--sample-size", type=int, default=256)
    parser.add_argument(
        "--repair-mode",
        choices=("point", "coordinate", "component"),
        default="point",
        help=(
            "point-wise min-conflicts, global best-target coordinate descent, "
            "or prime-power CRT-component descent with point fallback"
        ),
    )
    parser.add_argument(
        "--stream-targets",
        action="store_true",
        help=(
            "compute candidate target columns on demand instead of retaining "
            "the full points-by-candidates uint32 matrix; supports "
            "--repair-mode point and coordinate"
        ),
    )
    parser.add_argument(
        "--stream-cache-file",
        type=Path,
        help=(
            "validated NumPy checkpoint for streaming coordinate residues "
            "and cover counts; exact point prefixes remain reusable after "
            "lessons are appended"
        ),
    )
    parser.add_argument(
        "--required-coverage",
        type=int,
        default=1,
        help="minimum number of selected fibres covering every retained point",
    )
    parser.add_argument(
        "--repair-only",
        action="store_true",
        help=(
            "write a repaired finite-sample phase and exit before random or "
            "exact checking; this is heuristic preprocessing, never a cover "
            "certificate"
        ),
    )
    parser.add_argument("--random-checks", type=int, default=100000)
    parser.add_argument("--random-coordinate-bits", type=int, default=128)
    parser.add_argument("--random-batch", type=int, default=500)
    parser.add_argument(
        "--retain-random-tested",
        action="store_true",
        help="retain covered as well as missed random tests to prevent phase overfitting",
    )
    parser.add_argument("--exact-batch", type=int, default=200)
    parser.add_argument(
        "--expand-exact-top-components",
        default="",
        help=(
            "comma-separated prime components; replace each exact witness "
            "lesson by its complete Cartesian tile over the top digit of "
            "each selected component"
        ),
    )
    parser.add_argument(
        "--expand-exact-component-digits",
        default="",
        help=(
            "comma-separated prime:zero_based_digit pairs; replace each "
            "exact witness lesson by its Cartesian tile over those CRT "
            "digits, for example 2:0,3:0,5:0"
        ),
    )
    parser.add_argument(
        "--expand-exact-component-digit-groups",
        default="",
        help=(
            "semicolon-separated unions of component-digit tiles, with "
            "comma-separated prime:digit pairs inside each group; for "
            "example 2:0,3:0,5:0;7:0;11:0"
        ),
    )
    parser.add_argument(
        "--exact-when-random-misses-at-most",
        type=int,
        default=0,
        help=(
            "invoke the exact checker when the random audit has at most this "
            "many misses; exact witnesses and random misses are both retained"
        ),
    )
    parser.add_argument("--max-component", type=int, default=300000)
    parser.add_argument("--checker-solver", default="cadical195")
    parser.add_argument(
        "--retain-target-matrix-during-exact",
        action="store_true",
        help=(
            "keep the synthesis matrix while the exact checker runs and "
            "append only new lesson columns afterward; trades memory for "
            "faster repeated exact rounds"
        ),
    )
    parser.add_argument("--power", type=int, default=1)
    parser.add_argument(
        "--algebraic-primes",
        help=(
            "comma-separated exact origin sublattices for the checker; with "
            "--derived-targets, defaults to candidate-pool metadata"
        ),
    )
    parser.add_argument(
        "--through-origin-mod",
        type=int,
        default=1,
        help=(
            "restrict every fibre to target c=0 modulo gcd(h,Q); for Q=2 "
            "all even-order lines pass through the origin parity cell"
        ),
    )
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    source_payload = json.loads(args.candidate_pool.read_text())
    if args.derived_targets and not args.derived_pool:
        raise SystemExit("--derived-targets requires --derived-pool")
    derived_target_by_prime = {}
    if args.derived_targets:
        for row in source_payload["choices"]:
            prime = int(row["p"])
            if prime in derived_target_by_prime:
                raise RuntimeError(f"duplicate derived target prime {prime}")
            modulus = int(row["target_modulus"])
            residue = int(row["target_residue"]) % modulus
            if int(row["h"]) % modulus:
                raise RuntimeError(
                    f"derived target modulus does not divide h for prime {prime}"
                )
            derived_target_by_prime[prime] = (residue, modulus)

    candidates = exact_greedy.load_candidates(args.candidate_pool, args.derived_pool)
    if args.candidate_count is not None:
        candidates = candidates[: args.candidate_count]
    if args.period:
        candidates = [item for item in candidates if args.period % item[0] == 0]
    if args.min_parent_common:
        if not args.derived_pool:
            raise SystemExit("--min-parent-common requires --derived-pool")
        parent_common_by_prime = {
            int(row["p"]): int(row.get("parent_common", 1))
            for row in source_payload["choices"]
        }
        candidates = [
            item
            for item in candidates
            if parent_common_by_prime.get(item[1], 1)
            >= args.min_parent_common
        ]
    if args.power < 1:
        raise SystemExit("--power must be positive")
    if args.through_origin_mod < 1:
        raise SystemExit("--through-origin-mod must be positive")
    if args.random_coordinate_bits < 32:
        raise SystemExit("--random-coordinate-bits must be at least 32")
    if args.diversity_quota < 1:
        raise SystemExit("--diversity-quota must be positive")
    if args.diversity_target_cap < 0:
        raise SystemExit("--diversity-target-cap must be nonnegative")
    if args.required_coverage < 1:
        raise SystemExit("--required-coverage must be positive")
    if args.stream_targets and args.repair_mode == "component":
        raise SystemExit(
            "--stream-targets does not support --repair-mode component"
        )
    if args.stream_cache_file and not args.stream_targets:
        raise SystemExit("--stream-cache-file requires --stream-targets")
    if args.exact_when_random_misses_at_most < 0:
        raise SystemExit("--exact-when-random-misses-at-most must be nonnegative")
    raw_count = len(candidates)
    if args.derived_targets:
        candidates = [
            item for item in candidates if item[1] in derived_target_by_prime
        ]
    else:
        candidates = [
            item
            for item in candidates
            if ((item[1] - 1) // 2)
            % math.gcd(
                (item[1] - 1) // item[0],
                math.gcd(args.power, item[1] - 1),
            )
            == 0
        ]
    if not 0 <= args.fixed_prefix <= len(candidates):
        raise SystemExit("--fixed-prefix is outside the compatible candidate pool")
    by_prime = {p: index for index, (_h, p, _a, _b, _o2, _o3) in enumerate(candidates)}
    normalizers = tuple(
        int(value) for value in args.normalize_primes.split(",") if value
    )
    diversity_primes = tuple(
        int(value) for value in args.diversity_primes.split(",") if value
    )
    diversity_coordinate_moduli = tuple(
        int(value)
        for value in args.diversity_coordinate_moduli.split(",")
        if value
    )
    exact_top_components = tuple(
        int(value)
        for value in args.expand_exact_top_components.split(",")
        if value
    )
    if any(value < 2 for value in exact_top_components):
        raise SystemExit("--expand-exact-top-components must contain primes")
    exact_component_digits = []
    for raw_spec in args.expand_exact_component_digits.split(","):
        if not raw_spec:
            continue
        parts = raw_spec.split(":")
        if len(parts) != 2:
            raise SystemExit(
                "--expand-exact-component-digits entries must be prime:digit"
            )
        prime, digit = map(int, parts)
        if prime < 2 or digit < 0:
            raise SystemExit(
                "--expand-exact-component-digits contains an invalid pair"
            )
        exact_component_digits.append((prime, digit))
    exact_component_digits = tuple(exact_component_digits)
    exact_component_digit_groups = []
    for raw_group in args.expand_exact_component_digit_groups.split(";"):
        if not raw_group:
            continue
        group = []
        for raw_spec in raw_group.split(","):
            parts = raw_spec.split(":")
            if len(parts) != 2:
                raise SystemExit(
                    "--expand-exact-component-digit-groups entries "
                    "must be prime:digit"
                )
            prime, digit = map(int, parts)
            if prime < 2 or digit < 0:
                raise SystemExit(
                    "--expand-exact-component-digit-groups contains "
                    "an invalid pair"
                )
            group.append((prime, digit))
        exact_component_digit_groups.append(tuple(group))
    exact_component_digit_groups = tuple(
        exact_component_digit_groups
    )
    if sum(
        bool(mode)
        for mode in (
            exact_top_components,
            exact_component_digits,
            exact_component_digit_groups,
        )
    ) > 1:
        raise SystemExit(
            "choose only one exact component-tile expansion mode"
        )
    try:
        # Validate the components once before a long synthesis run.  The
        # returned one-point tile is deliberately discarded.
        expand_top_component_tiles(
            [(0, 0)],
            [item[0] for item in candidates],
            exact_top_components,
        )
        expand_component_digit_tiles(
            [(0, 0)],
            [item[0] for item in candidates],
            exact_component_digits,
        )
        expand_component_digit_tile_groups(
            [(0, 0)],
            [item[0] for item in candidates],
            exact_component_digit_groups,
        )
    except ValueError as error:
        raise SystemExit(str(error)) from error
    diversity_coordinate_schedule = []
    if args.diversity_coordinate_schedule:
        if args.diversity_stage_rounds < 1:
            raise SystemExit(
                "--diversity-stage-rounds must be positive with a schedule"
            )
        for raw_stage in args.diversity_coordinate_schedule.split(";"):
            raw_stage = raw_stage.strip()
            if raw_stage.lower() in ("", "-", "none"):
                diversity_coordinate_schedule.append(())
                continue
            stage = tuple(
                int(value) for value in raw_stage.split(",") if value
            )
            if not stage or any(value < 1 for value in stage):
                raise SystemExit("invalid diversity coordinate schedule stage")
            diversity_coordinate_schedule.append(stage)
    elif args.diversity_stage_rounds:
        raise SystemExit(
            "--diversity-stage-rounds requires a coordinate schedule"
        )
    if args.diversity_schedule_cycle and not diversity_coordinate_schedule:
        raise SystemExit(
            "--diversity-schedule-cycle requires a coordinate schedule"
        )
    freeze_target_coordinate_moduli = tuple(
        int(value)
        for value in args.freeze_target_coordinate_moduli.split(",")
        if value
    )
    if any(value < 1 for value in freeze_target_coordinate_moduli):
        raise SystemExit("--freeze-target-coordinate-moduli must be positive")
    freeze_target_period = math.lcm(
        1, *freeze_target_coordinate_moduli
    )
    checker_coordinate_cell = []
    for raw_spec in args.checker_coordinate_cell.split(","):
        if not raw_spec:
            continue
        parts = raw_spec.split(":")
        if len(parts) != 3:
            raise SystemExit(
                "--checker-coordinate-cell entries must be modulus:k:l"
            )
        modulus, kr, lr = map(int, parts)
        if modulus < 1:
            raise SystemExit("checker coordinate moduli must be positive")
        checker_coordinate_cell.append(
            (modulus, kr % modulus, lr % modulus)
        )
    checker_coordinate_cell = tuple(checker_coordinate_cell)
    if checker_coordinate_cell and args.random_checks:
        raise SystemExit(
            "--random-checks must be 0 with --checker-coordinate-cell"
        )
    fixed_primes = tuple(
        int(value) for value in args.fixed_primes.split(",") if value
    )
    if args.algebraic_primes is not None:
        algebraic_primes = tuple(
            int(value)
            for value in args.algebraic_primes.split(",")
            if value
        )
    elif args.derived_targets:
        algebraic_primes = tuple(
            int(value)
            for value in source_payload.get("algebraic_primes", ())
        )
    else:
        algebraic_primes = tuple(
            prime for prime in exact_uncovered.factor(args.power) if prime % 2
        )
    sophie_germain = (
        bool(source_payload.get("sophie_germain", False))
        if args.derived_targets
        else args.power % 4 == 0
    )
    if not all(prime in by_prime for prime in normalizers):
        raise SystemExit("normalization prime missing from selected candidates")
    if not all(prime in by_prime for prime in diversity_primes):
        raise SystemExit("diversity prime missing from selected candidates")
    if args.diversity_target_cap and not diversity_primes:
        raise SystemExit(
            "--diversity-target-cap requires --diversity-primes"
        )
    if not all(prime in by_prime for prime in fixed_primes):
        raise SystemExit("fixed prime missing from selected candidates")
    zero_fixed = {by_prime[prime] for prime in normalizers}
    saved_fixed = {by_prime[prime] for prime in fixed_primes}
    saved_fixed.update(range(args.fixed_prefix))
    fixed = zero_fixed | saved_fixed
    mutable = [index for index in range(len(candidates)) if index not in fixed]

    phase_source = (
        args.phase_file if args.phase_file.exists() else args.initial_phase_file
    )
    saved_phases = (
        json.loads(phase_source.read_text())
        if phase_source and phase_source.exists()
        else {}
    )
    if freeze_target_period > 1:
        missing_saved = [
            p for _h, p, _a, _b, _ord2, _ord3 in candidates
            if str(p) not in saved_phases
        ]
        if missing_saved:
            raise RuntimeError(
                "component freezing requires a saved phase for every "
                f"candidate; missing {missing_saved[:5]}"
            )
    valid_moduli = []
    valid_residues = []
    initial_values = []
    for h, p, _a, _b, _ord2, _ord3 in candidates:
        if args.derived_targets:
            residue, modulus = derived_target_by_prime[p]
        else:
            d = math.gcd(args.power, p - 1)
            step = (p - 1) // h
            divisor = math.gcd(step, d)
            half = (p - 1) // 2
            if half % divisor:
                raise RuntimeError(f"prime {p} has no power-compatible target")
            modulus = d // divisor
            residue = (
                0
                if modulus == 1
                else (half // divisor)
                * pow(step // divisor, -1, modulus)
                % modulus
            )
        origin_modulus = math.gcd(h, args.through_origin_mod)
        merged = merge_congruences(residue, modulus, 0, origin_modulus)
        if merged is None:
            raise RuntimeError(
                f"prime {p} has incompatible power and origin target constraints"
            )
        residue, modulus = merged
        frozen_modulus = math.gcd(h, freeze_target_period)
        if frozen_modulus > 1:
            frozen_residue = int(saved_phases[str(p)]) % frozen_modulus
            merged = merge_congruences(
                residue,
                modulus,
                frozen_residue,
                frozen_modulus,
            )
            if merged is None:
                raise RuntimeError(
                    f"saved frozen components for {p} are incompatible"
                )
            residue, modulus = merged
        if by_prime[p] in saved_fixed and str(p) not in saved_phases:
            raise RuntimeError(f"fixed prime {p} has no saved phase")
        value = int(saved_phases.get(str(p), residue)) % h
        if value % modulus != residue:
            if by_prime[p] in saved_fixed:
                raise RuntimeError(f"saved fixed phase for {p} is incompatible")
            value = residue
        valid_moduli.append(modulus)
        valid_residues.append(residue)
        initial_values.append(value)
    assignment = np.asarray(initial_values, dtype=np.uint32)
    valid_moduli = np.asarray(valid_moduli, dtype=np.uint32)
    valid_residues = np.asarray(valid_residues, dtype=np.uint32)
    candidate_moduli = np.asarray([item[0] for item in candidates], dtype=np.uint32)
    for index in zero_fixed:
        if valid_residues[index] != 0:
            raise RuntimeError("normalization target zero is not power-compatible")
        assignment[index] = 0

    points = []
    seen = set()

    def fixed_covered(point):
        k, l = point
        return any(
            (
                candidates[index][2] * (k % candidates[index][0])
                + candidates[index][3] * (l % candidates[index][0])
                - int(assignment[index])
            )
            % candidates[index][0]
            == 0
            for index in fixed
        )

    if args.points_file.exists():
        for raw_k, raw_l in json.loads(args.points_file.read_text()):
            point = (int(raw_k), int(raw_l))
            if any(
                point[0] % modulus != kr or point[1] % modulus != lr
                for modulus, kr, lr in checker_coordinate_cell
            ):
                continue
            if any(
                point[0] % prime == 0 and point[1] % prime == 0
                for prime in algebraic_primes
            ):
                continue
            if sophie_germain and point[0] % 4 == 2 and point[1] % 4 == 0:
                continue
            if args.drop_fixed_covered_points and fixed_covered(point):
                continue
            if point not in seen:
                seen.add(point)
                points.append(point)
    rng = random.Random(
        args.seed if args.seed is not None else 203_5000 + len(points)
    )
    if args.randomize_mutable_initial:
        for index in mutable:
            modulus = int(valid_moduli[index])
            residue = int(valid_residues[index])
            choices = (
                int(candidate_moduli[index]) - 1 - residue
            ) // modulus + 1
            assignment[index] = residue + modulus * rng.randrange(choices)
    print(
        f"candidates={len(candidates)}/{raw_count} points={len(points)} "
        f"density={sum(1 / item[0] for item in candidates):.12f} "
        f"randomized_initial={args.randomize_mutable_initial} "
        f"random_coordinate_bits={args.random_coordinate_bits}",
        flush=True,
    )

    targets = None
    stream_coordinate_state = None
    stream_cover = None
    if args.stream_targets:
        build_seconds = 0.0
        cache_prefix = 0
        cache_status = "disabled"
        if args.stream_cache_file:
            (
                stream_coordinate_state,
                stream_cover,
                cache_prefix,
                cache_status,
            ) = load_streaming_cache(
                args.stream_cache_file,
                candidates,
                assignment,
                points,
                np,
            )
            if stream_coordinate_state is not None and cache_prefix < len(
                points
            ):
                addition = streaming_coordinate_state(
                    points[cache_prefix:],
                    candidates,
                    np,
                    force_components=(
                        stream_coordinate_state["mode"] == "components"
                    ),
                )
                extended = append_streaming_coordinate_state(
                    stream_coordinate_state,
                    addition,
                    np,
                )
                if extended is None:
                    stream_coordinate_state = streaming_coordinate_state(
                        points,
                        candidates,
                        np,
                    )
                    stream_cover = None
                    cache_status = "mode-transition"
                else:
                    stream_coordinate_state = extended
                    if stream_cover is not None:
                        new_cover = streaming_cover_counts(
                            addition,
                            candidates,
                            assignment,
                            np,
                        )
                        stream_cover = np.concatenate(
                            (stream_cover, new_cover)
                        )
                    cache_status = f"{cache_status}+appended"
            if (
                stream_coordinate_state is not None
                and stream_cover is None
            ):
                stream_cover = streaming_cover_counts(
                    stream_coordinate_state,
                    candidates,
                    assignment,
                    np,
                )
                cache_status = f"{cache_status}+recounted"
            print(
                f"stream_cache={cache_status} prefix={cache_prefix} "
                f"points={len(points)} file={args.stream_cache_file}",
                flush=True,
            )
    else:
        targets, build_seconds = build_targets(points, candidates, np)
    for round_no in range(1, args.rounds + 1):
        repair_started = time.monotonic()
        if args.stream_targets:
            (
                solved,
                steps,
                before,
                after,
                stream_coordinate_state,
                stream_cover,
            ) = repair_streaming(
                points,
                candidates,
                assignment,
                mutable,
                rng,
                np,
                args.max_steps,
                args.sample_size,
                valid_moduli,
                valid_residues,
                candidate_moduli,
                args.required_coverage,
                coordinate_descent=args.repair_mode == "coordinate",
                coordinate_state=stream_coordinate_state,
                initial_cover=stream_cover,
                return_state=True,
            )
        else:
            repair_function = {
                "point": repair,
                "coordinate": repair_coordinate,
                "component": repair_component,
            }[args.repair_mode]
            solved, steps, before, after = repair_function(
                targets,
                assignment,
                mutable,
                rng,
                np,
                args.max_steps,
                args.sample_size,
                valid_moduli,
                valid_residues,
                candidate_moduli,
                args.required_coverage,
            )
        repair_seconds = time.monotonic() - repair_started
        print(
            f"round={round_no} points={len(points)} matrix_s={build_seconds:.3f} "
            f"repair={solved} steps={steps} misses={before}->{after} "
            f"repair_s={repair_seconds:.3f}",
            flush=True,
        )
        if not solved:
            args.phase_file.write_text(
                json.dumps(
                    {str(item[1]): int(assignment[index]) for index, item in enumerate(candidates)}
                )
                + "\n"
            )
            return 3
        args.phase_file.write_text(
            json.dumps(
                {str(item[1]): int(assignment[index]) for index, item in enumerate(candidates)}
            )
            + "\n"
        )
        if args.repair_only:
            if args.stream_cache_file:
                save_streaming_cache(
                    args.stream_cache_file,
                    stream_coordinate_state,
                    stream_cover,
                    candidates,
                    assignment,
                    points,
                    np,
                )
            print(
                f"round={round_no} repair_only=True checker=skipped",
                flush=True,
            )
            return 0

        random_miss_points, tested = random_misses(
            candidates,
            assignment,
            args.random_coordinate_bits,
            rng,
            args.random_checks,
            args.random_batch,
            algebraic_primes,
            sophie_germain,
        )
        misses = random_miss_points
        exact_misses = []
        engine = "random"
        meta = {"sat": True, "engine": engine}
        round_diversity_coordinate_moduli = diversity_coordinate_moduli
        if diversity_coordinate_schedule:
            stage_index = (
                (round_no - 1) // args.diversity_stage_rounds
            )
            if args.diversity_schedule_cycle:
                stage_index %= len(diversity_coordinate_schedule)
            else:
                stage_index = min(
                    stage_index,
                    len(diversity_coordinate_schedule) - 1,
                )
            round_diversity_coordinate_moduli = (
                diversity_coordinate_schedule[stage_index]
            )
        if len(random_miss_points) <= args.exact_when_random_misses_at_most:
            rows = [
                exact_greedy.as_row(item, int(assignment[index]))
                for index, item in enumerate(candidates)
            ]
            if (
                not args.retain_target_matrix_during_exact
                and targets is not None
            ):
                targets = None
                gc.collect()
            exact_misses, exact_meta = exact_uncovered.find_uncovered(
                rows,
                max_component=args.max_component,
                limit=args.exact_batch,
                solver_name=args.checker_solver,
                diversity_primes=(
                    diversity_primes
                    if diversity_primes
                    or round_diversity_coordinate_moduli
                    else normalizers
                ),
                diversity_coordinate_moduli=(
                    round_diversity_coordinate_moduli
                ),
                diversity_quota=args.diversity_quota,
                diversity_target_cap=args.diversity_target_cap,
                algebraic_primes=algebraic_primes,
                sophie_germain=sophie_germain,
                fixed_coordinate_residues=checker_coordinate_cell,
            )
            if not exact_misses and random_miss_points:
                raise AssertionError(
                    "exact checker proved coverage despite random misses"
                )
            misses = list(random_miss_points) + list(exact_misses)
            engine = "exact" if not random_miss_points else "random+exact"
            meta = {
                "sat": bool(exact_misses),
                "engine": engine,
                "random_misses": len(random_miss_points),
                "exact": exact_meta,
            }
            if not exact_misses:
                payload = {
                    "candidate_pool": str(args.candidate_pool),
                    "candidate_count": len(candidates),
                    "choices": rows,
                    "checker": meta,
                    "synthesis": "local_phase_cegis",
                    "power": args.power,
                    "derived_targets": args.derived_targets,
                    "algebraic_primes": algebraic_primes,
                    "sophie_germain": sophie_germain,
                    "checker_coordinate_cell": checker_coordinate_cell,
                }
                args.output.write_text(json.dumps(payload, indent=2) + "\n")
                result_label = (
                    "CELL_COVER" if checker_coordinate_cell else "COVER"
                )
                print(
                    f"{result_label} rows={len(rows)} output={args.output}",
                    flush=True,
                )
                return 0
        if args.checker_checkpoint_file:
            coordinate_fingerprints = []
            if round_diversity_coordinate_moduli:
                coordinate_fingerprints = sorted(
                    {
                        tuple(
                            coordinate
                            for modulus in round_diversity_coordinate_moduli
                            for coordinate in (
                                int(k) % modulus,
                                int(l) % modulus,
                            )
                        )
                        for k, l in exact_misses
                    }
                )
            args.checker_checkpoint_file.write_text(
                json.dumps(
                    {
                        "round": round_no,
                        "engine": engine,
                        "phase_file": str(args.phase_file),
                        "random_misses": random_miss_points,
                        "exact_misses": exact_misses,
                        "diversity_coordinate_moduli": list(
                            round_diversity_coordinate_moduli
                        ),
                        "exact_coordinate_fingerprints": (
                            coordinate_fingerprints
                        ),
                        "checker": meta,
                    },
                    indent=2,
                )
                + "\n"
            )
        expanded_exact_misses = []
        if exact_misses and (
            exact_top_components
            or exact_component_digits
            or exact_component_digit_groups
        ):
            if exact_component_digit_groups:
                expanded_exact_misses = (
                    expand_component_digit_tile_groups(
                        exact_misses,
                        candidate_moduli,
                        exact_component_digit_groups,
                    )
                )
            elif exact_component_digits:
                expanded_exact_misses = expand_component_digit_tiles(
                    exact_misses,
                    candidate_moduli,
                    exact_component_digits,
                )
            else:
                expanded_exact_misses = expand_top_component_tiles(
                    exact_misses,
                    candidate_moduli,
                    exact_top_components,
                )
            lesson_points = [
                *random_miss_points,
                *expanded_exact_misses,
            ]
        else:
            lesson_points = (
                tested
                if engine == "random" and args.retain_random_tested
                else misses
            )
        new_points = []
        for point in lesson_points:
            point = (int(point[0]), int(point[1]))
            if any(
                point[0] % modulus != kr or point[1] % modulus != lr
                for modulus, kr, lr in checker_coordinate_cell
            ):
                continue
            if any(
                point[0] % prime == 0 and point[1] % prime == 0
                for prime in algebraic_primes
            ):
                continue
            if (
                sophie_germain
                and point[0] % 4 == 2
                and point[1] % 4 == 0
            ):
                continue
            if args.drop_fixed_covered_points and fixed_covered(point):
                continue
            if point not in seen:
                seen.add(point)
                points.append(point)
                new_points.append(point)
        args.points_file.write_text(json.dumps(points) + "\n")
        print(
            f"round={round_no} checker={engine} "
            f"misses={len(misses)}/"
            f"{len(tested) if engine == 'random' else 'exact'} "
            f"random={len(random_miss_points)} "
            f"exact={len(exact_misses) if engine != 'random' else 0} "
            f"exact_tile={len(expanded_exact_misses) or '-'} "
            f"diversity={','.join(map(str, round_diversity_coordinate_moduli)) or '-'} "
            f"new={len(new_points)} total={len(points)}",
            flush=True,
        )
        if not new_points:
            raise RuntimeError("checker returned no new point")
        if round_no == args.rounds and not args.stream_cache_file:
            # The final witness batch is retained in the points checkpoint
            # for a later continuation, but no next learner round exists in
            # this invocation, so rebuilding its coordinate state is wasted.
            continue
        if args.stream_targets:
            update_started = time.monotonic()
            new_coordinate_state = streaming_coordinate_state(
                new_points,
                candidates,
                np,
                force_components=(
                    stream_coordinate_state["mode"] == "components"
                ),
            )
            new_cover = streaming_cover_counts(
                new_coordinate_state,
                candidates,
                assignment,
                np,
            )
            extended_coordinate_state = append_streaming_coordinate_state(
                stream_coordinate_state,
                new_coordinate_state,
                np,
            )
            if extended_coordinate_state is None:
                stream_coordinate_state = streaming_coordinate_state(
                    points,
                    candidates,
                    np,
                )
            else:
                stream_coordinate_state = extended_coordinate_state
            stream_cover = np.concatenate(
                (stream_cover, new_cover)
            )
            build_seconds = time.monotonic() - update_started
            if args.stream_cache_file:
                save_streaming_cache(
                    args.stream_cache_file,
                    stream_coordinate_state,
                    stream_cover,
                    candidates,
                    assignment,
                    points,
                    np,
                )
        elif engine != "random":
            if args.retain_target_matrix_during_exact:
                new_targets, build_seconds = build_targets(
                    new_points, candidates, np
                )
                targets = np.asfortranarray(
                    np.vstack((targets, new_targets))
                )
            else:
                targets, build_seconds = build_targets(points, candidates, np)
        else:
            new_targets, build_seconds = build_targets(new_points, candidates, np)
            targets = np.vstack((targets, new_targets))
    return 4


if __name__ == "__main__":
    raise SystemExit(main())
