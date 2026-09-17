#!/usr/bin/env python3
"""
hapboundary4.py

Pairwise estimation of a disease-haplotype boundary around a known anchor.

Key difference from the earlier any-state implementation:

1. --sample1 is a confirmed affected sample.
2. --anchor-allele specifies which phased haplotype of sample1 (1 or 2)
   carries the pathogenic allele at the anchor.
3. At the anchor window, the script compares only the two candidate states
   involving that fixed sample1 haplotype and infers the matching haplotype
   of sample2.
4. A strict boundary follows that fixed state without allowing label changes.
5. A switch-tolerant boundary keeps sample1's confirmed disease haplotype fixed,
   but permits sample2's B1/B2 label to switch when the alternative state is
   strongly and persistently supported across multiple consecutive windows.

Both strict and switch-tolerant boundaries are reported. This separates a true
loss of sharing from a likely statistical-phasing label switch in sample2.

The family file must contain three columns (header names are flexible):

    sample_id    family_id    affected

Accepted affected values include 1/0, affected/unaffected, case/control,
yes/no, true/false.

The anchor is still required even when the pathogenic deletion/variant was
excluded from the phased VCF. It identifies the genomic location where the
confirmed disease haplotype is anchored.

Only Python's standard library is required. Input must be VCF or VCF.GZ.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import sys
from bisect import bisect_left
from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Sequence, Tuple


STATE_NAMES = ("A1-B1", "A1-B2", "A2-B1", "A2-B2")
STATE_INDEX = {name: i for i, name in enumerate(STATE_NAMES)}
CONFIDENT = "CONFIDENT_SHARED"
POSSIBLE = "POSSIBLE_SHARED"
NON_SHARED = "NON_SHARED"


@dataclass(frozen=True)
class SampleMetadata:
    family_id: str
    affected: Optional[bool]


@dataclass(frozen=True)
class Marker:
    chrom: str
    pos: int
    cm: Optional[float]
    marker_id: str
    gt_a: str
    gt_b: str
    votes: Tuple[int, int, int, int]
    ibs0: int
    marker_type: str


@dataclass(frozen=True)
class WindowResult:
    index: int
    start_marker_index: int
    end_marker_index: int
    center_marker_index: int
    start_bp: int
    end_bp: int
    center_bp: int
    start_cm: Optional[float]
    end_cm: Optional[float]
    center_cm: Optional[float]
    span_cm: Optional[float]
    vote_counts: Tuple[int, int, int, int]
    best_count: int
    best_fraction: float
    best_states: Tuple[str, ...]
    selected_count: int
    selected_fraction: float
    classification: str
    ibs0_count: int
    tolerant_state: str = ""
    tolerant_count: int = 0
    tolerant_fraction: float = 0.0
    tolerant_classification: str = NON_SHARED
    switch_event: str = ""


@dataclass(frozen=True)
class SwitchEvent:
    direction: str
    window_index: int
    center_bp: int
    center_cm: Optional[float]
    from_state: str
    to_state: str
    from_count: int
    to_count: int
    margin: int
    supporting_windows: int


class GeneticMap:
    def __init__(self, positions: Sequence[int], cms: Sequence[float]) -> None:
        if len(positions) != len(cms) or len(positions) < 2:
            raise ValueError("Genetic map must contain at least two valid points.")

        paired = sorted(zip(positions, cms), key=lambda x: x[0])
        dedup_positions: List[int] = []
        dedup_cms: List[float] = []
        for pos, cm in paired:
            if dedup_positions and pos == dedup_positions[-1]:
                dedup_cms[-1] = cm
            else:
                dedup_positions.append(pos)
                dedup_cms.append(cm)

        if len(dedup_positions) < 2:
            raise ValueError("Genetic map has fewer than two unique positions.")
        self.positions = dedup_positions
        self.cms = dedup_cms

    def interpolate(self, pos: int) -> Optional[float]:
        idx = bisect_left(self.positions, pos)
        if idx < len(self.positions) and self.positions[idx] == pos:
            return self.cms[idx]
        if idx == 0 or idx == len(self.positions):
            return None

        left_pos = self.positions[idx - 1]
        right_pos = self.positions[idx]
        left_cm = self.cms[idx - 1]
        right_cm = self.cms[idx]
        if right_pos == left_pos:
            return left_cm
        fraction = (pos - left_pos) / (right_pos - left_pos)
        return left_cm + fraction * (right_cm - left_cm)


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom.upper()


def open_text(path: str):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def split_fields(line: str) -> List[str]:
    return re.split(r"\s+", line.strip())


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def clean_header(value: str) -> str:
    value = value.strip().lower().replace("#", "")
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def midpoint_int(a: int, b: int) -> int:
    return int(round((a + b) / 2.0))


def format_float(value: Optional[float], digits: int = 6) -> str:
    return "NA" if value is None else f"{value:.{digits}f}"


def bool_to_text(value: Optional[bool]) -> str:
    if value is None:
        return "NA"
    return "1" if value else "0"


def parse_affected(value: str) -> Optional[bool]:
    token = value.strip().lower()
    true_values = {"1", "affected", "case", "yes", "y", "true", "t", "patient"}
    false_values = {
        "0",
        "unaffected",
        "control",
        "no",
        "n",
        "false",
        "f",
        "healthy",
    }
    missing_values = {"", ".", "na", "nan", "unknown", "-9"}
    if token in true_values:
        return True
    if token in false_values:
        return False
    if token in missing_values:
        return None
    raise ValueError(f"Unrecognized affected value: {value!r}")


# ---------------------------------------------------------------------------
# Genetic map parsing
# ---------------------------------------------------------------------------


def detect_map_columns(
    rows: List[List[str]], target_chrom: str
) -> Tuple[int, Optional[int], int, int]:
    first = rows[0]
    cleaned = [clean_header(x) for x in first]

    chrom_aliases = {"chrom", "chromosome", "chr"}
    bp_aliases = {
        "position",
        "position_bp",
        "pos",
        "bp",
        "base_pair",
        "base_pair_position",
        "physical_position",
    }
    cm_aliases = {
        "cm",
        "position_cm",
        "genetic_map",
        "genetic_map_cm",
        "genetic_position",
        "map_cm",
        "cumulative_cm",
    }

    def find_alias(aliases: set[str]) -> Optional[int]:
        for i, name in enumerate(cleaned):
            if name in aliases:
                return i
        return None

    chrom_col = find_alias(chrom_aliases)
    bp_col = find_alias(bp_aliases)
    cm_col = find_alias(cm_aliases)
    if bp_col is not None and cm_col is not None:
        return 1, chrom_col, bp_col, cm_col

    n = len(first)
    if n >= 4:
        return 0, 0, 3, 2  # PLINK .map: CHR ID cM BP
    if n == 3:
        if normalize_chrom(first[0]) == target_chrom and is_number(first[1]):
            return 0, 0, 1, 2  # chrom bp cM
        return 0, None, 0, 2  # SHAPEIT: bp rate cM
    if n == 2:
        return 0, None, 0, 1
    raise ValueError("Could not infer genetic-map columns.")


def load_genetic_map(path: str, target_chrom: str) -> GeneticMap:
    rows: List[List[str]] = []
    with open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("##"):
                continue
            if stripped.startswith("#"):
                stripped = stripped[1:].strip()
            fields = split_fields(stripped)
            if fields:
                rows.append(fields)

    if not rows:
        raise ValueError(f"No usable rows found in genetic map: {path}")

    data_start, chrom_col, bp_col, cm_col = detect_map_columns(rows, target_chrom)
    positions: List[int] = []
    cms: List[float] = []

    for row in rows[data_start:]:
        required = [bp_col, cm_col] + ([chrom_col] if chrom_col is not None else [])
        if any(col is None or col >= len(row) for col in required):
            continue
        if chrom_col is not None and normalize_chrom(row[chrom_col]) != target_chrom:
            continue
        try:
            bp = int(float(row[bp_col]))
            cm = float(row[cm_col])
        except ValueError:
            continue
        positions.append(bp)
        cms.append(cm)

    if len(positions) < 2:
        raise ValueError(
            f"Fewer than two map points found for chromosome {target_chrom}."
        )
    return GeneticMap(positions, cms)


# ---------------------------------------------------------------------------
# Family/sample metadata
# ---------------------------------------------------------------------------


def load_family_metadata(path: str) -> Dict[str, SampleMetadata]:
    with open_text(path) as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("##")]
    if not lines:
        raise ValueError("Family file contains no usable rows.")

    delimiter = "\t" if "\t" in lines[0] else None
    header = (
        lines[0].rstrip("\n").split("\t")
        if delimiter == "\t"
        else split_fields(lines[0])
    )
    cleaned = [clean_header(x) for x in header]

    sample_aliases = ("sample_id", "sample", "iid", "individual_id", "id")
    family_aliases = ("family_id", "family", "fid")
    affected_aliases = ("affected", "phenotype", "status", "case_control")

    sample_col = next((cleaned.index(x) for x in sample_aliases if x in cleaned), None)
    family_col = next((cleaned.index(x) for x in family_aliases if x in cleaned), None)
    affected_col = next(
        (cleaned.index(x) for x in affected_aliases if x in cleaned), None
    )

    if sample_col is None or family_col is None or affected_col is None:
        raise ValueError(
            "Family file needs sample_id/IID, family_id/FID, and affected/status columns."
        )

    metadata: Dict[str, SampleMetadata] = {}
    for line_number, line in enumerate(lines[1:], start=2):
        fields = (
            line.rstrip("\n").split("\t")
            if delimiter == "\t"
            else split_fields(line)
        )
        if max(sample_col, family_col, affected_col) >= len(fields):
            continue
        sample = fields[sample_col]
        family = fields[family_col]
        try:
            affected = parse_affected(fields[affected_col])
        except ValueError as exc:
            raise ValueError(f"Family file line {line_number}: {exc}") from exc
        metadata[sample] = SampleMetadata(family_id=family, affected=affected)

    return metadata


# ---------------------------------------------------------------------------
# VCF parsing and marker voting
# ---------------------------------------------------------------------------


def parse_gt(sample_field: str, gt_index: int) -> Optional[Tuple[int, int, str]]:
    parts = sample_field.split(":")
    if gt_index >= len(parts):
        return None
    gt = parts[gt_index]
    if gt not in {"0|0", "0|1", "1|0", "1|1"}:
        return None
    return int(gt[0]), int(gt[2]), gt


def classify_marker_type(a: Tuple[int, int], b: Tuple[int, int]) -> str:
    a_hom = a[0] == a[1]
    b_hom = b[0] == b[1]
    if a_hom and b_hom:
        return "SAME_HOMO" if a[0] == b[0] else "IBS0"
    if a_hom and not b_hom:
        return "A_HOMO_B_HET"
    if not a_hom and b_hom:
        return "A_HET_B_HOMO"
    return "HET_HET"


def marker_votes(
    a: Tuple[int, int],
    b: Tuple[int, int],
    marker_mode: str,
) -> Optional[Tuple[Tuple[int, int, int, int], int, str]]:
    """
    Return (votes, ibs0, marker_type), or None when excluded from K.

    marker_mode='original':
        Exclude only same-homozygous sites. Include het-het and IBS0; IBS0 casts
        no votes and therefore lowers the selected vote fraction.

    marker_mode='hom-het':
        Include only homozygous-heterozygous sites. This follows the strict
        informativeness definition. IBS0 and het-het are excluded from K.
    """
    marker_type = classify_marker_type(a, b)
    if marker_type == "SAME_HOMO":
        return None
    if marker_mode == "hom-het" and marker_type not in {
        "A_HOMO_B_HET",
        "A_HET_B_HOMO",
    }:
        return None

    ibs0 = int(marker_type == "IBS0")
    votes = (
        int(a[0] == b[0]),
        int(a[0] == b[1]),
        int(a[1] == b[0]),
        int(a[1] == b[1]),
    )
    return votes, ibs0, marker_type


def load_pair_markers(
    vcf_path: str,
    target_chrom: str,
    sample_a: str,
    sample_b: str,
    genetic_map: GeneticMap,
    marker_mode: str,
) -> List[Marker]:
    markers: List[Marker] = []
    sample_a_col: Optional[int] = None
    sample_b_col: Optional[int] = None
    saw_header = False

    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                if sample_a not in samples:
                    raise ValueError(f"Sample not found in VCF: {sample_a}")
                if sample_b not in samples:
                    raise ValueError(f"Sample not found in VCF: {sample_b}")
                sample_a_col = 9 + samples.index(sample_a)
                sample_b_col = 9 + samples.index(sample_b)
                saw_header = True
                continue
            if line.startswith("#"):
                continue
            if not saw_header or sample_a_col is None or sample_b_col is None:
                raise ValueError("VCF #CHROM header was not found before data rows.")

            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(sample_a_col, sample_b_col):
                continue
            if normalize_chrom(fields[0]) != target_chrom:
                continue

            try:
                pos = int(fields[1])
            except ValueError:
                continue

            ref = fields[3]
            alt = fields[4]
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                continue

            format_keys = fields[8].split(":")
            if "GT" not in format_keys:
                continue
            gt_index = format_keys.index("GT")
            parsed_a = parse_gt(fields[sample_a_col], gt_index)
            parsed_b = parse_gt(fields[sample_b_col], gt_index)
            if parsed_a is None or parsed_b is None:
                continue

            vote_result = marker_votes(
                (parsed_a[0], parsed_a[1]),
                (parsed_b[0], parsed_b[1]),
                marker_mode,
            )
            if vote_result is None:
                continue
            votes, ibs0, marker_type = vote_result
            markers.append(
                Marker(
                    chrom=fields[0],
                    pos=pos,
                    cm=genetic_map.interpolate(pos),
                    marker_id=fields[2],
                    gt_a=parsed_a[2],
                    gt_b=parsed_b[2],
                    votes=votes,
                    ibs0=ibs0,
                    marker_type=marker_type,
                )
            )

    if not saw_header:
        raise ValueError("VCF #CHROM header not found.")
    markers.sort(key=lambda marker: marker.pos)
    return markers


# ---------------------------------------------------------------------------
# Anchor assignment and window classification
# ---------------------------------------------------------------------------


def raw_window_vote_counts(
    markers: Sequence[Marker], k: int
) -> List[Tuple[int, int, int, int]]:
    if k < 3 or k % 2 == 0:
        raise ValueError("K must be an odd integer >= 3.")
    if len(markers) < k:
        raise ValueError(
            f"Only {len(markers)} informative markers are available; K={k}."
        )

    n = len(markers)
    vote_prefix = [[0] * (n + 1) for _ in range(4)]
    for i, marker in enumerate(markers):
        for state in range(4):
            vote_prefix[state][i + 1] = vote_prefix[state][i] + marker.votes[state]

    counts_by_window: List[Tuple[int, int, int, int]] = []
    for start in range(0, n - k + 1):
        end = start + k
        counts_by_window.append(
            tuple(
                vote_prefix[state][end] - vote_prefix[state][start]
                for state in range(4)
            )
        )
    return counts_by_window


def choose_anchor_window_from_markers(
    markers: Sequence[Marker], k: int, anchor_bp: int
) -> int:
    """
    Prefer a K-marker window spanning the anchor. Otherwise choose the window
    whose center marker is closest to the anchor.
    """
    half = k // 2
    candidates: List[Tuple[int, int, int]] = []
    for start in range(0, len(markers) - k + 1):
        end = start + k - 1
        center = start + half
        spans_anchor = int(markers[start].pos <= anchor_bp <= markers[end].pos)
        distance = abs(markers[center].pos - anchor_bp)
        # Sort by: spanning anchor first, then nearest center, then narrower bp span.
        candidates.append((-spans_anchor, distance, markers[end].pos - markers[start].pos))
    return min(range(len(candidates)), key=lambda i: candidates[i])


def infer_sample2_anchor_allele(
    vote_counts: Tuple[int, int, int, int],
    sample1_anchor_allele: int,
    minimum_count: int,
    minimum_margin: int,
) -> Tuple[Optional[int], Optional[str], int, int, int, str]:
    """
    Infer sample2 haplotype at the anchor while fixing sample1's disease allele.

    Returns:
      (sample2_allele, selected_state, selected_count, competing_count,
       margin, status)
    """
    if sample1_anchor_allele == 1:
        candidates = ((1, "A1-B1", vote_counts[0]), (2, "A1-B2", vote_counts[1]))
    else:
        candidates = ((1, "A2-B1", vote_counts[2]), (2, "A2-B2", vote_counts[3]))

    ranked = sorted(candidates, key=lambda item: item[2], reverse=True)
    best_allele, best_state, best_count = ranked[0]
    _, _, second_count = ranked[1]
    margin = best_count - second_count

    if best_count < minimum_count:
        return None, None, best_count, second_count, margin, "INSUFFICIENT_ANCHOR_SUPPORT"
    if margin < minimum_margin:
        return None, None, best_count, second_count, margin, "AMBIGUOUS_ANCHOR_ASSIGNMENT"
    return best_allele, best_state, best_count, second_count, margin, "ASSIGNED"


def classify_windows_for_selected_state(
    markers: Sequence[Marker],
    k: int,
    selected_state: str,
    confident_fraction: float,
    possible_fraction: float,
) -> List[WindowResult]:
    if not (0.0 < possible_fraction <= confident_fraction <= 1.0):
        raise ValueError(
            "Require 0 < possible_fraction <= confident_fraction <= 1."
        )

    state_index = STATE_INDEX[selected_state]
    counts_by_window = raw_window_vote_counts(markers, k)
    confident_count = math.ceil(k * confident_fraction)
    possible_count = math.ceil(k * possible_fraction)
    half = k // 2

    ibs0_prefix = [0] * (len(markers) + 1)
    for i, marker in enumerate(markers):
        ibs0_prefix[i + 1] = ibs0_prefix[i] + marker.ibs0

    results: List[WindowResult] = []
    for start, counts in enumerate(counts_by_window):
        end = start + k - 1
        center = start + half
        best_count = max(counts)
        best_states = tuple(
            STATE_NAMES[i] for i, value in enumerate(counts) if value == best_count
        )
        selected_count = counts[state_index]

        if selected_count >= confident_count:
            classification = CONFIDENT
        elif selected_count >= possible_count:
            classification = POSSIBLE
        else:
            classification = NON_SHARED

        start_cm = markers[start].cm
        end_cm = markers[end].cm
        center_cm = markers[center].cm
        span_cm = (
            end_cm - start_cm
            if start_cm is not None and end_cm is not None
            else None
        )

        results.append(
            WindowResult(
                index=start,
                start_marker_index=start,
                end_marker_index=end,
                center_marker_index=center,
                start_bp=markers[start].pos,
                end_bp=markers[end].pos,
                center_bp=markers[center].pos,
                start_cm=start_cm,
                end_cm=end_cm,
                center_cm=center_cm,
                span_cm=span_cm,
                vote_counts=counts,
                best_count=best_count,
                best_fraction=best_count / k,
                best_states=best_states,
                selected_count=selected_count,
                selected_fraction=selected_count / k,
                classification=classification,
                ibs0_count=ibs0_prefix[end + 1] - ibs0_prefix[start],
            )
        )

    return results


def classify_vote_count(
    count: int,
    confident_count: int,
    possible_count: int,
) -> str:
    if count >= confident_count:
        return CONFIDENT
    if count >= possible_count:
        return POSSIBLE
    return NON_SHARED


def paired_sample2_state_indices(selected_state: str) -> Tuple[int, int]:
    """Return the two states that keep sample1's anchor allele fixed."""
    idx = STATE_INDEX[selected_state]
    if idx in (0, 1):
        return 0, 1  # A1-B1, A1-B2
    return 2, 3      # A2-B1, A2-B2


def apply_sample2_switch_tolerance(
    windows: Sequence[WindowResult],
    anchor_index: int,
    selected_state: str,
    k: int,
    confident_fraction: float,
    possible_fraction: float,
    switch_min_fraction: float,
    switch_min_margin: int,
    switch_min_windows: int,
) -> Tuple[List[WindowResult], List[SwitchEvent]]:
    """
    Track the confirmed sample1 haplotype while allowing only sample2 B1/B2
    label switches. A switch is accepted when the alternative state reaches
    switch_min_fraction, exceeds the current state by switch_min_margin, and
    does so for switch_min_windows consecutive sliding windows.

    The accepted switch is applied retroactively to the first window in the
    qualifying run. Left and right sides are traced independently from anchor.
    """
    if not (0.0 < switch_min_fraction <= 1.0):
        raise ValueError("--switch-min-fraction must be in (0, 1].")
    if switch_min_margin < 0:
        raise ValueError("--switch-min-margin must be >= 0.")
    if switch_min_windows < 1:
        raise ValueError("--switch-min-windows must be >= 1.")

    pair = paired_sample2_state_indices(selected_state)
    initial_idx = STATE_INDEX[selected_state]
    alternative_of = {pair[0]: pair[1], pair[1]: pair[0]}
    switch_count = math.ceil(k * switch_min_fraction)
    confident_count = math.ceil(k * confident_fraction)
    possible_count = math.ceil(k * possible_fraction)

    assigned: List[Optional[int]] = [None] * len(windows)
    assigned[anchor_index] = initial_idx
    events: List[SwitchEvent] = []

    def trace(direction: int) -> None:
        current = initial_idx
        # transition contains windows after the alternative first becomes both
        # possible and better than the current state. strong_run is the subset
        # satisfying the stricter switch-acceptance threshold.
        transition: List[int] = []
        strong_run: List[int] = []
        i = anchor_index + direction
        while 0 <= i < len(windows):
            alternative = alternative_of[current]
            current_count = windows[i].vote_counts[current]
            alternative_count = windows[i].vote_counts[alternative]

            transition_candidate = (
                alternative_count >= possible_count
                and alternative_count > current_count
            )
            strongly_qualifies = (
                alternative_count >= switch_count
                and alternative_count - current_count >= switch_min_margin
            )

            if transition_candidate:
                transition.append(i)
                if strongly_qualifies:
                    strong_run.append(i)
                else:
                    strong_run = []

                if len(strong_run) >= switch_min_windows:
                    first_index = transition[0]
                    old_state = current
                    current = alternative
                    # Once the switch is confirmed, bridge the mixed transition
                    # windows back to the first point where the alternative was
                    # possible and superior. This avoids a short artificial gap
                    # caused by overlapping K-marker windows.
                    for transition_index in transition:
                        assigned[transition_index] = current
                    first_window = windows[first_index]
                    events.append(
                        SwitchEvent(
                            direction="LEFT" if direction < 0 else "RIGHT",
                            window_index=first_index,
                            center_bp=first_window.center_bp,
                            center_cm=first_window.center_cm,
                            from_state=STATE_NAMES[old_state],
                            to_state=STATE_NAMES[current],
                            from_count=first_window.vote_counts[old_state],
                            to_count=first_window.vote_counts[current],
                            margin=(
                                first_window.vote_counts[current]
                                - first_window.vote_counts[old_state]
                            ),
                            supporting_windows=switch_min_windows,
                        )
                    )
                    transition = []
                    strong_run = []
            else:
                # A tentative crossover that was not followed by persistent
                # strong evidence is rejected and remains assigned to current.
                for transition_index in transition:
                    assigned[transition_index] = current
                transition = []
                strong_run = []
                assigned[i] = current
            i += direction

        for transition_index in transition:
            assigned[transition_index] = current

    trace(-1)
    trace(1)

    # Defensive fallback; all windows should have been assigned by the traces.
    for i, value in enumerate(assigned):
        if value is None:
            assigned[i] = initial_idx

    event_by_index: Dict[int, List[str]] = {}
    for event in events:
        event_by_index.setdefault(event.window_index, []).append(
            f"{event.direction}:{event.from_state}>{event.to_state}"
        )

    tolerant_windows: List[WindowResult] = []
    for i, window in enumerate(windows):
        state_idx = int(assigned[i])
        count = window.vote_counts[state_idx]
        tolerant_windows.append(
            replace(
                window,
                tolerant_state=STATE_NAMES[state_idx],
                tolerant_count=count,
                tolerant_fraction=count / k,
                tolerant_classification=classify_vote_count(
                    count, confident_count, possible_count
                ),
                switch_event=";".join(event_by_index.get(i, [])),
            )
        )

    events.sort(key=lambda event: event.center_bp)
    return tolerant_windows, events


def contiguous_region(
    windows: Sequence[WindowResult],
    anchor_index: int,
    allowed: set[str],
) -> Optional[Tuple[int, int]]:
    if windows[anchor_index].classification not in allowed:
        return None
    left = anchor_index
    right = anchor_index
    while left > 0 and windows[left - 1].classification in allowed:
        left -= 1
    while right + 1 < len(windows) and windows[right + 1].classification in allowed:
        right += 1
    return left, right


def region_boundary_bp(
    windows: Sequence[WindowResult], region: Optional[Tuple[int, int]]
) -> Tuple[Optional[int], Optional[int]]:
    if region is None:
        return None, None
    left, right = region

    if left == 0:
        start_bp = windows[left].start_bp
    else:
        start_bp = midpoint_int(windows[left - 1].center_bp, windows[left].center_bp)

    if right == len(windows) - 1:
        end_bp = windows[right].end_bp
    else:
        end_bp = midpoint_int(windows[right].center_bp, windows[right + 1].center_bp)

    return start_bp, end_bp


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write_windows(
    path: str,
    windows: Sequence[WindowResult],
    sample_a: str,
    sample_b: str,
    selected_state: str,
    k: int,
) -> None:
    header = [
        "sample1",
        "sample2",
        "window_index",
        "start_bp",
        "end_bp",
        "center_bp",
        "start_cM",
        "end_cM",
        "center_cM",
        "span_cM",
        "K",
        "strict_selected_state",
        "votes_A1_B1",
        "votes_A1_B2",
        "votes_A2_B1",
        "votes_A2_B2",
        "strict_selected_count",
        "strict_selected_fraction",
        "strict_classification",
        "tolerant_selected_state",
        "tolerant_selected_count",
        "tolerant_selected_fraction",
        "tolerant_classification",
        "switch_event",
        "any_state_best_count",
        "any_state_best_fraction",
        "any_state_best_state",
        "n_any_state_best_states",
        "IBS0_count",
    ]
    with open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for w in windows:
            writer.writerow(
                [
                    sample_a,
                    sample_b,
                    w.index,
                    w.start_bp,
                    w.end_bp,
                    w.center_bp,
                    format_float(w.start_cm),
                    format_float(w.end_cm),
                    format_float(w.center_cm),
                    format_float(w.span_cm),
                    k,
                    selected_state,
                    *w.vote_counts,
                    w.selected_count,
                    f"{w.selected_fraction:.6f}",
                    w.classification,
                    w.tolerant_state,
                    w.tolerant_count,
                    f"{w.tolerant_fraction:.6f}",
                    w.tolerant_classification,
                    w.switch_event,
                    w.best_count,
                    f"{w.best_fraction:.6f}",
                    ",".join(w.best_states),
                    len(w.best_states),
                    w.ibs0_count,
                ]
            )


def write_markers(
    path: str,
    markers: Sequence[Marker],
    sample_a: str,
    sample_b: str,
    selected_state: str,
) -> None:
    selected_index = STATE_INDEX[selected_state]
    header = [
        "sample1",
        "sample2",
        "chrom",
        "position_bp",
        "position_cM",
        "marker_id",
        "GT_sample1",
        "GT_sample2",
        "marker_type",
        "vote_A1_B1",
        "vote_A1_B2",
        "vote_A2_B1",
        "vote_A2_B2",
        "selected_state",
        "selected_state_vote",
        "IBS0",
    ]
    with open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for marker in markers:
            writer.writerow(
                [
                    sample_a,
                    sample_b,
                    marker.chrom,
                    marker.pos,
                    format_float(marker.cm),
                    marker.marker_id,
                    marker.gt_a,
                    marker.gt_b,
                    marker.marker_type,
                    *marker.votes,
                    selected_state,
                    marker.votes[selected_index],
                    marker.ibs0,
                ]
            )


def write_switches(path: str, events: Sequence[SwitchEvent]) -> None:
    header = [
        "direction",
        "window_index",
        "center_bp",
        "center_cM",
        "from_state",
        "to_state",
        "from_count",
        "to_count",
        "margin",
        "supporting_windows",
    ]
    with open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for event in events:
            writer.writerow(
                [
                    event.direction,
                    event.window_index,
                    event.center_bp,
                    format_float(event.center_cm),
                    event.from_state,
                    event.to_state,
                    event.from_count,
                    event.to_count,
                    event.margin,
                    event.supporting_windows,
                ]
            )


def region_metrics(
    bp: Tuple[Optional[int], Optional[int]],
    genetic_map: GeneticMap,
) -> Tuple[Optional[int], Optional[int], Optional[int], Optional[float], Optional[float], Optional[float]]:
    start_bp, end_bp = bp
    start_cm = genetic_map.interpolate(start_bp) if start_bp is not None else None
    end_cm = genetic_map.interpolate(end_bp) if end_bp is not None else None
    length_bp = (
        end_bp - start_bp
        if start_bp is not None and end_bp is not None
        else None
    )
    length_cm = (
        end_cm - start_cm
        if start_cm is not None and end_cm is not None
        else None
    )
    return start_bp, end_bp, length_bp, start_cm, end_cm, length_cm


def write_summary(
    path: str,
    sample_a: str,
    sample_b: str,
    metadata_a: SampleMetadata,
    metadata_b: SampleMetadata,
    target_chrom: str,
    anchor_bp: int,
    anchor_window: WindowResult,
    sample1_anchor_allele: int,
    sample2_anchor_allele: int,
    selected_state: str,
    anchor_selected_count: int,
    anchor_competing_count: int,
    anchor_margin: int,
    strict_confident_bp: Tuple[Optional[int], Optional[int]],
    strict_possible_bp: Tuple[Optional[int], Optional[int]],
    tolerant_confident_bp: Tuple[Optional[int], Optional[int]],
    tolerant_possible_bp: Tuple[Optional[int], Optional[int]],
    events: Sequence[SwitchEvent],
    genetic_map: GeneticMap,
    k: int,
    confident_fraction: float,
    possible_fraction: float,
    switch_min_fraction: float,
    switch_min_margin: int,
    switch_min_windows: int,
    marker_mode: str,
    n_markers: int,
) -> None:
    strict_conf = region_metrics(strict_confident_bp, genetic_map)
    strict_poss = region_metrics(strict_possible_bp, genetic_map)
    tolerant_conf = region_metrics(tolerant_confident_bp, genetic_map)
    tolerant_poss = region_metrics(tolerant_possible_bp, genetic_map)

    pair_type = (
        "WITHIN_FAMILY"
        if metadata_a.family_id == metadata_b.family_id
        else "BETWEEN_FAMILY"
    )
    comparison_type = (
        "AFFECTED_AFFECTED"
        if metadata_a.affected is True and metadata_b.affected is True
        else "AFFECTED_UNAFFECTED"
        if metadata_a.affected is True and metadata_b.affected is False
        else "OTHER_OR_UNKNOWN"
    )

    switch_positions = ",".join(str(event.center_bp) for event in events) or "NA"
    switch_descriptions = ";".join(
        f"{event.direction}:{event.from_state}>{event.to_state}@{event.center_bp}"
        for event in events
    ) or "NA"

    header = [
        "sample1", "sample2", "family1", "family2", "affected1", "affected2",
        "pair_type", "comparison_type", "chrom", "anchor_bp",
        "sample1_anchor_allele", "sample2_inferred_anchor_allele",
        "selected_anchor_state", "anchor_window_start_bp", "anchor_window_end_bp",
        "anchor_window_center_bp", "anchor_selected_count", "anchor_competing_count",
        "anchor_assignment_margin", "anchor_selected_fraction",
        "anchor_strict_classification", "anchor_any_state_best_count",
        "anchor_any_state_best_fraction", "anchor_any_state_best_state",
        "K", "confident_fraction", "possible_fraction", "marker_mode",
        "switch_min_fraction", "switch_min_margin", "switch_min_windows",
        "n_markers_used", "n_sample2_phase_switches", "phase_switch_positions_bp",
        "phase_switch_descriptions",
        "strict_confident_start_bp", "strict_confident_end_bp", "strict_confident_length_bp",
        "strict_confident_start_cM", "strict_confident_end_cM", "strict_confident_length_cM",
        "strict_possible_start_bp", "strict_possible_end_bp", "strict_possible_length_bp",
        "strict_possible_start_cM", "strict_possible_end_cM", "strict_possible_length_cM",
        "switch_tolerant_confident_start_bp", "switch_tolerant_confident_end_bp",
        "switch_tolerant_confident_length_bp", "switch_tolerant_confident_start_cM",
        "switch_tolerant_confident_end_cM", "switch_tolerant_confident_length_cM",
        "switch_tolerant_possible_start_bp", "switch_tolerant_possible_end_bp",
        "switch_tolerant_possible_length_bp", "switch_tolerant_possible_start_cM",
        "switch_tolerant_possible_end_cM", "switch_tolerant_possible_length_cM",
        "primary_boundary_mode", "status",
    ]

    def vals(metrics):
        sbp, ebp, lbp, scm, ecm, lcm = metrics
        return [
            sbp if sbp is not None else "NA",
            ebp if ebp is not None else "NA",
            lbp if lbp is not None else "NA",
            format_float(scm), format_float(ecm), format_float(lcm),
        ]

    status = (
        "SWITCH_TOLERANT_POSSIBLE_BOUNDARY"
        if tolerant_poss[0] is not None
        else "ANCHOR_HAPLOTYPE_NOT_SHARED"
    )

    row = [
        sample_a, sample_b, metadata_a.family_id, metadata_b.family_id,
        bool_to_text(metadata_a.affected), bool_to_text(metadata_b.affected),
        pair_type, comparison_type, target_chrom, anchor_bp,
        sample1_anchor_allele, sample2_anchor_allele, selected_state,
        anchor_window.start_bp, anchor_window.end_bp, anchor_window.center_bp,
        anchor_selected_count, anchor_competing_count, anchor_margin,
        f"{anchor_selected_count / k:.6f}", anchor_window.classification,
        anchor_window.best_count, f"{anchor_window.best_fraction:.6f}",
        ",".join(anchor_window.best_states), k,
        f"{confident_fraction:.6f}", f"{possible_fraction:.6f}", marker_mode,
        f"{switch_min_fraction:.6f}", switch_min_margin, switch_min_windows,
        n_markers, len(events), switch_positions, switch_descriptions,
        *vals(strict_conf), *vals(strict_poss), *vals(tolerant_conf), *vals(tolerant_poss),
        "SWITCH_TOLERANT", status,
    ]

    with open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerow(row)

# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def parse_anchor(value: str) -> Tuple[str, int]:
    if ":" not in value:
        raise argparse.ArgumentTypeError(
            "Anchor must be CHROM:POSITION, e.g. 2:47611613"
        )
    chrom, pos_text = value.rsplit(":", 1)
    try:
        pos = int(pos_text.replace(",", ""))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Anchor position must be an integer.") from exc
    return normalize_chrom(chrom), pos


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate strict and sample2-switch-tolerant pairwise disease-haplotype boundaries "
            "using a confirmed sample and its known disease-haplotype index."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="Phased VCF or VCF.GZ")
    parser.add_argument("-m", "--map", required=True, help="Genetic map text file")
    parser.add_argument(
        "-f",
        "--family",
        required=True,
        help="Metadata file with sample_id, family_id, and affected columns",
    )
    parser.add_argument(
        "-a",
        "--anchor",
        required=True,
        type=parse_anchor,
        metavar="CHROM:POS",
        help="Pathogenic-variant/deletion anchor coordinate",
    )
    parser.add_argument(
        "--sample1", required=True, help="Confirmed affected reference sample"
    )
    parser.add_argument("--sample2", required=True, help="Test sample")
    parser.add_argument(
        "--anchor-allele",
        "-anchorallele",
        required=True,
        type=int,
        choices=(1, 2),
        help="Disease-carrying phased haplotype of sample1 at the anchor (1 or 2)",
    )
    parser.add_argument(
        "-k",
        type=int,
        default=51,
        help="Odd number of informative markers per sliding window",
    )
    parser.add_argument(
        "--confident-fraction",
        type=float,
        default=2.0 / 3.0,
        help="Selected-state vote fraction for CONFIDENT_SHARED",
    )
    parser.add_argument(
        "--possible-fraction",
        type=float,
        default=0.5,
        help="Selected-state vote fraction for POSSIBLE_SHARED",
    )
    parser.add_argument(
        "--anchor-min-fraction",
        type=float,
        default=0.5,
        help="Minimum selected candidate fraction required to assign sample2 at anchor",
    )
    parser.add_argument(
        "--anchor-min-margin",
        type=int,
        default=1,
        help="Minimum vote-count difference between sample2 haplotype candidates",
    )
    parser.add_argument(
        "--switch-min-fraction",
        type=float,
        default=2.0 / 3.0,
        help=(
            "Minimum alternative-state fraction required for a sample2 B1/B2 "
            "phase-label switch"
        ),
    )
    parser.add_argument(
        "--switch-min-margin",
        type=int,
        default=5,
        help="Minimum alternative-minus-current vote margin for a phase-label switch",
    )
    parser.add_argument(
        "--switch-min-windows",
        type=int,
        default=5,
        help="Consecutive qualifying windows required to accept a phase-label switch",
    )
    parser.add_argument(
        "--marker-mode",
        choices=("original", "hom-het"),
        default="original",
        help=(
            "original: all informative/non-identical markers including het-het and IBS0; "
            "hom-het: strict homozygous-heterozygous markers only"
        ),
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help="Prefix for .summary.tsv and .windows.tsv",
    )
    parser.add_argument(
        "--write-markers",
        action="store_true",
        help="Also write per-marker votes to <prefix>.markers.tsv",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.sample1 == args.sample2:
        parser.error("--sample1 and --sample2 must be different samples.")
    if args.k < 3 or args.k % 2 == 0:
        parser.error("-k must be an odd integer >= 3.")
    if not (0.0 < args.anchor_min_fraction <= 1.0):
        parser.error("--anchor-min-fraction must be in (0, 1].")
    if args.anchor_min_margin < 0:
        parser.error("--anchor-min-margin must be >= 0.")
    if not (0.0 < args.switch_min_fraction <= 1.0):
        parser.error("--switch-min-fraction must be in (0, 1].")
    if args.switch_min_margin < 0:
        parser.error("--switch-min-margin must be >= 0.")
    if args.switch_min_windows < 1:
        parser.error("--switch-min-windows must be >= 1.")

    target_chrom, anchor_bp = args.anchor

    try:
        genetic_map = load_genetic_map(args.map, target_chrom)
        family_metadata = load_family_metadata(args.family)

        if args.sample1 not in family_metadata:
            raise ValueError(f"Sample1 is absent from family file: {args.sample1}")
        if args.sample2 not in family_metadata:
            raise ValueError(f"Sample2 is absent from family file: {args.sample2}")

        metadata_a = family_metadata[args.sample1]
        metadata_b = family_metadata[args.sample2]
        if metadata_a.affected is not True:
            raise ValueError(
                "--sample1 must be marked affected=1 in the family file because it "
                "is the confirmed disease-haplotype reference."
            )

        markers = load_pair_markers(
            args.input,
            target_chrom,
            args.sample1,
            args.sample2,
            genetic_map,
            args.marker_mode,
        )
        counts_by_window = raw_window_vote_counts(markers, args.k)
        anchor_index = choose_anchor_window_from_markers(markers, args.k, anchor_bp)
        anchor_counts = counts_by_window[anchor_index]
        anchor_min_count = math.ceil(args.k * args.anchor_min_fraction)

        (
            sample2_anchor_allele,
            selected_state,
            anchor_selected_count,
            anchor_competing_count,
            anchor_margin,
            anchor_assignment_status,
        ) = infer_sample2_anchor_allele(
            anchor_counts,
            args.anchor_allele,
            anchor_min_count,
            args.anchor_min_margin,
        )

        if anchor_assignment_status != "ASSIGNED":
            raise ValueError(
                f"Cannot assign sample2 haplotype at anchor: {anchor_assignment_status}; "
                f"fixed-sample1 candidate counts={anchor_counts}, "
                f"minimum_count={anchor_min_count}, minimum_margin={args.anchor_min_margin}."
            )
        assert sample2_anchor_allele is not None
        assert selected_state is not None

        strict_windows = classify_windows_for_selected_state(
            markers,
            args.k,
            selected_state,
            args.confident_fraction,
            args.possible_fraction,
        )
        windows, switch_events = apply_sample2_switch_tolerance(
            strict_windows,
            anchor_index,
            selected_state,
            args.k,
            args.confident_fraction,
            args.possible_fraction,
            args.switch_min_fraction,
            args.switch_min_margin,
            args.switch_min_windows,
        )
        anchor_window = windows[anchor_index]

        strict_confident_region = contiguous_region(
            windows, anchor_index, {CONFIDENT}
        )
        strict_possible_region = contiguous_region(
            windows, anchor_index, {CONFIDENT, POSSIBLE}
        )
        strict_confident_bp = region_boundary_bp(windows, strict_confident_region)
        strict_possible_bp = region_boundary_bp(windows, strict_possible_region)

        tolerant_proxy = [
            replace(window, classification=window.tolerant_classification)
            for window in windows
        ]
        tolerant_confident_region = contiguous_region(
            tolerant_proxy, anchor_index, {CONFIDENT}
        )
        tolerant_possible_region = contiguous_region(
            tolerant_proxy, anchor_index, {CONFIDENT, POSSIBLE}
        )
        tolerant_confident_bp = region_boundary_bp(
            tolerant_proxy, tolerant_confident_region
        )
        tolerant_possible_bp = region_boundary_bp(
            tolerant_proxy, tolerant_possible_region
        )

        output_dir = os.path.dirname(os.path.abspath(args.output_prefix))
        os.makedirs(output_dir, exist_ok=True)

        summary_path = args.output_prefix + ".summary.tsv"
        windows_path = args.output_prefix + ".windows.tsv"
        markers_path = args.output_prefix + ".markers.tsv"
        switches_path = args.output_prefix + ".switches.tsv"

        write_summary(
            summary_path,
            args.sample1,
            args.sample2,
            metadata_a,
            metadata_b,
            target_chrom,
            anchor_bp,
            anchor_window,
            args.anchor_allele,
            sample2_anchor_allele,
            selected_state,
            anchor_selected_count,
            anchor_competing_count,
            anchor_margin,
            strict_confident_bp,
            strict_possible_bp,
            tolerant_confident_bp,
            tolerant_possible_bp,
            switch_events,
            genetic_map,
            args.k,
            args.confident_fraction,
            args.possible_fraction,
            args.switch_min_fraction,
            args.switch_min_margin,
            args.switch_min_windows,
            args.marker_mode,
            len(markers),
        )
        write_windows(
            windows_path,
            windows,
            args.sample1,
            args.sample2,
            selected_state,
            args.k,
        )
        write_switches(switches_path, switch_events)
        if args.write_markers:
            write_markers(
                markers_path,
                markers,
                args.sample1,
                args.sample2,
                selected_state,
            )

        print(f"Markers used: {len(markers)}", file=sys.stderr)
        print(f"Windows: {len(windows)}", file=sys.stderr)
        print(
            f"Anchor window: {anchor_window.start_bp}-{anchor_window.end_bp} "
            f"(center={anchor_window.center_bp})",
            file=sys.stderr,
        )
        print(
            f"Selected disease-haplotype state: {selected_state} "
            f"({anchor_selected_count} vs {anchor_competing_count}; "
            f"margin={anchor_margin})",
            file=sys.stderr,
        )
        print(
            f"Sample2 status: affected={bool_to_text(metadata_b.affected)}; "
            f"inferred anchor haplotype={sample2_anchor_allele}",
            file=sys.stderr,
        )
        print(f"Wrote: {summary_path}", file=sys.stderr)
        print(f"Accepted sample2 phase-label switches: {len(switch_events)}", file=sys.stderr)
        for event in switch_events:
            print(
                f"  {event.direction}: {event.from_state}->{event.to_state} "
                f"near {event.center_bp} (margin={event.margin})",
                file=sys.stderr,
            )
        print(f"Wrote: {windows_path}", file=sys.stderr)
        print(f"Wrote: {switches_path}", file=sys.stderr)
        if args.write_markers:
            print(f"Wrote: {markers_path}", file=sys.stderr)
        return 0

    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
