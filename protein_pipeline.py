import random
from itertools import product

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

PKA = {"D": 3.9, "E": 4.1, "H": 6.0, "C": 8.3, "Y": 10.1, "K": 10.5, "R": 12.5, "Nterm": 8.0, "Cterm": 3.1}
VALID_AA = "ACDEFGHIKLMNPQRSTVWY"

# ── Conservative substitution groups ─────────────────────────────────────────
CONSERVATIVE_GROUPS = {
    "G": ["A"],
    "A": ["G", "V", "S"],
    "V": ["A", "I", "L"],
    "L": ["V", "I", "M"],
    "I": ["V", "L", "M"],
    "M": ["L", "I"],
    "F": ["Y", "W"],
    "Y": ["F", "W", "H"],
    "W": ["F", "Y"],
    "S": ["T", "A", "N"],
    "T": ["S", "V", "N"],
    "N": ["Q", "S", "D"],
    "Q": ["N", "E", "K"],
    "K": ["R", "Q"],
    "R": ["K", "H"],
    "H": ["R", "K", "Y"],
    "D": ["E", "N"],
    "E": ["D", "Q"],
    "C": ["S", "A"],
    "P": ["A", "G"],
}


def net_charge(seq, ph=7.0):
    c = 1.0 / (1.0 + 10 ** (ph - PKA["Nterm"]))
    c -= 1.0 / (1.0 + 10 ** (PKA["Cterm"] - ph))
    for aa in ("D", "E", "C", "Y"):
        c -= seq.count(aa) / (1.0 + 10 ** (PKA[aa] - ph))
    for aa in ("H", "K", "R"):
        c += seq.count(aa) / (1.0 + 10 ** (ph - PKA[aa]))
    return round(c, 3)


def analyze_sequence(seq):
    seq = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    bad = set(seq) - set(VALID_AA)
    if bad:
        return {"error": f"invalid chars: {bad}"}
    if len(seq) < 2:
        return {"error": "too short"}
    try:
        pa = ProteinAnalysis(seq)
        helix, turn, sheet = pa.secondary_structure_fraction()
        return {
            "length": len(seq),
            "molecular_weight": round(pa.molecular_weight(), 2),
            "isoelectric_point": round(pa.isoelectric_point(), 2),
            "gravy": round(pa.gravy(), 4),
            "aromaticity": round(pa.aromaticity(), 4),
            "instability_index": round(pa.instability_index(), 2),
            "stable": pa.instability_index() <= 40,
            "charge_at_pH7": net_charge(seq),
            "positive_res_RK": seq.count("R") + seq.count("K"),
            "negative_res_DE": seq.count("D") + seq.count("E"),
            "helix_fraction": round(helix, 4),
            "sheet_fraction": round(sheet, 4),
            "turn_fraction": round(turn, 4),
            "error": "",
        }
    except Exception as e:
        return {"error": str(e)}


def parse_mutation_regions(text):
    regions = []
    if not str(text).strip():
        return regions
    parts = [p.strip() for p in str(text).split(",") if p.strip()]
    for p in parts:
        if "-" in p:
            a, b = [x.strip() for x in p.split("-", 1)]
            if not a.isdigit() or not b.isdigit():
                raise ValueError(f"Invalid region: {p}")
            start, end = int(a), int(b)
        else:
            if not p.isdigit():
                raise ValueError(f"Invalid position: {p}")
            start = end = int(p)
        if start < 1 or end < 1:
            raise ValueError(f"Positions must be >= 1: {p}")
        if end < start:
            raise ValueError(f"Region end must be >= start: {p}")
        regions.append((start, end))
    return regions


def parse_copy_count(value):
    if pd.isna(value) or str(value).strip() == "":
        return 1
    s = str(value).strip()
    try:
        num = float(s)
    except ValueError as e:
        raise ValueError(f"Invalid copy count: {s}") from e
    if not num.is_integer():
        raise ValueError(f"Copy count must be an integer: {s}")
    n = int(num)
    if n < 0:
        raise ValueError(f"Copy count must be >= 0: {s}")
    return n


# ── Combo count helpers ───────────────────────────────────────────────────────
def max_unique_conservative_variants(seq, regions):
    """
    Returns the maximum number of unique sequences producible by conservative
    substitution across the given regions (product of substitute pool sizes,
    NOT including the original residue — since we always mutate away from it).
    Used to detect when num_copies > unique possibilities.
    Returns (max_unique, per_position) where per_position is a list of
    (pos, original_aa, n_substitutes).
    """
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    total = 1
    per_position = []
    seen = set()
    for start, end in regions:
        for pos in range(start, end + 1):
            if pos in seen:
                continue
            seen.add(pos)
            idx = pos - 1
            if idx < 0 or idx >= len(clean):
                continue
            aa = clean[idx]
            subs = CONSERVATIVE_GROUPS.get(aa, [])
            n = len(subs) if subs else 1  # if no subs, position is fixed → factor of 1
            per_position.append((pos, aa, n))
            total *= n
    return total, per_position


def count_conservative_combos(seq, regions):
    """
    Returns (total_combos, per_position_breakdown).
    total_combos = product of (n_substitutes + 1) across all positions in region,
    where +1 counts keeping the original residue as an option.
    per_position_breakdown = list of (position, original_aa, substitutes, n_choices).
    """
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    total = 1
    per_position = []
    seen = set()
    for start, end in regions:
        for pos in range(start, end + 1):
            if pos in seen:
                continue
            seen.add(pos)
            idx = pos - 1
            if idx < 0 or idx >= len(clean):
                continue
            aa = clean[idx]
            subs = CONSERVATIVE_GROUPS.get(aa, [])
            n_choices = len(subs) + 1  # include original
            per_position.append((pos, aa, subs, n_choices))
            total *= n_choices
    return total, per_position


def count_random_combos(seq, regions):
    """
    Returns (total_combos, n_positions).
    Each position can be any of 20 AAs, so total = 20^n_positions.
    """
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    seen = set()
    n_positions = 0
    for start, end in regions:
        for pos in range(start, end + 1):
            if pos in seen:
                continue
            seen.add(pos)
            idx = pos - 1
            if 0 <= idx < len(clean):
                n_positions += 1
    return 20 ** n_positions, n_positions


def count_weighted_combos(seq, regions, weighted_map):
    """
    Returns (total_combos, per_position_breakdown) for weighted mode.
    At each valid position:
      choices = original residue + all weighted substitutes for that residue.
    """
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    total = 1
    per_position = []
    seen = set()
    for start, end in regions:
        for pos in range(start, end + 1):
            if pos in seen:
                continue
            seen.add(pos)
            idx = pos - 1
            if idx < 0 or idx >= len(clean):
                continue
            aa = clean[idx]
            n_choices = len(weighted_map.get(aa, {}))
            if aa not in weighted_map.get(aa, {}):
                n_choices += 1
            per_position.append((pos, aa, n_choices))
            total *= n_choices
    return total, per_position


def parse_weighted_substitutions_df(weights_df):
    """
    Parse substitution weights CSV into:
      {FROM_AA: {TO_AA: probability, ...}, ...}

    Supported formats:
      1) Long format: columns from,to,weight (or source,target,probability etc)
      2) Matrix format: first column = source AA, remaining AA columns = weights
    """
    if weights_df is None or weights_df.empty:
        raise ValueError("weights CSV is empty")

    col_map = {str(c).strip().lower(): c for c in weights_df.columns}
    long_from_keys = ["from", "from_aa", "source", "src", "original"]
    long_to_keys = ["to", "to_aa", "target", "dst", "new"]
    long_w_keys = ["weight", "prob", "probability", "p"]

    from_col = next((col_map[k] for k in long_from_keys if k in col_map), None)
    to_col = next((col_map[k] for k in long_to_keys if k in col_map), None)
    w_col = next((col_map[k] for k in long_w_keys if k in col_map), None)

    weighted_map = {}
    if from_col and to_col and w_col:
        for _, row in weights_df.iterrows():
            src = str(row[from_col]).strip().upper()
            dst = str(row[to_col]).strip().upper()
            if not src or src == "NAN" or not dst or dst == "NAN":
                continue
            if src not in VALID_AA or dst not in VALID_AA:
                raise ValueError(f"Invalid amino acid in weights row: {src}->{dst}")
            try:
                w = float(row[w_col])
            except Exception as e:
                raise ValueError(f"Invalid weight for {src}->{dst}: {row[w_col]}") from e
            if w < 0:
                raise ValueError(f"Negative weight for {src}->{dst}")
            if w == 0:
                continue
            weighted_map.setdefault(src, {})
            weighted_map[src][dst] = weighted_map[src].get(dst, 0.0) + w
    else:
        # Matrix format: first column is source AA, remaining columns are to-AA.
        src_col = weights_df.columns[0]
        for _, row in weights_df.iterrows():
            src = str(row[src_col]).strip().upper()
            if not src or src == "NAN":
                continue
            if src not in VALID_AA:
                raise ValueError(f"Invalid source amino acid in matrix: {src}")
            for col in weights_df.columns[1:]:
                dst = str(col).strip().upper()
                if dst not in VALID_AA:
                    continue
                cell = row[col]
                if pd.isna(cell) or str(cell).strip() == "":
                    continue
                try:
                    w = float(cell)
                except Exception as e:
                    raise ValueError(f"Invalid matrix weight {src}->{dst}: {cell}") from e
                if w < 0:
                    raise ValueError(f"Negative matrix weight for {src}->{dst}")
                if w == 0:
                    continue
                weighted_map.setdefault(src, {})
                weighted_map[src][dst] = weighted_map[src].get(dst, 0.0) + w

    if not weighted_map:
        raise ValueError("No valid substitution weights found in CSV")
    return weighted_map


# ── Mode 1: random substitution ───────────────────────────────────────────────
def mutate_sequence_random(seq, regions):
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    bad = set(clean) - set(VALID_AA)
    if bad:
        raise ValueError(f"invalid chars: {bad}")
    if not clean:
        raise ValueError("empty sequence")
    arr = list(clean)
    changed_positions = []
    for start, end in regions:
        for pos in range(start, end + 1):
            idx = pos - 1
            if idx < 0 or idx >= len(arr):
                continue
            current = arr[idx]
            choices = [aa for aa in VALID_AA if aa != current]
            arr[idx] = random.choice(choices)
            changed_positions.append(pos)
    return "".join(arr), len(changed_positions)


# ── Mode 2: conservative substitution ────────────────────────────────────────
def mutate_sequence_conservative(seq, regions):
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    bad = set(clean) - set(VALID_AA)
    if bad:
        raise ValueError(f"invalid chars: {bad}")
    if not clean:
        raise ValueError("empty sequence")
    arr = list(clean)
    changed_positions = []
    skipped_positions = []
    for start, end in regions:
        for pos in range(start, end + 1):
            idx = pos - 1
            if idx < 0 or idx >= len(arr):
                continue
            current = arr[idx]
            choices = CONSERVATIVE_GROUPS.get(current, [])
            if choices:
                arr[idx] = random.choice(choices)
                changed_positions.append(pos)
            else:
                skipped_positions.append(pos)
    return "".join(arr), len(changed_positions), skipped_positions


def mutate_sequence_weighted(seq, regions, weighted_map):
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    bad = set(clean) - set(VALID_AA)
    if bad:
        raise ValueError(f"invalid chars: {bad}")
    if not clean:
        raise ValueError("empty sequence")

    arr = list(clean)
    changed_positions = []
    skipped_positions = []
    for start, end in regions:
        for pos in range(start, end + 1):
            idx = pos - 1
            if idx < 0 or idx >= len(arr):
                continue
            current = arr[idx]
            subs = weighted_map.get(current, {})
            if not subs:
                skipped_positions.append(pos)
                continue
            choices = list(subs.keys())
            weights = list(subs.values())
            chosen = random.choices(choices, weights=weights, k=1)[0]
            arr[idx] = chosen
            if chosen != current:
                changed_positions.append(pos)
    return "".join(arr), len(changed_positions), skipped_positions


# ── Mode 3: exhaustive single-position scanning ───────────────────────────────
def scan_sequence_single_position(seq, regions):
    """
    For each position in regions, generate one mutant per conservative substitute,
    mutating only that single position at a time.
    Returns list of {mutated_seq, mutated_position, original_aa, new_aa}.
    """
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    bad = set(clean) - set(VALID_AA)
    if bad:
        raise ValueError(f"invalid chars: {bad}")
    if not clean:
        raise ValueError("empty sequence")

    variants = []
    seen = set()
    for start, end in regions:
        for pos in range(start, end + 1):
            if pos in seen:
                continue
            seen.add(pos)
            idx = pos - 1
            if idx < 0 or idx >= len(clean):
                continue
            original_aa = clean[idx]
            for new_aa in CONSERVATIVE_GROUPS.get(original_aa, []):
                arr = list(clean)
                arr[idx] = new_aa
                variants.append({
                    "mutated_seq": "".join(arr),
                    "mutated_position": pos,
                    "original_aa": original_aa,
                    "new_aa": new_aa,
                })
    return variants


# ── Shared expansion logic ────────────────────────────────────────────────────
def expand_rows_with_mutations(df, seq_col, region_col, copies_col,
                               random_seed, mutation_mode="random", weighted_map=None):
    """
    mutation_mode: "random" | "conservative" | "scan" | "weighted"
    In scan mode, copies_col is ignored — all single-position variants are
    generated exhaustively.

    Returns (expanded_df, skipped_zero_copy_rows, duplicate_warnings)
    where duplicate_warnings is a list of dicts:
        {row: int, requested: int, max_unique: int, per_position: list}
    """
    random.seed(int(random_seed))
    expanded_rows = []
    skipped_zero_copy_rows = 0
    duplicate_warnings = []

    for idx, row in df.iterrows():
        raw_seq = row[seq_col]
        raw_region = row[region_col]
        raw_copies = row[copies_col]

        region_text = "" if pd.isna(raw_region) else str(raw_region).strip()
        row_error = ""

        try:
            regions = parse_mutation_regions(region_text)
        except ValueError as e:
            regions = []
            row_error = str(e)

        # ── Scan mode ─────────────────────────────────────────────────────────
        if mutation_mode == "scan":
            base_row = row.copy()
            base_row["source_row"] = idx + 1
            base_row["original_sequence"] = str(raw_seq)
            base_row["mutation_regions"] = region_text
            base_row["mutation_mode"] = mutation_mode

            if row_error:
                base_row["mutation_error"] = row_error
                base_row["mutated_position"] = ""
                base_row["original_aa"] = ""
                base_row["new_aa"] = ""
                expanded_rows.append(base_row)
                continue

            try:
                variants = scan_sequence_single_position(raw_seq, regions)
            except ValueError as e:
                base_row["mutation_error"] = str(e)
                base_row["mutated_position"] = ""
                base_row["original_aa"] = ""
                base_row["new_aa"] = ""
                expanded_rows.append(base_row)
                continue

            if not variants:
                base_row["mutation_error"] = "no conservative substitutes found in region"
                base_row["mutated_position"] = ""
                base_row["original_aa"] = ""
                base_row["new_aa"] = ""
                expanded_rows.append(base_row)
                continue

            for v in variants:
                new_row = base_row.copy()
                new_row[seq_col] = v["mutated_seq"]
                new_row["mutated_position"] = v["mutated_position"]
                new_row["original_aa"] = v["original_aa"]
                new_row["new_aa"] = v["new_aa"]
                new_row["mutation_error"] = ""
                expanded_rows.append(new_row)
            continue

        # ── Random / conservative / weighted modes ────────────────────────────
        try:
            num_copies = parse_copy_count(raw_copies)
        except ValueError as e:
            num_copies = 1
            row_error = str(e) if not row_error else f"{row_error}; {e}"

        if num_copies == 0:
            skipped_zero_copy_rows += 1
            continue

        # ── Duplicate check for conservative/random/weighted modes ────────────
        if not row_error and regions:
            try:
                if mutation_mode == "conservative":
                    max_unique, per_pos = max_unique_conservative_variants(str(raw_seq), regions)
                    if num_copies > max_unique:
                        duplicate_warnings.append({
                            "row": idx + 1,
                            "requested": num_copies,
                            "max_unique": max_unique,
                            "mode": "conservative",
                            "per_position": per_pos,
                        })
                elif mutation_mode == "random":
                    max_unique, n_positions = count_random_combos(str(raw_seq), regions)
                    if num_copies > max_unique:
                        duplicate_warnings.append({
                            "row": idx + 1,
                            "requested": num_copies,
                            "max_unique": max_unique,
                            "mode": "random",
                            "per_position": [(None, None, 19)] * n_positions,
                        })
                elif mutation_mode == "weighted" and weighted_map is not None:
                    max_unique, per_pos = count_weighted_combos(str(raw_seq), regions, weighted_map)
                    if num_copies > max_unique:
                        duplicate_warnings.append({
                            "row": idx + 1,
                            "requested": num_copies,
                            "max_unique": max_unique,
                            "mode": "weighted",
                            "per_position": per_pos,
                        })
            except Exception:
                pass

        for copy_i in range(1, num_copies + 1):
            new_row = row.copy()
            new_row["source_row"] = idx + 1
            new_row["copy_index"] = copy_i
            new_row["requested_copies"] = num_copies
            new_row["original_sequence"] = str(raw_seq)
            new_row["mutation_regions"] = region_text
            new_row["mutation_mode"] = mutation_mode

            if row_error:
                new_row[seq_col] = str(raw_seq)
                new_row["mutated_positions_count"] = 0
                new_row["conservative_skipped"] = ""
                new_row["mutation_error"] = row_error
            else:
                try:
                    if mutation_mode == "conservative":
                        mutated_seq, n_changed, skipped = mutate_sequence_conservative(raw_seq, regions)
                        new_row["conservative_skipped"] = ",".join(map(str, skipped)) if skipped else ""
                    elif mutation_mode == "weighted":
                        if weighted_map is None:
                            raise ValueError("weighted substitution map is missing")
                        mutated_seq, n_changed, skipped = mutate_sequence_weighted(raw_seq, regions, weighted_map)
                        new_row["conservative_skipped"] = ",".join(map(str, skipped)) if skipped else ""
                    else:
                        mutated_seq, n_changed = mutate_sequence_random(raw_seq, regions)
                        new_row["conservative_skipped"] = ""
                    new_row[seq_col] = mutated_seq
                    new_row["mutated_positions_count"] = n_changed
                    new_row["mutation_error"] = ""
                except ValueError as e:
                    new_row[seq_col] = str(raw_seq)
                    new_row["mutated_positions_count"] = 0
                    new_row["conservative_skipped"] = ""
                    new_row["mutation_error"] = str(e)

            expanded_rows.append(new_row)

    return pd.DataFrame(expanded_rows), skipped_zero_copy_rows, duplicate_warnings
