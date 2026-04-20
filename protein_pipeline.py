import random

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

PKA = {"D": 3.9, "E": 4.1, "H": 6.0, "C": 8.3, "Y": 10.1, "K": 10.5, "R": 12.5, "Nterm": 8.0, "Cterm": 3.1}
VALID_AA = "ACDEFGHIKLMNPQRSTVWY"


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


def mutate_sequence(seq, regions):
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


def expand_rows_with_mutations(df, seq_col, region_col, copies_col, random_seed):
    random.seed(int(random_seed))
    expanded_rows = []
    skipped_zero_copy_rows = 0

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

        try:
            num_copies = parse_copy_count(raw_copies)
        except ValueError as e:
            num_copies = 1
            row_error = str(e) if not row_error else f"{row_error}; {e}"

        if num_copies == 0:
            skipped_zero_copy_rows += 1
            continue

        for copy_i in range(1, num_copies + 1):
            new_row = row.copy()
            new_row["source_row"] = idx + 1
            new_row["copy_index"] = copy_i
            new_row["requested_copies"] = num_copies
            new_row["original_sequence"] = str(raw_seq)
            new_row["mutation_regions"] = region_text

            if row_error:
                new_row[seq_col] = str(raw_seq)
                new_row["mutated_positions_count"] = 0
                new_row["mutation_error"] = row_error
            else:
                try:
                    mutated_seq, n_changed = mutate_sequence(raw_seq, regions)
                    new_row[seq_col] = mutated_seq
                    new_row["mutated_positions_count"] = n_changed
                    new_row["mutation_error"] = ""
                except ValueError as e:
                    new_row[seq_col] = str(raw_seq)
                    new_row["mutated_positions_count"] = 0
                    new_row["mutation_error"] = str(e)
            expanded_rows.append(new_row)

    return pd.DataFrame(expanded_rows), skipped_zero_copy_rows
