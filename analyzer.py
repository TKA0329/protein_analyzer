"""
Batch Protein Property Analyzer — Streamlit App
Install: pip install biopython pandas streamlit
Run:     streamlit run analyzer.py
"""

import pandas as pd
import streamlit as st

from protein_pipeline import (
    CONSERVATIVE_GROUPS,
    analyze_sequence,
    count_conservative_combos,
    count_random_combos,
    count_weighted_combos,
    expand_rows_with_mutations,
    max_unique_conservative_variants,
    parse_weighted_substitutions_df,
    parse_mutation_regions,
)
from ui_constants import APP_CSS, EXPECTED_FORMAT_EXAMPLE


def render_header():
    st.title("🧬 Protein property analyzer")
    st.caption("Upload a CSV with sequence, mutation region, and copy count columns.")
    st.divider()


def render_upload_panel():
    col_up, col_hint = st.columns([2, 1])
    with col_up:
        uploaded = st.file_uploader("Upload CSV", type=["csv"], label_visibility="collapsed")
    with col_hint:
        st.markdown("**Expected format**")
        st.code(EXPECTED_FORMAT_EXAMPLE, language="text")
        st.caption(
            "First column = sequence, second = mutation region (e.g. `5-7` or `3,6,10`), "
            "third = number of copies (ignored in scan mode)."
        )
    return uploaded


def read_input_df(uploaded):
    try:
        return pd.read_csv(uploaded)
    except Exception as e:
        st.error(f"Could not read file: {e}")
        return None


def resolve_required_columns(df):
    if df.shape[1] < 3:
        st.error("CSV must contain at least 3 columns: sequence, mutation_region, num_copies.")
        return None
    seq_col, region_col, copies_col = df.columns[:3]
    st.success(
        f"Using columns: sequence=`{seq_col}`, "
        f"mutation region=`{region_col}`, copies=`{copies_col}`"
    )
    return seq_col, region_col, copies_col


def render_combo_counter(df, seq_col, region_col, mutation_mode, weighted_map=None):
    """Show the number of possible combinations for each row."""
    st.markdown("#### Possible combinations in your dataset")
    rows_info = []
    for i, (_, row) in enumerate(df.iterrows()):
        seq = str(row[seq_col]).strip().upper().replace(" ", "")
        region_text = "" if pd.isna(row[region_col]) else str(row[region_col]).strip()
        try:
            regions = parse_mutation_regions(region_text)
            if mutation_mode in ("conservative", "scan"):
                total, breakdown = count_conservative_combos(seq, regions)
                per_pos = ", ".join(
                    f"pos {p}: {aa}→[{','.join(subs)}] ({n} choices)"
                    for p, aa, subs, n in breakdown
                )
            elif mutation_mode == "weighted":
                if weighted_map is None:
                    total = "n/a"
                    per_pos = "Load weights CSV to compute weighted combinations"
                else:
                    total, breakdown = count_weighted_combos(seq, regions, weighted_map)
                    per_pos = ", ".join(
                        f"pos {p}: {aa} ({n} weighted choices incl original)"
                        for p, aa, n in breakdown
                    )
            else:
                total, n_pos = count_random_combos(seq, regions)
                per_pos = f"{n_pos} positions × 20 AAs each"
            rows_info.append({
                "Row": i + 1,
                "Sequence (truncated)": seq[:30] + ("…" if len(seq) > 30 else ""),
                "Region": region_text,
                "Total combinations": f"{total:,}" if isinstance(total, int) else str(total),
                "Breakdown": per_pos,
            })
        except Exception as e:
            rows_info.append({
                "Row": i + 1,
                "Sequence (truncated)": seq[:30],
                "Region": region_text,
                "Total combinations": "error",
                "Breakdown": str(e),
            })

    st.dataframe(pd.DataFrame(rows_info), use_container_width=True, hide_index=True)
    st.caption(
        "Conservative/scan combos count each position's substitute pool size + 1 (keeping original). "
        "Random combos = 20^n_positions. Weighted combos use your uploaded weight table."
    )


def render_mutation_controls():
    st.divider()
    st.subheader("Mutation settings")

    mode = st.radio(
        "Substitution mode",
        options=["Random", "Conservative", "Scan (single-position exhaustive)", "Weighted (semi-random)"],
        horizontal=True,
        help=(
            "**Random**: each position mutates to any other AA at random.\n\n"
            "**Conservative**: each position swaps to a physicochemically similar AA "
            "(polar→polar, nonpolar→nonpolar, etc.). Number of copies set per row.\n\n"
            "**Scan**: for each position in the region, generate every possible "
            "conservative substitute — one mutation at a time. "
            "Ignores the copies column. Best for identifying which residue matters.\n\n"
            "**Weighted (semi-random)**: substitutions are sampled from your uploaded "
            "weights/probabilities table."
        ),
    )

    if "Scan" in mode:
        mode_key = "scan"
    elif "Weighted" in mode:
        mode_key = "weighted"
    else:
        mode_key = mode.lower()

    if mode_key in ("conservative", "scan"):
        with st.expander("View conservative substitution groups"):
            rows = [
                {"Residue": aa, "Conservative substitutes": ", ".join(subs)}
                for aa, subs in sorted(CONSERVATIVE_GROUPS.items())
            ]
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    if mode_key == "scan":
        st.info(
            "Scan mode generates all single-position conservative variants. "
            "Each output row differs from the original by exactly one residue. "
            "The `copies` column in your CSV is ignored."
        )

    weights_uploaded = None
    if mode_key == "weighted":
        st.markdown("#### Weighted substitution table")
        st.caption(
            "Upload a CSV in either format:\n"
            "1) long: `from,to,weight`\n"
            "2) matrix: first column source AA, remaining AA columns as weights."
        )
        weights_uploaded = st.file_uploader(
            "Upload weighted substitutions CSV",
            type=["csv"],
            key="weighted_subs_csv",
        )

    random_seed = st.number_input(
        "Random seed (used in random/conservative/weighted modes)",
        min_value=0, max_value=999999, value=42, step=1,
        disabled=(mode_key == "scan"),
    )

    run_analysis = st.button("Generate variants and analyze")
    return mode_key, random_seed, run_analysis, weights_uploaded


def render_summary(results_df, total_count, mutation_mode):
    n_ok = results_df["error"].eq("").sum()
    n_err = total_count - n_ok
    ok = results_df[results_df["error"] == ""]

    st.subheader(f"Results — {n_ok} of {total_count} sequences processed")
    if n_err:
        st.warning(f"{n_err} sequence(s) had errors. Check the `error` column in the download.")

    if ok.empty:
        return

    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Variants", n_ok)
    m2.metric("Avg MW", f"{ok['molecular_weight'].mean() / 1000:.1f} kDa")
    m3.metric("Avg GRAVY", f"{ok['gravy'].mean():.3f}")
    m4.metric("Avg pI", f"{ok['isoelectric_point'].mean():.2f}")
    m5.metric("Stable", f"{ok['stable'].sum()} / {n_ok}")

    st.divider()
    st.subheader("Distributions")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**GRAVY score**")
        st.bar_chart(ok["gravy"], height=160, use_container_width=True)
        st.markdown("**Isoelectric point (pI)**")
        st.bar_chart(ok["isoelectric_point"], height=160, use_container_width=True)
    with c2:
        st.markdown("**Molecular weight (Da)**")
        st.bar_chart(ok["molecular_weight"], height=160, use_container_width=True)
        st.markdown("**Instability index**")
        st.bar_chart(ok["instability_index"], height=160, use_container_width=True)
    st.divider()


def render_outputs(out_df):
    st.subheader("Full results table")
    display_df = out_df.copy()
    if "stable" in display_df.columns:
        display_df["stable"] = display_df["stable"].map({True: "✓ stable", False: "✗ unstable"})
    st.dataframe(display_df, use_container_width=True, height=400)

    csv_bytes = out_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="⬇ Download results CSV",
        data=csv_bytes,
        file_name="protein_results.csv",
        mime="text/csv",
    )


def main():
    st.set_page_config(page_title="Protein Analyzer", page_icon="🧬", layout="wide")
    st.markdown(APP_CSS, unsafe_allow_html=True)
    render_header()

    uploaded = render_upload_panel()
    if not uploaded:
        st.info("Upload a CSV file above to get started.")
        return

    df = read_input_df(uploaded)
    if df is None:
        return

    cols = resolve_required_columns(df)
    if cols is None:
        return
    seq_col, region_col, copies_col = cols

    mutation_mode, random_seed, run_analysis, weights_uploaded = render_mutation_controls()

    weighted_map = None
    if mutation_mode == "weighted" and weights_uploaded:
        try:
            weights_df = pd.read_csv(weights_uploaded)
            weighted_map = parse_weighted_substitutions_df(weights_df)
        except Exception:
            weighted_map = None

    # ── Combo counter (always shown once file is loaded) ──────────────────────
    with st.expander("📊 How many combinations are possible for your sequences?", expanded=False):
        render_combo_counter(df, seq_col, region_col, mutation_mode, weighted_map=weighted_map)

    if not run_analysis:
        st.info("Configure settings above and click **Generate variants and analyze**.")
        return

    if mutation_mode == "weighted":
        if not weights_uploaded:
            st.error("Weighted mode requires a second CSV file with substitution weights.")
            return
        if weighted_map is None:
            try:
                weights_df = pd.read_csv(weights_uploaded)
                weighted_map = parse_weighted_substitutions_df(weights_df)
            except Exception as e:
                st.error(f"Could not parse weighted substitution CSV: {e}")
                return
        st.success("Loaded weighted substitution table.")

    working_df, skipped_zero_copy_rows, duplicate_warnings = expand_rows_with_mutations(
        df=df,
        seq_col=seq_col,
        region_col=region_col,
        copies_col=copies_col,
        random_seed=random_seed,
        mutation_mode=mutation_mode,
        weighted_map=weighted_map,
    )

    if working_df.empty:
        st.error("No sequences to analyze after expansion. Check your region and copies values.")
        return
    if skipped_zero_copy_rows:
        st.warning(f"{skipped_zero_copy_rows} input row(s) skipped (copy count was 0).")

    # ── Duplicate warnings ────────────────────────────────────────────────────
    if duplicate_warnings:
        for w in duplicate_warnings:
            if w["mode"] == "random":
                n_pos = len(w["per_position"])
                st.warning(
                    f"⚠️ Row {w['row']}: you requested **{w['requested']} copies** but only "
                    f"**{w['max_unique']:,} unique** random variants exist "
                    f"({n_pos} position{'s' if n_pos != 1 else ''} × 19 possible substitutes each). "
                    f"The extra copies will be duplicates."
                )
            elif w["mode"] == "weighted":
                pos_detail = ", ".join(
                    f"pos {p} ({aa}: {n} weighted choices incl original)"
                    for p, aa, n in w["per_position"]
                )
                st.warning(
                    f"⚠️ Row {w['row']}: you requested **{w['requested']} copies** but only "
                    f"**{w['max_unique']} unique** weighted variants exist for this region "
                    f"({pos_detail}). "
                    f"The extra copies will be duplicates."
                )
            else:
                pos_detail = ", ".join(
                    f"pos {p} ({aa}: {n} substitute{'s' if n != 1 else ''})"
                    for p, aa, n in w["per_position"]
                )
                st.warning(
                    f"⚠️ Row {w['row']}: you requested **{w['requested']} copies** but only "
                    f"**{w['max_unique']} unique** conservative variants exist for this region "
                    f"({pos_detail}). "
                    f"The extra copies will be duplicates."
                )

    with st.spinner(f"Analyzing {len(working_df)} sequences…"):
        results = [analyze_sequence(s) for s in working_df[seq_col].astype(str)]

    results_df = pd.DataFrame(results)
    out_df = pd.concat([working_df.reset_index(drop=True), results_df], axis=1)

    render_summary(results_df, len(working_df), mutation_mode)
    render_outputs(out_df)


if __name__ == "__main__":
    main()
