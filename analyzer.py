"""
Batch Protein Property Analyzer — Streamlit App
Install: pip install biopython pandas streamlit
Run:     streamlit run analyzer.py
"""

import pandas as pd
import streamlit as st

from protein_pipeline import analyze_sequence, expand_rows_with_mutations
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
        st.caption("First column = sequence, second = mutation region (optional), third = copies to generate.")
    return uploaded


def read_input_df(uploaded):
    try:
        return pd.read_csv(uploaded)
    except Exception as e:
        st.error(f"Could not read file: {e}")
        return None


def resolve_required_columns(df):
    if df.shape[1] < 3:
        st.error("CSV must contain at least 3 columns: sequence, mutation_region, num_random_copies.")
        return None
    seq_col, region_col, copies_col = df.columns[:3]
    st.success(f"Using columns: sequence=`{seq_col}`, mutation region=`{region_col}`, copies=`{copies_col}`")
    return seq_col, region_col, copies_col


def render_mutation_controls():
    st.divider()
    st.subheader("Row-wise random substitution module")
    st.caption("Each row can have its own mutation region and number of random copies.")
    random_seed = st.number_input(
        "Random seed",
        min_value=0,
        max_value=999999,
        value=42,
        step=1,
    )
    run_analysis = st.button("Generate row-wise random copies and analyze")
    return random_seed, run_analysis


def render_summary(results_df, total_count):
    n_ok = results_df["error"].eq("").sum()
    n_err = total_count - n_ok
    ok = results_df[results_df["error"] == ""]

    st.subheader(f"Results — {n_ok} of {total_count} sequences processed")
    if n_err:
        st.warning(f"{n_err} sequence(s) had errors. Check the `error` column in the download.")

    if ok.empty:
        return

    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Sequences", n_ok)
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

    random_seed, run_analysis = render_mutation_controls()
    if not run_analysis:
        st.info("Set the CSV and click the button above to run mutation + analysis.")
        return

    working_df, skipped_zero_copy_rows = expand_rows_with_mutations(
        df=df,
        seq_col=seq_col,
        region_col=region_col,
        copies_col=copies_col,
        random_seed=random_seed,
    )
    if working_df.empty:
        st.error("No sequences to analyze after expansion. Check `num_random_copies` values.")
        return
    if skipped_zero_copy_rows:
        st.warning(f"{skipped_zero_copy_rows} input row(s) were skipped because copy count was 0.")

    with st.spinner(f"Analyzing {len(working_df)} sequences…"):
        results = [analyze_sequence(s) for s in working_df[seq_col].astype(str)]
    results_df = pd.DataFrame(results)
    out_df = pd.concat([working_df.reset_index(drop=True), results_df], axis=1)

    render_summary(results_df, len(working_df))
    render_outputs(out_df)


if __name__ == "__main__":
    main()
