"""
Batch Protein Property Analyzer — Streamlit App
Install: pip install biopython pandas streamlit
Run:     streamlit run streamlit_protein_analyzer.py
"""

import io
import pandas as pd
import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Protein Analyzer",
    page_icon="🧬",
    layout="wide",
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Mono:wght@400;500&family=DM+Sans:wght@300;400;500&display=swap');

html, body, [class*="css"] { font-family: 'DM Sans', sans-serif; }
code, .stCode, pre { font-family: 'DM Mono', monospace !important; }

/* Top bar */
.block-container { padding-top: 2rem; padding-bottom: 2rem; }

/* Metric cards */
[data-testid="metric-container"] {
    background: #f8f7f4;
    border: 1px solid #e8e5de;
    border-radius: 10px;
    padding: 1rem;
}
[data-testid="stMetricLabel"] { font-size: 0.72rem; color: #888; letter-spacing: 0.06em; text-transform: uppercase; }
[data-testid="stMetricValue"] { font-size: 1.5rem; font-weight: 500; }

/* Table */
[data-testid="stDataFrame"] { border: 1px solid #e8e5de; border-radius: 10px; overflow: hidden; }

/* Upload box */
[data-testid="stFileUploadDropzone"] {
    border: 1.5px dashed #c8c5be !important;
    border-radius: 10px !important;
    background: #fafaf8 !important;
}

/* Buttons */
.stDownloadButton > button {
    background: #1a1a1a !important;
    color: #fff !important;
    border: none !important;
    border-radius: 8px !important;
    font-size: 0.85rem !important;
    padding: 0.5rem 1.2rem !important;
    font-family: 'DM Sans', sans-serif !important;
}
.stDownloadButton > button:hover { opacity: 0.85 !important; }

/* Status badges */
.badge-stable   { background:#e6f4ec; color:#1a7a3e; padding:2px 10px; border-radius:20px; font-size:0.78rem; }
.badge-unstable { background:#fdecea; color:#b91c1c; padding:2px 10px; border-radius:20px; font-size:0.78rem; }

h1 { font-weight: 500 !important; letter-spacing: -0.02em; }
h3 { font-weight: 500 !important; color: #444; }
</style>
""", unsafe_allow_html=True)

# ── Helpers ───────────────────────────────────────────────────────────────────
PKA = {"D":3.9,"E":4.1,"H":6.0,"C":8.3,"Y":10.1,"K":10.5,"R":12.5,"Nterm":8.0,"Cterm":3.1}

def net_charge(seq, ph=7.0):
    c  = 1.0 / (1.0 + 10**(ph - PKA["Nterm"]))
    c -= 1.0 / (1.0 + 10**(PKA["Cterm"] - ph))
    for aa in ("D","E","C","Y"): c -= seq.count(aa) / (1.0 + 10**(PKA[aa] - ph))
    for aa in ("H","K","R"):     c += seq.count(aa) / (1.0 + 10**(ph - PKA[aa]))
    return round(c, 3)

def analyze_sequence(seq):
    seq = str(seq).strip().upper().replace(" ","").replace("\n","")
    bad = set(seq) - set("ACDEFGHIKLMNPQRSTVWY")
    if bad:   return {"error": f"invalid chars: {bad}"}
    if len(seq) < 2: return {"error": "too short"}
    try:
        pa = ProteinAnalysis(seq)
        helix, turn, sheet = pa.secondary_structure_fraction()
        return {
            "length":            len(seq),
            "molecular_weight":  round(pa.molecular_weight(), 2),
            "isoelectric_point": round(pa.isoelectric_point(), 2),
            "gravy":             round(pa.gravy(), 4),
            "aromaticity":       round(pa.aromaticity(), 4),
            "instability_index": round(pa.instability_index(), 2),
            "stable":            pa.instability_index() <= 40,
            "charge_at_pH7":     net_charge(seq),
            "positive_res_RK":   seq.count("R") + seq.count("K"),
            "negative_res_DE":   seq.count("D") + seq.count("E"),
            "helix_fraction":    round(helix, 4),
            "sheet_fraction":    round(sheet, 4),
            "turn_fraction":     round(turn, 4),
            "error":             "",
        }
    except Exception as e:
        return {"error": str(e)}

def find_seq_column(df):
    candidates = [c for c in df.columns
                  if c.strip().lower() in
                  ("sequence","seq","protein","aa_sequence",
                   "protein_sequence","amino_acid_sequence")]
    return candidates[0] if candidates else None

# ── UI ────────────────────────────────────────────────────────────────────────
st.title("🧬 Protein property analyzer")
st.caption("Upload a CSV with an amino acid sequence column — get physicochemical properties back instantly.")

st.divider()

col_up, col_hint = st.columns([2, 1])

with col_up:
    uploaded = st.file_uploader("Upload CSV", type=["csv"], label_visibility="collapsed")

with col_hint:
    st.markdown("**Expected format**")
    st.code("sequence\nMPYEKHVEQ...\nRQGGGAPAG...", language="text")
    st.caption("One sequence per row. Column can be named `sequence`, `seq`, `protein`, etc.")

if uploaded:
    try:
        df = pd.read_csv(uploaded)
    except Exception as e:
        st.error(f"Could not read file: {e}")
        st.stop()

    # Auto-detect or let user pick column
    seq_col = find_seq_column(df)
    if not seq_col:
        seq_col = st.selectbox(
            "Could not auto-detect sequence column. Please select it:",
            options=df.columns.tolist()
        )
    else:
        st.success(f"Detected sequence column: **{seq_col}**")

    st.divider()

    # ── Run analysis ──────────────────────────────────────────────────────────
    with st.spinner(f"Analyzing {len(df)} sequences…"):
        results = [analyze_sequence(s) for s in df[seq_col].astype(str)]
    results_df = pd.DataFrame(results)
    out_df = pd.concat([df.reset_index(drop=True), results_df], axis=1)

    n_ok  = results_df["error"].eq("").sum()
    n_err = len(df) - n_ok
    ok    = results_df[results_df["error"] == ""]

    # ── Summary metrics ───────────────────────────────────────────────────────
    st.subheader(f"Results — {n_ok} of {len(df)} sequences processed")
    if n_err:
        st.warning(f"{n_err} sequence(s) had errors. Check the `error` column in the download.")

    if not ok.empty:
        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("Sequences",        n_ok)
        m2.metric("Avg MW",           f"{ok['molecular_weight'].mean()/1000:.1f} kDa")
        m3.metric("Avg GRAVY",        f"{ok['gravy'].mean():.3f}")
        m4.metric("Avg pI",           f"{ok['isoelectric_point'].mean():.2f}")
        m5.metric("Stable",           f"{ok['stable'].sum()} / {n_ok}")

        st.divider()

        # ── Per-property distributions ────────────────────────────────────────
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

    # ── Full results table ────────────────────────────────────────────────────
    st.subheader("Full results table")
    display_df = out_df.copy()
    if "stable" in display_df.columns:
        display_df["stable"] = display_df["stable"].map({True: "✓ stable", False: "✗ unstable"})
    st.dataframe(display_df, use_container_width=True, height=400)

    # ── Download ──────────────────────────────────────────────────────────────
    csv_bytes = out_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="⬇ Download results CSV",
        data=csv_bytes,
        file_name="protein_results.csv",
        mime="text/csv",
    )

else:
    st.info("Upload a CSV file above to get started.")