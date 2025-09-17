# streamlit is used for the user interface and final visual
# pandas is for handling tables
# subprocess used to call outside blast tools
# tempfile save for inputs and outputs
# os and shutil to check files
import streamlit as st
import pandas as pd
import subprocess, tempfile, os, shutil

#set tab and title and page width and show the main title 
st.set_page_config(page_title="Allergen BLASTP", layout="wide")
st.title("Crop Protein vs Allergen Similarity")

# shows the instructions for this app
st.markdown("""
Paste a **crop protein sequence** in FASTA format *or* upload a `.faa/.fa/.fasta` file.
This app will compare your protein(s) against the **allergen database**  built from `data/allergen_proteins.faa`.
""")

# database location for the blast to go find
DEFAULT_DB_PREFIX = "data/allergen_blastdb" 

# quick helper
def run_blastp(query_faa: str, db_prefix: str, out_tsv: str):
    cmd = [
        "blastp",
        "-query", query_faa,
        "-db", db_prefix,
        "-out", out_tsv,
        "-outfmt", "6 qseqid aseqid pident length evalue bitscore qlen alen"
    ]
    subprocess.run(cmd, check=True)

# Input widgets
col1, col2 = st.columns(2)
with col1:
    seq_text = st.text_area("Paste FASTA (one or more protein sequences)", height=200, placeholder=">MyProtein\nMSEQUENCE...")
with col2:
    file_up = st.file_uploader("...or upload a FASTA file", type=["fa","faa","fasta","txt"])

db_prefix = DEFAULT_DB_PREFIX
st.caption(f"Using allergen BLAST database prefix: `{db_prefix}`")

# Threshold controls
st.subheader("Result filters (adjustable)")
pident_min = st.slider("Minimum % identity", 0, 100, 60)
alen_min   = st.slider("Minimum alignment length (aa)", 0, 1000, 100, step=10)
evalue_max = st.text_input("Maximum e-value (scientific notation ok)", "1e-10")

if st.button("Run BLASTP"):
    # Prep query FASTA
    if not file_up and not seq_text.strip():
        st.error("Please paste a FASTA or upload a file.")
        st.stop()

    # Write query to a temp file
    qpath = None
    try:
        if file_up:
            suffix = ".faa"
            qpath = tempfile.NamedTemporaryFile(delete=False, suffix=suffix).name
            with open(qpath, "wb") as f:
                f.write(file_up.read())
        else:
            if not seq_text.strip().startswith(">"):
                st.error("Pasted input must be valid FASTA (start with '>').")
                st.stop()
            qpath = tempfile.NamedTemporaryFile(delete=False, suffix=".faa", mode="w")
            qpath.write(seq_text.strip())
            qpath.close()
            qpath = qpath.name

        # Temp output
        out_tsv = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv").name

        # Ensure DB files exist
        missing = [ext for ext in (".pin", ".psq", ".phr") if not os.path.exists(db_prefix + ext)]
        if missing:
            st.error(
                "BLAST database files not found. Please build them first:\n"
                f"`makeblastdb -in data/allergen_proteins.faa -dbtype prot -out {db_prefix}`"
            )
            st.stop()

        # Run BLASTP
        try:
            run_blastp(qpath, db_prefix, out_tsv)
        except subprocess.CalledProcessError as e:
            st.error(f"BLASTP failed: {e}")
            st.stop()

        # Load results
        cols = ["qseqid","sseqid","pident","length","evalue","bitscore","qlen","slen"]
        if os.path.getsize(out_tsv) == 0:
            st.info("No hits at all.")
        else:
            df = pd.read_csv(out_tsv, sep="\t", header=None, names=cols)
            # Coerce types
            df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
            df["length"] = pd.to_numeric(df["length"], errors="coerce")
            df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")

            # Label
            df["status"] = "Screen"
            df.loc[(df["pident"]>=80) & (df["length"]>=150) & (df["evalue"]<=1e-20), "status"] = "High"
            df.loc[(df["pident"]>=60) & (df["length"]>=100) & (df["evalue"]<=1e-10) & (df["status"]!="High"), "status"] = "Medium"

            # Apply filters from UI
            try:
                emax = float(evalue_max)
            except ValueError:
                st.warning("Invalid e-value; using 1e-10 default.")
                emax = 1e-10

            mask = (
                (df["pident"] >= pident_min) &
                (df["length"] >= alen_min) &
                (df["evalue"] <= emax)
            )
            filt = df[mask].copy()

            st.subheader("Results")
            st.write(f"Total hits: {len(df)} | After filters: {len(filt)}")
            st.dataframe(filt, use_container_width=True)

            st.download_button(
                "Download filtered results (CSV)",
                filt.to_csv(index=False),
                file_name="blastp_results_filtered.csv",
                mime="text/csv"
            )

            # Quick summary by status
            st.subheader("Summary by status")
            st.write(filt["status"].value_counts().rename_axis("status").reset_index(name="count"))

    finally:
        # Cleanup temp files
        for p in (qpath, locals().get("out_tsv", None)):
            if p and os.path.exists(p):
                try: os.remove(p)
                except: pass

st.markdown("---")
st.caption("Tip: Keep informative FASTA headers so results are easy to read.")

