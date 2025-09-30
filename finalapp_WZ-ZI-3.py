import streamlit as st
import pandas as pd
import subprocess, tempfile, os, re

st.set_page_config(page_title="BLAST: Crop vs Allergen (BYO references)", layout="wide")
st.title("BLAST • Crop vs Allergen (User-provided sequences)")

tab_p, tab_x = st.tabs(["BLASTP (protein→protein)", "BLASTX (DNA/RNA→protein)"])

# ---------------- validate AAs and NTs ----------------
AA_RE = re.compile(r'^[ACDEFGHIKLMNPQRSTVWYXBZJUO\-\*\s]+$', re.IGNORECASE)
NT_RE = re.compile(r'^[ACGTURYKMSWBDHVN\-\s]+$', re.IGNORECASE)

def looks_like_fasta(text: str, aa: bool=True) -> bool:
    text = text.strip()
    if not text.startswith(">"):
        return False
    seq_lines = [ln.strip() for ln in text.splitlines() if ln and not ln.startswith(">")]
    if not seq_lines:
        return False
    check = AA_RE if aa else NT_RE
    return any(bool(check.match(ln)) for ln in seq_lines)

# ------------- Temp file helpers -------------
def write_temp_text(text: str, suffix: str):
    fp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix, mode="w")
    fp.write(text.strip() + "\n")
    fp.close()
    return fp.name

def write_temp_upload(uploaded, suffix: str):

    # Read and decode the uploaded file
    raw_text = uploaded.read().decode("utf-8")
    lines = raw_text.strip().splitlines()

    # Fill headers white space with - to retain entire header in output
    edited_lines = []
    for line in lines:
        if line.startswith(">"):
            header = line[1:].strip().replace(' ','-')
            edited_lines.append(f'>{header}')
        else:
            edited_lines.append(line)

    fp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix, encoding="utf-8", mode="w")
    fp.write("\n".join(edited_lines) + "\n")
    fp.close()

    return fp.name

def write_temp_upload_all(uploaded, suffix: str):
    fp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    fp.write(uploaded.read())
    fp.close()
    return fp.name

# ------------- BLAST -----------------
def build_db(faa_path: str, prefix: str):
    subprocess.run(["makeblastdb", "-in", faa_path, "-dbtype", "prot", "-out", prefix], check=True)

def run_blast(tool: str, query_path: str, db_prefix: str, out_tsv: str):
    # tool is 'blastp' or 'blastx'
    subprocess.run([
        tool,
        "-query", query_path,
        "-db", db_prefix,
        "-out", out_tsv,
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
        "-evalue", "1e-10"
    ], check=True)

# ------------- Header parsing ----------------
BRACKETS_SPECIES = re.compile(r'\[([^\[\]]+)\]')  # captures " [Species name] "
PIPE_SPLIT = re.compile(r'\|')
ACC1_FORMAT = re.compile(r'[A-Z]{2}_[0-9]+\.[0-9]')
ACC2_FORMAT = re.compile(r'[A-Z]{3}[0-9]+\.[0-9]')

def parse_id_rich_protein(header: str, default_label: str = ""):
    """
    Try to parse accession, species, and name/gene from diverse FASTA headers.
    Expected best case: Accession|Species|GeneOrDesc
    Fallbacks:
      - If square-bracket species exists at end: use it.
      - Use entire header as desc if nothing else.
    """
    raw = header
    header = header.replace('-',' ')

    # Try pipe-separated first
    ACC1 = ACC1_FORMAT.search(header)
    ACC2 = ACC2_FORMAT.search(header)

    if ACC1:
        acc = ACC1.group(0)
        name_species = header.split(acc)[1].strip()
        name_species_parts = name_species.split('[')
        if len(name_species_parts) > 1:
            name = name_species_parts[0].strip().replace('-',' ')
            species = name_species_parts[1].replace('[','').replace(']','').replace('-',' ').strip()
        else:
            name = "Error"
            species = "Error"
    
    elif ACC2: 
        acc = ACC2.group(0)
        name_species = header.split(acc)[1].strip()
        name_species_parts = name_species.split('[')
        if len(name_species_parts) > 1:
            name = name_species_parts[0].strip().replace('-',' ')
            species = name_species_parts[1].replace('[','').replace(']','').replace('-',' ').strip()
        else:
            name = "Error"
            species = "Error"
    
    else:
        acc = raw
        species = ""
        name = ""

    return acc, species, name

def parse_id_rich_cds(header: str, default_label: str = ""):
    """
    Try to parse accession, species, and name/gene from diverse FASTA headers.
    Expected best case: Accession|Species|GeneOrDesc
    Fallbacks:
      - If square-bracket species exists at end: use it.
      - Use entire header as desc if nothing else.
    """
    raw = header
    header = header.replace('-',' ')

    # Try pipe-separated first
    match_acc = re.search(r'protein_id=([^\]\s]+)', header)
    match_gene = re.search(r'gene=([^\]]+)', header)
    match_name = re.search(r'protein=([^\]]+)', header)
    acc = ""
    gene = ""
    name = ""
    if match_acc:
        acc = match_acc.group(1).strip()
    else:
        acc = raw

    if match_gene:
        gene = match_gene.group(1).strip()
    
    if match_name:
        name = match_name.group(1).strip()
    
    return acc, gene, name

    """
    #elif len(parts) == 2:
        #acc, species = parts[0].strip(), parts[1].strip()
    #elif len(parts) == 1:
        #acc = parts[0].strip()

    # If species is still empty, try bracketed species at end
    if not species:
        m = BRACKETS_SPECIES.search(header)
        if m:
            species = m.group(1).strip().replace(" ", "_")

    # If name is empty, try removing acc/species and take remaining words
    if not name:
        # remove bracketed species
        tmp = BRACKETS_SPECIES.sub("", header)
        # remove accession if it appears at start
        if acc and tmp.startswith(acc):
            tmp = tmp[len(acc):].lstrip(" |:-")
        # remove species if it appears after a pipe
        if species and species in tmp:
            tmp = tmp.replace(species, "")
        name = tmp.strip()
        # collapse spaces -> use underscores for compactness
        name = re.sub(r'\s+', '_', name).strip("_")

    # If everything is still empty, use defaults
    if not acc:
        acc = raw[1:] if raw.startswith(">") else raw
    if default_label and not name:
        name = default_label

    # to be safe avoid commas which can make CSV messy
    name = name.replace(",", " ")
    species = species.replace(",", " ")

    return acc or "", species or (default_label.replace(" ", "_") if default_label else ""), name or (default_label if default_label else "")

    # Try pipe-separated first
    parts = PIPE_SPLIT.split(header)
    if len(parts) >= 3:
        acc, species, name = parts[0].strip(), parts[1].strip(), parts[2].strip()
    elif len(parts) == 2:
        acc, species = parts[0].strip(), parts[1].strip()
    elif len(parts) == 1:
        acc = parts[0].strip()

    # If species is still empty, try bracketed species at end
    if not species:
        m = BRACKETS_SPECIES.search(header)
        if m:
            species = m.group(1).strip().replace(" ", "_")

    # If name is empty, try removing acc/species and take remaining words
    if not name:
        # remove bracketed species
        tmp = BRACKETS_SPECIES.sub("", header)
        # remove accession if it appears at start
        if acc and tmp.startswith(acc):
            tmp = tmp[len(acc):].lstrip(" |:-")
        # remove species if it appears after a pipe
        if species and species in tmp:
            tmp = tmp.replace(species, "")
        name = tmp.strip()
        # collapse spaces -> use underscores for compactness
        name = re.sub(r'\s+', '_', name).strip("_")

    # If everything is still empty, use defaults
    if not acc:
        acc = raw[1:] if raw.startswith(">") else raw
    if default_label and not name:
        name = default_label

    # to be safe avoid commas which can make CSV messy
    name = name.replace(",", " ")
    species = species.replace(",", " ")

    return acc or "", species or (default_label.replace(" ", "_") if default_label else ""), name or (default_label if default_label else "")
    """
# ------------- Results rendering -------------
def render_results(out_path: str, crop_label: str, crop_is_protein: bool, pident_min: int, type: str):
    if not os.path.exists(out_path) or os.path.getsize(out_path) == 0:
        st.info("No hits found.")
        return

    # EXACTLY the 6 columns we asked for above
    cols = ["qseqid","sseqid","pident","length","evalue","bitscore"]
    df = pd.read_csv(out_path, sep="\t", header=None, names=cols, engine="python", quoting=3)

    # Safety check: if the file didn't parse to 6 cols, show the first raw line
    if df.shape[1] != len(cols):
        st.error(f"Unexpected BLAST output shape: got {df.shape[1]} columns, expected {len(cols)}.")
        with open(out_path, "r") as fh:
            first = fh.readline().strip()
        st.write("First raw line from BLAST:")
        st.code(first or "(empty)", language="text")
        return

    # Coerce numerics
    for c in ["pident","length","evalue","bitscore"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Triage labels (optional visual hint)
    df["status"] = "Hits"
    df.loc[(df["pident"]>=80), "status"] = "High"
    df.loc[(df["pident"]>=50) & (df["status"]!="High"), "status"] = "Medium"

    # Parse richer IDs from headers to show crop/allergen name & species
    if type == "protein":
        query_parsed = df["qseqid"].apply(lambda s: parse_id_rich_protein(str(s), default_label=crop_label))
    elif type == "cds":
        query_parsed = df["qseqid"].apply(lambda s: parse_id_rich_cds(str(s), default_label=crop_label))
        
    all_parsed  = df["sseqid"].apply(lambda s: parse_id_rich_protein(str(s), default_label="Allergen"))

    c_acc, c_sp, c_name = zip(*query_parsed)
    a_acc, a_sp, a_gene = zip(*all_parsed)

    if type == "protein":
        df.insert(0, "crop_acc", c_acc)
        df.insert(1, "crop_species", c_sp)
        df.insert(2, "crop_protein", c_name)
        df.insert(3, "allergen_acc", a_acc)
        df.insert(4, "allergen_species", a_sp)
        df.insert(5, "allergen_protein", a_gene)
    elif type == "cds":
        df.insert(0, "crop_acc", c_acc)
        df.insert(1, "crop_gene", c_sp)
        df.insert(2, "crop_protein", c_name)
        df.insert(3, "allergen_acc", a_acc)
        df.insert(4, "allergen_species", a_sp)
        df.insert(5, "allergen_protein", a_gene)

    # E-value as scientific notation string
    df["evalue_sci"] = df["evalue"].apply(lambda x: f"{x:.2e}" if pd.notna(x) else "")

    # Only filter: % identity
    fdf = df[df["pident"] >= pident_min].copy()

    # Sort for readability
    fdf = fdf.sort_values(by=["evalue","pident"], ascending=[True, False])
    
    if type == "protein":
        display_cols = [
            "crop_acc","crop_species","crop_protein",
            "allergen_acc","allergen_species","allergen_protein",
            "pident","length","evalue_sci","bitscore","status"
        ]
    elif type == "cds":
        display_cols = [
            "crop_acc","crop_gene","crop_protein",
            "allergen_acc","allergen_species","allergen_protein",
            "pident","length","evalue_sci","bitscore","status"
            ]

    st.subheader("Results")
    st.write(f"Total hits: {len(df)} | After %identity filter (≥{pident_min}%): {len(fdf)}")
    st.dataframe(fdf[display_cols], use_container_width=True)

    st.download_button(
        "Download filtered results (CSV)",
        fdf[display_cols].to_csv(index=False),
        file_name="blast_results_filtered.csv",
        mime="text/csv"
    )

    st.subheader("Summary by status")
    st.write(fdf["status"].value_counts().rename_axis("status").reset_index(name="count"))

# -------------------- BLASTP --------------------
with tab_p:
    st.header("BLASTP — protein queries vs protein allergens")

    # Controls FIRST (before run)
    c_pident = st.slider("Minimum % identity to display", 0, 100, 50, key="pident_p")
    crop_label = st.text_input("Optional crop label (used if headers lack name/species)", "")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Crop protein FASTA (query)** — paste or upload")
        crop_text = st.text_area("Paste crop protein FASTA", height=180, key="p_crop_text")
        crop_upload = st.file_uploader("...or upload crop protein FASTA", type=["fa","faa","fasta","txt"], key="p_crop_up")
    with c2:
        st.markdown("**Allergen protein FASTA (reference)** — one or many allergens in a single FASTA file")
        all_text = st.text_area("Paste allergen(s) protein FASTA", height=180, key="p_all_text")
        all_upload = st.file_uploader("...or upload allergen(s) protein FASTA", type=["fa","faa","fasta","txt"], key="p_all_up")

    if st.button("Run BLASTP"):
        if not (crop_text.strip() or crop_upload) or not (all_text.strip() or all_upload):
            st.error("Provide BOTH crop protein FASTA and allergen protein FASTA (paste or upload).")
            st.stop()
        # validate if file was uploaded or text was pasted in input box
        if crop_text and not looks_like_fasta(crop_text, aa=True) and not crop_upload:
            st.error("Crop FASTA doesn't look like protein FASTA.")
            st.stop()
        if all_text and not looks_like_fasta(all_text, aa=True) and not all_upload:
            st.error("Allergen FASTA doesn't look like protein FASTA.")
            st.stop()

        try:
            q_faa = write_temp_upload(crop_upload, ".faa") if crop_upload else write_temp_text(crop_text, ".faa")
            a_faa = write_temp_upload(all_upload, ".faa") if all_upload else write_temp_text(all_text, ".faa")
            db_prefix = tempfile.NamedTemporaryFile(delete=False).name

            with st.spinner("Building allergen BLAST database..."):
                build_db(a_faa, db_prefix)

            out_tsv = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv").name
            with st.spinner("Running BLASTP..."):
                run_blast("blastp", q_faa, db_prefix, out_tsv)

            render_results(out_tsv, crop_label=crop_label, crop_is_protein=True, pident_min=c_pident, type="protein")

        except subprocess.CalledProcessError as e:
            st.error(f"BLAST error: {e}")
        finally:
            # cleanup temp files and DB
            for pth in [locals().get("q_faa"), locals().get("a_faa"),
                        db_prefix + ".pin" if 'db_prefix' in locals() else None,
                        db_prefix + ".psq" if 'db_prefix' in locals() else None,
                        db_prefix + ".phr" if 'db_prefix' in locals() else None,
                        locals().get("out_tsv")]:
                if pth and os.path.exists(pth):
                    try: os.remove(pth)
                    except: pass

# -------------------- BLASTX --------------------
with tab_x:
    st.header("BLASTX — nucleotide queries vs protein allergens")

    # Controls FIRST (before run)
    x_pident = st.slider("Minimum % identity to display", 0, 100, 50, key="pident_x")
    crop_label_x = st.text_input("Optional crop label (used if headers lack name/species) ", "", key="crop_label_x")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Crop DNA/RNA FASTA (query)** — paste or upload")
        crop_text = st.text_area("Paste crop nucleotide FASTA", height=180, key="x_crop_text")
        crop_upload = st.file_uploader("...or upload crop nucleotide FASTA", type=["fa","fasta","fna","ffn","txt"], key="x_crop_up")
    with c2:
        st.markdown("**Allergen protein FASTA (reference)** — one or many allergens in a single FASTA file")
        all_text = st.text_area("Paste allergen(s) protein FASTA", height=180, key="x_all_text")
        all_upload = st.file_uploader("...or upload allergen(s) protein FASTA", type=["fa","faa","fasta","txt"], key="x_all_up")

    if st.button("Run BLASTX"):
        if not (crop_text.strip() or crop_upload) or not (all_text.strip() or all_upload):
            st.error("Provide BOTH crop nucleotide FASTA and allergen protein FASTA (paste or upload).")
            st.stop()
        # validate
        if crop_text and not looks_like_fasta(crop_text, aa=False) and not crop_upload:
            st.error("Crop FASTA doesn't look like nucleotide FASTA.")
            st.stop()
        if all_text and not looks_like_fasta(all_text, aa=True) and not all_upload:
            st.error("Allergen FASTA doesn't look like protein FASTA.")
            st.stop()

        try:
            q_fna = write_temp_upload(crop_upload, ".fna") if crop_upload else write_temp_text(crop_text, ".fna")
            a_faa = write_temp_upload(all_upload, ".faa") if all_upload else write_temp_text(all_text, ".faa")
            db_prefix = tempfile.NamedTemporaryFile(delete=False).name

            with st.spinner("Building allergen BLAST database..."):
                build_db(a_faa, db_prefix)

            out_tsv = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv").name
            with st.spinner("Running BLASTX..."):
                run_blast("blastx", q_fna, db_prefix, out_tsv)

            render_results(out_tsv, crop_label=crop_label_x, crop_is_protein=False, pident_min=x_pident, type="cds")

        except subprocess.CalledProcessError as e:
            st.error(f"BLAST error: {e}")
        finally:
            for pth in [locals().get("q_fna"), locals().get("a_faa"),
                        db_prefix + ".pin" if 'db_prefix' in locals() else None,
                        db_prefix + ".psq" if 'db_prefix' in locals() else None,
                        db_prefix + ".phr" if 'db_prefix' in locals() else None,
                        locals().get("out_tsv")]:
                if pth and os.path.exists(pth):
                    try: os.remove(pth)
                    except: pass

st.markdown("---")
st.caption("Tip: Use headers like `>Accession|Species|GeneOrDescription`. If your query headers lack name/species, fill the optional crop label above.")
