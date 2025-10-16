# in your main code, add st.markdown("[View Simplified Results](simple_results)") at the end, then create a new python file called simple_results and add this:

import streamlit as st
import pandas as pd
import plotly.express as px

st.title("Simplified Results Viewer")

uploaded = st.file_uploader("Upload your filtered results (CSV)", type="csv")

if uploaded:
    df = pd.read_csv(uploaded)

    # --- Rename columns for clarity (if they exist) ---
    rename_dict = {
        "crop_acc": "Crop Protein Accession Number",#both
        "crop_gene": "Crop Gene",#x
        "crop_protein": "Crop Protein",#both
        "allergen_acc": "Allergen Accession Number",#both
        "allergen_species": "Allergen Species",#both
        "allergen_protein": "Allergen Protein",#both
        "pident": "Percent Identical",#both
        "length": "Length",#both
        "evalue": "E-Value Score",#both
        "bitscore": "Bit Score",#both
        "status": "Status",#both
    }

    # Only rename existing columns
    df.rename(columns={k: v for k, v in rename_dict.items() if k in df.columns}, inplace=True)

    # --- Add filters ---
    st.subheader("Filter Results")
    allergen_species_options = df["Allergen Species"].dropna().unique() if "Allergen Species" in df.columns else []
    allergen_protein_options = df["Allergen Accession Number"].dropna().unique() if "Allergen Accession Number" in df.columns else []

    selected_species = st.selectbox("Select Allergen Species (optional):", ["All"] + list(allergen_species_options))
    selected_allergen = st.selectbox("Select Allergen (optional):", ["All"] + list(allergen_protein_options))

    filtered_df = df.copy()
    if selected_species != "All":
        filtered_df = filtered_df[filtered_df["Allergen Species"] == selected_species]
    if selected_allergen != "All":
        filtered_df = filtered_df[filtered_df["Allergen Accession Number"] == selected_allergen]

    # --- Summary by status ---
    if "Status" in filtered_df.columns:
        st.subheader("Simplified Summary")
        summary = filtered_df.groupby("Status").size().reset_index(name="Count")
        st.table(summary)

    # --- Top 5 hits ---
    st.subheader("Top 5 BLAST Hits")
    sort_col = None
    for col in ["Percent Identical", "Bit Score", "Score"]:
        if col in filtered_df.columns:
            sort_col = col
            break

    if sort_col:
        top_hits = filtered_df.sort_values(by=sort_col, ascending=False).head(5)
    else:
        top_hits = filtered_df.head(5)
    
    def acc_link(acc):
        return f'<a href="https://www.ncbi.nlm.nih.gov/protein/{acc}" target="_blank" rel="noopener noreferrer">{acc}</a>' if acc else ""

    top_hits = top_hits.copy()
    top_hits["Crop Protein Accession Number"] = top_hits["Crop Protein Accession Number"].apply(acc_link)
    top_hits["Allergen Accession Number"] = top_hits["Allergen Accession Number"].apply(acc_link)
    
    st.markdown(
        f"""
        <div style="display: flex; justify-content: center;">
            {top_hits.to_html(escape=False, index=False)}
        </div>
        """,
        unsafe_allow_html=True
    )
    
    # --- Proteins with high similarity ---
    if "Status" in filtered_df.columns and "High" in filtered_df["Status"].values:
        st.subheader("Proteins Detected with High Similarity to Allergens")
        cols_to_show = [
            "Crop Protein Accession Number", "Crop Gene", "Crop Protein",
            "Allergen Accession Number", "Allergen Species", "Allergen Protein",
            "Percent Identical", "Bit Score", "Status"
        ]
        available_cols = [c for c in cols_to_show if c in filtered_df.columns]
        st.dataframe(filtered_df[filtered_df["Status"] == "High"][available_cols])

    # --- Bubble Chart Explanation ---
    st.subheader("BLAST Hit Quality Overview")
    st.markdown("""
    **How to read this chart:**
    - Points **farther to the right** have **higher sequence similarity** (more identical amino acids or bases).
    - Points **higher up** have **stronger alignments** (higher bit scores mean better BLAST matches).
    - **Larger bubbles** represent **longer alignments** — meaning more of the sequence was compared.
    - **Colors** indicate different crop species.
    """)

    # --- Bubble Chart ---
    required_cols = {"Percent Identical", "Bit Score"}
    if required_cols.issubset(filtered_df.columns):
        if "Length" not in filtered_df.columns:
            filtered_df["Length"] = 1  # fallback if missing
        if "Crop Accession Number" not in filtered_df.columns:
            filtered_df["Crop Protein Accession Number"] = "Unknown"

        fig = px.scatter(
            filtered_df,
            x="Percent Identical",
            y="Bit Score",
            size="Length",
            color="Allergen Accession Number",
            hover_data=["Allergen Species", "Crop Protein", "Allergen Protein"],
            title="Allergen Match Strength by Crop Species",
            size_max=30,
        )
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("Bubble chart unavailable — missing required columns (Percent Identical and Bit Score).")

else:
    st.info("Please upload a CSV file to begin.")
