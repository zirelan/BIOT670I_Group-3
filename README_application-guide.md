# Application user guide

## BLASTp and BLASTx - Allergen Detection

After launching the application as described on the [main](README.md) welcome page, the interface below will be loaded. Both the BLASTp and BLASTx interfaces are functionally the same.

### Tools - BLASTp and BLASTx

Toggle between BLASTp and BLASTx by selecting the tool underneath the page title

<img width="1885" height="194" alt="image" src="https://github.com/user-attachments/assets/a768baee-dc9e-46e2-92cb-26797ebe8206" />

### Search Criteria

The following minimums can be set to filter out results. For high-probability allergen hits, the following criteria are recommended:
- Percent Sequence Identity = 70
- Alignment Length = 100
- e-Value = 1e-10

<img width="1885" height="250" alt="image" src="https://github.com/user-attachments/assets/254a70a6-a4b7-46bd-b20f-69f40380b888" />

### Sequences

The crop header can be used if a sequence lacks detail in it's current header.

Query or allergen sequences can be pasted or uploaded in FASTA format.

BLASTp query and allergen headers and BLASTx allergen headers should have the following format:

`ACCESSION NAME [ORGANISM]`

`NP_001291328.1 11S globulin subunit beta like precursor [Sesamum indicum]`

BLASTx query headers should have the following tags:

`protein_id=` `gene=` `protein=`

If these formats are not detected, the sequence identifier will be all characters until the first whitespace of the header.

<img width="1885" height="400" alt="image" src="https://github.com/user-attachments/assets/831adffd-353f-4d6f-901a-d612c2605985" />

### Running

Once all the above is completed, select the run button and wait for your results!

<img width="1885" height="814" alt="image" src="https://github.com/user-attachments/assets/829f2c31-3826-4372-b4c5-f7a332aae47c" />

### Results

Once the BLAST search is complete, hits will be recorded and displayed in a scrollable table. If accessions are in the NCBI database, links will populate to direct to the NCBI record when selected.

Hits are classified as High, Medium, and Low confidence matches.

The results can be exported as a csv and imported into the results viewer.

<img width="1885" height="917" alt="image" src="https://github.com/user-attachments/assets/d010379a-1594-44ad-84a2-48bbbad004dd" />

## Results Viewer

### Upload results file

Navigate to the results viewer as described in [main](README.md) and upload the csv export of your blast search.

<img width="1885" height="272" alt="image" src="https://github.com/user-attachments/assets/bb8a98ed-97a5-4f54-b4ce-ebe91403cdde" />

### Filtering and Top Results

The first section of the results viewer is focused on high-impact results and displays a pie chart to visualize the allergenic organism to each query sequence hit and the top 5 blast hits by percent sequence identity. These sections allow for quick visual inspection and understanding of which allergen organism the query is most similar to. If accessions are present in the NCBI database, you can navigate to their page by selecting the accession in the top 5 hits table.

<img width="1885" height="916" alt="image" src="https://github.com/user-attachments/assets/a44c3f03-ec69-41ea-a8ee-f889342736ed" />

### High Classified Results Table

All hits with a Status classification of High will be displayed in a scrollable table

<img width="1885" height="283" alt="image" src="https://github.com/user-attachments/assets/334599f7-3c33-4213-8274-223cf5e34b7b" />

### BLAST Hits Scoring Visualization

BLAST hits within the bubble chart are plotted against bit score and percent identity. As points move up and right (having higher bit score and percent identity), the hit becomes more significant. Clusters of hits in the top right strongly suggest sequence similarity to an allergenic protein. The legend on the right has selectable allergen accessions to filter out/in query sequence matches to that allergen. The counts table to the right depicts how many matches each allergenic accession had to the query sequences.

<img width="1394" height="692" alt="image" src="https://github.com/user-attachments/assets/9f2d74f5-760b-41f9-bc41-c0f01d76c949" />
