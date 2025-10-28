# Welcome
AllergenScan is a bioinformatics tool designed to identify potential allergens from protein or coding sequence data using BLASTp or BLASTx. The tool includes an application providing a user-interface to conduct BLASTp or BLASTx searches against a custom allergenic database. Results can be exported in csv format and uploaded into the results viewer for a more detailed analysis.


## Set-up
To clone this repository, verify `git` is installed on your machine and run:

`git clone https://github.com/zirelan/BIOT670I_Group-3.git`

### Required tools
It is recommended to use a version of conda to install needed tools in an isolated environment. Required tools and packages are:
  - python3.12
  - blast
  - streamlit
  - pandas
  - plotly

If you choose to use conda, navigate to the cloned repository directory `BIOT670I_Group-3` and run the following to create an environment for the application:

`conda env create -f allergen-scan-env.yml`

Once the environment is created successfully, you now have all tools required to run the detection and visualization application.


## Using the Application
### Launch the application
To start the application, navigate to the cloned repository directory `BIOT670I_Group-3` and run the following (if needed: activate conda environment):

`conda activate allergen-scan-env`

`streamlit run allergen-scan.py`


### Navigation
The following icon in the top left can be selected to navigate between the detection and visualization pages.

<img width="38" height="38" alt="image" src="https://github.com/user-attachments/assets/1f9c7275-0706-439b-9d87-ca2f6dfb1ecc" />


### Required Files
Both BLASTp and BLASTx require an allergenic database to search query sequences against. A FASTA formatted file containing allergenic protein sequences is required and can be used for both methods.

For BLASTp, users can paste FASTA formatted sequence(s) to run against an allergenic database or use the file uploader to import a FASTA formatted file containing protein sequences.

For BLASTx, users can paste FASTA formatted sequence(s) to run against an allergenic database or use the file uploader to import a FASTA formatted file containing coding sequences.

The results viewer requires a csv file produced by the detection application. This should be downloaded and uploaded on the results viewer page.

### Detailed User Guide
Navigate to the [application guide](README_application-guide.md) for a detailed walkthrough of how to use this application.

## Appendix

### Crop Data Collection Tools
Visit the following for instructions on how to install NCBI Genome datasets using the command line. If browsing the NCBI Genome database, you will need to identify cultivars with annotations available to source protein and coding sequence data.

[Data Collection](README_data-collection.md)


### Crop protein and cds FASTAs
Visit the following dropbox to download peanut, soybean, or sesame query/crop FASTAs

[Crop Protein and CDS FASTAs](https://www.dropbox.com/scl/fo/bwcl29pyo41evatbbkpgm/ALsCv_aGQWCRlqEq1wBAut8?rlkey=sclcidm77ryiwx2y14eh2v1n1&st=31ti2znf&dl=0)


### Protein database
Visit the following dropbox to download organsim specific allergenic protein FASTAs or compiled peanut, soybean, wheat, and sesame allergens.

[Protein Databases](https://www.dropbox.com/scl/fo/xylt6z1774zg0zedh2ns0/AE4oXtg07LBp8Ac35QVxt7Y?rlkey=qu3g93llnur4yg3um20r6shua&st=7oqbgv1y&dl=0)

## Authors
Zach Irelan | Wahaj Zuberi | Amelia Bradford | Laurie Orell | Mustafa Talha Sigindere
