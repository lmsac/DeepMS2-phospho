# Spectral Library Generation Using Prosit

### 1. Prepare Peptide Lists
Follow the instruction described in *README* to prepare peptide list files.

Set the peptide list directory as working directory and run `scripts/prosit/peptide_to_prosit_csv.R`.
```R
setwd("{PATH_TO_DATA}")
source("{PATH_TO_CODE}/scripts/prosit/peptide_to_prosit_csv.R")
```

Collision energy (CE) parameter can be set in the script.

Peptide list files are converted to Prosit peptide CSV format and output in the same directory. 
- 130418_03_SynthMix1_Sphospho1_charge2_ce{CE}.prosit_peptide.csv
- 130418_03_SynthMix1_Sphospho1_charge3_ce{CE}.prosit_peptide.csv
- 130418_03_SynthMix1_Tphospho1_charge2_ce{CE}.prosit_peptide.csv
- 130418_03_SynthMix1_Tphospho1_charge3_ce{CE}.prosit_peptide.csv
- 130418_03_SynthMix1_Yphospho1_charge2_ce{CE}.prosit_peptide.csv
- 130418_03_SynthMix1_Yphospho1_charge3_ce{CE}.prosit_peptide.csv


### 2. Predict MS/MS Peak Intensities
Use Prosit (https://www.proteomicsdb.org/prosit/) to predicted MS/MS peak intensities from peptide lists. Select Spectronaut CSV as the output format.

Unzip the downloaded archive and rename the prediction result file:
- 130418_03_SynthMix1_Sphospho1_charge2_ce{CE}.prosit.csv
- 130418_03_SynthMix1_Sphospho1_charge3_ce{CE}.prosit.csv
- 130418_03_SynthMix1_Tphospho1_charge2_ce{CE}.prosit.csv
- 130418_03_SynthMix1_Tphospho1_charge3_ce{CE}.prosit.csv
- 130418_03_SynthMix1_Yphospho1_charge2_ce{CE}.prosit.csv
- 130418_03_SynthMix1_Yphospho1_charge3_ce{CE}.prosit.csv


### 3. Generate Spectral Libraries
Set the prediction directory as working directory and run `scripts/prosit/prosit_to_ions.R`.
```R
setwd("{PATH_TO_DATA}")
source("{PATH_TO_CODE}/scripts/prosit/prosit_to_ions.R")
```

The MS2PIP CSV are converted to JSON files.
- 130418_03_SynthMix1_Sphospho1_charge2_ce{CE}.prosit.ions.json
- 130418_03_SynthMix1_Sphospho1_charge3_ce{CE}.prosit.ions.json
- 130418_03_SynthMix1_Tphospho1_charge2_ce{CE}.prosit.ions.json
- 130418_03_SynthMix1_Tphospho1_charge3_ce{CE}.prosit.ions.json
- 130418_03_SynthMix1_Yphospho1_charge2_ce{CE}.prosit.ions.json
- 130418_03_SynthMix1_Yphospho1_charge3_ce{CE}.prosit.ions.json

Follow the instruction described in *README* to generate spectral libraries.
