# Spectral Library Generation Using MS2PIP

### 1. Prepare Peptide Lists
Follow the instruction described in *README* to prepare peptide list files.

Set the peptide list directory as working directory and run `scripts/ms2pip/peptide_to_ms2pip_peprec.R`.
```R
setwd("{PATH_TO_DATA}")
source("{PATH_TO_CODE}/scripts/ms2pip/peptide_to_ms2pip_peprec.R")
```

Peptide list files are converted to MS2PIP peptide format and output in the same directory. 
- 130418_03_SynthMix1_Sphospho1_charge2.prediction.peprec.txt
- 130418_03_SynthMix1_Sphospho1_charge3.prediction.peprec.txt
- 130418_03_SynthMix1_Tphospho1_charge2.prediction.peprec.txt
- 130418_03_SynthMix1_Tphospho1_charge3.prediction.peprec.txt
- 130418_03_SynthMix1_Yphospho1_charge2.prediction.peprec.txt
- 130418_03_SynthMix1_Yphospho1_charge3.prediction.peprec.txt


### 2. Predict MS/MS Peak Intensities
Use MS2PIP (https://iomics.ugent.be/ms2pip/) to predicted MS/MS peak intensities from peptide lists. Select CSV as the output format.

Unzip the downloaded archive and rename the prediction result file:
- 130418_03_SynthMix1_Sphospho1_charge2.ms2pip.csv
- 130418_03_SynthMix1_Sphospho1_charge3.ms2pip.csv
- 130418_03_SynthMix1_Tphospho1_charge2.ms2pip.csv
- 130418_03_SynthMix1_Tphospho1_charge3.ms2pip.csv
- 130418_03_SynthMix1_Yphospho1_charge2.ms2pip.csv
- 130418_03_SynthMix1_Yphospho1_charge3.ms2pip.csv


### 3. Generate Spectral Libraries
Set the prediction directory as working directory and run `scripts/ms2pip/ms2pip_csv_to_ions.R`.
```R
setwd("{PATH_TO_DATA}")
source("{PATH_TO_CODE}/scripts/ms2pip/ms2pip_csv_to_ions.R")
```

The MS2PIP CSV are converted to JSON files.
- 130418_03_SynthMix1_Sphospho1_charge2.prediction.ions.json
- 130418_03_SynthMix1_Sphospho1_charge3.prediction.ions.json
- 130418_03_SynthMix1_Tphospho1_charge2.prediction.ions.json
- 130418_03_SynthMix1_Tphospho1_charge3.prediction.ions.json
- 130418_03_SynthMix1_Yphospho1_charge2.prediction.ions.json
- 130418_03_SynthMix1_Yphospho1_charge3.prediction.ions.json


Follow the instruction described in *README* to generate spectral libraries.
