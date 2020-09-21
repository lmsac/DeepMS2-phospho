# Spectral Library Generation Using pDeep2

### 1. Prepare Peptide Lists
Follow the instruction described in *README* to prepare peptide list files.

Set the peptide list directory as working directory and run `scripts/pdeep2/peptide_to_pdeep2_txt.R`.
```R
setwd("{PATH_TO_DATA}")
source("{PATH_TO_CODE}/scripts/pdeep2/peptide_to_pdeep2_txt.R")
```

Peptide list files are converted to pDeep2 peptide format and output in the same directory. 
- 130418_03_SynthMix1_Sphospho1_charge2.peptide.txt
- 130418_03_SynthMix1_Sphospho1_charge3.peptide.txt
- 130418_03_SynthMix1_Tphospho1_charge2.peptide.txt
- 130418_03_SynthMix1_Tphospho1_charge3.peptide.txt
- 130418_03_SynthMix1_Yphospho1_charge2.peptide.txt
- 130418_03_SynthMix1_Yphospho1_charge3.peptide.txt


### 2. Predict MS/MS Peak Intensities
Use pDeep2 (https://github.com/pFindStudio/pDeep) to predicted MS/MS peak intensities from peptide lists.

Rename the prediction result file:
- 130418_03_SynthMix1_Sphospho1_charge2.pdeep2.mgf
- 130418_03_SynthMix1_Sphospho1_charge3.pdeep2.mgf
- 130418_03_SynthMix1_Tphospho1_charge2.pdeep2.mgf
- 130418_03_SynthMix1_Tphospho1_charge3.pdeep2.mgf
- 130418_03_SynthMix1_Yphospho1_charge2.pdeep2.mgf
- 130418_03_SynthMix1_Yphospho1_charge3.pdeep2.mgf


### 3. Generate Spectral Libraries
Set the prediction directory as working directory and run `scripts/pdeep2/pdeep2_mgf_to_ions.R`.
```R
setwd("{PATH_TO_DATA}")
source("{PATH_TO_CODE}/scripts/pdeep2/pdeep2_mgf_to_ions.R")
```

The MS2PIP CSV are converted to JSON files.
- 130418_03_SynthMix1_Sphospho1_charge2.pdeep2.ions.json
- 130418_03_SynthMix1_Sphospho1_charge3.pdeep2.ions.json
- 130418_03_SynthMix1_Tphospho1_charge2.pdeep2.ions.json
- 130418_03_SynthMix1_Tphospho1_charge3.pdeep2.ions.json
- 130418_03_SynthMix1_Yphospho1_charge2.pdeep2.ions.json
- 130418_03_SynthMix1_Yphospho1_charge3.pdeep2.ions.json

Follow the instruction described in *README* to generate spectral libraries.

If the used pDeep2 model was trained with phosphopeptide data and neutral loss prediction is enabled, there is no need to perform neutral loss optimization when generating spectral libraries.
