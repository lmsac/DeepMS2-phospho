# DeepMS2-phospho
Using deep learning to generate in silico spectral libraries for site localization of phosphopeptides. 

## System Requirements
- Python >= 3.5. Anaconda is recommended.
- Keras with TensorFlow backend.
- R. As an alternative, the latest version of Microsoft R Open should be fine.
- RStudio is recommended but optional.

DeepMS2-phospho has been tested on a workstation with Intel Xeon E5-2690 v3 CPU, 16 GB RAM, and Microsoft Windows Server 2016 Version 1607 (OS Build 14393.2430) 64-bit operating system.

## Installation

### 1. Install Python (Anaconda)
Download Anaconda Installer form https://www.anaconda.com/distribution/.

DeepMS2-phospho has been tested using Anaconda 4.2.0 (Python 3.5.2).

### 2. Install TensorFlow and Keras
Install TensorFlow using `pip`:
```bash
pip install --upgrade tensorflow
```
For GPU-supported version,
```bash
pip install --upgrade tensorflow-gpu
```
TensorFlow documentation is available at https://www.tensorflow.org/.

Install Keras:
```bash
pip install keras
```
Keras documentation is available at https://keras.io/.

DeepMS2-phospho has been tested using Keras 2.2.4 and TensorFlow 1.11.

### 3. Install R and RStudio
R is available at https://www.r-project.org/. As an alternative, Microsoft R Open is available at https://mran.microsoft.com/open/. 

RStudio is available at https://www.rstudio.com/.

DeepMS2-phospho has been tested using Microsoft R Open 3.5.1 and RStudio 1.1.447.

Start R, ensure packages `readr` and `rjson` have been installed.
```R
install.packages("readr")
install.packages("rjson")
```

## Input Data Formats

### 1. Prepare a PSM List
A peptide-spectrum match (PSM) list is stored in a comma-separated values (CSV) file including the following columns: `file`, `scan`, `charge` and `peptide`.

An example file can be found at [`data/example/130418_03_SynthMix1.psm.csv`](data/example/130418_03_SynthMix1.psm.csv).
```csv
"file","scan","charge","peptide"
"130418_03_SynthMix1_1pmol_10MS2_20min_NoDynExc",1708,2,"GTPSQS*PVVGR"
"130418_03_SynthMix1_1pmol_10MS2_20min_NoDynExc",1709,2,"ESLKEEDES*DDDNM@"
"130418_03_SynthMix1_1pmol_10MS2_20min_NoDynExc",1720,2,"ESLKEEDES*DDDNM@"
"130418_03_SynthMix1_1pmol_10MS2_20min_NoDynExc",1727,2,"GTPSQS*PVVGR"
"130418_03_SynthMix1_1pmol_10MS2_20min_NoDynExc",1734,2,"ESLKEEDES*DDDNM@"
...
```

In the `peptide` column, `*` represents phosphorylation (STY), and `@` represents oxidation (M). Other modifications are not supported.

Users can generate PSM list files from database search tools or write them using Excel manually.

### 2. Prepare MS/MS Spectra
DeepMS2-phospho supports Mascot Generic Format (MGF) files as MS/MS spectra format. The file names of MGF files should match with those in the `file` column of the PSM list file (without extension).

Users can convert raw LC-MS/MS data to MGF files using MSConvert from ProteoWizard (https://github.com/ProteoWizard/pwiz).

An example MGF file can be found at [`data/example/130418_03_SynthMix1_1pmol_10MS2_20min_NoDynExc.mgf`](data/example/130418_03_SynthMix1_1pmol_10MS2_20min_NoDynExc.mgf).


## Spectral Library Generation

### 1. Prepare Peptide Lists
A peptide list is stored in a comma-separated values (CSV) file including the following columns: `sequence` and `modification`.  
```csv
"sequence","modification"
"GTPSQSPVVGR","S4(Phospho)"
"GTPSQSPVVGR","S6(Phospho)"
"ESLKEEDESDDDNM","S2(Phospho)"
"ESLKEEDESDDDNM","S9(Phospho)"
"SPSPYYSR","S1(Phospho)"
```

Peptide list files can be generated from a PSM list file.

Start R and run `init.R` to load the functions.
```R
source("{PATH_TO_CODE}/init.R")
```

Set the PSM list directory as working directory and run `scripts/generate_phosphopeptides.R`.
```R
setwd("{PATH_TO_DATA}")
source("{PATH_TO_CODE}/scripts/generate_phosphopeptides.R")
```

Three peptide list files are output in the same directory. 
- 130418_03_SynthMix1_Sphospho1.peptide.csv
- 130418_03_SynthMix1_Tphospho1.peptide.csv
- 130418_03_SynthMix1_Yphospho1.peptide.csv


### 2. Predict MS/MS Peak Intensities
Use DeepDIA to predicted MS/MS peak intensities from peptide lists.

For detailed manual, see https://github.com/lmsac/DeepDIA.

DeepDIA outputs the predicted JSON files.
- 130418_03_SynthMix1_Sphospho1_charge2.prediction.ions.json
- 130418_03_SynthMix1_Sphospho1_charge3.prediction.ions.json
- 130418_03_SynthMix1_Tphospho1_charge2.prediction.ions.json
- 130418_03_SynthMix1_Tphospho1_charge3.prediction.ions.json
- 130418_03_SynthMix1_Yphospho1_charge2.prediction.ions.json
- 130418_03_SynthMix1_Yphospho1_charge3.prediction.ions.json

User can also use MS2PIP, Prosit and pDeep2 for MS/MS peak intensity prediction. See [`docs`](docs) for detailed manuals.


### 3. Generate Spectral Libraries
Set the prediction directory as working directory and run one of the following scripts:
- `scripts/ions_to_msp.R`: Convert the JSON files to spectral library files without any optimization.
- `scripts/ions_to_msp_lossST.R`: Convert the JSON files to spectral library files, in which phosphoric acid neutral loss peaks are generated for phosphorylated S/T.
- `scripts/ions_to_msp_budding.R`: Convert the JSON files to spectral library files with "budding" peaks.
- `scripts/ions_to_msp_lossST_budding.R`: Convert the JSON files to spectral library files with neutral loss and "budding" peaks.
```R
setwd("{PATH_TO_DATA}")
source("{PATH_TO_CODE}/scripts/ions_to_msp[_lossST][_budding].R")
```
Spectral library (MSP) files are output in the same directory. 
- 130418_03_SynthMix1_Sphospho1_charge2.prediction[_loss50][_bud100byp2].msp
- 130418_03_SynthMix1_Sphospho1_charge3.prediction[_loss50][_bud100byp2].msp
- 130418_03_SynthMix1_Tphospho1_charge2.prediction[_loss50][_bud100byp2].msp
- 130418_03_SynthMix1_Tphospho1_charge3.prediction[_loss50][_bud100byp2].msp
- 130418_03_SynthMix1_Yphospho1_charge2.prediction[_bud100by2].msp
- 130418_03_SynthMix1_Yphospho1_charge3.prediction[_bud100by2].msp


## Phosphorylation Site Localization

### 1. Perform Spectral Matching
Open `scripts/spectral_matching.R` and change the parameters:
```R
psmFile = '{PATH_TO_PSM_FILE}' # Path to the PSM list file
referencePath = '{PATH_TO_MSP}' # Path to the spectral library (MSP) directory
spectraPath = '{PATH_TO_MGF}' # Path to the MS/MS spectra (MGF) directory
resultFile = 'spectral_matching.result.csv' # Output file name

fragmentTolerance = 50 # m/z tolerance
```

Save and run the script.
```R
source("{PATH_TO_CODE}/scripts/spectral_matching.R")
```

The spectral matching result file is output in the same directory. 
- spectral_matching.result.csv


### 2. Output Result File
Place the PSM file and spectral matching result file in the working directory.

Run `scripts/write_result.R`.
```R
source("{PATH_TO_CODE}/scripts/write_result.R")
```


## Combining Spectral Matching with phosphoRS Results
Spectral matching is a sensitive method that supplements the probability-based approaches (e.g. phosphoRS) for phosphorylation site identification. Users can combine the spectral matching results with phosphoRS results. See [`docs`](docs) for detailed manuals.


## Publications
Yang, Y., Horvatovich, P., Qiao, L. Fragment mass spectrum prediction facilitates site localization of phosphorylation. *J Proteome Res* **20**, 634–644 (2021). https://doi.org/10.1021/acs.jproteome.0c00580

## License
DeepMS2-phospho functions (in [`code/functions`](code/functions)) are distributed under a BSD license. See the LICENSE file for details.
Scripts (in [`code/scripts`](code/scripts)) that require `rjson` and `readr` package are under a GPL-2 licence.

The example data (in [`data/example`](data/example)) are derived from a data set at ProteomeXchange (http://proteomecentral.proteomexchange.org/) with the accession key `PXD000474`. (Suni, V. et al. J. Proteome Res. 2015, 14 (5), 2348-2359.)

## Contacts
Please report any problems directly to the github issue tracker. Also, you can send feedback to liang_qiao@fudan.edu.cn.
