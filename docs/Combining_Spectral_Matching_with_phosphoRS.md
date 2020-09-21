# Combining Spectral Matching with phosphoRS Results

### 0. Denpendency
Download phosphoRS-cli, a command line interface of phosphoRS (version 3.1): 
https://github.com/lmsac/phosphoRS-cli

### 1. Prepare a PSM List
Follow the instruction described in *README* to prepare a PSM list file.

### 2. Generate phosphoRS Input File
Place the PSM list file in the same directory as the MGF files to be searched.
Set the PSM list and MGF directory as working directory and run `scripts/phosphoRS/generate_phosphoRS_input_xml`.

You can change the parameters in the script according to your experimental settings.
```
activationTypes = 'HCD'
massTolerance = 0.02
scoreNeutralLoss = FALSE
```

### 3. Run phosphoRS-cli
Follow the instruction described at https://github.com/lmsac/phosphoRS-cli.

### 4. Run Spectral Matching
Follow the instruction described in *README* to generate in silico spectral library and perform spectral matching.

### 5. Output Result File
Place the PSM file, the spectral matching result file, and the phosphoRS result file in the working directory.

Rename the files as 
- *.psm.csv
- *.result.csv
- *.phosphoRS.csv

Run `scripts/phosphoRS/write_result_dot_phosphoRS.R.

