
# DBMLFF: Linear Scaling Machine Learning Force Fields via Electron Density Decomposition for Molecular Electrolytes
﻿
## File Description
﻿
### example_EC
We provide a complete working example for building DBMLFF for EC (ethylene carbonate) molecules and performing MD simulations using DBMLFF. The operational workflow can be found in the user guide: `Tutorial.md`
﻿
This directory contains:
- Step-by-step DBMLFF construction process
- Molecular dynamics simulation examples
﻿
### model
We provide pre-trained DBMLFF models for other molecules that can be directly used for DBMLFF simulations. These include molecules that have appeared in the paper:
﻿
**Available Pre-trained Models:**
- **DMC** (Dimethyl Carbonate)
- **DME** (Methyl ether) 
- **EA** (Ethyl Acetate)
- **EC** (Ethylene Carbonate)
- **EGDME** (Ethylene Glycol Dimethyl Ether)
- **FSI** (Bis(fluorosulfonyl)imide)
- **Li** (Lithium ion)
- **PF6** (Hexafluorophosphate)
- **TTE** (1,1,2,2-Tetrafluoroethyl-2,2,3,3-tetrafluoropropyl ether)
﻿
## Usage of Pre-trained Models
﻿
### Direct Application
```bash
# Copy pre-trained model to your simulation directory
cp model/mol.EC  your_simulation_directory/mol.1
```
﻿
### MSD_Ion_conductor - Ionic Conductor Mean Squared Displacement Toolkit

## Overview
This toolkit provides a complete workflow for loading DBMLFF trajectory files and calculating the Mean Squared Displacement (MSD) of ionic conductors, including correction factors for ion association effects.

## Usage Instructions

### 0. Unwrapped Coordinates
```bash
# Submit job
qsub qsub_unwrapped.sh

# Or run directly
./run_unwrapped.sh
```
**Main file**: `pwmatMOVEMENT_MSD.py`

### 1. Mean Squared Displacement Calculation
```bash
# Submit job
qsub qsub_msd.sh

# Or run directly
./run_msd.sh
```
**Main files**:
- `msd.py` - Calculate MSD
- `linear_msd_slop.py` - Analyze MSD slope

### 2. Ion Association Correction Factor (Alpha)
```bash
# Submit job
qsub qsub_compose.sh
```
**Main files**:
- `msd_compose.py` - Calculate ion association correction factor
- `plot_compose_slop_fit.py` - Plot correction factor fitting results

## Workflow
1. **Trajectory Preprocessing**: Use step 0 to convert DBMLFF trajectories to unwrapped coordinates
2. **MSD Calculation**: Use step 1 to calculate mean squared displacement and obtain diffusion coefficients
3. **Association Effect Correction**: Use step 2 to calculate the correction factor α for ion association effects

## Output Files
Each step will generate corresponding analysis results and chart files. Refer to the documentation of each script for specific output formats.

## Notes
- Ensure correct DBMLFF trajectory file format
- There are dependencies between steps, please execute in order
- Adjust parameters in scripts according to your actual system requirements

## Citation
If you use these pre-trained models in your research, please cite our paper on DBMLFF methodology.
﻿
## Support
For technical support and additional molecular requests, please contact the development team.